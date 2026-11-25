#!/usr/bin/env python3
"""Assign the reads of one run to lineages (PROTOCOL.md 4.1, read assignment).

Every read is aligned on its own (mates separately) against an index made by make_index.py, with
minimap2 keeping up to 50 secondary alignments. Each read then gets one outcome:

  counted for a lineage   its best alignments (those with the highest score) all pass the filter,
                          at least 95% identity over at least 80% of the read; all fall in that
                          lineage; none overlaps a shared stretch listed in references/overlaps.tsv
                          (from v0.2 also the organelle stretches in bacterial genomes, PROTOCOL.md
                          section 10, entry 10; these count under the same reason); and there are
                          fewer than 51 of them
  counted for none        it has an alignment that passes the filter, but one of the conditions
                          above fails (each reason is counted)
  routed                  no alignment passes the filter: the read goes on to SingleM and SSUplex

The fungal read count for E2 is the number of reads whose best alignments all pass the filter,
avoid the shared stretches and fall in host fungi (any of them, ties among host fungi included).
Settings the protocol leaves open are in PROTOCOL.md section 10, entry 4.

Writes in --out:
  counts.tsv          total reads, reads per lineage, the fungal read count for E2, and the reads
                      counted for none or routed, by reason
  routed.fastq.gz     the routed reads as they were in the input (routed.fasta.gz for FASTA input)
  spread.tsv          the spread check (v0.2; section 10, entry 10): for each bacterial genome with reads
                      counted for its lineage in the first input file (mate 1 of a paired run), the
                      reads, their distinct start positions S, the 1-kb windows they hit, the genome's
                      windows G, the windows S starts would hit if spread evenly, E = G(1 - e^(-S/G)),
                      and the ratio of windows hit to E. A read whose best alignments fall on several
                      genomes of its lineage counts for each.
  spread_bits.json.gz the start positions and windows behind spread.tsv, as bit sets (v0.3; section 10,
                      entry 13), so that the spread of two read sets taken together can be worked out
                      exactly (combined_spread), as T3 does for its spikes
  weighted_counts.tsv only with --weights (v0.4; section 10, entry 14): every read's header carries a
                      weight (w:f:<number>), and this file has the rows of counts.tsv with the weights
                      added up instead of the reads counted; T4's pieces of Logan unitigs use it
  assign.json         versions, command, index, input files and counts

Usage (paired runs: give both files; each mate is counted on its own):
  conda activate phase0-reads
  python scripts/assign_reads.py --index ~/phase0-data/index/full --platform short \\
      --reads SRR123_1.fastq.gz SRR123_2.fastq.gz --out ~/phase0-data/runs/SRR123/full --threads 16
"""

import argparse
import base64
import csv
import datetime as dt
import gzip
import hashlib
import json
import math
import os
import subprocess
import sys
import threading
from array import array
from collections import Counter
from pathlib import Path

__version__ = "0.4"

MINIMAP2_ARGS = {   # PROTOCOL.md 4.1; -c, --frag=no and --paf-no-hit: section 10, entry 4
    "short": ["-x", "sr", "--secondary=yes", "-N", "50", "-c", "--frag=no", "--paf-no-hit"],
    "hifi": ["-x", "map-hifi", "--secondary=yes", "-N", "50", "-c", "--paf-no-hit"],
}
INDEX_FILE = {"short": "refs.sr.mmi", "hifi": "refs.hifi.mmi"}
MAX_EQUAL = 51                    # with -N 50, minimap2 reports at most 51 alignments per read
MIN_ID_PCT, MIN_COV_PCT = 95, 80
NONE_REASONS = ["best alignment fails the filter", "51 or more equally best alignments",
                "overlaps a shared stretch", "tied across lineages"]
ROUTED_REASONS = ["no alignment", "no alignment passes the filter"]
BACTERIAL_KINDS = ("endosymbiont", "contaminant")
WINDOW = 1000                     # the spread check (section 10, entry 10)
SPREAD_COLS = ["lineage", "genome", "reads", "starts", "windows", "genome_windows", "expected", "ratio"]
ORGANELLE_ID_PCT, ORGANELLE_COV_PCT = 90, 50      # the organelle check on routed reads (entry 10)


# ------------------------------------------------------------------ reads ---

def open_reads(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rb") if magic == b"\x1f\x8b" else open(path, "rb")


def records(path):
    """Yield (header, sequence, quality) as bytes without newlines; quality is None for FASTA."""
    with open_reads(path) as fh:
        it = iter(fh)
        first = next((line for line in it if line.strip()), None)
        if first is None:
            return
        if first.startswith(b"@"):
            h = first
            while h is not None:
                try:
                    s, p, q = next(it), next(it), next(it)
                except StopIteration:
                    raise ValueError(f"{path}: the last FASTQ record is cut short") from None
                s, q = s.rstrip(b"\r\n"), q.rstrip(b"\r\n")
                if not h.startswith(b"@") or not p.startswith(b"+") or len(s) != len(q):
                    raise ValueError(f"{path}: not a 4-line FASTQ record at {h[:60]!r}")
                yield h.rstrip(b"\r\n"), s, q
                h = next((line for line in it if line.strip()), None)
        elif first.startswith(b">"):
            h, chunks = first.rstrip(b"\r\n"), []
            for line in it:
                if line.startswith(b">"):
                    yield h, b"".join(chunks), None
                    h, chunks = line.rstrip(b"\r\n"), []
                else:
                    chunks.append(line.strip())
            yield h, b"".join(chunks), None
        else:
            raise ValueError(f"{path}: neither FASTQ nor FASTA")


def header_weight(h):
    """The weight a read's header carries as w:f:<number> (--weights, entry 14)."""
    i = h.find(b" w:f:")
    if i < 0:
        raise ValueError(f"no w:f: weight in the header {h[:60]!r}")
    w = float(h[i + 5:].split(None, 1)[0])
    if not w >= 0:
        raise ValueError(f"a weight must be 0 or more: {h[:60]!r}")
    return w


def feed(path, sink, state):
    """Write every read to minimap2, renamed to its position in the file (0, 1, 2, ...). With
    state["weights"] a list, each read's weight is added to it before the read is written."""
    n, kind, weights = 0, None, state.get("weights")
    try:
        for h, s, q in records(path):
            if weights is not None:
                weights.append(header_weight(h))
            if kind is None:
                kind = "fasta" if q is None else "fastq"
            if q is None:
                sink.write(b">%d\n%s\n" % (n, s))
            else:
                sink.write(b"@%d\n%s\n+\n%s\n" % (n, s, q))
            n += 1
    except BaseException as e:   # reported by the main thread
        state["error"] = e
    finally:
        state["reads"], state["format"] = n, kind
        try:
            sink.close()
        except OSError:
            pass


class Bitset:
    def __init__(self):
        self.b = bytearray()

    def add(self, i):
        j = i >> 3
        if j >= len(self.b):
            self.b.extend(bytes(max(j + 1 - len(self.b), 1 << 16)))
        self.b[j] |= 1 << (i & 7)

    def __contains__(self, i):
        j = i >> 3
        return j < len(self.b) and bool(self.b[j] >> (i & 7) & 1)

    def count(self):
        return bin(int.from_bytes(self.b, "little")).count("1")


# --------------------------------------------------------------- the rule ---

def parse_paf(stream):
    """Yield (read position, alignments) in output order. An alignment is
    (score, matches, block, qstart, qend, qlen, target, tstart, tend)."""
    current, alns = None, []
    for line in stream:
        f = line.rstrip(b"\n").split(b"\t")
        qi = int(f[0])
        if qi != current:
            if current is not None:
                yield current, alns
            current, alns = qi, []
        if f[5] == b"*":
            continue
        score = None
        for tag in f[12:]:
            if tag.startswith(b"AS:i:"):
                score = int(tag[5:])
                break
        if score is None:
            raise ValueError("a PAF line has no AS tag (minimap2 must run with -c)")
        alns.append((score, int(f[9]), int(f[10]), int(f[2]), int(f[3]), int(f[1]), f[5], int(f[7]), int(f[8])))
    if current is not None:
        yield current, alns


def passes(a):
    """At least 95% identity (PAF col 10 / col 11) over at least 80% of the read (integers, no rounding)."""
    return a[1] * 100 >= MIN_ID_PCT * a[2] and (a[4] - a[3]) * 100 >= MIN_COV_PCT * a[5]


def in_mask(a, mask):
    return any(a[7] < e and s < a[8] for s, e in mask.get(a[6], ()))


def evaluate(alns, lineage_of, host_lineages, mask):
    """Return (outcome, detail, fungal); outcome is 'lineage', 'none' or 'routed'."""
    if not alns:
        return "routed", ROUTED_REASONS[0], False
    if not any(passes(a) for a in alns):
        return "routed", ROUTED_REASONS[1], False
    best = max(a[0] for a in alns)
    top = [a for a in alns if a[0] == best]
    all_pass = all(passes(a) for a in top)
    masked = any(in_mask(a, mask) for a in top)
    lineages = {lineage_of[a[6]] for a in top}
    fungal = all_pass and not masked and lineages <= host_lineages
    if not all_pass:
        return "none", NONE_REASONS[0], fungal
    if len(top) >= MAX_EQUAL:
        return "none", NONE_REASONS[1], fungal
    if masked:
        return "none", NONE_REASONS[2], fungal
    if len(lineages) > 1:
        return "none", NONE_REASONS[3], fungal
    return "lineage", lineages.pop(), fungal


# ------------------------------------------------------------------- main ---

def load_index(index_dir):
    info = json.loads((index_dir / "index.json").read_text())
    lineage_of, order, hosts = {}, [], set()
    with open(index_dir / "sequences.tsv", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lineage_of[r["sequence"].encode("ascii")] = r["lineage"]
            if r["lineage"] not in order:
                order.append(r["lineage"])
            if r["kind"] == "host fungus":
                hosts.add(r["lineage"])
    return info, lineage_of, order, hosts


def load_sequences(index_dir):
    with open(index_dir / "sequences.tsv", newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_mask(path):
    """Stretches whose reads count for none: the shared stretches (both sides) and, from v0.2, the organelle
    stretches in bacterial genomes (the bacterial side only; entry 10)."""
    mask = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["kind"] == "shared_stretch":
                mask.setdefault(r["host_contig"].encode("ascii"), []).append(
                    (int(r["host_start"]), int(r["host_end"])))
            elif r["kind"] != "organelle_stretch":
                continue
            mask.setdefault(r["bacterial_contig"].encode("ascii"), []).append(
                (int(r["bacterial_start"]), int(r["bacterial_end"])))
    return mask


class Spread:
    """The spread check (entry 10): where the reads counted for each bacterial lineage start, per genome."""

    def __init__(self, seq_rows):
        self.where, self.lineage = {}, {}
        self.bp, self.windows = Counter(), Counter()
        for r in seq_rows:
            if r["kind"] not in BACTERIAL_KINDS:
                continue
            acc, n = r["accession"], int(r["length"])
            self.where[r["sequence"].encode("ascii")] = (acc, self.bp[acc], self.windows[acc])
            self.lineage[acc] = r["lineage"]
            self.bp[acc] += n
            self.windows[acc] += -(-n // WINDOW)
        self.lineages = set(self.lineage.values())
        self.starts, self.hit, self.reads = {}, {}, Counter()

    def add(self, alns):
        """One read counted for a bacterial lineage: its best alignments."""
        best = max(a[0] for a in alns)
        seen = set()
        for a in alns:
            if a[0] != best:
                continue
            acc, bp_off, w_off = self.where[a[6]]
            self.starts.setdefault(acc, Bitset()).add(bp_off + a[7])
            self.hit.setdefault(acc, Bitset()).add(w_off + a[7] // WINDOW)
            seen.add(acc)
        for acc in seen:
            self.reads[acc] += 1

    def bits(self):
        """The bit sets behind rows(), per genome, for spread_bits.json.gz (entry 13)."""
        return [{"genome": acc, "lineage": self.lineage[acc], "reads": self.reads[acc],
                 "genome_windows": self.windows[acc],
                 "starts": base64.b64encode(bytes(self.starts[acc].b).rstrip(b"\0")).decode("ascii"),
                 "windows": base64.b64encode(bytes(self.hit[acc].b).rstrip(b"\0")).decode("ascii")}
                for acc in sorted(self.reads, key=lambda x: (self.lineage[x], x))]

    def rows(self):
        out = []
        for acc in sorted(self.reads, key=lambda x: (self.lineage[x], x)):
            s, w, g = self.starts[acc].count(), self.hit[acc].count(), self.windows[acc]
            e = expected_windows(s, g)
            out.append({"lineage": self.lineage[acc], "genome": acc, "reads": self.reads[acc], "starts": s,
                        "windows": w, "genome_windows": g, "expected": f"{e:.1f}",
                        "ratio": f"{w / e:.4f}" if e else ""})
        return out


def expected_windows(starts, genome_windows):
    """Windows that this many distinct start positions would hit if spread evenly: G(1 - e^(-S/G))."""
    g = genome_windows
    return g * (1 - math.exp(-starts / g)) if g else 0.0


def read_bits(path):
    """(index refs_md5, rows) from a spread_bits.json.gz, with the bit sets as integers."""
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        doc = json.load(fh)
    for r in doc["genomes"]:
        r["starts"] = int.from_bytes(base64.b64decode(r["starts"]), "little")
        r["windows"] = int.from_bytes(base64.b64decode(r["windows"]), "little")
    return doc["refs_md5"], doc["genomes"]


def combined_spread(parts):
    """The spread check (entry 10) for read sets taken together (entry 13): {lineage: highest ratio among its
    genomes}, as read_spread gives for one set. parts: read_bits() results, all from the same index."""
    if len({md5 for md5, _ in parts}) > 1:
        raise ValueError("the spread bit sets come from different indexes")
    merged = {}
    for _, rows in parts:
        for r in rows:
            m = merged.setdefault(r["genome"], {"lineage": r["lineage"], "g": r["genome_windows"], "s": 0, "w": 0})
            if (m["lineage"], m["g"]) != (r["lineage"], r["genome_windows"]):
                raise ValueError(f"{r['genome']} differs between the spread bit sets")
            m["s"] |= r["starts"]
            m["w"] |= r["windows"]
    best = {}
    for m in merged.values():
        e = expected_windows(m["s"].bit_count(), m["g"])
        ratio = float(f"{m['w'].bit_count() / e:.4f}") if e else 0.0
        best[m["lineage"]] = max(ratio, best.get(m["lineage"], 0.0))
    return best


def read_spread(path):
    """{lineage: highest ratio among its genomes} from a spread.tsv."""
    best = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            v = float(r["ratio"]) if r["ratio"] else 0.0
            best[r["lineage"]] = max(v, best.get(r["lineage"], 0.0))
    return best


def organelle_names(paf_lines):
    """The reads (names, as str) whose alignments at 90% identity or more (PAF col 10 / col 11) together cover
    at least half of the read (entry 10). paf_lines: bytes lines of minimap2 PAF output."""
    covered = {}
    for line in paf_lines:
        f = line.rstrip(b"\n").split(b"\t")
        if len(f) < 12 or f[5] == b"*" or int(f[9]) * 100 < ORGANELLE_ID_PCT * int(f[10]):
            continue
        covered.setdefault(f[0].decode("ascii"), (int(f[1]), []))[1].append((int(f[2]), int(f[3])))
    out = set()
    for name, (qlen, iv) in covered.items():
        total, end = 0, 0
        for s, e in sorted(iv):
            if e > end:
                total += e - max(s, end)
                end = e
        if total * 100 >= ORGANELLE_COV_PCT * qlen:
            out.add(name)
    return out


def organelle_reads(reads, organelles_fna, platform, minimap2, threads, log_path):
    """The organelle check (entry 10): the names of the reads in the FASTQ or FASTA file `reads` that are
    organelle DNA, by organelle_names, when aligned to the organelle genomes alone with the platform's
    minimap2 settings."""
    cmd = [minimap2, *MINIMAP2_ARGS[platform], "-t", str(threads), str(organelles_fna), str(reads)]
    with open(log_path, "a") as err:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=err)
    if p.returncode != 0:
        raise RuntimeError(f"minimap2 failed on {reads} against the organelle genomes (exit {p.returncode}); "
                           f"see {log_path}")
    return organelle_names(p.stdout.splitlines())


def assign_file(path, cmd, log_path, lineage_of, hosts, mask, tally, spread=None, wtally=None):
    """Map one file and add its reads to tally (and to spread, if given; and their weights to wtally, if
    given). Returns (reads, format, routed positions)."""
    with open(log_path, "a") as err:
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=err, bufsize=1 << 20)
    state = {"weights": array("d")} if wtally is not None else {}       # 8 bytes a read
    feeder = threading.Thread(target=feed, args=(path, proc.stdin, state), daemon=True)
    feeder.start()
    routed, expected = Bitset(), 0

    def count(pos, outcome, detail, fungal):
        if outcome == "lineage":
            tally["lineage"][detail] += 1
        else:
            tally[outcome][detail] += 1
            if outcome == "routed":
                routed.add(pos)
        if fungal:
            tally["fungal"] += 1
        if wtally is not None:
            try:
                w = state["weights"][pos]
            except IndexError:
                raise RuntimeError(f"no weight for read {pos} of {path}: {state.get('error')}") from None
            wtally[outcome][detail] += w
            wtally["total"] += w
            if fungal:
                wtally["fungal"] += w

    try:
        for pos, alns in parse_paf(proc.stdout):
            if pos < expected:
                raise RuntimeError("minimap2 output is not in input order")
            while expected < pos:                 # a read with no output line at all
                count(expected, "routed", ROUTED_REASONS[0], False)
                expected += 1
            try:
                outcome, detail, fungal = evaluate(alns, lineage_of, hosts, mask)
                count(pos, outcome, detail, fungal)
                if spread is not None and outcome == "lineage" and detail in spread.lineages:
                    spread.add(alns)
            except KeyError as e:
                raise RuntimeError(f"minimap2 aligned a read to {e}, which is not in sequences.tsv: "
                                   "the index and its table do not match") from None
            expected = pos + 1
    except BaseException:
        proc.kill()
        feeder.join(timeout=60)
        raise
    feeder.join()
    rc = proc.wait()
    if state.get("error"):
        raise state["error"]
    if rc != 0:
        raise RuntimeError(f"minimap2 failed on {path} (exit {rc}); see {log_path}")
    n = state["reads"]
    if expected > n:
        raise RuntimeError(f"minimap2 reported more reads than {path} holds")
    while expected < n:
        count(expected, "routed", ROUTED_REASONS[0], False)
        expected += 1
    return n, state["format"], routed


def write_routed(inputs, dest):
    """Second pass: copy the routed reads, with their original headers."""
    tmp = dest.with_name(dest.name + ".part")
    with gzip.open(tmp, "wb", compresslevel=1) as out:
        for path, routed in inputs:
            for i, (h, s, q) in enumerate(records(path)):
                if i in routed:
                    out.write(h + b"\n" + s + b"\n" + (b"+\n" + q + b"\n" if q is not None else b""))
    return tmp


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", required=True, help="an index set made by make_index.py (index/full or index/t2)")
    ap.add_argument("--platform", required=True, choices=sorted(MINIMAP2_ARGS))
    ap.add_argument("--reads", required=True, nargs="+", help="FASTQ or FASTA files, gzipped or not")
    ap.add_argument("--out", required=True)
    ap.add_argument("--refdir", default="references", help="holds overlaps.tsv")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--weights", action="store_true",
                    help="every read's header carries a weight (w:f:<number>); also write weighted_counts.tsv")
    args = ap.parse_args(argv)
    index_dir, out = Path(args.index).expanduser(), Path(args.out).expanduser()
    started = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    try:
        info, lineage_of, order, hosts = load_index(index_dir)
        spread = Spread(load_sequences(index_dir))
        overlaps = Path(args.refdir) / "overlaps.tsv"
        mask = load_mask(overlaps)
        mmi = index_dir / INDEX_FILE[args.platform]
        if not mmi.exists():
            raise FileNotFoundError(f"{mmi} is missing")
        for p in args.reads:
            if not Path(p).is_file():
                raise FileNotFoundError(f"{p} is missing")
        version = subprocess.run([args.minimap2, "--version"], capture_output=True, text=True,
                                 check=True).stdout.strip()
    except (OSError, ValueError, KeyError, json.JSONDecodeError, subprocess.CalledProcessError) as e:
        print(f"Cannot start: {e}", file=sys.stderr)
        return 1
    if version != info["minimap2"]:
        print(f"minimap2 here is {version}, but the index was built with {info['minimap2']}. "
              "Use the phase0-reads environment.", file=sys.stderr)
        return 1

    out.mkdir(parents=True, exist_ok=True)
    cmd = [args.minimap2, *MINIMAP2_ARGS[args.platform], "-t", str(args.threads), str(mmi), "-"]
    tally = {"lineage": Counter(), "none": Counter(), "routed": Counter(), "fungal": 0}
    wtally = {"lineage": Counter(), "none": Counter(), "routed": Counter(), "fungal": 0.0, "total": 0.0} \
        if args.weights else None
    inputs, routed_sets, formats = [], [], set()
    try:
        for k, p in enumerate(args.reads):
            print(f"aligning {p}", flush=True)
            n, fmt, routed = assign_file(p, cmd, out / "minimap2.log", lineage_of, hosts, mask, tally,
                                         spread if k == 0 else None, wtally)
            inputs.append({"path": str(Path(p).resolve()), "bytes": Path(p).stat().st_size, "reads": n})
            routed_sets.append((p, routed))
            if fmt:
                formats.add(fmt)
        if len(formats) > 1:
            raise ValueError("the inputs mix FASTQ and FASTA")
        routed_name = "routed.fasta.gz" if formats == {"fasta"} else "routed.fastq.gz"
        print("writing the routed reads", flush=True)
        routed_tmp = write_routed(routed_sets, out / routed_name)
    except (OSError, ValueError, RuntimeError) as e:
        print(f"Read assignment failed, nothing was written to {out}/: {e}", file=sys.stderr)
        return 1

    total = sum(i["reads"] for i in inputs)
    accounted = sum(tally["lineage"].values()) + sum(tally["none"].values()) + sum(tally["routed"].values())
    if accounted != total:
        print(f"Internal error: {accounted} outcomes for {total} reads", file=sys.stderr)
        return 1

    rows = [("total", "reads in input", total)]
    rows += [("lineage", lin, tally["lineage"][lin]) for lin in order]
    rows += [("fungal", "fungal reads for E2", tally["fungal"])]
    rows += [("none", r, tally["none"][r]) for r in NONE_REASONS]
    rows += [("routed", r, tally["routed"][r]) for r in ROUTED_REASONS]
    tmp = out / "counts.tsv.part"
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["section", "name", "reads"])
        w.writerows(rows)
    tmp_w = None
    if wtally is not None:
        wsum = sum(sum(wtally[k].values()) for k in ("lineage", "none", "routed"))
        if abs(wsum - wtally["total"]) > 1e-6 * max(1.0, wtally["total"]):
            print(f"Internal error: weights {wsum} by outcome, {wtally['total']} in all", file=sys.stderr)
            return 1
        wrows = [("total", "reads in input", wtally["total"])]
        wrows += [("lineage", lin, wtally["lineage"][lin]) for lin in order]
        wrows += [("fungal", "fungal reads for E2", wtally["fungal"])]
        wrows += [("none", r, wtally["none"][r]) for r in NONE_REASONS]
        wrows += [("routed", r, wtally["routed"][r]) for r in ROUTED_REASONS]
        tmp_w = out / "weighted_counts.tsv.part"
        with open(tmp_w, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["section", "name", "weight"])
            w.writerows((a, b, f"{v:.4f}") for a, b, v in wrows)
    meta = {"made_by": f"assign_reads.py v{__version__}", "started": started,
            "finished": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
            "platform": args.platform, "minimap2": version,
            "command": " ".join(cmd[:-2] + [f"index/{info['set']}/{INDEX_FILE[args.platform]}", "<reads, renamed>"]),
            "index": {"dir": str(index_dir.resolve()), **{k: info[k] for k in ("set", "made_by", "refs_md5", "bp")}},
            "overlaps_md5": md5_of(overlaps), "inputs": inputs, "weights": bool(args.weights),
            "counts": {f"{s}: {n}": r for s, n, r in rows}}
    meta["spread_from"] = inputs[0]["path"] if inputs else None
    tmp_spread = out / "spread.tsv.part"
    with open(tmp_spread, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=SPREAD_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(spread.rows())
    tmp_bits = out / "spread_bits.json.gz.part"
    with gzip.open(tmp_bits, "wt", encoding="utf-8") as fh:
        json.dump({"made_by": meta["made_by"], "refs_md5": info["refs_md5"], "window": WINDOW,
                   "genomes": spread.bits()}, fh)
    tmp_json = out / "assign.json.part"
    tmp_json.write_text(json.dumps(meta, indent=1) + "\n")
    os.replace(tmp_spread, out / "spread.tsv")
    os.replace(tmp_bits, out / "spread_bits.json.gz")
    os.replace(tmp, out / "counts.tsv")
    os.replace(routed_tmp, out / routed_name)
    if tmp_w is not None:
        os.replace(tmp_w, out / "weighted_counts.tsv")
    os.replace(tmp_json, out / "assign.json")

    print(f"{total} reads: {sum(tally['lineage'].values())} counted for a lineage, "
          f"{sum(tally['none'].values())} for none, {sum(tally['routed'].values())} routed to SingleM and SSUplex")
    for lin in order:
        if tally["lineage"][lin]:
            print(f"  {lin}: {tally['lineage'][lin]}")
    print(f"-> {out}/counts.tsv, {routed_name}, spread.tsv, spread_bits.json.gz, assign.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
