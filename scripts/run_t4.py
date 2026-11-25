#!/usr/bin/env python3
"""Test T4, substrate: T1 and T2 again on Logan unitigs (PROTOCOL.md sections 5 and 6; section 10, entry 14).

For every short-read run in t1/runs.tsv and t2/runs.tsv, Logan's unitigs (v1.2) are downloaded from
https://s3.amazonaws.com/logan-pub/u/<run>/<run>.unitigs.fa.zst (no account needed). A run Logan does not
hold is listed and left out; HiFi runs are left out too. Each unitig is cut into pieces as long as the run's
reads (r: the mean length of the first 100,000 reads of its raw reads). Pieces start r - 30 bases apart, so
each 31-mer is in one piece, and the last piece ends where the unitig ends. Each piece stands for reads: the
unitig's average 31-mer count (ka:f:) x the 31-mers it adds / the 31-mers in a read (r - 30). The pieces then
go through T1 and T2 as reads do, with these weights added up where reads are counted (assign_reads.py v0.4,
--weights): the same t, m and rules. The spread check counts piece starts as it counts read starts. sylph,
which only reports, is not run.

  T1  every T1 run Logan holds, full index; K1 as it would read on Logan (reported)
  T2  every T2 run Logan holds, in both arms: kept in (full index, all of SILVA) and held out (T2 index,
      SILVA without the removal list, placements derived from all of SILVA and checked as in T2); K2 as it
      would read on Logan (reported)

Every level and call with signal on either substrate is compared with raw reads (t1/levels_entry10.tsv,
t1/calls_entry10.tsv, t2/levels.tsv, t2/calls.tsv). Calls of a line or isolate with runs left out are listed
but not compared. The switch (section 6): of the positive runs of stage 0a (stage0a/runs.tsv) whose symbiont
is present in raw reads (T1, or T2's kept-in arm) and that Logan holds, the share whose symbiont is not
present on Logan. Above 10%, the survey runs on raw reads.

Writes to t4/ (commit these): logan.tsv (each run: in Logan or not, read length, unitigs, pieces, weight),
levels.tsv and calls.tsv (on Logan, with test and arm columns), compare.tsv and T4.md. Work files go to
WORKDIR/t4/; a run's download is removed once its pieces are made. A rerun skips what is done.

Needs T1 (entry 10) and T2 done, zstd (only to unpack the downloads), and what T2 needs.

Usage:
  conda activate phase0-reads
  python scripts/run_t4.py --workdir ~/phase0-data --threads 30 --jobs 4
"""

import argparse
import gzip
import http.client
import json
import os
import random
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
import assign_reads as AR  # noqa: E402
import call_lineages as C  # noqa: E402
import get_marker_databases as DB  # noqa: E402
import get_silva as GS  # noqa: E402
import rescore_stage0a_m as R  # noqa: E402
import run_stage0a_m as M2  # noqa: E402
import run_t1 as T  # noqa: E402
import run_t2 as RT  # noqa: E402

__version__ = "0.1"

LOGAN_URL = "https://s3.amazonaws.com/logan-pub/u"          # Logan v1.2 (entry 14)
K = 31                                                        # Logan's k
READ_SAMPLE = 100_000                                         # reads that set a run's read length (entry 14)
SWITCH = 0.10                                                 # section 6
RETRY_WAIT = 30                                               # seconds before a download is tried again
MRE = RT.MRE
KA = re.compile(rb"\sk[am]:f:(\S+)")                          # ka:f:, or km:f: in some large runs (a Logan bug)
HEADERS = {"User-Agent": f"amf-endosym-phase0 run_t4.py/{__version__}"}
EARLIER = [("t1", "place", "batches"), ("t2", "place_full", "batches"), ("t3", "place", "batches")] + \
    list(RT.FULL_CACHES)                                      # placements made with entry 8's command, all of SILVA
ARMS = [("T1", "full"), ("T2", "kept_in"), ("T2", "held_out")]
INDEX_OF = {"full": "full", "kept_in": "full", "held_out": "t2"}
ARM_NAME = {"full": "full references", "kept_in": "kept in", "held_out": "held out"}
LOGAN_COLS = ["run", "test", "unit", "in_logan", "note", "bytes", "last_modified", "unitigs_md5", "read_length",
              "unitigs", "unitig_bp", "pieces", "piece_bp", "weight"]
COMPARE_COLS = ["what", "test", "arm", "run_or_unit", "lineage", "raw", "logan", "raw_n", "logan_n", "same", "note"]


def now():
    return T.now()


# -------------------------------------------------------------------- Logan ---

def logan_check(run, base, cache):
    """Whether Logan holds the run's unitigs (an HTTP HEAD request); cached."""
    if cache.exists():
        return json.loads(cache.read_text())
    url = f"{base}/{run}/{run}.unitigs.fa.zst"
    try:
        req = urllib.request.Request(url, method="HEAD", headers=HEADERS)
        with urllib.request.urlopen(req, timeout=60) as resp:
            info = {"run": run, "url": url, "in_logan": True, "bytes": int(resp.headers.get("Content-Length") or -1),
                    "last_modified": resp.headers.get("Last-Modified", ""),
                    "etag": resp.headers.get("ETag", "").strip('"')}
    except urllib.error.HTTPError as e:
        if e.code not in (403, 404):
            raise RuntimeError(f"{url}: HTTP {e.code}") from None
        info = {"run": run, "url": url, "in_logan": False, "http": e.code}
    except (urllib.error.URLError, http.client.HTTPException, OSError) as e:
        raise RuntimeError(f"{url}: {e}") from None
    info["checked"] = now()
    M2.write_atomic(cache, json.dumps(info, indent=1) + "\n")
    return info


def download(info, dest):
    """The run's unitigs, checked against the size Logan gave; three tries."""
    part, err = dest.with_name(dest.name + ".part"), ""
    for attempt in range(3):
        if attempt:
            time.sleep(RETRY_WAIT)
        try:
            with urllib.request.urlopen(urllib.request.Request(info["url"], headers=HEADERS), timeout=300) as resp, \
                    open(part, "wb") as out:
                shutil.copyfileobj(resp, out, 1 << 22)
            n = part.stat().st_size
            if info["bytes"] < 0 or n == info["bytes"]:
                os.replace(part, dest)
                return
            err = f"{n:,} bytes arrived, {info['bytes']:,} expected"
        except (urllib.error.URLError, http.client.HTTPException, OSError) as e:
            err = str(e) or type(e).__name__
        part.unlink(missing_ok=True)
    raise RuntimeError(f"{info['url']}: the download failed 3 times ({err})")


def raw_reads(run, work):
    """The first file of the run's raw reads, as T1 or stage 0a fetched them."""
    for d in (work / "t1" / "reads" / run, work / "stage0a" / "reads" / run):
        if (d / "fetch.json").exists():
            return d / json.loads((d / "fetch.json").read_text())["files"][0]
    raise RuntimeError("T1 and stage 0a have no raw reads for it, and T4 takes the read length from them")


def read_length(path):
    """The mean length of the first READ_SAMPLE reads, rounded (entry 14)."""
    n = total = 0
    gen = AR.records(path)
    for _, s, _ in gen:
        n, total = n + 1, total + len(s)
        if n >= READ_SAMPLE:
            break
    gen.close()
    if n == 0 or round(total / n) <= K:
        raise RuntimeError(f"{path}: no reads, or reads too short for pieces of 31-mers")
    return round(total / n)


def unitig_pieces(seq, ka, r, k=K):
    """A unitig cut into pieces of r bases (entry 14). Pieces start r - k + 1 apart, so each k-mer is in one of
    them; the last piece ends where the unitig ends and counts only the k-mers it adds. A unitig shorter than r
    is one piece. Each piece's weight: ka x the k-mers it counts / (r - k + 1). Returns [(piece, weight)]."""
    n_k, step = len(seq) - k + 1, r - k + 1
    if n_k <= 0:
        return []
    if len(seq) <= r:
        return [(seq, ka * n_k / step)]
    out, s = [], 0
    while s + r <= len(seq):
        out.append((seq[s:s + r], ka))
        s += step
    if s < n_k:
        out.append((seq[len(seq) - r:], ka * (n_k - s) / step))
    return out


def unitigs(stream):
    """(header, sequence) from a FASTA byte stream."""
    h, chunks = None, []
    for line in stream:
        if line.startswith(b">"):
            if h is not None:
                yield h, b"".join(chunks)
            h, chunks = line.rstrip(b"\r\n"), []
        else:
            chunks.append(line.rstrip(b"\r\n"))
    if h is not None:
        yield h, b"".join(chunks)


def make_pieces(run, zst, pdir, r, args):
    """The run's pieces as FASTQ, each weight in its header (w:f:), unless done."""
    done = pdir / "pieces.json"
    if done.exists():
        return json.loads(done.read_text())
    pdir.mkdir(parents=True, exist_ok=True)
    out = pdir / f"{run}_pieces.fastq.gz"
    part = out.with_name(out.name + ".part")
    st, wsum = Counter(), 0.0
    with open(pdir / "zstd.log", "w") as log:
        proc = subprocess.Popen([args.zstd, "-dc", str(zst)], stdout=subprocess.PIPE, stderr=log)
        try:
            with gzip.open(part, "wb", compresslevel=1) as fh:
                for h, s in unitigs(proc.stdout):
                    hit = KA.search(h)
                    ka = float(hit.group(1)) if hit else -1.0
                    if not ka > 0:
                        raise RuntimeError(f"{run}: a unitig header has no ka:f: count above 0: {h[:80]!r}")
                    uid = h[1:].split(None, 1)[0]
                    st["unitigs"] += 1
                    st["unitig_bp"] += len(s)
                    ps = unitig_pieces(s, ka, r)
                    st["unitigs_under_31_bp"] += not ps
                    for i, (p, w) in enumerate(ps):
                        ws = b"%.6g" % w
                        fh.write(b"@%s.%d w:f:%s\n%s\n+\n%s\n" % (uid, i, ws, p, b"I" * len(p)))
                        st["pieces"] += 1
                        st["piece_bp"] += len(p)
                        wsum += float(ws)
        except BaseException:
            proc.kill()
            proc.wait()
            part.unlink(missing_ok=True)
            raise
        rc = proc.wait()
    if rc != 0:
        part.unlink(missing_ok=True)
        zst.unlink(missing_ok=True)
        raise RuntimeError(f"zstd could not unpack the download (exit {rc}; see {pdir / 'zstd.log'}). The download "
                           "was removed, so a rerun gets it again")
    if st["pieces"] == 0:
        part.unlink(missing_ok=True)
        raise RuntimeError("the download holds no unitigs")
    os.replace(part, out)
    meta = {"made_by": f"run_t4.py v{__version__}", "made": now(), "run": run, "read_length": r, "k": K,
            "unitigs_md5": AR.md5_of(zst), "unitigs_bytes": zst.stat().st_size, **st, "weight": round(wsum, 4),
            "file": out.name}
    M2.write_atomic(done, json.dumps(meta, indent=1) + "\n")
    return meta


# ------------------------------------------------------------ the pipeline ---

def assign_w(run, piece_file, out, index, args):
    """assign_reads.py --weights on the pieces, unless done."""
    if (out / "weighted_counts.tsv").exists() and (out / "assign.json").exists():
        return
    cmd = [sys.executable, str(SCRIPTS / "assign_reads.py"), "--index", str(index), "--platform", "short",
           "--reads", str(piece_file), "--out", str(out), "--refdir", args.refdir, "--threads", str(args.threads),
           "--minimap2", args.minimap2, "--weights"]
    out.mkdir(parents=True, exist_ok=True)
    log = out.parent / f"{run}.assign.log"
    with open(log, "w") as fh:
        p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    if p.returncode != 0 or not (out / "weighted_counts.tsv").exists():
        raise RuntimeError(f"assign_reads.py failed; see {log}")


def weighted_counts(adir):
    return {(x["section"], x["name"]): float(x["weight"]) for x in M2.read_tsv(adir / "weighted_counts.tsv")}


def found_weights(rdir, found):
    """{(mate, position): weight} of the found pieces, from the headers in the route's found-read files."""
    want = defaultdict(set)
    for (mate, _), (pos, _, _) in found.items():
        want[mate].add(pos)
    out = {}
    for mate, positions in want.items():
        for h, _, _ in AR.records(rdir / f"{mate}.found.fastq.gz"):
            pos = int(M2.first_word(h))
            if pos in positions:
                out[(mate, pos)] = AR.header_weight(h)
    return out


def run_data(r, wc, adir, ctx, route=None, mre=None):
    """One run on Logan, as T1 and T2 make it (run_t2.arm_data), with each piece's weight where a read counts
    1. route: (found, weights, SSU pieces, SingleM taxonomies, hits, taxonomy, organelle DNA, split reads), or
    None for a run without the route."""
    d = {"reads": wc[("total", "reads in input")], "fungal": wc[("fungal", "fungal reads for E2")],
         "mapped": {lin: wc.get(("lineage", lin), 0.0) for lin in ctx["bacterial"]},
         "spread": AR.read_spread(adir / "spread.tsv") if r["kind"] == "DNA" else None,
         "route": None, "split": Counter(), "placements": defaultdict(Counter), "stats": Counter(), "sylph": {},
         "mre": Counter(), "mre_stats": Counter()}
    if route is None:
        return d
    found, wts, pieces, sm_tax, hits, tax, organelle, split = route
    d["route"], st, ms = Counter(), d["stats"], d["mre_stats"]
    for key, (pos, a, b) in found.items():
        w, is_mre = wts[(key[0], pos)], mre is not None and key in mre
        st["found"] += w
        ms["found"] += w * is_mre
        piece = pieces[key] if a else None
        ssu = C.ssu_label(hits.get(R.seq_id(piece), []), tax, len(piece)) if a else None
        sml = C.singlem_label(sm_tax[key]) if b else None
        label = C.read_label(ssu, sml)
        if label is None:
            continue
        if (key[0], pos) in organelle:
            st["organelle_dna"] += w
            ms["organelle_dna"] += w * is_mre
            continue
        if (key[0], pos) in split:
            lin = split[(key[0], pos)]
            st["split_bact" if lin else "split_host"] += w
            ms["split"] += w * is_mre
            if lin:
                d["split"][lin] += w
            continue
        d["route"][label] += w
        if is_mre:
            d["mre"][label] += w
            ms["labelled"] += w
        if label[0] == "novel":
            lname = C.NOVEL_SSU if label[1] == "SSU" else C.NOVEL_SINGLEM
            d["placements"][lname][(ssu if ssu and ssu[0] == "kept" else sml)[2]] += w
    for c in d["placements"].values():           # shown in calls.tsv only
        for p, n in c.items():
            c[p] = int(n) if n == int(n) else round(n, 1)
    return d


# ---------------------------------------------------------------- compare ---

def num(v):
    return "" if v in ("", None) else f"{float(v):.1f}"


def fmt(v, spec=".4g"):
    return "" if v in ("", None) else format(v, spec)


def compare(test, arm, raw, logan, have, left, left_projects):
    """compare.tsv rows for one test and arm: levels (runs Logan holds) and calls with signal on either
    substrate. Calls of units with runs left out are not compared."""
    (raw_levels, raw_calls), (levels, calls) = raw, logan
    rows = []
    lv = {(x["run"], x["lineage"]): x for x in levels}
    raw_lv = {(x["run"], x["lineage"]): x for x in raw_levels if x["run"] in have}
    for key in sorted(set(lv) | set(raw_lv)):
        a, b = raw_lv.get(key, {}), lv.get(key, {})
        la, lb = a.get("level", "none"), b.get("level", "none")
        if not {la, lb} & {"present", "low"}:
            continue
        bp = (b or a)["bioproject"]
        rows.append({"what": "level", "test": test, "arm": arm, "run_or_unit": key[0], "lineage": key[1],
                     "raw": la, "logan": lb, "raw_n": num(a.get("n")), "logan_n": num(b.get("n")),
                     "same": "yes" if la == lb else "no",
                     "note": "its BioProject has runs left out (cross-talk)" if bp in left_projects else ""})
    cl = {(c["unit"], c["lineage"]): c for c in calls}
    raw_cl = {(c["unit"], c["lineage"]): c for c in raw_calls}
    for key in sorted(set(cl) | set(raw_cl)):
        a, b = raw_cl.get(key, {}), cl.get(key, {})
        if not (a.get("signal") == "yes" or b.get("signal")):
            continue
        ca, cb, out = a.get("call", "none"), b.get("call", "none"), left.get(key[0])
        rows.append({"what": "call", "test": test, "arm": arm, "run_or_unit": key[0], "lineage": key[1],
                     "raw": ca, "logan": cb, "raw_n": "", "logan_n": "",
                     "same": "n/a" if out else ("yes" if ca == cb else "no"),
                     "note": "not compared, runs left out: " + ", ".join(f"{x} ({why})" for x, why in out)
                     if out else ""})
    return rows


def switch(positives, by_run, raw_by, logan_by, have, why):
    """The Logan switch (section 6, entry 14). Returns (rows, share missed or None, positives left out)."""
    rows, left = [], []
    for p in positives:
        run, lin = p["run"], p["symbiont_lineage"]
        src = ("T1", "full") if by_run[run]["test"] == "T1" else ("T2", "kept_in")
        raw = raw_by[src].get((run, lin), {}).get("level", "none")
        if run not in have:
            left.append((run, why[run], raw))
            continue
        logan = logan_by[src].get((run, lin), {}).get("level", "none")
        rows.append({"run": run, "test": by_run[run]["test"], "lineage": lin, "raw": raw, "logan": logan,
                     "counted": raw == "present", "missed": raw == "present" and logan != "present"})
    counted = [x for x in rows if x["counted"]]
    return rows, (sum(x["missed"] for x in counted) / len(counted) if counted else None), left


# ------------------------------------------------------------------ outputs ---

def table(cols, rows):
    return "\t".join(cols) + "\n" + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows)


def write_outputs(t4_dir, logan_rows, built, cmp_rows, sw, k1_reasons, k2, meta):
    def level_row(test, arm, x):
        v = {**x, "test": test, "arm": arm, "n": num(x["n"]), "reads": num(x["reads"]),
             "fungal_reads": num(x["fungal"]), "N": fmt(x["N"]), "fraction": fmt(x["fraction"]),
             "informative": T.yn(x["informative"]), "cross_talk": T.yn(x["cross_talk"]),
             "spread": fmt(x["spread"], ".2f"), "spread_cap": T.yn(x["spread_cap"]),
             "split_reads": num(x["split_reads"]), "sylph_ani": ""}
        return "\t".join(str(v[c]) for c in ["test", "arm"] + T.LEVEL_COLS) + "\n"

    def call_row(test, arm, c):
        v = {**c, "test": test, "arm": arm}
        return "\t".join(T.yn(v[k]) if k in ("signal", "e1", "e2", "e5") else str(v[k])
                         for k in ["test", "arm"] + T.CALL_COLS) + "\n"

    out = {"logan.tsv": table(LOGAN_COLS, logan_rows),
           "levels.tsv": "\t".join(["test", "arm"] + T.LEVEL_COLS) + "\n" +
           "".join(level_row(t, a, x) for (t, a), b in built.items() for x in b[0]),
           "calls.tsv": "\t".join(["test", "arm"] + T.CALL_COLS) + "\n" +
           "".join(call_row(t, a, c) for (t, a), b in built.items() for c in b[1]),
           "compare.tsv": table(COMPARE_COLS, cmp_rows)}
    rows, share, left = sw
    n_in = sum(x["in_logan"] == "yes" for x in logan_rows)
    counted, missed = sum(x["counted"] for x in rows), sum(x["missed"] for x in rows)
    if share is None:
        verdict = "No positive run is both present in raw reads and in Logan, so the switch cannot be worked out."
    else:
        verdict = (f"Logan misses {missed} of the {counted} positive runs that are present in raw reads and that "
                   f"Logan holds ({share:.0%}): " + ("more than 10%, so **the survey runs on raw reads**."
                                                    if share > SWITCH else
                                                    "10% or less, so **the survey can run on Logan unitigs**."))
    md = ["# Test T4: Logan unitigs against raw reads", "",
          f"Made by run_t4.py v{__version__} (run_t1.py v{T.__version__}, run_t2.py v{RT.__version__}, "
          f"call_lineages.py v{C.__version__}) on {meta['made']}, with {meta['assign']} (--weights), minimap2 "
          f"{meta['minimap2']}, {meta['tools']}, VSEARCH {meta['vsearch']} and {meta['zstd']}. Unitigs: Logan v1.2 "
          f"({meta['logan_url']}).",
          f"Rule: PROTOCOL.md 4.1, section 5 (T4) and section 6, with section 10, entries 4 to 14. t = {meta['t']:.4g} "
          f"and m = {meta['m']:.4g} (entry 10). Each unitig is cut into pieces as long as the run's reads; each piece "
          "stands for its unitig's average 31-mer count x the 31-mers it adds / the 31-mers in a read, and these "
          "weights are added up where reads are counted.", "",
          "## The switch (section 6)", "", verdict, "",
          "| Run | Test | Symbiont | Raw reads | Logan | Counted | Missed |", "|---|---|---|---|---|---|---|"]
    md += [f"| {x['run']} | {x['test']} | {x['lineage']} | {x['raw']} | {x['logan']} | {T.yn(x['counted'])} | "
           f"{T.yn(x['missed'])} |" for x in rows]
    md += ["", "Counted: present in raw reads. Missed: counted, and not present on Logan."]
    if left:
        md += ["", "Positive runs left out: " + ", ".join(f"{r} ({why}; {raw} in raw reads)" for r, why, raw in left)
               + "."]
    md += ["", "## Runs", "", f"Logan holds {n_in} of the {len(logan_rows)} runs (`logan.tsv`).", ""]
    out_rows = [x for x in logan_rows if x["in_logan"] != "yes" or x["note"]]
    if out_rows:
        md += ["| Run | Test | Line or isolate | Left out because |", "|---|---|---|---|"]
        md += [f"| {x['run']} | {x['test']} | {x['unit']} | {x['note']} |" for x in out_rows]
    for test, arm in [a for a in ARMS if a in built]:
        mine = [x for x in cmp_rows if (x["test"], x["arm"]) == (test, arm)]
        md += ["", f"## {test}, {ARM_NAME[arm]}", ""]
        lv = [x for x in mine if x["what"] == "level"]
        cl = [x for x in mine if x["what"] == "call" and x["same"] != "n/a"]
        nc = [x for x in mine if x["what"] == "call" and x["same"] == "n/a"]
        md.append(f"Levels: {sum(x['same'] == 'yes' for x in lv)} of {len(lv)} the same (runs Logan holds, "
                  "lineages present or low on either substrate). "
                  f"Calls: {sum(x['same'] == 'yes' for x in cl)} of {len(cl)} the same (calls with signal on either "
                  "substrate)" + (f"; {len(nc)} more not compared, as their line or isolate has runs left out."
                                  if nc else "."))
        diff = [x for x in lv + cl if x["same"] == "no"]
        if diff:
            md += ["", "| Level or call | Run or unit | Lineage | Raw reads | Logan | Note |",
                   "|---|---|---|---|---|---|"]
            md += [f"| {x['what']} | {x['run_or_unit']} | {x['lineage']} | {x['raw']}" +
                   (f" ({x['raw_n']})" if x["raw_n"] else "") + f" | {x['logan']}" +
                   (f" ({x['logan_n']})" if x["logan_n"] else "") + f" | {x['note']} |" for x in diff[:100]]
            if len(diff) > 100:
                md.append(f"| ... | {len(diff) - 100} more in `compare.tsv` | | | | |")
    md += ["", "## K1 and K2 on Logan (reported; the kill criteria are judged on raw reads)", "",
           "K1: " + ("would pass." if not k1_reasons else "would fail: " + "; ".join(k1_reasons) + ".")]
    if k2 is not None:
        k2_rows, reasons2 = k2
        md += ["", "K2: " + ("would pass." if not reasons2 else "would fail: " + "; ".join(reasons2) + "."), "",
               "| Isolate | MRE (full index) | Found | Labelled | Lineage that holds most | Its MRE | Call (held out) "
               "| K2 |", "|---|---|---|---|---|---|---|---|"]
        md += [f"| {x['unit']} | {num(x['step1'])} | {num(x['found'])} | {num(x['labelled'])} | "
               f"{x['lineage'] or 'none'} | {num(x['lineage_mre'])} | {x['call']} | "
               f"{'passes' if x['passed'] else 'FAILS'} |" for x in k2_rows]
    md += ["", "On Logan, n and the read counts are weights: the reads the pieces stand for, so they are on the "
               "scale of reads. A run's unitigs hold every 31-mer seen more than twice in all its reads (Logan's "
               "documentation); the raw reads are the first 20,000,000 spots. sylph is not run on pieces. "
               "`compare.tsv` has every comparison; `levels.tsv` and `calls.tsv` have the results on Logan."]
    out["T4.md"] = "\n".join(md) + "\n"
    for name, text in out.items():
        (t4_dir / (name + ".part")).write_text(text, encoding="utf-8")
    for name in out:
        os.replace(t4_dir / (name + ".part"), t4_dir / name)


# --------------------------------------------------------------------- main ---

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--refdir", default="references")
    ap.add_argument("--t1-dir", default="t1")
    ap.add_argument("--t2-dir", default="t2")
    ap.add_argument("--t4-dir", default="t4", help="the outputs to commit go here")
    ap.add_argument("--stage0a-dir", default="stage0a")
    ap.add_argument("--work-name", default="t4", help="work files go to WORKDIR/<this> [t4]")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--jobs", type=int, default=3, help="runs downloaded and cut at the same time [3]")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--ssuplex", default="ssuplex")
    ap.add_argument("--conda", default=os.environ.get("CONDA_EXE") or shutil.which("conda"))
    ap.add_argument("--singlem-env", default="phase0-markers")
    ap.add_argument("--vsearch", default="vsearch")
    ap.add_argument("--zstd", default="zstd")
    ap.add_argument("--logan-url", default=LOGAN_URL, help=f"where Logan's unitigs are [{LOGAN_URL}]")
    args = ap.parse_args(argv)
    work = Path(args.workdir).expanduser()
    wd = work / args.work_name
    t1_dir, t2_dir, t4_dir, stage_dir = Path(args.t1_dir), Path(args.t2_dir), Path(args.t4_dir), Path(args.stage0a_dir)

    try:
        for d, names, script in ((t1_dir, ("levels_entry10.tsv", "calls_entry10.tsv"), "run_t1.py"),
                                 (t2_dir, ("levels.tsv", "calls.tsv"), "run_t2.py")):
            for name in names:
                if not (d / name).exists():
                    raise RuntimeError(f"{d / name} is missing: run {script} first (and commit its outputs)")
        t1_runs, t2_runs = T.read_runs(t1_dir / "runs.tsv"), T.read_runs(t2_dir / "runs.tsv")
        runs = [dict(r, test="T1") for r in t1_runs] + [dict(r, test="T2") for r in t2_runs]
        by_run = {r["run"]: r for r in runs}
        if len(by_run) != len(runs):
            raise RuntimeError("a run is in both t1/runs.tsv and t2/runs.tsv")
        raw = {("T1", "full"): (M2.read_tsv(t1_dir / "levels_entry10.tsv"), M2.read_tsv(t1_dir / "calls_entry10.tsv"))}
        for arm in ("kept_in", "held_out"):
            raw[("T2", arm)] = ([x for x in M2.read_tsv(t2_dir / "levels.tsv") if x["arm"] == arm],
                                [x for x in M2.read_tsv(t2_dir / "calls.tsv") if x["arm"] == arm])
        positives = M2.read_tsv(stage_dir / "runs.tsv")
        wrong = [p["run"] for p in positives if p["run"] not in by_run or by_run[p["run"]]["test"] != p["test"]]
        if wrong:
            raise RuntimeError(f"positive runs of {stage_dir}/runs.tsv not in the T1 or T2 run list they name: "
                               + ", ".join(wrong))
        t, m = T.stage_values(stage_dir)
        overlaps = Path(args.refdir) / "overlaps.tsv"
        mm2 = subprocess.run([args.minimap2, "--version"], capture_output=True, text=True, check=True).stdout.strip()
        lin_rows = M2.read_tsv(Path(args.refdir) / "lineages.tsv")
        ctx = {ix: RT.arm_ctx(work / "index" / ix, ix, overlaps, mm2, lin_rows) for ix in ("full", "t2")}
        if MRE in ctx["t2"]["bacterial"] or MRE not in ctx["full"]["bacterial"]:
            raise RuntimeError("the T2 index must leave out every CaMg/MRE genome, and the full index hold them")
        silva = json.loads((work / "db" / "silva" / "manifest.json").read_text())
        silva_t2 = json.loads((work / "db" / "silva_t2" / "manifest.json").read_text())
        vsearch = GS.vsearch_version(args.vsearch)
        if vsearch != silva["vsearch"] or vsearch != silva_t2["vsearch"]:
            raise RuntimeError(f"VSEARCH {vsearch} is not the version get_silva.py and make_t2_silva.py used")
        removed_text = (t2_dir / "silva_removed.tsv").read_text(encoding="utf-8")
        if RT.md5_text(removed_text) != silva_t2["removed_list_md5"] or silva_t2["source_md5"] != silva["md5"]:
            raise RuntimeError(f"{work / 'db' / 'silva_t2'} was not made from {t2_dir}/silva_removed.tsv and this "
                               "SILVA; run make_t2_silva.py again")
        if AR.md5_of(silva_t2["fasta"]) != silva_t2["fasta_md5"]:
            raise RuntimeError(f"{silva_t2['fasta']} is not the FASTA make_t2_silva.py wrote; run it again")
        removed = {line.split("\t", 1)[0] for line in removed_text.splitlines()[1:] if line}
        db = json.loads((work / "db" / "manifest.json").read_text())
        versions = DB.check_versions(args)
        if versions != db["versions"]:
            raise RuntimeError(f"the tools here ({versions}) are not those get_marker_databases.py checked "
                               f"({db['versions']})")
        zstd_v = DB.tool_output([args.zstd, "--version"]).strip()
        if not re.search(r"zstandard|zstd", zstd_v, re.I):
            raise RuntimeError(f"'{args.zstd} --version' said {zstd_v!r}")
        zstd_v = zstd_v.strip("* \n").split(", by")[0]
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.CalledProcessError) as e:
        msg = str(e)
        if "cannot run" in msg and args.zstd in msg:
            msg += " (T4 needs zstd to unpack Logan's files; see README.md, Test T4)"
        print(f"Cannot start: {msg}", file=sys.stderr)
        return 1
    for sub in ("logan", "unitigs", "pieces", "assign", "assign_t2", "route", "route_t2", "organelle", "split"):
        (wd / sub).mkdir(parents=True, exist_ok=True)
    failed, why = [], {r["run"]: "HiFi" for r in runs if r["platform"] != "short"}

    # 1. Logan: which runs it holds; download; pieces
    short = [r for r in runs if r["platform"] == "short"]

    def one(r):
        run = r["run"]
        try:
            info = logan_check(run, args.logan_url, wd / "logan" / f"{run}.json")
            pdir, zst = wd / "pieces" / run, wd / "unitigs" / f"{run}.unitigs.fa.zst"
            if info["in_logan"] and not (pdir / "pieces.json").exists():
                rl = read_length(raw_reads(run, work))
                if not zst.exists():
                    print(f"downloading {run} ({info['bytes']:,} bytes)", flush=True)
                    download(info, zst)
                print(f"cutting the unitigs of {run} into pieces of {rl} bases", flush=True)
                make_pieces(run, zst, pdir, rl, args)
            zst.unlink(missing_ok=True)
            return None
        except Exception as e:          # noqa: BLE001  reported, and nothing is written
            return f"{run}: {e}" if isinstance(e, (OSError, RuntimeError)) else f"{run}: {type(e).__name__}: {e}"

    print(f"Logan: {len(short)} short-read runs, {args.jobs} at a time", flush=True)
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        for err in pool.map(one, short):
            if err:
                failed.append(err)
                print(f"FAILED: {err}", flush=True)
    have = set()
    for r in short:
        info_p = wd / "logan" / f"{r['run']}.json"
        if (wd / "pieces" / r["run"] / "pieces.json").exists():
            have.add(r["run"])
        elif info_p.exists() and not json.loads(info_p.read_text())["in_logan"]:
            why[r["run"]] = "not in Logan"

    # 2. assignment and route: the full index for every run (T1, and T2 kept in), the T2 index for T2 runs
    wcs, adirs, rdirs = {}, {}, {}
    for r in short:
        run = r["run"]
        if run not in have:
            continue
        pm = json.loads((wd / "pieces" / run / "pieces.json").read_text())
        pf = wd / "pieces" / run / pm["file"]
        for ix in ["full"] + (["t2"] if r["test"] == "T2" else []):
            key = (run, ix)
            adir = wd / ("assign" if ix == "full" else "assign_t2") / run
            rdir = wd / ("route" if ix == "full" else "route_t2") / run
            try:
                if not (adir / "weighted_counts.tsv").exists():
                    print(f"assigning the pieces of {run} ({ix} index)", flush=True)
                assign_w(run, pf, adir, ctx[ix]["index"], args)
                files, _ = T.assigned(adir, "short", ctx[ix])
                wcs[key], adirs[key] = weighted_counts(adir), adir
                if r["route"] == "none":
                    continue
                if not (rdir / "route.json").exists():
                    print(f"route: {run} ({ix} index)", flush=True)
                if ix == "full":
                    T.route_new(r, files, rdir, args, db, ctx[ix])
                else:
                    T.route_view(r, files, rdirs[(run, "full")], rdir, args, db, ctx[ix])
                rdirs[key] = rdir
            except (OSError, KeyError, ValueError, RuntimeError, subprocess.CalledProcessError) as e:
                failed.append(f"{run} ({ix} index): {e}")
                print(f"FAILED: {run} ({ix} index): {e}", flush=True)
                break
    if failed:
        print(f"\n{len(failed)} run(s) not done, so nothing was written to {t4_dir}/. Run the same command again to "
              "finish them:\n" + "\n".join(failed), file=sys.stderr)
        return 1

    # 3. placement, data, levels and calls
    try:
        found, pieces, sm_tax, wts = {}, {}, {}, {}
        for key, rdir in rdirs.items():
            r = by_run[key[0]]
            found[key] = T.routed_found(r, rdir)
            pieces[key] = R.ssu_sequences(rdir, {k for k, v in found[key].items() if v[1]})
            sm_tax[key] = T.singlem_taxonomies(rdir, {k for k, v in found[key].items() if v[2]})
            wts[key] = found_weights(rdir, found[key])
        q_full, q_t2 = {}, {}
        for key, ps in pieces.items():
            for s in ps.values():
                (q_t2 if key[1] == "t2" else q_full).setdefault(R.seq_id(s), s)
        both = {**q_full, **q_t2}
        earlier = [work.joinpath(*p) for p in EARLIER]
        print(f"SSU pieces against all of SILVA: {len(both):,}", flush=True)
        full = RT.placements(both, wd / "place_full" / "batches", silva["fasta"], args,
                             wd / "place_full" / "vsearch.log", earlier)
        t2_hits, todo = RT.derive({q: full[q] for q in q_t2}, removed)
        searched = RT.placements({q: both[q] for q in todo}, wd / "place_t2" / "batches", silva_t2["fasta"], args,
                                 wd / "place_t2" / "vsearch.log")
        sample = sorted(t2_hits)
        sample = sorted(random.Random(1).sample(sample, RT.CHECK_SAMPLE)) if len(sample) > RT.CHECK_SAMPLE else sample
        direct = RT.placements({q: both[q] for q in sample}, wd / "place_check" / "batches", silva_t2["fasta"], args,
                               wd / "place_check" / "vsearch.log")
        differ = [q for q in sample if not RT.same(direct[q], t2_hits[q])]
        if differ:
            raise RuntimeError(f"{len(differ)} of {len(sample)} derived placements differ from a direct search against "
                               f"the T2 database (for example {differ[0]}); the derivation cannot be used")
        t2_hits.update(searched)
        hits = {"full": {q: full[q] for q in q_full}, "t2": t2_hits}
        tax = {}
        for ix, manifest, qs in (("full", silva, q_full), ("t2", silva_t2, q_t2)):
            targets = {t_ for q in qs for t_, _ in hits[ix].get(q, [])}
            if ix == "t2" and targets & removed:
                raise RuntimeError("a held-out placement hits a removed sequence")
            tax[ix] = {k: [x.strip() for x in v.split(";")] for k, v in
                       R.load_taxonomy(manifest["taxonomy"], targets).items()}

        built, data = {}, {}
        for test, arm in ARMS:
            ix = INDEX_OF[arm]
            mine = [r for r in (t1_runs if test == "T1" else t2_runs) if r["run"] in have]
            if not mine:
                continue
            d = {}
            for r in mine:
                key, route = (r["run"], ix), None
                if key in rdirs:
                    print(f"{test}, {ARM_NAME[arm]}: route pieces of {r['run']}", flush=True)
                    org = T.organelle_check(found[key], rdirs[key], wd / "organelle" / f"{r['run']}_{ix}.json", args,
                                            ctx[ix])
                    spl = T.split_check(found[key], rdirs[key], wd / "split" / f"{r['run']}_{ix}.json", args, ctx[ix])
                    route = (found[key], wts[key], pieces[key], sm_tax[key], hits[ix], tax[ix], org, spl)
                mre = RT.mre_reads(rdirs[(r["run"], "full")]) if arm == "held_out" else None
                d[r["run"]] = run_data(r, wcs[key], adirs[key], ctx[ix], route, mre)
            data[(test, arm)] = d
            built[(test, arm)] = (mine, *T.build(mine, d, t, m, ctx[ix]))
        if ("T1", "full") in built:
            mine, levels1, calls1, _ = built[("T1", "full")]
            k1_reasons = T.k1(calls1, levels1, mine)
        else:
            k1_reasons = ["no T1 run is in Logan"]
        k2 = None
        if ("T2", "held_out") in built:
            mine, _, calls_h, joiner_h = built[("T2", "held_out")]
            k2 = RT.k2(mine, data[("T2", "held_out")], calls_h, joiner_h,
                       {r["run"]: wcs[(r["run"], "full")] for r in mine})
    except (OSError, KeyError, ValueError, RuntimeError, subprocess.CalledProcessError) as e:
        print(f"Nothing was written to {t4_dir}/: {e}", file=sys.stderr)
        return 1

    # 4. compared with raw reads; the switch
    left, left_projects = defaultdict(list), set()
    for r in runs:
        if r["run"] not in have:
            left[r["unit"]].append((r["run"], why.get(r["run"], "not done")))
            left_projects.add(r["bioproject"])
    cmp_rows = []
    for (test, arm), (_, levels, calls, _) in built.items():
        cmp_rows += compare(test, arm, raw[(test, arm)], (levels, calls), have, left, left_projects)
    raw_by = {k: {(x["run"], x["lineage"]): x for x in v[0]} for k, v in raw.items()}
    logan_by = {k: {(x["run"], x["lineage"]): x for x in built[k][1]} if k in built else {} for k in raw}
    sw = switch(positives, by_run, raw_by, logan_by, have, why)
    logan_rows = []
    for r in runs:
        info_p, pm_p = wd / "logan" / f"{r['run']}.json", wd / "pieces" / r["run"] / "pieces.json"
        info = json.loads(info_p.read_text()) if info_p.exists() else {}
        pm = json.loads(pm_p.read_text()) if r["run"] in have else {}
        note = {"HiFi": "HiFi: not repeated on Logan (entry 14)",
                "not in Logan": f"not in Logan v1.2 (HTTP {info.get('http', '')})"}.get(why.get(r["run"]), "")
        logan_rows.append({"run": r["run"], "test": r["test"], "unit": r["unit"],
                           "in_logan": {True: "yes", False: "no"}.get(info.get("in_logan"), ""), "note": note,
                           "bytes": info.get("bytes", ""), "last_modified": info.get("last_modified", ""),
                           **{k: pm.get(k, "") for k in ("unitigs_md5", "read_length", "unitigs", "unitig_bp",
                                                         "pieces", "piece_bp", "weight")}})
    a1 = json.loads((next(iter(adirs.values())) / "assign.json").read_text()) if adirs else {"made_by": ""}
    meta = {"made": now(), "assign": a1["made_by"], "minimap2": mm2, "vsearch": vsearch, "zstd": zstd_v,
            "tools": ", ".join(f"{k} {v}" for k, v in sorted(db["versions"].items())), "t": t, "m": m,
            "logan_url": args.logan_url}
    t4_dir.mkdir(parents=True, exist_ok=True)
    write_outputs(t4_dir, logan_rows, {k: v[1:] for k, v in built.items()}, cmp_rows, sw, k1_reasons, k2, meta)
    rows, share, _ = sw
    print(f"\nLogan holds {len(have)} of the {len(short)} short-read runs. The switch: " +
          ("cannot be worked out" if share is None else
           f"{share:.0%} of the positive runs missed, so the survey runs on " +
           ("raw reads" if share > SWITCH else "Logan unitigs")))
    print(f"-> {t4_dir}/logan.tsv, levels.tsv, calls.tsv, compare.tsv, T4.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
