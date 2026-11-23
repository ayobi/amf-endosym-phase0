#!/usr/bin/env python3
"""Test T1, the paired positive and negative (PROTOCOL.md sections 5 and 6; section 10, entry 9).

For every run in t1/runs.tsv:
  1. reads and read assignment: the first 20,000,000 spots (the HiFi run: 1,000,000 reads), fetched
     with fastq-dump and assigned with assign_reads.py against the full index. Runs marked
     reuse=stage0a use stage 0a's reads and counts (same reads, index and rule).
  2. the route (route = ssu+singlem or ssu): SSUplex, and SingleM for DNA runs, on every read, and the
     reads they find aligned again as in step 1 (run_stage0a_m.py's functions). Stage 0a runs reuse
     step 2's files. Runs with route = none (the reported HiFi and JGI runs) are mapped only.
  3. the routed reads the tools found: SSU pieces placed with entry 8's VSEARCH command and labelled,
     SingleM reads labelled (DNA runs only); reads that the references explain (alignments of 95%
     identity or more that together cover 80% of the read) leave the route, as split reads of a
     bacterial lineage or as host reads.
  4. sylph against the index's bacteria.syldb (a cross-check: it does not change calls).
Then levels and calls per line (call_lineages.py), and, only when every run is done, writes to t1/ in
the repository (commit these):
  levels.tsv   per run and lineage: n, N, the level, whether the run is informative, cross-talk, sylph
  calls.tsv    per line (unit) and lineage: the call, its basis and the evidence
  T1.md        the calls, the confusion table against controls.tsv, and K1 (section 6)
Work files go to WORKDIR/t1/. Reruns skip what is done.

Needs entry 8's m (stage0a/marker_share_entry8.tsv), the phase0-reads environment (sra-tools,
minimap2, SSUplex, VSEARCH, sylph), SingleM in phase0-markers, db/manifest.json and db/silva/.

Usage:
  conda activate phase0-reads
  python scripts/run_t1.py --workdir ~/phase0-data --threads 30
"""

import argparse
import csv
import datetime as dt
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
import assign_reads as AR          # noqa: E402  step 1's reader, minimap2 command and rule
import call_lineages as C          # noqa: E402  labels, levels and calls
import get_marker_databases as DB  # noqa: E402  tool versions
import get_silva as GS             # noqa: E402  VSEARCH version
import rescore_stage0a_m as R      # noqa: E402  entry 8's VSEARCH command, SSU pieces
import run_stage0a_m as M2         # noqa: E402  step 2's route functions

__version__ = "0.1"

SPOTS = {"short": 20_000_000, "hifi": 1_000_000}            # entries 5 and 9
FASTQ_DUMP_ARGS = ["--split-files", "--skip-technical", "--gzip"]
SPLIT_ID, SPLIT_COV = 95, 80                                # entry 9, as the read filter of 4.1
RUN_COLS = ["run", "unit", "control_id", "scored", "kind", "platform", "prep", "layout", "bioproject", "year",
            "lab", "isolate", "route", "reuse", "note"]
PREPS = ("dna", "poly-A", "total RNA", "unstated")
B_PLUS, B_MINUS, PLANT = "BEG34 B+", "BEG34 B-", "Lotus plant-only roots"
K1_PRESENT = ("GMAR_BEG34_HIC_Bplus", "GMAR_BEG34_GENOME")   # section 6: these runs must be at least present
CAGG_GENERA = {"Glomeribacter", "Mycoavidus"}                 # entry 9: K1 counts route lineages in either
LEVEL_COLS = ["run", "unit", "kind", "prep", "bioproject", "lineage", "source", "n", "N", "fraction", "level",
              "informative", "cross_talk", "reads", "fungal_reads", "split_reads", "sylph_ani"]
CALL_COLS = ["unit", "scored", "lineage", "source", "category", "call", "basis", "signal", "informative_runs",
             "present_runs", "low_runs", "e1", "e2", "e5", "members", "placements"]


def now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")


def read_runs(path):
    rows, problems, seen = M2.read_tsv(path), [], set()
    for r in rows:
        if any(c not in r for c in RUN_COLS):
            raise ValueError(f"{path} needs the columns {', '.join(RUN_COLS)}")
        run = r["run"]
        if run in seen:
            problems.append(f"{run} is listed twice")
        seen.add(run)
        want_route = {("DNA", "short"): ("ssu+singlem", "none"), ("RNA", "short"): ("ssu", "none"),
                      ("DNA", "hifi"): ("none",)}.get((r["kind"], r["platform"]))
        if want_route is None:
            problems.append(f"{run}: kind must be DNA or RNA and platform short or hifi (HiFi runs are DNA)")
        elif r["route"] not in want_route:
            problems.append(f"{run}: route must be {' or '.join(want_route)} for a {r['kind']} {r['platform']} run")
        if r["prep"] not in PREPS or (r["kind"] == "DNA") != (r["prep"] == "dna"):
            problems.append(f"{run}: prep must be dna for DNA runs, and poly-A, total RNA or unstated for RNA runs")
        if r["scored"] not in ("yes", "no") or r["layout"] not in ("PAIRED", "SINGLE"):
            problems.append(f"{run}: scored must be yes or no, layout PAIRED or SINGLE")
        if r["reuse"] not in ("", "stage0a") or (r["reuse"] and r["platform"] != "short"):
            problems.append(f"{run}: reuse must be empty or stage0a (short-read runs only)")
        if not r["unit"] or not r["bioproject"] or not r["lab"] or not r["isolate"]:
            problems.append(f"{run}: unit, bioproject, lab and isolate must be given")
    for unit in {r["unit"] for r in rows}:
        mine = [r for r in rows if r["unit"] == unit]
        if len({(r["scored"], r["isolate"]) for r in mine}) != 1:
            problems.append(f"{unit}: its runs disagree on scored or isolate")
    if problems:
        raise ValueError("\n".join(problems))
    return rows


# ------------------------------------------------------------ reads, mapping ---

def fetch(run, reads_dir, spots, args):
    """The first `spots` spots of one run (as run_stage0a.py fetches them). Returns (info, fetched now?)."""
    final = reads_dir / run
    if (final / "fetch.json").exists():
        return json.loads((final / "fetch.json").read_text()), False
    tmp, log = reads_dir / f"{run}.part", reads_dir / f"{run}.fastq-dump.log"
    started = time.time()
    cmd = [args.fastq_dump, *FASTQ_DUMP_ARGS, "-X", str(spots), "-O", str(tmp), run]
    for attempt in range(1, 4):
        shutil.rmtree(tmp, ignore_errors=True)
        tmp.mkdir(parents=True)
        p = subprocess.run(cmd, capture_output=True, text=True)
        with open(log, "a") as fh:
            fh.write(f"attempt {attempt}: {' '.join(cmd)}\nexit {p.returncode}\n{p.stdout}{p.stderr}\n")
        written = re.findall(r"Written (\d+) spots", p.stdout + p.stderr)
        files = sorted(f.name for f in tmp.glob("*.fastq.gz"))
        if p.returncode == 0 and written and files:
            break
    else:
        shutil.rmtree(tmp, ignore_errors=True)
        raise RuntimeError(f"{run}: fastq-dump failed 3 times; see {log}")
    n = int(written[-1])
    if not 0 < n <= spots:
        shutil.rmtree(tmp, ignore_errors=True)
        raise RuntimeError(f"{run}: fastq-dump wrote {n} spots (expected 1 to {spots}); see {log}")
    info = {"run": run, "spots": n, "files": files, "bytes": {f: (tmp / f).stat().st_size for f in files},
            "command": " ".join(["fastq-dump", *FASTQ_DUMP_ARGS, "-X", str(spots), "-O", f"t1/reads/{run}", run]),
            "fetched": now()}
    (tmp / "fetch.json").write_text(json.dumps(info, indent=1) + "\n")
    shutil.rmtree(final, ignore_errors=True)
    os.replace(tmp, final)
    for c in (Path.home() / "ncbi" / "public" / "sra").glob(f"{run}.sra.cache"):
        if c.stat().st_mtime >= started:
            c.unlink()
    return info, True


def assign(run, files, out, platform, args):
    """assign_reads.py on one run unless it is done."""
    if (out / "counts.tsv").exists():
        return
    cmd = [sys.executable, str(SCRIPTS / "assign_reads.py"), "--index", str(args.index), "--platform", platform,
           "--reads", *map(str, files), "--out", str(out), "--refdir", args.refdir, "--threads", str(args.threads),
           "--minimap2", args.minimap2]
    out.mkdir(parents=True, exist_ok=True)
    with open(out.parent / f"{run}.assign.log", "w") as log:
        p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if p.returncode != 0 or not (out / "counts.tsv").exists():
        raise RuntimeError(f"{run}: assign_reads.py failed; see {out.parent / (run + '.assign.log')}")


def assigned(adir, platform, ctx):
    """The mate files and counts of an assigned run, once it is clear they match this index and rule."""
    a = json.loads((adir / "assign.json").read_text())
    if (a["platform"], a["minimap2"], a["index"]["refs_md5"], a["overlaps_md5"]) != \
            (platform, ctx["minimap2"], ctx["refs_md5"], ctx["overlaps_md5"]):
        raise RuntimeError(f"{adir} was assigned with another platform, minimap2, index or overlaps.tsv")
    files = []
    for inp in a["inputs"]:
        p = Path(inp["path"])
        if not p.is_file() or p.stat().st_size != inp["bytes"]:
            raise RuntimeError(f"{p} is missing or is not the file that was assigned")
        files.append((p, inp["reads"]))
    return files, M2.read_counts(adir / "counts.tsv")


# ------------------------------------------------------------------ the route ---

def route_new(r, files, rdir, args, db, ctx):
    """Step 2's route on a run stage 0a did not have: the same files in rdir. SingleM runs on DNA only."""
    if (rdir / "route.json").exists():
        return
    rdir.mkdir(parents=True, exist_ok=True)
    log, dna = rdir / "tools.log", r["route"] == "ssu+singlem"
    table, per_origin = [], {}
    for fq, expected in files:
        stem = fq.name.split(".")[0]
        fa, n_fa = rdir / f"{stem}.reads.fasta", None
        done = (rdir / f"{stem}.ssuplex.json").exists() and (not dna or (rdir / f"{stem}.singlem_names.txt").exists())
        if not done:
            n_fa = M2.mate_fasta(fq, fa, expected)
        ssu, info = M2.ssuplex_names(fa, n_fa, stem, rdir, args, db, log)
        sm = M2.singlem_names(fa, stem, rdir, args, db, log) if dna else []
        fa.unlink(missing_ok=True)
        rows = M2.collect_found(fq, stem, rdir, ssu, sm, expected)
        outcomes = M2.align_found(rdir / f"{stem}.found.fastq.gz", args, ctx, log) if rows else {}
        for pos, name, a, b in rows:
            table.append(f"{stem}\t{pos}\t{name}\t{a}\t{b}\t" + "\t".join(outcomes[pos]) + "\n")
        per_origin[stem] = info["per_origin"]
    part = rdir / "found.tsv.gz.part"
    with gzip.open(part, "wt", encoding="utf-8") as fh:
        fh.write("mate\tposition\tread\tssuplex\tsinglem\toutcome\tdetail\n" + "".join(table))
    os.replace(part, rdir / "found.tsv.gz")
    M2.write_atomic(rdir / "route.json", json.dumps({"made_by": f"run_t1.py v{__version__}", "finished": now(),
                                                     "singlem": dna, "ssuplex_reads_per_origin": per_origin},
                                                    indent=1) + "\n")
    shutil.rmtree(rdir / "tmp", ignore_errors=True)


def routed_found(r, rdir):
    """{(mate, read): (position, found by SSUplex, found by SingleM)} for the routed reads the tools found.
    SingleM counts in DNA runs only (entry 7)."""
    use_sm, out = r["route"] == "ssu+singlem", {}
    with gzip.open(rdir / "found.tsv.gz", "rt", encoding="utf-8") as fh:
        for x in csv.DictReader(fh, delimiter="\t"):
            a, b = x["ssuplex"] == "1", use_sm and x["singlem"] == "1"
            if x["outcome"] == "routed" and (a or b):
                out[(x["mate"], x["read"])] = (int(x["position"]), a, b)
    return out


def singlem_taxonomies(rdir, wanted):
    """(mate, read) -> the taxonomies of the SingleM OTUs the read is in, for the wanted reads."""
    tx = defaultdict(list)
    for f in sorted(rdir.glob("*.singlem.otu.tsv.gz")):
        stem = f.name.split(".")[0]
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            i, t = header.index("read_names"), header.index("taxonomy")
            for line in fh:
                row = line.rstrip("\n").split("\t")
                if len(row) > max(i, t):
                    for name in row[i].split():
                        if (stem, name) in wanted:
                            tx[(stem, name)].append(row[t])
    missing = wanted - tx.keys()
    if missing:
        raise RuntimeError(f"{len(missing)} SingleM reads are not in its OTU tables in {rdir}")
    return tx


def place(queries, bdir, args, log):
    """VSEARCH with entry 8's command on every query not placed yet, in batches of R.BATCH; returns
    id -> [(target, identity)], the top hits."""
    bdir.mkdir(parents=True, exist_ok=True)
    for fa in bdir.glob("batch_*.fa"):                       # a batch cut short: do it again
        if not fa.with_suffix(".tsv").exists():
            fa.unlink()
    done, hits = set(), defaultdict(list)
    for fa in sorted(bdir.glob("batch_*.fa")):
        done.update(line[1:].strip() for line in open(fa) if line.startswith(">"))
        for line in open(fa.with_suffix(".tsv")):
            q, t, i = line.rstrip("\n").split("\t")
            hits[q].append((t, float(i)))
    todo = [q for q in queries if q not in done]
    k = 1 + max((int(fa.stem.split("_")[1]) for fa in bdir.glob("batch_*.fa")), default=0)
    n_batches = (len(todo) + R.BATCH - 1) // R.BATCH
    print(f"placing {len(todo)} SSU sequences ({len(queries) - len(todo)} placed before)", flush=True)
    for j, start in enumerate(range(0, len(todo), R.BATCH), 1):
        chunk = todo[start:start + R.BATCH]
        fa, out = bdir / f"batch_{k:04d}.fa", bdir / f"batch_{k:04d}.tsv"
        fa.write_text("".join(f">{q}\n{queries[q]}\n" for q in chunk))
        part = out.with_name(out.name + ".part")
        cmd = [args.vsearch, "--usearch_global", str(fa), "--db", args.silva_fasta, *R.VSEARCH_ARGS,
               "--userout", str(part), "--userfields", "query+target+id", "--threads", str(args.threads), "--quiet"]
        t0 = time.time()
        with open(log, "a") as fh:
            fh.write("$ " + " ".join(cmd) + "\n")
            fh.flush()
            p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
        if p.returncode != 0 or not part.exists():
            fa.unlink(missing_ok=True)
            part.unlink(missing_ok=True)
            raise RuntimeError(f"VSEARCH failed on batch {k} (exit {p.returncode}); see {log}")
        os.replace(part, out)
        for line in open(out):
            q, t, i = line.rstrip("\n").split("\t")
            hits[q].append((t, float(i)))
        secs = time.time() - t0
        print(f"  batch {j}/{n_batches}: {len(chunk)} sequences in {secs / 60:.1f} min; about "
              f"{secs * (n_batches - j) / 3600:.1f} h left", flush=True)
        k += 1
    return hits


def explained_reads(found_fq, positions, args, ctx, log):
    """Align the given reads of a found.fastq.gz again (step 1's command, named by position) and return
    {position: bacterial lineage or ''} for those the references explain: their alignments of 95%
    identity or more together cover at least 80% of the read (entry 9). '' means hosts or organelles alone."""
    cmd = [args.minimap2, *AR.MINIMAP2_ARGS["short"], "-t", str(args.threads), str(ctx["mmi"]), "-"]
    with open(log, "a") as err:
        err.write("$ " + " ".join(cmd) + f" < {found_fq.name} (routed reads, named by position)\n")
        err.flush()
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=err, bufsize=1 << 20)
    state = {}

    def feed():
        try:
            for h, s, q in AR.records(found_fq):
                pos = int(M2.first_word(h))
                if pos in positions:
                    proc.stdin.write(b"@%d\n%s\n+\n%s\n" % (pos, s, q))
        except BaseException as e:   # reported by the main thread
            state["error"] = e
        finally:
            try:
                proc.stdin.close()
            except OSError:
                pass

    feeder = threading.Thread(target=feed, daemon=True)
    feeder.start()
    out, seen = {}, set()
    try:
        for pos, alns in AR.parse_paf(proc.stdout):
            seen.add(pos)
            if alns and AR.evaluate(alns, ctx["lineage_of"], ctx["hosts"], ctx["mask"])[0] != "routed":
                raise RuntimeError(f"read {pos} of {found_fq.name} is not routed when aligned again")
            ok = sorted((a[3], a[4]) for a in alns if a[1] * 100 >= SPLIT_ID * a[2])
            covered, end = 0, -1
            for s, e in ok:
                if e > end:
                    covered += e - max(s, end)
                    end = e
            if ok and covered * 100 >= SPLIT_COV * alns[0][5]:
                bact = [(a[4] - a[3], ctx["lineage_of"][a[6]]) for a in alns
                        if a[1] * 100 >= SPLIT_ID * a[2] and ctx["lineage_of"][a[6]] in ctx["bacterial"]]
                out[pos] = max(bact)[1] if bact else ""
    except BaseException:
        proc.kill()
        feeder.join(timeout=60)
        raise
    feeder.join()
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"minimap2 failed (exit {rc}); see {log}")
    if state.get("error"):
        raise state["error"]
    if seen != set(positions):
        raise RuntimeError(f"minimap2 did not report every read it was given from {found_fq.name}")
    return out


def split_check(found, rdir, cache, args, ctx):
    """{(mate, position): lineage or ''} for the routed found reads the references explain; cached."""
    if cache.exists():
        res = {}
        for k, v in json.loads(cache.read_text()).items():
            mate, pos = k.rsplit(":", 1)
            res[(mate, int(pos))] = v
        return res
    by_mate = defaultdict(set)
    for (mate, _), (pos, _, _) in found.items():
        by_mate[mate].add(pos)
    res = {}
    for mate, positions in sorted(by_mate.items()):
        for pos, lin in explained_reads(rdir / f"{mate}.found.fastq.gz", positions, args, ctx,
                                        cache.with_suffix(".log")).items():
            res[(mate, pos)] = lin
    M2.write_atomic(cache, json.dumps({f"{m}:{p}": v for (m, p), v in sorted(res.items())}) + "\n")
    return res


# ---------------------------------------------------------------------- sylph ---

def sylph(files, out, args, ctx):
    """sylph profile on the run's reads, against the index's bacteria.syldb; returns {lineage: best ANI}."""
    if not out.exists():
        paths = [str(p) for p, _ in files]
        part = out.with_name(out.name + ".part")
        cmd = [args.sylph, "profile", str(ctx["syldb"]), "-t", str(args.threads), "-o", str(part)]
        cmd += ["-1", paths[0], "-2", paths[1]] if len(paths) == 2 else paths
        with open(out.with_suffix(".log"), "a") as fh:
            fh.write("$ " + " ".join(cmd) + "\n")
            fh.flush()
            p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
        if p.returncode != 0 or not part.exists():
            part.unlink(missing_ok=True)
            raise RuntimeError(f"sylph failed (exit {p.returncode}); see {out.with_suffix('.log')}")
        os.replace(part, out)
    best = {}
    for row in M2.read_tsv(out):
        acc = re.sub(r"\.(fna|fa|fasta)(\.gz)?$", "", Path(row["Genome_file"]).name)
        lin = ctx["acc_lineage"].get(acc)
        if lin is None:
            raise RuntimeError(f"sylph reports {row['Genome_file']}, which is not in references/lineages.tsv")
        try:
            ani = float(row["Adjusted_ANI"])
        except ValueError:
            continue
        best[lin] = max(best.get(lin, 0.0), ani)
    return best


# --------------------------------------------------------------------- calls ---

def build(runs, data, t, m, ctx):
    """Levels for every lineage and run, then calls for every unit. Returns (levels, calls, joiner)."""
    genera = sorted({k[1] for d in data.values() if d["route"] for k in d["route"] if k[0] == "genus"})
    joiner = C.Joiner(genera)
    route_counts = {}
    for run, d in data.items():
        if d["route"] is not None:
            c = Counter()
            for k, n in d["route"].items():
                c[C.lineage_name(k, joiner)] += n
            route_counts[run] = c
    route_lins = sorted(joiner.members) + [x for x in (C.NOVEL_SSU, C.NOVEL_SINGLEM)
                                            if any(x in c for c in route_counts.values())]
    lineages = [(x, "mapped") for x in ctx["bacterial_order"]] + [(x, "route") for x in route_lins]
    levels = []
    for lin, source in lineages:
        recs = []
        for r in runs:
            d = data[r["run"]]
            if source == "route" and (d["route"] is None or (lin == C.NOVEL_SINGLEM and r["route"] != "ssu+singlem")):
                continue
            n = d["mapped"].get(lin, 0) if source == "mapped" else route_counts[r["run"]].get(lin, 0)
            N = t * d["reads"] * (1 if source == "mapped" else m)
            recs.append({"run": r["run"], "unit": r["unit"], "kind": r["kind"], "prep": r["prep"],
                         "bioproject": r["bioproject"], "year": r["year"], "lab": r["lab"], "isolate": r["isolate"],
                         "lineage": lin, "source": source, "n": n, "N": N, "reads": d["reads"],
                         "fungal": d["fungal"], "absent_ok": r["kind"] == "DNA" or r["prep"] == "total RNA",
                         "split_reads": d["split"].get(lin, 0) if source == "mapped" else "",
                         "sylph_ani": d["sylph"].get(lin, "") if source == "mapped" else ""})
        levels += C.run_levels(recs)
    units = []
    for r in runs:
        if r["unit"] not in [u for u, _ in units]:
            units.append((r["unit"], r))
    calls = []
    for lin, source in lineages:
        recs = [x for x in levels if x["lineage"] == lin]
        if source == "mapped":
            category, members = C.category_of_mapped(lin), ""
        elif lin in joiner.members:
            category, members = C.category_of_route(joiner.members[lin]), ", ".join(joiner.members[lin])
        else:
            category, members = "novel", ""
        for unit, first in units:
            mine = [x for x in recs if x["unit"] == unit]
            if not mine:
                continue
            labs = {x["lab"] for x in mine if x["level"] == "present"}
            background = defaultdict(list)
            for y in recs:
                if y["lab"] in labs and y["isolate"] not in (first["isolate"], "none") and y["level"] == "present":
                    background[y["isolate"]].append(y["fraction"])
            res = C.call(mine, category, dict(background))
            places = ""
            if category == "novel":
                pc = Counter()
                for x in mine:
                    pc.update(data[x["run"]]["placements"].get(lin, {}))
                places = "; ".join(f"{p} ({n})" for p, n in pc.most_common(5))
            calls.append({"unit": unit, "scored": first["scored"], "lineage": lin, "source": source,
                          "members": members, "placements": places, **res})
    return levels, calls, joiner


def k1(calls, levels, runs):
    """Section 6, K1: the reasons T1 fails (an empty list means it passes)."""
    by = {(c["unit"], c["lineage"]): c for c in calls}
    reasons = []
    cagg_like = [c for c in calls if c["unit"] == B_MINUS and c["call"] == "resident" and
                 (c["lineage"] == "CaGg" or (c["source"] == "route" and
                                              {p for g in c["members"].split(", ") if g for p in C.genus_parts(g)}
                                              & CAGG_GENERA))]
    reasons += [f"the B- line is called resident for {c['lineage']}" for c in cagg_like]
    reasons += [f"the plant-only roots are called resident for {c['lineage']}" for c in calls
                if c["unit"] == PLANT and c["call"] == "resident"]
    b = by.get((B_PLUS, "CaGg"))
    if b is None or b["call"] != "resident":
        reasons.append(f"the B+ line is not called resident for CaGg ({b['call'] if b else 'no call made'})")
    want = {r["run"] for r in runs if r["control_id"] in K1_PRESENT}
    for x in levels:
        if x["lineage"] == "CaGg" and x["run"] in want and x["level"] != "present":
            reasons.append(f"CaGg is {x['level']}, not present, in {x['run']}")
    return reasons


# -------------------------------------------------------------------- outputs ---

def fmt(v, spec=".4g"):
    return "" if v in ("", None) else format(v, spec)


def write_outputs(t1_dir, runs, data, levels, calls, reasons, controls, meta):
    tmp = {}

    def level_row(x):
        v = {**x, "fungal_reads": x["fungal"], "informative": yn(x["informative"]), "cross_talk": yn(x["cross_talk"]),
             "N": fmt(x["N"]), "fraction": fmt(x["fraction"]), "sylph_ani": fmt(x["sylph_ani"], ".2f")}
        return "\t".join(str(v[c]) for c in LEVEL_COLS) + "\n"

    tmp["levels.tsv"] = "\t".join(LEVEL_COLS) + "\n" + "".join(level_row(x) for x in levels)
    tmp["calls.tsv"] = "\t".join(CALL_COLS) + "\n" + "".join(
        "\t".join(("yes" if c[k] else "no") if k in ("signal", "e1", "e2", "e5") else str(c[k]) for k in CALL_COLS)
        + "\n" for c in calls)
    units = []
    for r in runs:
        if r["unit"] not in units:
            units.append(r["unit"])
    by = {(c["unit"], c["lineage"]): c for c in calls}
    md = ["# Test T1: the paired positive and negative", "",
          f"Made by run_t1.py v{__version__} (call_lineages.py v{C.__version__}) on {meta['made']}, with "
          f"{meta['assign']}, minimap2 {meta['minimap2']} and the full index (refs.fna md5 {meta['refs_md5']}); "
          f"SSUplex, SingleM and HMMER {meta['tools']}; VSEARCH {meta['vsearch']} against SILVA `{meta['silva']}`; "
          f"{meta['sylph']}.",
          f"Rule: PROTOCOL.md 4.1, section 5 (T1) and section 6 (K1), with section 10, entries 4 to 9. "
          f"t = {meta['t']:.4g} (entry 5), m = {meta['m']:.4g} (entry 8).", "",
          f"## K1: {'passed' if not reasons else 'FAILED'}", ""]
    md += [f"- {x}" for x in reasons] or ["None of the K1 conditions holds: the B- line is not called resident "
                                          "for CaGg, the plant-only roots are called resident for nothing, and the "
                                          "B+ line is resident for CaGg, present in its Hi-C and genome runs."]
    md += ["", "## Against controls.tsv", "", "| Line | Scored | Expected | CaGg | Basis | Lineages called resident |",
           "|---|---|---|---|---|---|"]
    for u in units:
        cids = sorted({r["control_id"] for r in runs if r["unit"] == u})
        exp = "; ".join(sorted({controls.get(c, "?") for c in cids}))
        c = by.get((u, "CaGg"), {})
        res = ", ".join(x["lineage"] for x in calls if x["unit"] == u and x["call"] == "resident") or "none"
        md.append(f"| {u} | {runs_of(runs, u)[0]['scored']} | {exp} | {c.get('call', '')} | {c.get('basis', '')} "
                  f"| {res} |")
    md += ["", "## Calls other than absent and no call", "",
           "| Line | Lineage | Category | Call | Basis | Present / low runs | E1 | E2 | E5 |",
           "|---|---|---|---|---|---|---|---|---|"]
    for c in calls:
        if c["call"] not in ("absent", "no call"):
            md.append(f"| {c['unit']} | {c['lineage']} | {c['category']} | {c['call']} | {c['basis']} | "
                      f"{c['present_runs']} / {c['low_runs']} | {yn(c['e1'])} | {yn(c['e2'])} | {yn(c['e5'])} |")
    md += ["", "## CaGg per run", "", "| Run | Line | Kind | Reads | n | N | Level | Informative | sylph ANI |",
           "|---|---|---|---|---|---|---|---|---|"]
    for x in levels:
        if x["lineage"] == "CaGg":
            md.append(f"| {x['run']} | {x['unit']} | {x['kind']} {x['prep'] if x['kind'] == 'RNA' else ''} | "
                      f"{x['reads']} | {x['n']} | {x['N']:.4g} | {x['level']}"
                      f"{' (cross-talk)' if x['cross_talk'] else ''} | {yn(x['informative'])} | "
                      f"{fmt(x['sylph_ani'], '.2f') or '-'} |")
    md += ["", "## The route per run", "",
           "| Run | Routed reads found | SSU kept | SSU dropped: outside Bacteria / organelle | SingleM kept | "
           "Split reads: bacterial / host |", "|---|---|---|---|---|---|"]
    for r in runs:
        d = data[r["run"]]
        if d["route"] is None:
            continue
        s = d["stats"]
        md.append(f"| {r['run']} | {s['found']} | {s['ssu_kept']} | {s['ssu_outside']} / {s['ssu_organelle']} | "
                  f"{s['sm_kept'] if r['route'] == 'ssu+singlem' else 'n/a'} | {s['split_bact']} / {s['split_host']} |")
    md += ["", "Reads are counted per read (mates separately). `levels.tsv` gives every lineage in every run; "
               "`calls.tsv` gives every call with its evidence; novel lineages list their most common placements. "
               "sylph is a cross-check and does not change calls (4.2)."]
    tmp["T1.md"] = "\n".join(md) + "\n"
    for name, text in tmp.items():
        (t1_dir / (name + ".part")).write_text(text, encoding="utf-8")
    for name in tmp:
        os.replace(t1_dir / (name + ".part"), t1_dir / name)


def runs_of(runs, unit):
    return [r for r in runs if r["unit"] == unit]


def yn(v):
    return "yes" if v else "no"


# ----------------------------------------------------------------------- main ---

def stage_values(stage_dir):
    """t (entry 5) and m (entry 8), from the committed stage 0a outputs."""
    fr = [x for x in M2.read_tsv(stage_dir / "fractions.tsv") if x["sets_t"] == "yes"]
    e8 = stage_dir / "marker_share_entry8.tsv"
    if not e8.exists():
        raise RuntimeError(f"{e8} is missing: entry 8's rescore must finish (and be committed) first")
    ms = [float(x["share"]) for x in M2.read_tsv(e8) if x["sets_m"] == "yes" and x["share"]]
    if not fr or not ms:
        raise RuntimeError("t or m cannot be read from stage0a/")
    return min(float(x["fraction"]) for x in fr) / 10, min(ms)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--index", help="the full index set [WORKDIR/index/full]")
    ap.add_argument("--refdir", default="references")
    ap.add_argument("--t1-dir", default="t1", help="holds runs.tsv; the outputs to commit go here too")
    ap.add_argument("--stage0a-dir", default="stage0a")
    ap.add_argument("--controls", default="controls.tsv")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--jobs", type=int, default=3, help="runs fetched at the same time [3]")
    ap.add_argument("--fastq-dump", default="fastq-dump")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--ssuplex", default="ssuplex")
    ap.add_argument("--conda", default=os.environ.get("CONDA_EXE") or shutil.which("conda"))
    ap.add_argument("--singlem-env", default="phase0-markers")
    ap.add_argument("--vsearch", default="vsearch")
    ap.add_argument("--sylph", default="sylph")
    args = ap.parse_args(argv)
    work, t1_dir, stage_dir = Path(args.workdir).expanduser(), Path(args.t1_dir), Path(args.stage0a_dir)
    args.index = Path(args.index).expanduser() if args.index else work / "index" / "full"

    try:
        runs = read_runs(t1_dir / "runs.tsv")
        controls = {x["control_id"]: x["expected_call"] for x in M2.read_tsv(Path(args.controls))}
        missing = sorted({r["control_id"] for r in runs} - controls.keys())
        if missing:
            raise RuntimeError(f"control_id not in {args.controls}: {', '.join(missing)}")
        t, m = stage_values(stage_dir)
        info, lineage_of, _, hosts = AR.load_index(args.index)
        if info.get("set") != "full":
            raise RuntimeError(f"{args.index} is the {info.get('set')} index; T1 uses the full references")
        overlaps = Path(args.refdir) / "overlaps.tsv"
        mm2 = subprocess.run([args.minimap2, "--version"], capture_output=True, text=True, check=True).stdout.strip()
        if mm2 != info["minimap2"]:
            raise RuntimeError(f"minimap2 here is {mm2}, but the index was built with {info['minimap2']}")
        sylph_v = DB.tool_output([args.sylph, "--version"]).strip().splitlines()[0]
        if sylph_v != info["sylph"]:
            raise RuntimeError(f"sylph here is '{sylph_v}', but the index was built with '{info['sylph']}'")
        silva = json.loads((work / "db" / "silva" / "manifest.json").read_text())
        vsearch = GS.vsearch_version(args.vsearch)
        if vsearch != silva["vsearch"]:
            raise RuntimeError(f"VSEARCH {vsearch} is not the {silva['vsearch']} get_silva.py checked")
        args.silva_fasta = silva["fasta"]
        db = json.loads((work / "db" / "manifest.json").read_text())
        if any(r["route"] != "none" and not r["reuse"] and
               not (work / "t1" / "route" / r["run"] / "route.json").exists() for r in runs):
            versions = DB.check_versions(args)
            if versions != db["versions"]:
                raise RuntimeError(f"the tools here ({versions}) are not those get_marker_databases.py checked "
                                   f"({db['versions']})")
        lin_rows = M2.read_tsv(Path(args.refdir) / "lineages.tsv")
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.CalledProcessError) as e:
        print(f"Cannot start: {e}", file=sys.stderr)
        return 1
    bact_kinds = {x["lineage"]: x["kind"] for x in lin_rows if x["kind"] in ("endosymbiont", "contaminant")}
    order = [x for x in ("CaGg", "CaMg/MRE", "Mycoavidus") if x in bact_kinds] + \
        sorted(x for x in bact_kinds if x not in ("CaGg", "CaMg/MRE", "Mycoavidus"))
    ctx = {"mmi": args.index / AR.INDEX_FILE["short"], "lineage_of": lineage_of, "hosts": hosts,
           "mask": AR.load_mask(overlaps), "minimap2": mm2, "refs_md5": info["refs_md5"],
           "overlaps_md5": AR.md5_of(overlaps), "bacterial": set(bact_kinds), "bacterial_order": order,
           "syldb": args.index / "bacteria.syldb", "acc_lineage": {x["accession"]: x["lineage"] for x in lin_rows}}

    tw = work / "t1"
    for sub in ("reads", "assign", "route", "split", "sylph"):
        (tw / sub).mkdir(parents=True, exist_ok=True)
    failed, files, counts = [], {}, {}

    def one_fetch(r):
        try:
            info_, now_ = fetch(r["run"], tw / "reads", SPOTS[r["platform"]], args)
            print(f"{'fetched' if now_ else 'already fetched'}: {r['run']} ({info_['spots']:,} spots)", flush=True)
            return r["run"], None
        except (OSError, RuntimeError) as e:
            return r["run"], str(e)

    new = [r for r in runs if not r["reuse"]]
    print(f"fetching {len(new)} runs, {args.jobs} at a time ({len(runs) - len(new)} reuse stage 0a)", flush=True)
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        for run, err in pool.map(one_fetch, new):
            if err:
                failed.append(err)
                print(f"FAILED: {err}", flush=True)
    for r in runs:
        run = r["run"]
        try:
            if r["reuse"]:
                files[run], counts[run] = assigned(work / "stage0a" / "assign" / run, "short", ctx)
                continue
            if not (tw / "reads" / run / "fetch.json").exists():
                continue
            fi = json.loads((tw / "reads" / run / "fetch.json").read_text())
            if not (tw / "assign" / run / "counts.tsv").exists():
                print(f"assigning {run}", flush=True)
            assign(run, [tw / "reads" / run / f for f in fi["files"]], tw / "assign" / run, r["platform"], args)
            files[run], counts[run] = assigned(tw / "assign" / run, r["platform"], ctx)
        except (OSError, KeyError, RuntimeError, json.JSONDecodeError) as e:
            failed.append(f"{run}: {e}")
            print(f"FAILED: {run}: {e}", flush=True)

    rdirs = {}
    for r in runs:
        run = r["run"]
        if r["route"] == "none" or run not in files:
            continue
        rdirs[run] = work / "stage0a" / "route" / run if r["reuse"] else tw / "route" / run
        try:
            if r["reuse"]:
                if not (rdirs[run] / "route.json").exists():
                    raise RuntimeError(f"stage 0a step 2 did not finish {run}")
            else:
                if not (rdirs[run] / "route.json").exists():
                    print(f"route: {run}", flush=True)
                route_new(r, files[run], rdirs[run], args, db, ctx)
        except (OSError, KeyError, ValueError, RuntimeError) as e:
            failed.append(f"{run}: {e}")
            print(f"FAILED: {run}: {e}", flush=True)
            del rdirs[run]
    if failed:
        print(f"\n{len(failed)} run(s) not done, so nothing was written to {t1_dir}/. Run the same command again "
              "to finish them:\n" + "\n".join(failed), file=sys.stderr)
        return 1

    try:
        found, pieces, sm_tax = {}, {}, {}
        for r in runs:
            run = r["run"]
            if run not in rdirs:
                continue
            found[run] = routed_found(r, rdirs[run])
            pieces[run] = R.ssu_sequences(rdirs[run], {k for k, v in found[run].items() if v[1]})
            sm_tax[run] = singlem_taxonomies(rdirs[run], {k for k, v in found[run].items() if v[2]})
        queries = {}
        for run in pieces:
            for s in pieces[run].values():
                queries.setdefault(R.seq_id(s), s)
        hits = place(queries, tw / "place" / "batches", args, tw / "place" / "vsearch.log")
        tax = {k: [x.strip() for x in v.split(";")]
               for k, v in R.load_taxonomy(silva["taxonomy"], {t for q in queries for t, _ in hits.get(q, [])}).items()}
        data = {}
        for r in runs:
            run, c = r["run"], counts[r["run"]]
            d = data[run] = {"reads": c[("total", "reads in input")], "fungal": c[("fungal", "fungal reads for E2")],
                             "mapped": {lin: c.get(("lineage", lin), 0) for lin in ctx["bacterial"]},
                             "route": None, "split": Counter(), "placements": defaultdict(Counter), "stats": Counter()}
            if run in rdirs:
                print(f"route reads: {run}", flush=True)
                split = split_check(found[run], rdirs[run], tw / "split" / f"{run}.json", args, ctx)
                d["route"], st = Counter(), d["stats"]
                for key, (pos, a, b) in found[run].items():
                    st["found"] += 1
                    ssu = C.ssu_label(hits.get(R.seq_id(pieces[run][key]), []), tax) if a else None
                    sml = C.singlem_label(sm_tax[run][key]) if b else None
                    if ssu:
                        st["ssu_kept" if ssu[0] == "kept" else "ssu_" + ssu[0].split()[0]] += 1
                    if sml and sml[0] == "kept":
                        st["sm_kept"] += 1
                    label = C.read_label(ssu, sml)
                    if label is None:
                        continue
                    if (key[0], pos) in split:
                        lin = split[(key[0], pos)]
                        st["split_bact" if lin else "split_host"] += 1
                        if lin:
                            d["split"][lin] += 1
                        continue
                    d["route"][label] += 1
                    if label[0] == "novel":
                        lname = C.NOVEL_SSU if label[1] == "SSU" else C.NOVEL_SINGLEM
                        d["placements"][lname][(ssu if ssu and ssu[0] == "kept" else sml)[2]] += 1
            print(f"sylph: {run}", flush=True)
            d["sylph"] = sylph(files[run], tw / "sylph" / f"{run}.tsv", args, ctx)
        levels, calls, _ = build(runs, data, t, m, ctx)
        reasons = k1(calls, levels, runs)
    except (OSError, KeyError, ValueError, RuntimeError) as e:
        print(f"Nothing was written to {t1_dir}/: {e}", file=sys.stderr)
        return 1
    a1 = json.loads((tw / "assign" / new[0]["run"] / "assign.json").read_text()) if new else {"made_by": "?"}
    meta = {"made": now(), "assign": a1["made_by"], "minimap2": mm2, "refs_md5": info["refs_md5"],
            "tools": ", ".join(f"{k} {v}" for k, v in sorted(db["versions"].items())), "vsearch": vsearch,
            "silva": silva["file"], "sylph": sylph_v, "t": t, "m": m}
    write_outputs(t1_dir, runs, data, levels, calls, reasons, controls, meta)
    print(f"\nK1: {'passed' if not reasons else 'FAILED: ' + '; '.join(reasons)}")
    for u in (B_PLUS, B_MINUS, PLANT):
        c = next((x for x in calls if x["unit"] == u and x["lineage"] == "CaGg"), None)
        if c:
            print(f"  {u}: CaGg {c['call']} ({c['basis']})")
    print(f"-> {t1_dir}/levels.tsv, {t1_dir}/calls.tsv, {t1_dir}/T1.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
