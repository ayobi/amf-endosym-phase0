#!/usr/bin/env python3
"""Stage 0a, step 1: measure the reference titer t (PROTOCOL.md 4.1 and 5; section 10, entry 5).

For every run in stage0a/runs.tsv (the T1 and T2 positives):
  1. fetch the first 20,000,000 spots (read pairs, or reads for single-end runs) from SRA with
     fastq-dump, gzipped, into --workdir/stage0a/reads/<run>/
  2. assign the reads with assign_reads.py against the full index, into --workdir/stage0a/assign/<run>/
Then, only when every run is done, write to stage0a/ in the repository (commit these):
  fractions.tsv   per run: reads, reads of its known symbiont, the fraction, fungal and organelle reads
  STAGE0A.md      t = the lowest fraction among the DNA positives (sets_t = yes), divided by 10

The lineage read fraction is the reads counted for the run's known symbiont (CaGg for T1, CaMg/MRE
for T2) divided by all reads, both counted per read (mates separately). RNA positives are measured
and reported, but do not set t.

Reruns skip what is done: a run with reads already fetched is not fetched again, and a run already
assigned is not assigned again. The reads are kept for step 2 (the marker share m).

Needs the phase0-reads environment (sra-tools, minimap2) and the index from make_index.py.

Usage:
  conda activate phase0-reads
  python scripts/run_stage0a.py --workdir ~/phase0-data --threads 16
"""

import argparse
import csv
import datetime as dt
import json
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

__version__ = "0.1"

SPOTS = 20_000_000                                           # section 10, entry 5
FASTQ_DUMP_ARGS = ["--split-files", "--skip-technical", "--gzip"]
SCRIPTS = Path(__file__).resolve().parent
RUN_COLS = ["run", "test", "kind", "line", "symbiont_lineage", "sets_t", "bioproject", "layout", "note"]
FRACTION_COLS = ["run", "test", "kind", "line", "symbiont_lineage", "sets_t", "spots", "reads", "symbiont_reads",
                 "fraction", "fungal_reads", "organelle_reads", "note"]
LINEAGES = ("CaGg", "CaMg/MRE")


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_runs(path):
    rows, problems, seen = read_tsv(path), [], set()
    for r in rows:
        if any(c not in r for c in RUN_COLS):
            raise ValueError(f"{path} needs the columns {', '.join(RUN_COLS)}")
        if r["run"] in seen:
            problems.append(f"{r['run']} is listed twice")
        seen.add(r["run"])
        if r["kind"] not in ("DNA", "RNA") or r["sets_t"] not in ("yes", "no"):
            problems.append(f"{r['run']}: kind must be DNA or RNA, and sets_t yes or no")
        if r["kind"] == "RNA" and r["sets_t"] == "yes":
            problems.append(f"{r['run']}: RNA runs do not set t (section 10, entry 5)")
        if r["symbiont_lineage"] not in LINEAGES:
            problems.append(f"{r['run']}: symbiont_lineage must be one of {', '.join(LINEAGES)}")
    if not any(r["sets_t"] == "yes" for r in rows):
        problems.append("no run sets t")
    if problems:
        raise ValueError("\n".join(problems))
    return rows


def sra_cache_dir():
    return Path.home() / "ncbi" / "public" / "sra"


def fetch(run, reads_dir, args):
    """Fetch the first SPOTS spots of one run. Returns (info, fetched now?)."""
    final = reads_dir / run
    if (final / "fetch.json").exists():
        return json.loads((final / "fetch.json").read_text()), False
    tmp = reads_dir / f"{run}.part"
    log = reads_dir / f"{run}.fastq-dump.log"
    started = time.time()
    cmd = [args.fastq_dump, *FASTQ_DUMP_ARGS, "-X", str(SPOTS), "-O", str(tmp), run]
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
    spots = int(written[-1])
    if not 0 < spots <= SPOTS:
        shutil.rmtree(tmp, ignore_errors=True)
        raise RuntimeError(f"{run}: fastq-dump wrote {spots} spots (expected 1 to {SPOTS}); see {log}")
    info = {"run": run, "spots": spots, "files": files,
            "bytes": {f: (tmp / f).stat().st_size for f in files},
            "command": " ".join(["fastq-dump", *FASTQ_DUMP_ARGS, "-X", str(SPOTS), "-O", f"stage0a/reads/{run}", run]),
            "fetched": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")}
    (tmp / "fetch.json").write_text(json.dumps(info, indent=1) + "\n")
    shutil.rmtree(final, ignore_errors=True)
    os.replace(tmp, final)
    # A remote fastq-dump can leave a partial cache of the run; remove only one made by this fetch.
    for c in sra_cache_dir().glob(f"{run}.sra.cache"):
        if c.stat().st_mtime >= started:
            c.unlink()
    return info, True


def assign(run, files, out, args):
    """Run assign_reads.py on one run unless it is done. Returns True if it ran now."""
    if (out / "counts.tsv").exists():
        return False
    cmd = [sys.executable, str(SCRIPTS / "assign_reads.py"), "--index", str(args.index), "--platform", "short",
           "--reads", *map(str, files), "--out", str(out), "--refdir", args.refdir,
           "--threads", str(args.threads), "--minimap2", args.minimap2]
    out.mkdir(parents=True, exist_ok=True)
    with open(out.parent / f"{run}.assign.log", "w") as log:
        p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if p.returncode != 0 or not (out / "counts.tsv").exists():
        raise RuntimeError(f"{run}: assign_reads.py failed; see {out.parent / (run + '.assign.log')}")
    return True


def read_counts(path):
    return {(r["section"], r["name"]): int(r["reads"]) for r in read_tsv(path)}


def summarize(runs, infos, assign_dir):
    rows = []
    for r in runs:
        c = read_counts(assign_dir / r["run"] / "counts.tsv")
        total = c[("total", "reads in input")]
        n = c.get(("lineage", r["symbiont_lineage"]))
        if n is None:
            raise RuntimeError(f"{r['run']}: {r['symbiont_lineage']} is not a lineage in its counts (wrong index?)")
        if total == 0:
            raise RuntimeError(f"{r['run']}: no reads")
        rows.append({**{k: r[k] for k in ("run", "test", "kind", "line", "symbiont_lineage", "sets_t", "note")},
                     "spots": infos[r["run"]]["spots"], "reads": total, "symbiont_reads": n,
                     "fraction": n / total, "fungal_reads": c[("fungal", "fungal reads for E2")],
                     "organelle_reads": c.get(("lineage", "plant organelles"), 0)})
    setting = [x for x in rows if x["sets_t"] == "yes"]
    zero = [x["run"] for x in setting if x["symbiont_reads"] == 0]
    if zero:
        raise RuntimeError("t would be 0: no symbiont reads in the DNA positive(s) " + ", ".join(zero) +
                           ". Stop here; this needs a deviations entry before stage 0a can finish.")
    lowest = min(setting, key=lambda x: x["fraction"])
    return rows, lowest, lowest["fraction"] / 10


def write_outputs(stage_dir, rows, lowest, t, meta):
    tmp = stage_dir / "fractions.tsv.part"
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FRACTION_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for x in rows:
            w.writerow({**x, "fraction": f"{x['fraction']:.10g}"})
    md = ["# Stage 0a, step 1: reference titer t", "",
          f"Made by run_stage0a.py v{__version__} on {meta['made']}, with {meta['assign']} and minimap2 "
          f"{meta['minimap2']}, against the full index (refs.fna md5 {meta['refs_md5']}; `references/INDEX.md`). "
          f"Reads: the first {SPOTS:,} spots of each run, fetched with {meta['fastq_dump']}.",
          "Rule: PROTOCOL.md 4.1 (reference titer) with section 10, entry 5: t is the lowest lineage read fraction "
          "among the DNA positives, divided by 10. RNA positives are measured and reported; they do not set t.", "",
          f"**t = {t:.4g}**: the lowest DNA-positive fraction is {lowest['fraction']:.4g}, in {lowest['run']} "
          f"({lowest['line']}, {lowest['note']}).", "",
          "| Run | Test | Kind | Line | Symbiont | Sets t | Reads | Symbiont reads | Fraction |",
          "|---|---|---|---|---|---|---|---|---|"]
    for x in rows:
        md.append(f"| {x['run']} | {x['test']} | {x['kind']} | {x['line']} | {x['symbiont_lineage']} | {x['sets_t']} | "
                  f"{x['reads']} | {x['symbiont_reads']} | {x['fraction']:.4g} |")
    md += ["", "Reads are counted per read (mates separately). Details per run, with fungal and organelle reads: "
               "`fractions.tsv`. Step 2 measures the marker share m (section 10, entry 6)."]
    tmp_md = stage_dir / "STAGE0A.md.part"
    tmp_md.write_text("\n".join(md) + "\n", encoding="utf-8")
    os.replace(tmp, stage_dir / "fractions.tsv")
    os.replace(tmp_md, stage_dir / "STAGE0A.md")


def version_of(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True)
    return (p.stdout or p.stderr).strip().splitlines()[0] if (p.stdout or p.stderr).strip() else "unknown"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--index", help="the full index set [WORKDIR/index/full]")
    ap.add_argument("--refdir", default="references")
    ap.add_argument("--stage0a-dir", default="stage0a", help="holds runs.tsv; the outputs to commit go here too")
    ap.add_argument("--threads", type=int, default=8, help="threads for minimap2")
    ap.add_argument("--jobs", type=int, default=3, help="runs fetched at the same time [3]")
    ap.add_argument("--fastq-dump", default="fastq-dump")
    ap.add_argument("--minimap2", default="minimap2")
    args = ap.parse_args(argv)
    work, stage_dir = Path(args.workdir).expanduser(), Path(args.stage0a_dir)
    args.index = Path(args.index).expanduser() if args.index else work / "index" / "full"

    try:
        runs = read_runs(stage_dir / "runs.tsv")
        index_info = json.loads((args.index / "index.json").read_text())
    except (OSError, ValueError, json.JSONDecodeError) as e:
        print(f"Cannot start: {e}", file=sys.stderr)
        return 1
    if index_info.get("set") != "full":
        print(f"{args.index} is the {index_info.get('set')} index; stage 0a uses the full references.", file=sys.stderr)
        return 1

    reads_dir, assign_dir = work / "stage0a" / "reads", work / "stage0a" / "assign"
    reads_dir.mkdir(parents=True, exist_ok=True)
    assign_dir.mkdir(parents=True, exist_ok=True)

    infos, failed = {}, []

    def one(r):
        try:
            info, now = fetch(r["run"], reads_dir, args)
            print(f"{'fetched' if now else 'already fetched'}: {r['run']} ({info['spots']:,} spots)", flush=True)
            return r["run"], info, None
        except (OSError, RuntimeError) as e:
            return r["run"], None, str(e)

    print(f"fetching {len(runs)} runs, {args.jobs} at a time", flush=True)
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        for run, info, err in pool.map(one, runs):
            if err:
                failed.append(err)
                print(f"FAILED: {err}", flush=True)
            else:
                infos[run] = info

    for r in runs:
        if r["run"] not in infos:
            continue
        files = [reads_dir / r["run"] / f for f in infos[r["run"]]["files"]]
        try:
            print(f"assigning {r['run']}", flush=True)
            now = assign(r["run"], files, assign_dir / r["run"], args)
            if not now:
                print(f"  already assigned: {r['run']}", flush=True)
        except (OSError, RuntimeError) as e:
            failed.append(str(e))
            print(f"FAILED: {e}", flush=True)

    if failed:
        print(f"\n{len(failed)} run(s) not done, so nothing was written to {stage_dir}/. Run the same command again "
              "to finish them:\n" + "\n".join(failed), file=sys.stderr)
        return 1
    try:
        rows, lowest, t = summarize(runs, infos, assign_dir)
    except (OSError, KeyError, RuntimeError) as e:
        print(f"Nothing was written to {stage_dir}/: {e}", file=sys.stderr)
        return 1
    first = json.loads((assign_dir / runs[0]["run"] / "assign.json").read_text())
    meta = {"made": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC"), "assign": first["made_by"],
            "minimap2": first["minimap2"], "refs_md5": index_info["refs_md5"],
            "fastq_dump": version_of([args.fastq_dump, "--version"]).replace(" : ", " ")}
    write_outputs(stage_dir, rows, lowest, t, meta)
    print(f"\nt = {t:.4g} (lowest DNA-positive fraction {lowest['fraction']:.4g}, in {lowest['run']})")
    print(f"-> {stage_dir}/fractions.tsv, {stage_dir}/STAGE0A.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
