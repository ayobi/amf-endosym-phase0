#!/usr/bin/env python3
"""Stage 0a, step 2: measure the marker share m (PROTOCOL.md 4.1 and 4.2; section 10, entry 6).

For every run in stage0a/runs.tsv, on the same reads as step 1, each mate file on its own:
  1. SSUplex sorts the SSU rRNA reads by origin. The reads it calls bacterial are found, since only
     those go on in the route (4.2).
  2. SingleM looks for marker-gene reads. Every read in its OTU table is found.
     Both tools get the reads as FASTA; neither uses base qualities.
  3. The found reads are aligned again exactly as in step 1 (same minimap2 command, index and read
     names), so each gets the outcome it had in step 1.
Per run, the route's share for the run's known symbiont is
    symbiont reads found by SSUplex or SingleM  /  symbiont reads found by mapping (step 1).

When every run is done, writes to stage0a/ in the repository (commit these):
  marker_share.tsv   per run: symbiont reads, those found by SSUplex, by SingleM and by either, the share
  STAGE0A_M.md       m, the lowest share among the DNA positives, and the two checks in 4.2: SSUplex finds
                     CaGg and MRE SSU (if not, phyloFlash replaces SSUplex and m is not set), and SingleM
                     sees MRE (if not, MRE rest on SSU)
Work files go to WORKDIR/stage0a/route/<run>/. Reruns skip what is done.

Needs step 1 finished, the databases from get_marker_databases.py, the phase0-reads environment
(SSUplex, HMMER, minimap2) and SingleM in the phase0-markers environment.

Usage:
  conda activate phase0-reads
  python scripts/run_stage0a_m.py --workdir ~/phase0-data --threads 16
"""

import argparse
import csv
import datetime as dt
import gzip
import json
import os
import shutil
import subprocess
import sys
import threading
from collections import Counter
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
import assign_reads as AR          # noqa: E402  step 1's reader, minimap2 command and rule
import get_marker_databases as DB  # noqa: E402  tool versions, SSUplex's summary

__version__ = "0.1"

FOUND_MIN = 10   # section 10, entry 6: a tool finds a symbiont when at least 10 of its reads are that symbiont's
SHARE_COLS = ["run", "test", "kind", "line", "symbiont_lineage", "sets_m", "symbiont_reads", "ssu_reads",
              "ssu_symbiont", "singlem_reads", "singlem_symbiont", "both_symbiont", "route_symbiont", "share", "note"]


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_counts(path):
    return {(r["section"], r["name"]): int(r["reads"]) for r in read_tsv(path)}


def first_word(header):
    return header[1:].split(None, 1)[0]


def write_atomic(path, text):
    tmp = path.with_name(path.name + ".part")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def read_names(path):
    return [n for n in path.read_text().split("\n") if n]


def run_logged(cmd, log, what, tmpdir):
    tmpdir.mkdir(parents=True, exist_ok=True)
    with open(log, "a") as fh:
        fh.write("$ " + " ".join(map(str, cmd)) + "\n")
        fh.flush()
        p = subprocess.run(list(map(str, cmd)), stdout=fh, stderr=subprocess.STDOUT,
                           env={**os.environ, "TMPDIR": str(tmpdir)})
    if p.returncode != 0:
        raise RuntimeError(f"{what} failed (exit {p.returncode}); see {log}")


# ------------------------------------------------------------- one mate file ---

def step1_files(run, work, ctx):
    """Step 1's mate files and counts, once it is clear that aligning again here repeats step 1."""
    a1 = json.loads((work / "stage0a" / "assign" / run / "assign.json").read_text())
    if (a1["platform"], a1["minimap2"], a1["index"]["refs_md5"], a1["overlaps_md5"]) != \
            ("short", ctx["minimap2"], ctx["refs_md5"], ctx["overlaps_md5"]):
        raise RuntimeError("step 1 used another platform, minimap2, index or overlaps.tsv than this step would")
    files = []
    for inp in a1["inputs"]:
        p = Path(inp["path"])
        if not p.is_file() or p.stat().st_size != inp["bytes"]:
            raise RuntimeError(f"{p} is missing or is not the file step 1 counted")
        files.append((p, inp["reads"]))
    return files, read_counts(work / "stage0a" / "assign" / run / "counts.tsv")


def mate_fasta(fq, fa, expected):
    """Write one mate file as FASTA for SSUplex and SingleM (name = first word of the header). Reads of
    length 0, if any, are left out: there is nothing in them to find. Returns the reads written."""
    n = written = 0
    tmp = fa.with_name(fa.name + ".part")
    with open(tmp, "wb") as out:
        for h, s, _ in AR.records(fq):
            n += 1
            if s:
                out.write(b">" + first_word(h) + b"\n" + s + b"\n")
                written += 1
    if n != expected:
        tmp.unlink()
        raise RuntimeError(f"{fq.name} holds {n} reads, but step 1 counted {expected}")
    os.replace(tmp, fa)
    return written


def ssuplex_names(fa, n_fa, stem, rdir, args, db, log):
    """Names of the reads SSUplex calls bacterial, and its read counts per origin."""
    done = rdir / f"{stem}.ssuplex.json"
    if done.exists():
        return read_names(rdir / f"{stem}.ssu_names.txt"), json.loads(done.read_text())
    print(f"  {stem}: SSUplex", flush=True)
    run_logged([args.ssuplex, "-i", fa, "-o", rdir / f"{stem}.ssuplex", "--hmm-dir", db["ssuplex_hmm_dir"],
                "--threads", args.threads, "--force"], log, f"SSUplex on {stem}", rdir / "tmp")
    summary = DB.parse_ssuplex_summary(rdir / f"{stem}.ssuplex.summary.txt")
    if summary["total"] != n_fa:
        raise RuntimeError(f"SSUplex read {summary['total']} of the {n_fa} reads of {stem}")
    bacteria, names = rdir / f"{stem}.ssuplex.bacteria.fasta", []
    if bacteria.exists():
        with open(bacteria, "rb") as fh:
            names = [first_word(line).decode() for line in fh if line.startswith(b">")]
    if len(names) != summary["per_origin"].get("bacteria", 0) or len(set(names)) != len(names):
        raise RuntimeError(f"SSUplex's bacteria.fasta for {stem} does not match its summary")
    for other in ("unclassified", "archaea", "eukaryota", "mitochondria", "chloroplast"):   # not needed; the
        (rdir / f"{stem}.ssuplex.{other}.fasta").unlink(missing_ok=True)                   # first is nearly all reads
    (rdir / f"{stem}.ssuplex.extraction.tsv").unlink(missing_ok=True)
    write_atomic(rdir / f"{stem}.ssu_names.txt", "".join(x + "\n" for x in names))
    info = {"reads": summary["total"], "per_origin": summary["per_origin"]}
    write_atomic(done, json.dumps(info) + "\n")
    return names, info


def singlem_names(fa, stem, rdir, args, db, log):
    """Names of the reads in SingleM's OTU table (its read_names column)."""
    done = rdir / f"{stem}.singlem_names.txt"
    if done.exists():
        return read_names(done)
    print(f"  {stem}: SingleM", flush=True)
    otu = rdir / f"{stem}.singlem.otu.tsv"
    otu.unlink(missing_ok=True)                  # a partial table from an interrupted run
    run_logged([args.conda, "run", "-n", args.singlem_env, "--no-capture-output", "singlem", "pipe", "--forward", fa,
                "--otu-table", otu, "--output-extras", "--metapackage", db["singlem_metapackage"],
                "--threads", args.threads], log, f"SingleM on {stem}", rdir / "tmp")
    if not otu.exists():
        raise RuntimeError(f"SingleM wrote no OTU table for {stem}; see {log}")
    names = set()
    with open(otu, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if "read_names" not in header:
            raise RuntimeError(f"{otu.name} has no read_names column")
        col = header.index("read_names")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > col:
                names.update(f[col].split())
    with open(otu, "rb") as src, gzip.open(str(otu) + ".gz", "wb") as dst:
        shutil.copyfileobj(src, dst)
    otu.unlink()
    names = sorted(names)
    write_atomic(done, "".join(x + "\n" for x in names))
    return names


def collect_found(fq, stem, rdir, ssu, sm, expected):
    """One pass over a mate file: copy out the reads either tool found, named by their position in the
    file (as step 1 named them for minimap2), and note which tool found each."""
    done = rdir / f"{stem}.found.tsv"
    if done.exists():
        return [(int(r["position"]), r["read"], int(r["ssuplex"]), int(r["singlem"])) for r in read_tsv(done)]
    ssu, sm = set(ssu), set(sm)
    rows, seen, n = [], set(), 0
    tmp = rdir / f"{stem}.found.fastq.gz.part"
    with gzip.open(tmp, "wb", compresslevel=1) as out:
        for pos, (h, s, q) in enumerate(AR.records(fq)):
            n += 1
            if q is None:
                raise RuntimeError(f"{fq.name} is not FASTQ")
            name = first_word(h).decode("ascii", "replace")
            a, b = name in ssu, name in sm
            if a or b:
                seen.add(name)
                out.write(b"@%d %s\n%s\n+\n%s\n" % (pos, h[1:], s, q))
                rows.append((pos, name, int(a), int(b)))
    missing = sorted((ssu | sm) - seen)
    if n != expected or missing:
        tmp.unlink()
        raise RuntimeError(f"{fq.name}: {n} reads (step 1 counted {expected}); {len(missing)} names from the tools "
                           f"not in it{': ' + ', '.join(missing[:3]) if missing else ''}")
    os.replace(tmp, rdir / f"{stem}.found.fastq.gz")
    write_atomic(done, "position\tread\tssuplex\tsinglem\n" + "".join(f"{p}\t{x}\t{a}\t{b}\n" for p, x, a, b in rows))
    return rows


def align_found(found_fq, args, ctx, log):
    """Align the found reads exactly as step 1 did; return {position: (outcome, detail)}."""
    cmd = [args.minimap2, *AR.MINIMAP2_ARGS["short"], "-t", str(args.threads), str(ctx["mmi"]), "-"]
    with open(log, "a") as err:
        err.write("$ " + " ".join(cmd) + f" < {found_fq.name} (reads named by position)\n")
        err.flush()
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=err, bufsize=1 << 20)
    given, state = [], {}

    def feed():
        try:
            for h, s, q in AR.records(found_fq):
                pos = int(first_word(h))
                given.append(pos)
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
    outcomes, last = {}, -1
    try:
        for pos, alns in AR.parse_paf(proc.stdout):
            if pos <= last:
                raise RuntimeError("minimap2 output is not in input order")
            last = pos
            try:
                outcomes[pos] = AR.evaluate(alns, ctx["lineage_of"], ctx["hosts"], ctx["mask"])[:2]
            except KeyError as e:
                raise RuntimeError(f"minimap2 aligned a read to {e}, which is not in sequences.tsv") from None
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
    if not set(outcomes) <= set(given):
        raise RuntimeError("minimap2 reported reads it was not given")
    return {pos: outcomes.get(pos, ("routed", AR.ROUTED_REASONS[0])) for pos in given}


# ------------------------------------------------------------------ one run ---

def one_run(r, work, args, db, ctx):
    run, lin = r["run"], r["symbiont_lineage"]
    rdir = work / "stage0a" / "route" / run
    if (rdir / "route.json").exists():
        return json.loads((rdir / "route.json").read_text())["row"]
    files, step1 = step1_files(run, work, ctx)
    if ("lineage", lin) not in step1:
        raise RuntimeError(f"step 1 has no count for {lin}")
    rdir.mkdir(parents=True, exist_ok=True)
    log = rdir / "tools.log"
    found, per_origin = {}, {}
    for fq, expected in files:
        stem = fq.name.split(".")[0]
        fa, n_fa = rdir / f"{stem}.reads.fasta", None
        if not ((rdir / f"{stem}.ssuplex.json").exists() and (rdir / f"{stem}.singlem_names.txt").exists()):
            n_fa = mate_fasta(fq, fa, expected)
        ssu, info = ssuplex_names(fa, n_fa, stem, rdir, args, db, log)
        sm = singlem_names(fa, stem, rdir, args, db, log)
        fa.unlink(missing_ok=True)
        found[stem] = collect_found(fq, stem, rdir, ssu, sm, expected)
        per_origin[stem] = info["per_origin"]

    tally, c, table = Counter(), Counter(), []
    for stem, rows in found.items():
        outcomes = align_found(rdir / f"{stem}.found.fastq.gz", args, ctx, log) if rows else {}
        for pos, name, a, b in rows:
            outcome, detail = outcomes[pos]
            tally[(outcome, detail)] += 1
            table.append(f"{stem}\t{pos}\t{name}\t{a}\t{b}\t{outcome}\t{detail}\n")
            c["ssu_reads"] += a
            c["singlem_reads"] += b
            if outcome == "lineage" and detail == lin:
                c["ssu_symbiont"] += a
                c["singlem_symbiont"] += b
                c["both_symbiont"] += a & b
                c["route_symbiont"] += 1
    over = [f"{s} '{d}': {k} found, {step1.get((s, d), 0)} in step 1" for (s, d), k in sorted(tally.items())
            if k > step1.get((s, d), 0)]
    if over:
        raise RuntimeError("the found reads do not agree with step 1 (" + "; ".join(over) + ")")

    n = step1[("lineage", lin)]
    row = {**{k: r[k] for k in ("run", "test", "kind", "line", "symbiont_lineage", "note")}, "sets_m": r["sets_t"],
           "symbiont_reads": n, **{k: c[k] for k in ("ssu_reads", "ssu_symbiont", "singlem_reads", "singlem_symbiont",
                                                      "both_symbiont", "route_symbiont")},
           "share": c["route_symbiont"] / n if n else None}
    with gzip.open(rdir / "found.tsv.gz", "wt", encoding="utf-8") as fh:
        fh.write("mate\tposition\tread\tssuplex\tsinglem\toutcome\tdetail\n" + "".join(table))
    write_atomic(rdir / "route.json", json.dumps({
        "made_by": f"run_stage0a_m.py v{__version__}",
        "finished": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "row": row, "ssuplex_reads_per_origin": per_origin,
        "found_reads_by_outcome": {f"{s}: {d}": k for (s, d), k in sorted(tally.items())},
        "mate_files": [str(fq) for fq, _ in files]}, indent=1) + "\n")
    shutil.rmtree(rdir / "tmp", ignore_errors=True)
    return row


# ------------------------------------------------------------------ summary ---

def summarize(rows):
    setting = [x for x in rows if x["sets_m"] == "yes"]
    lineages = {x["symbiont_lineage"] for x in setting}
    if not {"CaGg", "CaMg/MRE"} <= lineages:
        raise RuntimeError("the checks need DNA positives for both CaGg and CaMg/MRE")
    lowest_ssu = min(setting, key=lambda x: x["ssu_symbiont"])
    mre = [x for x in setting if x["symbiont_lineage"] == "CaMg/MRE"]
    lowest_sm = min(mre, key=lambda x: x["singlem_symbiont"])
    lowest = min(setting, key=lambda x: x["share"])
    ssuplex_ok = lowest_ssu["ssu_symbiont"] >= FOUND_MIN
    return {"ssuplex_ok": ssuplex_ok, "lowest_ssu": lowest_ssu, "lowest_sm": lowest_sm, "lowest": lowest,
            "singlem_sees_mre": lowest_sm["singlem_symbiont"] >= FOUND_MIN,
            "m": lowest["share"] if ssuplex_ok else None}


def fmt_share(v, spec=".4g"):
    return "n/a" if v is None else format(v, spec)


def write_outputs(stage_dir, rows, s, meta, db):
    tmp = stage_dir / "marker_share.tsv.part"
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=SHARE_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for x in rows:
            w.writerow({**x, "share": "" if x["share"] is None else f"{x['share']:.10g}"})
    v, lo, lssu, lsm = db["versions"], s["lowest"], s["lowest_ssu"], s["lowest_sm"]
    md = ["# Stage 0a, step 2: marker share m", "",
          f"Made by run_stage0a_m.py v{__version__} on {meta['made']}, with SSUplex {v['ssuplex']} (HMMER {v['hmmer']}, "
          f"profiles from Metaxa2 2.2.3), SingleM {v['singlem']} with metapackage `{db['singlem_metapackage_name']}`, and "
          f"minimap2 {meta['minimap2']} against the full index (refs.fna md5 {meta['refs_md5']}). Reads: the same as "
          "step 1, each mate file on its own.", "",
          "Rule: PROTOCOL.md 4.1 (marker share) and 4.2 (checks), with section 10, entry 6. For each run, the share of the "
          "known symbiont's reads that the route finds: the reads SSUplex calls bacterial SSU, plus the reads in "
          "SingleM's OTU table. The found reads were aligned again exactly as in step 1. m is the lowest share among "
          "the DNA positives.", ""]
    if s["m"] is not None:
        md.append(f"**m = {s['m']:.4g}**: the lowest DNA-positive share, in {lo['run']} ({lo['line']}, {lo['note']}).")
    else:
        md.append("**m is not set**, because SSUplex failed its check (below). Under section 4.2, phyloFlash replaces "
                  "SSUplex; m has to be measured again with it, after a new deviations entry.")
    md += ["", "Checks (section 4.2). A tool finds a symbiont when at least 10 of the reads it found map to that "
               "symbiont, in every DNA positive (section 10, entry 6):",
           f"- SSUplex finds CaGg and MRE SSU: **{'passed' if s['ssuplex_ok'] else 'failed'}**. Fewest symbiont reads "
           f"among the bacterial SSU reads of a DNA positive: {lssu['ssu_symbiont']}, in {lssu['run']}.",
           f"- SingleM sees MRE: **{'yes' if s['singlem_sees_mre'] else 'no, so MRE rest on SSU'}**. Fewest MRE reads "
           f"in SingleM's table for a T2 positive: {lsm['singlem_symbiont']}, in {lsm['run']}.", "",
           "| Run | Test | Kind | Symbiont | Sets m | Symbiont reads | By SSUplex | By SingleM | By either | Share |",
           "|---|---|---|---|---|---|---|---|---|---|"]
    for x in rows:
        md.append(f"| {x['run']} | {x['test']} | {x['kind']} | {x['symbiont_lineage']} | {x['sets_m']} | "
                  f"{x['symbiont_reads']} | {x['ssu_symbiont']} | {x['singlem_symbiont']} | {x['route_symbiont']} | "
                  f"{fmt_share(x['share'])} |")
    hm, qs = db["ssuplex_hmms"], db["ssuplex_quickstart"]
    md += ["", "Symbiont reads: those mapping counted for the run's symbiont in step 1. By SSUplex, by SingleM: the "
               "symbiont reads among the reads each tool found. By either: a read found by both counts once. All reads "
               "found by each tool, symbiont or not: `marker_share.tsv`.", "",
           "Databases (`get_marker_databases.py`):",
           f"- SingleM metapackage `{db['singlem_metapackage_name']}` from doi:{db['singlem_metapackage_doi']}, md5 "
           f"{db['singlem_metapackage_md5']} as in `references/DATABASES.md`, checked with `singlem data --verify-only`.",
           f"- Metaxa2 2.2.3 ({db['metaxa2_url']}, md5 {db['metaxa2_md5']}), made into SSUplex's profile files: " +
           "; ".join(f"{o} ({hm[o]['from']}, {hm[o]['profiles']} profile{'s' * (hm[o]['profiles'] != 1)}, "
                     f"md5 {hm[o]['md5']})" for o in hm) + ".",
           f"- SSUplex quick start on Metaxa2's test.fasta: {qs['total'] - qs['per_origin']['unclassified']} of "
           f"{qs['total']} classified, as in the SSUplex README."]
    tmp_md = stage_dir / "STAGE0A_M.md.part"
    tmp_md.write_text("\n".join(md) + "\n", encoding="utf-8")
    os.replace(tmp, stage_dir / "marker_share.tsv")
    os.replace(tmp_md, stage_dir / "STAGE0A_M.md")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--index", help="the full index set [WORKDIR/index/full]")
    ap.add_argument("--refdir", default="references", help="holds overlaps.tsv")
    ap.add_argument("--stage0a-dir", default="stage0a")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--ssuplex", default="ssuplex")
    ap.add_argument("--conda", default=os.environ.get("CONDA_EXE") or shutil.which("conda"))
    ap.add_argument("--singlem-env", default="phase0-markers")
    args = ap.parse_args(argv)
    work, stage_dir = Path(args.workdir).expanduser(), Path(args.stage0a_dir)
    index_dir = Path(args.index).expanduser() if args.index else work / "index" / "full"

    try:
        runs = read_tsv(stage_dir / "runs.tsv")
        if not (stage_dir / "fractions.tsv").exists():
            raise RuntimeError(f"{stage_dir}/fractions.tsv is missing: finish step 1 first")
        manifest = work / "db" / "manifest.json"
        if not manifest.exists():
            raise RuntimeError(f"{manifest} is missing: run get_marker_databases.py first")
        db = json.loads(manifest.read_text())
        info, lineage_of, _, hosts = AR.load_index(index_dir)
        if info.get("set") != "full":
            raise RuntimeError(f"{index_dir} is the {info.get('set')} index; stage 0a uses the full references")
        overlaps = Path(args.refdir) / "overlaps.tsv"
        mask = AR.load_mask(overlaps)
        mm2 = subprocess.run([args.minimap2, "--version"], capture_output=True, text=True, check=True).stdout.strip()
        if mm2 != info["minimap2"]:
            raise RuntimeError(f"minimap2 here is {mm2}, but the index was built with {info['minimap2']}")
        versions = DB.check_versions(args)
        if versions != db["versions"]:
            raise RuntimeError(f"the tools here ({versions}) are not those get_marker_databases.py checked "
                               f"({db['versions']})")
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.CalledProcessError) as e:
        print(f"Cannot start: {e}", file=sys.stderr)
        return 1
    ctx = {"mmi": index_dir / AR.INDEX_FILE["short"], "lineage_of": lineage_of, "hosts": hosts, "mask": mask,
           "minimap2": mm2, "refs_md5": info["refs_md5"], "overlaps_md5": AR.md5_of(overlaps)}

    order = [r for r in runs if r["sets_t"] == "yes"] + [r for r in runs if r["sets_t"] != "yes"]
    done, failed = {}, []
    for i, r in enumerate(order, 1):
        print(f"[{i}/{len(order)}] {r['run']} ({r['kind']}, {r['symbiont_lineage']})", flush=True)
        try:
            x = done[r["run"]] = one_run(r, work, args, db, ctx)
            print(f"  symbiont reads {x['symbiont_reads']}: found by SSUplex {x['ssu_symbiont']}, by SingleM "
                  f"{x['singlem_symbiont']}, by either {x['route_symbiont']}; share {fmt_share(x['share'])}", flush=True)
        except (OSError, RuntimeError, KeyError, ValueError) as e:
            failed.append(f"{r['run']}: {e}")
            print(f"  FAILED: {e}", flush=True)
    if failed:
        print(f"\n{len(failed)} run(s) not done, so nothing was written to {stage_dir}/. Run the same command again to "
              "finish them:\n" + "\n".join(failed), file=sys.stderr)
        return 1
    rows = [done[r["run"]] for r in runs]
    try:
        s = summarize(rows)
    except RuntimeError as e:
        print(f"Nothing was written to {stage_dir}/: {e}", file=sys.stderr)
        return 1
    meta = {"made": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC"), "minimap2": mm2,
            "refs_md5": info["refs_md5"]}
    write_outputs(stage_dir, rows, s, meta, db)
    print()
    if s["m"] is not None:
        print(f"m = {s['m']:.4g} (lowest DNA-positive share, in {s['lowest']['run']})")
    else:
        print("SSUplex failed its check, so m is not set (section 4.2: phyloFlash replaces SSUplex).")
    print(f"SSUplex finds CaGg and MRE SSU: {'passed' if s['ssuplex_ok'] else 'failed'}; SingleM sees MRE: "
          f"{'yes' if s['singlem_sees_mre'] else 'no (MRE rest on SSU)'}")
    print(f"-> {stage_dir}/marker_share.tsv, {stage_dir}/STAGE0A_M.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
