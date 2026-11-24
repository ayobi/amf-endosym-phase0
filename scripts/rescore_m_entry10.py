#!/usr/bin/env python3
"""Stage 0a, step 2 under entry 10: the marker share m again (PROTOCOL.md section 10, entry 10).

Entry 10 changes the index and what the route keeps, so m (entry 8's) is measured again, from step 2's
work files and step 1 run again with the entry 10 index (run_stage0a.py --work-name stage0a_entry10):
  - the symbiont reads: step 2's found reads are aligned again with the entry 10 index, exactly as step 1
    aligns reads; the symbiont's reads are those counted for the run's symbiont. The share's denominator
    is the run's symbiont reads in step 1 with that index.
  - organelle DNA: a symbiont read does not count if its alignments to the organelle genomes alone, at
    90% identity or more, together cover at least half of it, as such reads leave the route in T1.
  - SSU: a read SSUplex calls bacterial counts if T1 would keep its SSU piece (call_lineages.py's
    ssu_label, on entry 8's VSEARCH hits): not outside Bacteria or in an organelle (entry 7), and not
    unplaced (best hits at 94.5% or more in several genera, or a novel piece under 60 bp). A piece that
    matches a genus counts whatever its length, as in T1.
  - SingleM: as in entry 7 (a read counts if SingleM places it in Bacteria; DNA runs only).
m is the lowest share among the DNA positives. RNA runs are reported only: there the SSU reads entry 8
placed (a random sample of at most --rna-sample) are used again if they are still symbiont reads, with
any new symbiont SSU reads, and the share they keep is applied to all of the run's symbiont SSU reads.

Writes to stage0a/ in the repository (commit these): marker_share_entry10.tsv and STAGE0A_M_ENTRY10.md.
Work files go to WORKDIR/stage0a/rescore_entry10/. Entry 8's VSEARCH results
(WORKDIR/stage0a/rescore_entry8/batches/) are used again, since the command is the same: only SSU pieces
entry 8 did not place are placed now. A rerun skips finished batches and runs.

Usage (phase0-reads environment; after run_stage0a.py --work-name stage0a_entry10 --suffix _entry10):
  conda activate phase0-reads
  python scripts/rescore_m_entry10.py --workdir ~/phase0-data --threads 30
"""

import argparse
import csv
import datetime as dt
import gzip
import json
import os
import random
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
import assign_reads as AR          # noqa: E402  the index, the mask and the organelle check
import call_lineages as C          # noqa: E402  the SSU labels of T1
import get_silva as GS             # noqa: E402  VSEARCH version
import rescore_stage0a_m as R      # noqa: E402  entry 8's VSEARCH command, SSU pieces, SingleM domains
import run_stage0a_m as M2         # noqa: E402  step 2's alignment of the found reads

__version__ = "0.1"

FOUND_MIN = 10                     # section 10, entry 6 (the checks in 4.2)
WORK_NAME, ENTRY8_WORK = "rescore_entry10", "rescore_entry8"
OUT_TSV, OUT_MD = "marker_share_entry10.tsv", "STAGE0A_M_ENTRY10.md"
COLS = ["run", "test", "kind", "line", "symbiont_lineage", "sets_m", "symbiont_reads", "found_symbiont",
        "organelle_dna", "ssu_found", "ssu_placed", "ssu_outside_bacteria", "ssu_organelle", "ssu_several_genera",
        "ssu_too_short", "ssu_kept", "singlem_found", "singlem_kept", "route_kept", "share_entry8", "share", "note"]


def now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")


# ------------------------------------------------------------- one run's reads ---

def step1(run, args, ctx):
    """Step 1 with the entry 10 index: the run's counts, once it is clear they match this index and rule."""
    adir = args.assign_dir / run
    a = json.loads((adir / "assign.json").read_text())
    if (a["platform"], a["minimap2"], a["index"]["refs_md5"], a["overlaps_md5"]) != \
            ("short", ctx["minimap2"], ctx["refs_md5"], ctx["overlaps_md5"]):
        raise RuntimeError(f"{adir} was assigned with another platform, minimap2, index or overlaps.tsv")
    return M2.read_counts(adir / "counts.tsv")


def symbiont_reads(r, work, args, ctx, wdir):
    """Step 2's found reads aligned again with the entry 10 index: the symbiont's reads (as (mate, position,
    read, by SSUplex, by SingleM)), the ones that are organelle DNA, and the outcome counts. Cached per run."""
    run, lin = r["run"], r["symbiont_lineage"]
    cache = wdir / "runs" / f"{run}.json"
    if cache.exists():
        c = json.loads(cache.read_text())
        if (c["refs_md5"], c["overlaps_md5"]) == (ctx["refs_md5"], ctx["overlaps_md5"]):
            return [tuple(x) for x in c["symbiont"]], {tuple(x) for x in c["organelle"]}, c["outcomes"]
    rdir = work / "stage0a" / "route" / run
    rows = []
    with gzip.open(rdir / "found.tsv.gz", "rt", encoding="utf-8") as fh:
        for x in csv.DictReader(fh, delimiter="\t"):
            rows.append((x["mate"], int(x["position"]), x["read"], x["ssuplex"] == "1", x["singlem"] == "1"))
    log = wdir / "runs" / f"{run}.minimap2.log"
    outcomes = {}
    for mate in sorted({x[0] for x in rows}):
        for pos, oc in M2.align_found(rdir / f"{mate}.found.fastq.gz", args, ctx, log).items():
            outcomes[(mate, pos)] = oc
    missing = [x for x in rows if (x[0], x[1]) not in outcomes]
    if missing:
        raise RuntimeError(f"{run}: {len(missing)} found reads are not in the found.fastq.gz files")
    sym = [x for x in rows if outcomes[(x[0], x[1])] == ("lineage", lin)]
    tally = Counter(f"{o}: {d}" for o, d in (outcomes[(x[0], x[1])] for x in rows))

    fq = wdir / "runs" / f"{run}.symbiont.fq"
    wanted = defaultdict(set)
    for x in sym:
        wanted[x[0]].add(x[1])
    with open(fq, "wb") as out:
        for mate, positions in sorted(wanted.items()):
            for h, s, q in AR.records(rdir / f"{mate}.found.fastq.gz"):
                pos = int(M2.first_word(h))
                if pos in positions:
                    out.write(b"@%s:%d\n%s\n+\n%s\n" % (mate.encode(), pos, s, q))
    names = AR.organelle_reads(fq, ctx["organelles"], "short", args.minimap2, args.threads, log) if sym else set()
    fq.unlink()
    organelle = {(x[0], x[1]) for x in sym if f"{x[0]}:{x[1]}" in names}
    M2.write_atomic(cache, json.dumps({"refs_md5": ctx["refs_md5"], "overlaps_md5": ctx["overlaps_md5"],
                                       "symbiont": sym, "organelle": sorted(organelle),
                                       "outcomes": dict(sorted(tally.items()))}) + "\n")
    return sym, organelle, dict(tally)


def entry8_placed(r, work, rna_sample):
    """The SSU reads entry 8 placed for an RNA run: its symbiont SSU reads, or a random sample of them."""
    ssu_all = []
    with gzip.open(work / "stage0a" / "route" / r["run"] / "found.tsv.gz", "rt", encoding="utf-8") as fh:
        for x in csv.DictReader(fh, delimiter="\t"):
            if x["outcome"] == "lineage" and x["detail"] == r["symbiont_lineage"] and x["ssuplex"] == "1":
                ssu_all.append((x["mate"], x["read"]))
    ssu_all.sort()
    if rna_sample and len(ssu_all) > rna_sample:
        return set(random.Random(1).sample(ssu_all, rna_sample)), set(ssu_all)
    return set(ssu_all), set(ssu_all)


def gather(r, work, args, ctx, wdir):
    run, lin, dna = r["run"], r["symbiont_lineage"], r["kind"] == "DNA"
    counts = step1(run, args, ctx)
    if ("lineage", lin) not in counts:
        raise RuntimeError(f"{run}: step 1 has no count for {lin}")
    sym, organelle, outcomes = symbiont_reads(r, work, args, ctx, wdir)
    over = [f"{k}: {v} found, {counts.get(tuple(k.split(': ', 1)), 0)} in step 1" for k, v in sorted(outcomes.items())
            if v > counts.get(tuple(k.split(": ", 1)), 0)]
    if over:
        raise RuntimeError(f"{run}: the found reads do not agree with step 1 ({'; '.join(over)})")
    rdir = work / "stage0a" / "route" / run
    ssu = sorted((m, name) for m, p, name, a, _ in sym if a and (m, p) not in organelle)
    sm = sorted((m, name) for m, p, name, _, b in sym if b and (m, p) not in organelle) if dna else []
    placed = ssu
    if not dna:
        sample, before = entry8_placed(r, work, args.rna_sample)
        placed = [k for k in ssu if k in sample or k not in before]
    return {"r": r, "dna": dna, "n": counts[("lineage", lin)], "found": len(sym), "organelle": len(organelle),
            "ssu": ssu, "sm": sm, "placed": placed, "seqs": R.ssu_sequences(rdir, set(placed)),
            "sm_domains": R.singlem_domains(rdir, set(sm)) if sm else {}}


# ------------------------------------------------------------------ VSEARCH ---

def read_batches(bdir):
    """(ids placed, id -> [(target, identity)]) from a folder of finished batches."""
    done, hits = set(), defaultdict(list)
    for fa in sorted(bdir.glob("batch_*.fa")):
        tsv = fa.with_suffix(".tsv")
        if not tsv.exists():
            continue
        done.update(line[1:].strip() for line in open(fa) if line.startswith(">"))
        for line in open(tsv):
            q, t, i = line.rstrip("\n").split("\t")
            hits[q].append((t, float(i)))
    return done, hits


def place(queries, bdir, entry8_dir, args, log):
    """Entry 8's VSEARCH command on every query neither entry 8 nor an earlier run of this script placed;
    returns id -> [(target, identity)], the top hits."""
    bdir.mkdir(parents=True, exist_ok=True)
    for fa in bdir.glob("batch_*.fa"):                       # a batch cut short: do it again
        if not fa.with_suffix(".tsv").exists():
            fa.unlink()
    done8, hits8 = read_batches(entry8_dir) if entry8_dir.exists() else (set(), {})
    done, hits = read_batches(bdir)
    for q in done8:
        if q not in done:
            hits[q] = hits8.get(q, [])
    done |= done8
    todo = [q for q in queries if q not in done]
    k = 1 + max((int(fa.stem.split("_")[1]) for fa in bdir.glob("batch_*.fa")), default=0)
    n_batches = (len(todo) + R.BATCH - 1) // R.BATCH
    print(f"  {len(queries)} distinct SSU pieces: {len(queries) - len(todo)} placed before, {len(todo)} to place",
          flush=True)
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


# ------------------------------------------------------------------ results ---

def score(g, hits, tax):
    r = g["r"]
    v, kept = Counter(), set()
    for k in g["placed"]:
        seq = g["seqs"][k]
        status, _, why = C.ssu_label(hits.get(R.seq_id(seq), []), tax, len(seq))
        v[why if status == C.UNPLACED else status] += 1
        if status == "kept":
            kept.add(k)
    if g["dna"]:
        sm_kept = {k for k in g["sm"] if "d__Bacteria" in g["sm_domains"][k]}
        ssu_kept, route, sm_n = len(kept), len(kept | sm_kept), len(sm_kept)
    else:                                                     # SSU only; a sample may stand for the rest
        frac = len(kept) / len(g["placed"]) if g["placed"] else 0.0
        ssu_kept = round(len(g["ssu"]) * frac)
        sm_n, route = "", ssu_kept
    n = g["n"]
    return {**{k: r[k] for k in ("run", "test", "kind", "line", "symbiont_lineage", "note")}, "sets_m": r["sets_t"],
            "symbiont_reads": n, "found_symbiont": g["found"], "organelle_dna": g["organelle"],
            "ssu_found": len(g["ssu"]), "ssu_placed": len(g["placed"]),
            "ssu_outside_bacteria": v["outside Bacteria"], "ssu_organelle": v["organelle"],
            "ssu_several_genera": v[C.SEVERAL_GENERA], "ssu_too_short": v[C.TOO_SHORT], "ssu_kept": ssu_kept,
            "singlem_found": len(g["sm"]) if g["dna"] else "", "singlem_kept": sm_n, "route_kept": route,
            "share_entry8": g["share_entry8"], "share": route / n if n else None}


def summarize(rows):
    setting = [x for x in rows if x["sets_m"] == "yes"]
    lineages = {x["symbiont_lineage"] for x in setting}
    if not {"CaGg", "CaMg/MRE"} <= lineages:
        raise RuntimeError("the checks need DNA positives for both CaGg and CaMg/MRE")
    if any(x["share"] is None for x in setting):
        raise RuntimeError("a DNA positive has no symbiont reads in step 1 with the entry 10 index: "
                           + ", ".join(x["run"] for x in setting if x["share"] is None))
    lowest_ssu = min(setting, key=lambda x: x["ssu_kept"])
    lowest_sm = min((x for x in setting if x["symbiont_lineage"] == "CaMg/MRE"), key=lambda x: x["singlem_kept"])
    lowest = min(setting, key=lambda x: x["share"])
    ssuplex_ok = lowest_ssu["ssu_kept"] >= FOUND_MIN
    return {"ssuplex_ok": ssuplex_ok, "lowest_ssu": lowest_ssu, "lowest_sm": lowest_sm, "lowest": lowest,
            "singlem_sees_mre": lowest_sm["singlem_kept"] >= FOUND_MIN, "m": lowest["share"] if ssuplex_ok else None}


def fmt(v, spec=".4g"):
    return "n/a" if v is None or v == "" else format(v, spec)


def write_outputs(stage_dir, rows, s, meta):
    tmp = stage_dir / (OUT_TSV + ".part")
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for x in rows:
            w.writerow({**x, "share": "" if x["share"] is None else f"{x['share']:.10g}",
                        "share_entry8": f"{x['share_entry8']:.10g}"})
    lo, lssu, lsm, sil = s["lowest"], s["lowest_ssu"], s["lowest_sm"], meta["silva"]
    md = ["# Stage 0a, step 2 under entry 10: marker share m", "",
          f"Made by rescore_m_entry10.py v{__version__} on {meta['made']} (call_lineages.py v{C.__version__}), from "
          f"step 2's work files, with minimap2 {meta['minimap2']} and the entry 10 index (refs.fna md5 "
          f"{meta['refs_md5']}; `references/INDEX.md`), step 1 again with that index (`fractions_entry10.tsv`), and "
          f"VSEARCH {meta['vsearch']} against SILVA `{sil['file']}` (md5 {sil['md5']}).", "",
          "Rule: PROTOCOL.md 4.1 and 4.2, with section 10, entries 6, 7, 8 and 10. Step 2's found reads were aligned "
          "again with the entry 10 index; a read counts for the symbiont if that alignment counts it so. It is "
          "found by the route if SingleM places it in Bacteria (DNA runs only), or if SSUplex calls it bacterial "
          "and T1 would keep its SSU piece: not outside Bacteria or in an organelle, not in several genera at 94.5% "
          "or more, and at least 60 bp unless it matches a genus. A read whose alignments to the organelle genomes "
          "alone, at 90% identity or more, cover half of it is organelle DNA and is not found. m is the lowest share "
          "among the DNA positives, and it replaces entry 8's m.", ""]
    if s["m"] is not None:
        md.append(f"**m = {s['m']:.4g}** (entry 8: {meta['m_entry8']:.4g}): the lowest DNA-positive share, in "
                  f"{lo['run']} ({lo['line']}, {lo['note']}).")
    else:
        md.append("**m is not set**, because SSUplex failed its check (below). Under section 4.2, phyloFlash replaces "
                  "SSUplex; m has to be measured again with it, after a new deviations entry.")
    md += ["", "Checks (section 4.2), under the same rules. A tool finds a symbiont when at least 10 of the reads it "
               "found count for that symbiont, in every DNA positive (section 10, entry 6):",
           f"- SSUplex finds CaGg and MRE SSU: **{'passed' if s['ssuplex_ok'] else 'failed'}**. Fewest kept symbiont "
           f"SSU reads in a DNA positive: {lssu['ssu_kept']}, in {lssu['run']}.",
           f"- SingleM sees MRE: **{'yes' if s['singlem_sees_mre'] else 'no, so MRE rest on SSU'}**. Fewest MRE reads "
           f"SingleM places in Bacteria, in a T2 positive: {lsm['singlem_kept']}, in {lsm['run']}.", "",
           "| Run | Kind | Symbiont | Sets m | Symbiont reads | Found | Organelle DNA | SSU kept | "
           "SSU unplaced: several genera / under 60 bp | SingleM kept | Share (entry 8) | Share |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for x in rows:
        md.append(f"| {x['run']} | {x['kind']} | {x['symbiont_lineage']} | {x['sets_m']} | {x['symbiont_reads']} | "
                  f"{x['found_symbiont']} | {x['organelle_dna']} | {x['ssu_kept']} | {x['ssu_several_genera']} / "
                  f"{x['ssu_too_short']} | {fmt(x['singlem_kept'], 'd')} | {x['share_entry8']:.4g} | "
                  f"{fmt(x['share'])} |")
    sampled = [x["run"] for x in rows if x["kind"] != "DNA" and x["ssu_placed"] < x["ssu_found"]]
    md += ["", "Symbiont reads: step 1 with the entry 10 index. Found: step 2's found reads (by SSUplex or SingleM) "
               "that the entry 10 index counts for the symbiont. Organelle DNA: those left out by the organelle "
               f"check. Per run, with the SSU pieces dropped by reason: `{OUT_TSV}`. SingleM is not used in RNA runs "
               "(n/a)."]
    if sampled:
        md.append(f"RNA runs where the reads entry 8 placed stand for the rest (--rna-sample {meta['rna_sample']}): "
                  f"{', '.join(sampled)}. There SSU kept applies the kept share of the placed reads to all of the "
                  "run's symbiont SSU reads; the other SSU columns count the placed reads. RNA runs are reported "
                  "only.")
    tmp_md = stage_dir / (OUT_MD + ".part")
    tmp_md.write_text("\n".join(md) + "\n", encoding="utf-8")
    os.replace(tmp, stage_dir / OUT_TSV)
    os.replace(tmp_md, stage_dir / OUT_MD)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--index", help="the entry 10 full index [WORKDIR/index/full]")
    ap.add_argument("--assign-dir", help="step 1 with that index [WORKDIR/stage0a_entry10/assign]")
    ap.add_argument("--refdir", default="references", help="holds overlaps.tsv")
    ap.add_argument("--stage0a-dir", default="stage0a")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--vsearch", default="vsearch")
    ap.add_argument("--rna-sample", type=int, default=500, help="entry 8's value [500]; 0 = all")
    args = ap.parse_args(argv)
    work, stage_dir = Path(args.workdir).expanduser(), Path(args.stage0a_dir)
    index_dir = Path(args.index).expanduser() if args.index else work / "index" / "full"
    args.assign_dir = Path(args.assign_dir).expanduser() if args.assign_dir else work / "stage0a_entry10" / "assign"
    try:
        runs = R.read_tsv(stage_dir / "runs.tsv")
        e8 = {x["run"]: x for x in R.read_tsv(stage_dir / "marker_share_entry8.tsv")}
        if set(e8) != {r["run"] for r in runs}:
            raise RuntimeError("stage0a/marker_share_entry8.tsv and stage0a/runs.tsv list different runs")
        if not (stage_dir / "fractions_entry10.tsv").exists():
            raise RuntimeError(f"{stage_dir}/fractions_entry10.tsv is missing: run run_stage0a.py with --work-name "
                               "stage0a_entry10 --suffix _entry10 first")
        info, lineage_of, _, hosts = AR.load_index(index_dir)
        if info.get("set") != "full":
            raise RuntimeError(f"{index_dir} is the {info.get('set')} index; stage 0a uses the full references")
        if not (index_dir / "organelles.fna").exists():
            raise RuntimeError(f"{index_dir} has no organelles.fna: it was built before entry 10 (make_index.py v0.2)")
        overlaps = Path(args.refdir) / "overlaps.tsv"
        mm2 = subprocess.run([args.minimap2, "--version"], capture_output=True, text=True, check=True).stdout.strip()
        if mm2 != info["minimap2"]:
            raise RuntimeError(f"minimap2 here is {mm2}, but the index was built with {info['minimap2']}")
        silva = json.loads((work / "db" / "silva" / "manifest.json").read_text())
        vsearch = GS.vsearch_version(args.vsearch)
        if vsearch != silva["vsearch"]:
            raise RuntimeError(f"VSEARCH {vsearch} is not the {silva['vsearch']} get_silva.py checked")
        args.silva_fasta = silva["fasta"]
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.CalledProcessError) as e:
        print(f"Cannot start: {e} (get_silva.py writes db/silva/manifest.json)", file=sys.stderr)
        return 1
    ctx = {"mmi": index_dir / AR.INDEX_FILE["short"], "lineage_of": lineage_of, "hosts": hosts,
           "mask": AR.load_mask(overlaps), "minimap2": mm2, "refs_md5": info["refs_md5"],
           "overlaps_md5": AR.md5_of(overlaps), "organelles": index_dir / "organelles.fna"}

    wdir = work / "stage0a" / WORK_NAME
    (wdir / "runs").mkdir(parents=True, exist_ok=True)
    order = [r for r in runs if r["sets_t"] == "yes"] + [r for r in runs if r["sets_t"] != "yes"]
    groups, rows = {}, {}
    try:
        for i, r in enumerate(order, 1):
            print(f"[{i}/{len(order)}] {r['run']} ({r['kind']}, {r['symbiont_lineage']}): aligning the found reads "
                  "again", flush=True)
            g = groups[r["run"]] = gather(r, work, args, ctx, wdir)
            g["share_entry8"] = float(e8[r["run"]]["share"])
            print(f"  symbiont reads {g['n']}; found {g['found']}, organelle DNA {g['organelle']}", flush=True)
        for phase, members in (("DNA", [g for g in groups.values() if g["dna"]]),
                               ("RNA", [g for g in groups.values() if not g["dna"]])):
            queries = {}
            for g in members:
                for k in g["placed"]:
                    queries.setdefault(R.seq_id(g["seqs"][k]), g["seqs"][k])
            print(f"{phase} runs: placing SSU pieces", flush=True)
            hits = place(queries, wdir / "batches", work / "stage0a" / ENTRY8_WORK / "batches", args,
                         wdir / "vsearch.log")
            tax = {k: [x.strip() for x in v.split(";")] for k, v in
                   R.load_taxonomy(silva["taxonomy"], {t for q in queries for t, _ in hits.get(q, [])}).items()}
            for g in members:
                x = rows[g["r"]["run"]] = score(g, hits, tax)
                print(f"  {x['run']}: SSU kept {x['ssu_kept']} of {x['ssu_found']}, SingleM kept "
                      f"{fmt(x['singlem_kept'], 'd')}; share {x['share_entry8']:.4g} -> {fmt(x['share'])}", flush=True)
            if phase == "DNA":
                low = min((x for x in rows.values() if x["sets_m"] == "yes"),
                          key=lambda x: float("inf") if x["share"] is None else x["share"])
                print(f"m from the DNA positives: {fmt(low['share'])} (in {low['run']}); the RNA runs next", flush=True)
        ordered = [rows[r["run"]] for r in runs]
        s = summarize(ordered)
    except (OSError, RuntimeError, KeyError, ValueError) as e:
        print(f"\nStopped, so nothing was written to {stage_dir}/: {e}\nRun the same command again: finished runs "
              "and VSEARCH batches are kept.", file=sys.stderr)
        return 1

    m8 = min(float(x["share"]) for x in e8.values() if x["sets_m"] == "yes" and x["share"])
    meta = {"made": now(), "minimap2": mm2, "refs_md5": info["refs_md5"], "vsearch": vsearch, "silva": silva,
            "rna_sample": args.rna_sample, "m_entry8": m8}
    write_outputs(stage_dir, ordered, s, meta)
    print()
    print(f"m = {s['m']:.4g} (entry 8: {m8:.4g})" if s["m"] is not None else
          "SSUplex failed its check, so m is not set (section 4.2: phyloFlash replaces SSUplex).")
    print(f"SSUplex finds CaGg and MRE SSU: {'passed' if s['ssuplex_ok'] else 'failed'}; SingleM sees MRE: "
          f"{'yes' if s['singlem_sees_mre'] else 'no (MRE rest on SSU)'}")
    print(f"-> {stage_dir}/{OUT_TSV}, {stage_dir}/{OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
