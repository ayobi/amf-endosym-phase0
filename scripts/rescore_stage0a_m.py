#!/usr/bin/env python3
"""Stage 0a, step 2 rescored under entry 7's route rules: the marker share m again, from step 2's work
files (PROTOCOL.md section 10, entry 7).

Step 2 counted a read as found by the route if SSUplex called it bacterial or it was in SingleM's OTU
table. Entry 7 narrows that:
  - SingleM: a read counts only if SingleM places it in Bacteria. In RNA runs SingleM is not used.
  - SSU: the reads SSUplex calls bacterial are placed with VSEARCH against SILVA 144 SSU Ref NR99
    (PROTOCOL.md 4.1, plus --strand both). A read is dropped if one of its top hits lies outside
    Bacteria or in chloroplasts or mitochondria. Reads with no hit at 75% or more stay.
The reads and their mapping outcomes are step 2's (found.tsv.gz); only the found sets get smaller. Only
the symbiont's reads are placed, since m counts only those. In RNA runs, which are reported only, at
most --rna-sample of them are placed (a random sample, seed 1), and the share they keep is applied.

DNA runs are placed first; the log gives m as soon as they are done. When everything is done, writes
to stage0a/ in the repository (commit these): marker_share_entry7.tsv and STAGE0A_M_ENTRY7.md.
Work files go to WORKDIR/stage0a/rescore/; a rerun skips finished VSEARCH batches.

Usage (phase0-reads environment; needs step 2's work files and get_silva.py's database):
  conda activate phase0-reads
  python scripts/rescore_stage0a_m.py --workdir ~/phase0-data --threads 16
"""

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import random
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
import get_silva as GS  # noqa: E402  vsearch_version

__version__ = "0.1"

FOUND_MIN = 10                     # section 10, entry 6 (the checks in 4.2)
VSEARCH_ARGS = ["--id", "0.75", "--maxaccepts", "0", "--maxrejects", "0", "--top_hits_only",  # 4.1
                "--strand", "both"]                                                            # entry 7
BATCH = 500
ORGANELLE = re.compile(r"^(chloroplast|mitochondri)", re.IGNORECASE)
KEEP = ("Bacteria", "no hit")
COLS = ["run", "test", "kind", "line", "symbiont_lineage", "sets_m", "symbiont_reads", "ssu_step2", "ssu_placed",
        "ssu_outside_bacteria", "ssu_organelle", "ssu_no_hit", "ssu_kept", "singlem_step2", "singlem_kept", "route_kept",
        "share_step2", "share", "note"]


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def seq_id(seq):
    return hashlib.sha1(seq.upper().encode()).hexdigest()[:20]


# ------------------------------------------------------------ step 2's files ---

def singlem_domains(rdir, wanted):
    """(mate, read) -> the domains SingleM gives the read ('' when it gives none), for the wanted reads."""
    domains = defaultdict(set)
    for f in sorted(rdir.glob("*.singlem.otu.tsv.gz")):
        stem = f.name.split(".")[0]
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            i, t = header.index("read_names"), header.index("taxonomy")
            for line in fh:
                row = line.rstrip("\n").split("\t")
                parts = row[t].split("; ")
                domain = parts[1] if len(parts) > 1 else ""
                for name in row[i].split():
                    if (stem, name) in wanted:
                        domains[(stem, name)].add(domain)
    missing = wanted - domains.keys()
    if missing:
        raise RuntimeError(f"{len(missing)} SingleM reads are not in its OTU tables in {rdir}")
    return domains


def ssu_sequences(rdir, wanted):
    """(mate, read) -> the SSU segment SSUplex cut from the read, for the wanted reads."""
    seqs = {}
    for mate in sorted({m for m, _ in wanted}):
        name, buf = None, []
        with open(rdir / f"{mate}.ssuplex.bacteria.fasta", encoding="ascii") as fh:
            for line in list(fh) + [">"]:
                if line.startswith(">"):
                    if name is not None and (mate, name) in wanted:
                        seqs[(mate, name)] = "".join(buf)
                    name, buf = (line[1:].split(None, 1) or [None])[0], []
                else:
                    buf.append(line.strip())
    missing = wanted - seqs.keys()
    if missing:
        raise RuntimeError(f"{len(missing)} SSU reads are not in SSUplex's bacteria.fasta in {rdir}")
    return seqs


def gather(r, work, rna_sample):
    """The run's symbiont reads from step 2, with what each tool said about them."""
    rdir = work / "stage0a" / "route" / r["run"]
    step2 = json.loads((rdir / "route.json").read_text())["row"]
    sym = []
    with gzip.open(rdir / "found.tsv.gz", "rt", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["outcome"] == "lineage" and row["detail"] == r["symbiont_lineage"]:
                sym.append(((row["mate"], row["read"]), row["ssuplex"] == "1", row["singlem"] == "1"))
    ssu_all = sorted(k for k, a, _ in sym if a)
    sm_all = {k for k, _, b in sym if b}
    if (len(sym), len(ssu_all), len(sm_all)) != (step2["route_symbiont"], step2["ssu_symbiont"], step2["singlem_symbiont"]):
        raise RuntimeError(f"{r['run']}: found.tsv.gz does not match step 2's route.json")
    dna = r["kind"] == "DNA"
    placed = ssu_all
    if not dna and rna_sample and len(ssu_all) > rna_sample:
        placed = sorted(random.Random(1).sample(ssu_all, rna_sample))
    return {"r": r, "step2": step2, "dna": dna, "ssu_all": ssu_all, "placed": placed,
            "seqs": ssu_sequences(rdir, set(placed)),
            "sm_domains": singlem_domains(rdir, sm_all) if dna and sm_all else {}, "sm_all": sm_all}


# ------------------------------------------------------------------ VSEARCH ---

def place(queries, bdir, args, log):
    """VSEARCH on every query sequence not placed yet, in batches; returns id -> [target, ...] (top hits)."""
    bdir.mkdir(parents=True, exist_ok=True)
    for fa in bdir.glob("batch_*.fa"):                        # a batch cut short: do it again
        if not fa.with_suffix(".tsv").exists():
            fa.unlink()
    done, hits = set(), defaultdict(list)
    for fa in sorted(bdir.glob("batch_*.fa")):
        done.update(line[1:].strip() for line in open(fa) if line.startswith(">"))
        for line in open(fa.with_suffix(".tsv")):
            q, t, _ = line.rstrip("\n").split("\t")
            hits[q].append(t)
    todo = [q for q in queries if q not in done]
    k = 1 + max((int(fa.stem.split("_")[1]) for fa in bdir.glob("batch_*.fa")), default=0)
    n_batches = (len(todo) + BATCH - 1) // BATCH
    for i, start in enumerate(range(0, len(todo), BATCH), 1):
        chunk = todo[start:start + BATCH]
        fa, out = bdir / f"batch_{k:04d}.fa", bdir / f"batch_{k:04d}.tsv"
        fa.write_text("".join(f">{q}\n{queries[q]}\n" for q in chunk))
        part = out.with_name(out.name + ".part")
        cmd = [args.vsearch, "--usearch_global", str(fa), "--db", args.silva_fasta, *VSEARCH_ARGS,
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
            q, t, _ = line.rstrip("\n").split("\t")
            hits[q].append(t)
        secs = time.time() - t0
        print(f"  batch {i}/{n_batches}: {len(chunk)} sequences in {secs / 60:.1f} min "
              f"({secs / len(chunk):.2f} s each); about {secs * (n_batches - i) / 3600:.1f} h left", flush=True)
        k += 1
    return hits


def load_taxonomy(path, targets):
    tax = {}
    with open(path, encoding="ascii", errors="replace") as fh:
        next(fh)
        for line in fh:
            acc, _, path_ = line.rstrip("\n").partition("\t")
            if acc in targets:
                tax[acc] = path_
    missing = targets - tax.keys()
    if missing:
        raise RuntimeError(f"{len(missing)} VSEARCH targets are not in taxonomy.tsv")
    return tax


def verdict(targets, tax):
    """'Bacteria', 'no hit', 'outside Bacteria' or 'organelle', from a read's top hits."""
    if not targets:
        return "no hit"
    paths = [[p.strip() for p in tax[t].split(";")] for t in targets]
    if any(p[0] != "Bacteria" for p in paths):
        return "outside Bacteria"
    if any(ORGANELLE.match(x) for p in paths for x in p):
        return "organelle"
    return "Bacteria"


# ------------------------------------------------------------------ results ---

def score(g, verdicts):
    r, s2 = g["r"], g["step2"]
    v = Counter(verdicts[seq_id(g["seqs"][k])] for k in g["placed"])
    kept_placed = {k for k in g["placed"] if verdicts[seq_id(g["seqs"][k])] in KEEP}
    if g["dna"]:
        sm_kept = {k for k in g["sm_all"] if "d__Bacteria" in g["sm_domains"][k]}
        ssu_kept, route = len(kept_placed), len(kept_placed | sm_kept)
        sm_kept = len(sm_kept)
    else:                                                     # SSU only; a sample may stand for the rest
        frac = len(kept_placed) / len(g["placed"]) if g["placed"] else 0.0
        ssu_kept = round(len(g["ssu_all"]) * frac)
        sm_kept, route = "", ssu_kept
    n = s2["symbiont_reads"]
    return {**{k: r[k] for k in ("run", "test", "kind", "line", "symbiont_lineage", "note")}, "sets_m": r["sets_t"],
            "symbiont_reads": n, "ssu_step2": s2["ssu_symbiont"], "ssu_placed": len(g["placed"]),
            "ssu_outside_bacteria": v["outside Bacteria"], "ssu_organelle": v["organelle"], "ssu_no_hit": v["no hit"],
            "ssu_kept": ssu_kept,
            "singlem_step2": s2["singlem_symbiont"], "singlem_kept": sm_kept, "route_kept": route,
            "share_step2": s2["share"], "share": route / n if n else None}


def summarize(rows):
    setting = [x for x in rows if x["sets_m"] == "yes"]
    lowest_ssu = min(setting, key=lambda x: x["ssu_kept"])
    lowest_sm = min((x for x in setting if x["symbiont_lineage"] == "CaMg/MRE"), key=lambda x: x["singlem_kept"])
    lowest = min(setting, key=lambda x: x["share"])
    ssuplex_ok = lowest_ssu["ssu_kept"] >= FOUND_MIN
    return {"ssuplex_ok": ssuplex_ok, "lowest_ssu": lowest_ssu, "lowest_sm": lowest_sm, "lowest": lowest,
            "singlem_sees_mre": lowest_sm["singlem_kept"] >= FOUND_MIN, "m": lowest["share"] if ssuplex_ok else None}


def fmt(v, spec=".4g"):
    return "n/a" if v is None or v == "" else format(v, spec)


def write_outputs(stage_dir, rows, s, meta):
    tmp = stage_dir / "marker_share_entry7.tsv.part"
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for x in rows:
            w.writerow({**x, "share": "" if x["share"] is None else f"{x['share']:.10g}",
                        "share_step2": f"{x['share_step2']:.10g}"})
    lo, lssu, lsm, sil = s["lowest"], s["lowest_ssu"], s["lowest_sm"], meta["silva"]
    sampled = [x["run"] for x in rows if x["kind"] != "DNA" and x["ssu_placed"] < x["ssu_step2"]]
    md = ["# Stage 0a, step 2 rescored under entry 7: marker share m", "",
          f"Made by rescore_stage0a_m.py v{__version__} on {meta['made']}, from step 2's work files (run_stage0a_m.py "
          f"v0.1), with VSEARCH {meta['vsearch']} against SILVA `{sil['file']}` (md5 {sil['md5']}, "
          f"{sil['sequences']:,} sequences).", "",
          "Rule: PROTOCOL.md 4.1 and 4.2, with section 10, entries 6 and 7. A symbiont read counts as found by the route "
          "if SingleM places it in Bacteria (DNA runs only), or if SSUplex calls it bacterial and none of its top SILVA "
          "hits lies outside Bacteria or in chloroplasts or mitochondria (reads with no hit at 75% or more stay). The "
          "reads and their mapping outcomes are step 2's. m is the lowest share among the DNA positives, and it "
          "replaces step 2's m.", ""]
    if s["m"] is not None:
        md.append(f"**m = {s['m']:.4g}** (step 2: {meta['m_step2']:.4g}): the lowest DNA-positive share, in "
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
           "| Run | Kind | Symbiont | Sets m | Symbiont reads | SSU (step 2) | SSU kept | SingleM (step 2) | "
           "SingleM kept | Share (step 2) | Share |",
           "|---|---|---|---|---|---|---|---|---|---|---|"]
    for x in rows:
        md.append(f"| {x['run']} | {x['kind']} | {x['symbiont_lineage']} | {x['sets_m']} | {x['symbiont_reads']} | "
                  f"{x['ssu_step2']} | {x['ssu_kept']} | {x['singlem_step2']} | {fmt(x['singlem_kept'], 'd')} | "
                  f"{x['share_step2']:.4g} | {fmt(x['share'])} |")
    out_b = sum(x["ssu_outside_bacteria"] for x in rows)
    org = sum(x["ssu_organelle"] for x in rows)
    no_hit = sum(x["ssu_no_hit"] for x in rows)
    md += ["", f"Of the symbiont SSU reads placed, {out_b} had a top hit outside Bacteria and {org} a top hit in "
               f"chloroplasts or mitochondria (dropped), and {no_hit} had no hit at 75% or more (kept); "
               "`marker_share_entry7.tsv` gives them per run. SingleM is not used in RNA runs (n/a)."]
    if sampled:
        md.append(f"RNA runs with more symbiont SSU reads than --rna-sample ({meta['rna_sample']}): "
                  f"{', '.join(sampled)}. There a random sample was placed (seed 1), and SSU kept applies its kept share "
                  "to all of the run's symbiont SSU reads. RNA runs are reported only.")
    tmp_md = stage_dir / "STAGE0A_M_ENTRY7.md.part"
    tmp_md.write_text("\n".join(md) + "\n", encoding="utf-8")
    os.replace(tmp, stage_dir / "marker_share_entry7.tsv")
    os.replace(tmp_md, stage_dir / "STAGE0A_M_ENTRY7.md")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--stage0a-dir", default="stage0a")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--vsearch", default="vsearch")
    ap.add_argument("--rna-sample", type=int, default=500, help="most symbiont SSU reads placed per RNA run; 0 = all")
    args = ap.parse_args(argv)
    work, stage_dir = Path(args.workdir).expanduser(), Path(args.stage0a_dir)
    try:
        runs = read_tsv(stage_dir / "runs.tsv")
        step2 = {x["run"]: x for x in read_tsv(stage_dir / "marker_share.tsv")}
        silva = json.loads((work / "db" / "silva" / "manifest.json").read_text())
        vsearch = GS.vsearch_version(args.vsearch)
        if vsearch != silva["vsearch"]:
            raise RuntimeError(f"VSEARCH {vsearch} is not the {silva['vsearch']} get_silva.py checked")
        args.silva_fasta = silva["fasta"]
        print("reading step 2's work files", flush=True)
        order = [r for r in runs if r["sets_t"] == "yes"] + [r for r in runs if r["sets_t"] != "yes"]
        groups = {r["run"]: gather(r, work, args.rna_sample) for r in order}
        if set(groups) != set(step2):
            raise RuntimeError("stage0a/marker_share.tsv and stage0a/runs.tsv list different runs")
    except (OSError, ValueError, KeyError, RuntimeError) as e:
        print(f"Cannot start: {e} (get_silva.py writes db/silva/manifest.json)", file=sys.stderr)
        return 1

    rdir = work / "stage0a" / "rescore"
    rdir.mkdir(parents=True, exist_ok=True)
    log = rdir / "vsearch.log"
    verdicts, rows = {}, {}
    try:
        for phase, members in (("DNA", [g for g in groups.values() if g["dna"]]),
                               ("RNA", [g for g in groups.values() if not g["dna"]])):
            queries = {}
            for g in members:
                for k in g["placed"]:
                    queries.setdefault(seq_id(g["seqs"][k]), g["seqs"][k])
            print(f"{phase} runs: {sum(len(g['placed']) for g in members)} symbiont SSU reads, "
                  f"{len(queries)} distinct sequences to place", flush=True)
            hits = place(queries, rdir / "batches", args, log)
            tax = load_taxonomy(silva["taxonomy"], {t for q in queries for t in hits.get(q, [])})
            verdicts.update({q: verdict(hits.get(q, []), tax) for q in queries})
            for g in members:
                x = rows[g["r"]["run"]] = score(g, verdicts)
                print(f"  {x['run']}: SSU kept {x['ssu_kept']} of {x['ssu_step2']}, SingleM kept "
                      f"{fmt(x['singlem_kept'], 'd')} of {x['singlem_step2']}; share {x['share_step2']:.4g} -> "
                      f"{fmt(x['share'])}", flush=True)
            if phase == "DNA":
                low = min((x for x in rows.values() if x["sets_m"] == "yes"), key=lambda x: x["share"])
                print(f"m from the DNA positives: {low['share']:.4g} (in {low['run']}); placing the RNA runs next",
                      flush=True)
    except (OSError, RuntimeError, KeyError, ValueError) as e:
        print(f"\nStopped, so nothing was written to {stage_dir}/: {e}\nRun the same command again: finished VSEARCH "
              "batches are kept.", file=sys.stderr)
        return 1

    ordered = [rows[r["run"]] for r in runs]
    s = summarize(ordered)
    low2 = min((x for x in step2.values() if x["sets_m"] == "yes"), key=lambda x: float(x["share"]))
    meta = {"made": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC"), "vsearch": vsearch,
            "silva": silva, "rna_sample": args.rna_sample, "m_step2": float(low2["share"])}
    write_outputs(stage_dir, ordered, s, meta)
    print()
    print(f"m = {s['m']:.4g} (step 2: {meta['m_step2']:.4g})" if s["m"] is not None else
          "SSUplex failed its check, so m is not set (section 4.2: phyloFlash replaces SSUplex).")
    print(f"SSUplex finds CaGg and MRE SSU: {'passed' if s['ssuplex_ok'] else 'failed'}; SingleM sees MRE: "
          f"{'yes' if s['singlem_sees_mre'] else 'no (MRE rest on SSU)'}")
    print(f"-> {stage_dir}/marker_share_entry7.tsv, {stage_dir}/STAGE0A_M_ENTRY7.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
