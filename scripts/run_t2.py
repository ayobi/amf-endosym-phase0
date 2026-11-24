#!/usr/bin/env python3
"""Test T2, unknown lineages (PROTOCOL.md sections 5 and 6; section 10, entries 9 to 11).

The 13 MRE positives in t2/runs.tsv (D. epigaea S1, R. irregularis A2, Rhizophagus sp. LPA8) are called
twice, each time as T1 calls its runs (run_t1.py v0.2, with entry 10's rules), on stage 0a's reads (the
first 20,000,000 spots) and with the outputs SSUplex and SingleM made in stage 0a's step 2:
  held out (scored)    the T2 index, which leaves out the 9 MRE genomes, and SILVA without the sequences
                       in t2/silva_removed.tsv (make_t2_silva.py). The runs are assigned again with the
                       T2 index.
  kept in (reported)   the full index and all of SILVA, with stage 0a's assignment (entry 10 index).
In each arm the found reads are aligned again with that arm's index, the routed ones are labelled,
organelle DNA and split reads leave the route, and sylph runs against that arm's database.

SSU placement. Every routed SSU piece of either arm is placed against all of SILVA with entry 8's command
(placements made before with that command, in stage 0a and in T1, are used again). For the held-out arm,
a piece's top hits against SILVA without the removed sequences follow from those: VSEARCH compares the
piece with every sequence (--maxaccepts 0 --maxrejects 0) and reports only the best ones
(--top_hits_only), so they are its full-SILVA top hits that were not removed, if any are left, and none
if it had none. Only pieces whose top hits were all removed are searched again, against the T2 database.
A sample of up to 300 of the other pieces is searched again too; any difference stops the run.

K2 (section 6), per held-out isolate (entry 11). The held-out MRE's reads are the found reads that the
full index counts for CaMg/MRE. In the held-out arm, the lineage where the MRE should be is the route
lineage that holds the most of them in the isolate's runs. K2 fails for the isolate if no route lineage
holds any of them, or if that lineage's call is "absent" or "no call".

Writes to t2/ (commit these): levels.tsv and calls.tsv (both arms, with an arm column) and T2.md (K2,
the held-out MRE per isolate, the calls). Work files go to WORKDIR/t2/. Reruns skip what is done.

Needs t and m from entry 10 (stage0a/fractions_entry10.tsv, stage0a/marker_share_entry10.tsv), stage 0a's
reads and step 2 files, the full and T2 indexes (make_index.py v0.2), db/silva/ and db/silva_t2/.

Usage:
  conda activate phase0-reads
  python scripts/run_t2.py --workdir ~/phase0-data --threads 30
"""

import argparse
import copy
import csv
import gzip
import hashlib
import json
import os
import random
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
import assign_reads as AR          # noqa: E402  index, mask, spread
import call_lineages as C          # noqa: E402  labels, levels, calls
import get_marker_databases as DB  # noqa: E402  tool versions
import get_silva as GS             # noqa: E402  VSEARCH version
import rescore_stage0a_m as R      # noqa: E402  SSU pieces, sequence ids, taxonomy
import run_stage0a_m as M2         # noqa: E402  tables
import run_t1 as T                 # noqa: E402  the T1 pipeline, run here per arm

__version__ = "0.1"

ARMS = ("held_out", "kept_in")
ARM_NAME = {"held_out": "held out (scored)", "kept_in": "kept in (reported)"}
INDEX_SET = {"held_out": "t2", "kept_in": "full"}
MRE = "CaMg/MRE"
CHECK_SAMPLE = 300                                       # entry 11
FULL_CACHES = (("stage0a", "rescore_entry8", "batches"), ("stage0a", "rescore_entry10", "batches"),
               ("t1", "place", "batches"))               # earlier placements with entry 8's command, all of SILVA
LEVEL_COLS = ["arm"] + T.LEVEL_COLS
CALL_COLS = ["arm"] + T.CALL_COLS


def now():
    return T.now()


def md5_text(text):
    return hashlib.md5(text.encode("utf-8")).hexdigest()


# ------------------------------------------------------------------ the arms ---

def arm_ctx(index_dir, want_set, overlaps, mm2, lin_rows):
    """Everything an arm needs about its index. Bacterial lineages come from the index itself, so the T2
    index has no CaMg/MRE lineage at all."""
    info, lineage_of, _, hosts = AR.load_index(index_dir)
    if info.get("set") != want_set:
        raise RuntimeError(f"{index_dir} is the {info.get('set')} index, not {want_set}")
    if not (index_dir / "organelles.fna").exists():
        raise RuntimeError(f"{index_dir} has no organelles.fna: it was built before entry 10 (make_index.py v0.2)")
    if mm2 != info["minimap2"]:
        raise RuntimeError(f"minimap2 here is {mm2}, but {index_dir} was built with {info['minimap2']}")
    seq_rows = AR.load_sequences(index_dir)
    bact = {}
    for x in seq_rows:
        if x["kind"] in AR.BACTERIAL_KINDS:
            bact.setdefault(x["lineage"], x["kind"])
    first = [x for x in ("CaGg", MRE, "Mycoavidus") if x in bact]
    order = first + sorted(x for x in bact if x not in first)
    return {"mmi": index_dir / AR.INDEX_FILE["short"], "lineage_of": lineage_of, "hosts": hosts,
            "mask": AR.load_mask(overlaps), "minimap2": mm2, "refs_md5": info["refs_md5"],
            "overlaps_md5": AR.md5_of(overlaps), "bacterial": set(bact), "bacterial_order": order,
            "syldb": index_dir / "bacteria.syldb", "acc_lineage": {x["accession"]: x["lineage"] for x in lin_rows},
            "organelles": index_dir / "organelles.fna", "index": index_dir, "info": info}


def mre_reads(view):
    """(mate, read) of the found reads that this view's index counts for CaMg/MRE."""
    out = set()
    with gzip.open(view / "found.tsv.gz", "rt", encoding="utf-8") as fh:
        for x in csv.DictReader(fh, delimiter="\t"):
            if x["outcome"] == "lineage" and x["detail"] == MRE:
                out.add((x["mate"], x["read"]))
    return out


# ------------------------------------------------------------------- placement ---

def read_batches(dirs):
    """(ids placed, id -> [(target, identity)]) from folders of finished VSEARCH batches."""
    done, hits = set(), defaultdict(list)
    for d in dirs:
        if not d.exists():
            continue
        for fa in sorted(d.glob("batch_*.fa")):
            tsv = fa.with_suffix(".tsv")
            if not tsv.exists():
                continue
            ids = [line[1:].strip() for line in open(fa) if line.startswith(">")]
            new = [q for q in ids if q not in done]
            done.update(new)
            new = set(new)
            for line in open(tsv):
                q, t, i = line.rstrip("\n").split("\t")
                if q in new:
                    hits[q].append((t, float(i)))
    return done, hits


def placements(queries, bdir, fasta, args, log, also=()):
    """Entry 8's VSEARCH command against `fasta` for every query that neither the folders in `also` nor
    bdir hold yet; returns id -> [(target, identity)] for every query."""
    done, before = read_batches(also)
    a = copy.copy(args)
    a.silva_fasta = str(fasta)
    new = T.place({q: s for q, s in queries.items() if q not in done}, bdir, a, log)
    return {q: (before.get(q, []) if q in done else new.get(q, [])) for q in queries}


def derive(full, removed):
    """The T2 top hits that follow from the full-SILVA ones, and the queries that need a search."""
    t2, todo = {}, []
    for q, hits in full.items():
        kept = [h for h in hits if h[0] not in removed]
        if kept or not hits:
            t2[q] = kept
        else:
            todo.append(q)
    return t2, todo


def same(a, b):
    return sorted((t, round(i, 1)) for t, i in a) == sorted((t, round(i, 1)) for t, i in b)


# --------------------------------------------------------------- one arm's data ---

def arm_data(arm, runs, files, counts, adirs, views, ctx, hits, tax, pieces, found, sm_tax, wd, args, mre=None):
    """Per run: counts, spread, route counts, placements and route statistics, as T1 makes them. In the
    held-out arm, `mre` gives the held-out MRE's reads per run, and their route labels are counted."""
    data = {}
    for r in runs:
        run, c = r["run"], counts[r["run"]]
        d = data[run] = {"reads": c[("total", "reads in input")], "fungal": c[("fungal", "fungal reads for E2")],
                         "mapped": {lin: c.get(("lineage", lin), 0) for lin in ctx["bacterial"]},
                         "spread": AR.read_spread(adirs[run] / "spread.tsv") if r["kind"] == "DNA" else None,
                         "route": Counter(), "split": Counter(), "placements": defaultdict(Counter),
                         "stats": Counter(), "mre": Counter(), "mre_stats": Counter()}
        print(f"{ARM_NAME[arm]}: route reads of {run}", flush=True)
        organelle = T.organelle_check(found[run], views[run], wd / f"organelle_{arm}" / f"{run}.json", args, ctx)
        split = T.split_check(found[run], views[run], wd / f"split_{arm}" / f"{run}.json", args, ctx)
        st, mine = d["stats"], (mre or {}).get(run, set())
        for key, (pos, a, b) in found[run].items():
            st["found"] += 1
            is_mre = key in mine
            if is_mre:
                d["mre_stats"]["found"] += 1
            piece = pieces[run][key] if a else None
            ssu = C.ssu_label(hits.get(R.seq_id(piece), []), tax, len(piece)) if a else None
            sml = C.singlem_label(sm_tax[run][key]) if b else None
            if ssu:
                if ssu[0] == C.UNPLACED:
                    st["ssu_several" if ssu[2] == C.SEVERAL_GENERA else "ssu_short"] += 1
                else:
                    st["ssu_kept" if ssu[0] == "kept" else "ssu_" + ssu[0].split()[0]] += 1
            if sml and sml[0] == "kept":
                st["sm_kept"] += 1
            label = C.read_label(ssu, sml)
            if label is None:
                continue
            if (key[0], pos) in organelle:
                st["organelle_dna"] += 1
                d["mre_stats"]["organelle_dna"] += is_mre
                continue
            if (key[0], pos) in split:
                lin = split[(key[0], pos)]
                st["split_bact" if lin else "split_host"] += 1
                d["mre_stats"]["split"] += is_mre
                if lin:
                    d["split"][lin] += 1
                continue
            d["route"][label] += 1
            if is_mre:
                d["mre"][label] += 1
                d["mre_stats"]["labelled"] += 1
            if label[0] == "novel":
                lname = C.NOVEL_SSU if label[1] == "SSU" else C.NOVEL_SINGLEM
                d["placements"][lname][(ssu if ssu and ssu[0] == "kept" else sml)[2]] += 1
        print(f"{ARM_NAME[arm]}: sylph on {run}", flush=True)
        d["sylph"] = T.sylph(files[run], wd / f"sylph_{arm}" / f"{run}.tsv", args, ctx)
    return data


# ---------------------------------------------------------------------- K2 ---

def k2(runs, data, calls, joiner, kept_counts):
    """Per held-out isolate (unit): the held-out MRE's reads, the lineage that holds the most of them,
    its call, and whether K2 passes there. Returns (rows, reasons)."""
    by = {(c["unit"], c["lineage"]): c for c in calls}
    rows, reasons = [], []
    for unit in dict.fromkeys(r["unit"] for r in runs):
        mine = [r for r in runs if r["unit"] == unit]
        held, st = Counter(), Counter()
        for r in mine:
            for label, n in data[r["run"]]["mre"].items():
                held[C.lineage_name(label, joiner)] += n
            st.update(data[r["run"]]["mre_stats"])
        step1 = sum(kept_counts[r["run"]].get(("lineage", MRE), 0) for r in mine)
        main = max(held, key=lambda x: (held[x], x)) if held else None
        call = by[(unit, main)]["call"] if main else "none"
        ok = main is not None and call not in ("absent", "no call")
        if main is None:
            reasons.append(f"{unit}: no route lineage holds the held-out MRE's reads")
        elif not ok:
            reasons.append(f"{unit}: {main}, where the held-out MRE is, is called {call}")
        rows.append({"unit": unit, "step1": step1, "found": st["found"], "labelled": st["labelled"],
                     "organelle_dna": st["organelle_dna"], "split": st["split"], "lineage": main or "",
                     "lineage_mre": held[main] if main else 0, "call": call,
                     "basis": by[(unit, main)]["basis"] if main else "", "passed": ok,
                     "others": [(x, n) for x, n in held.most_common() if x != main]})
    return rows, reasons


# -------------------------------------------------------------------- outputs ---

def yn(v):
    return "yes" if v else "no"


def fmt(v, spec=".4g"):
    return "" if v in ("", None) else format(v, spec)


def calls_table(calls, unit_order):
    rows = ["| Isolate | Lineage | Source | Call | Basis | Present / low / informative runs | Placements |",
            "|---|---|---|---|---|---|---|"]
    for unit in unit_order:
        for c in calls:
            if c["unit"] == unit and c["signal"]:
                rows.append(f"| {unit} | {c['lineage']} | {c['source']} | {c['call']} | {c['basis']} | "
                            f"{c['present_runs']} / {c['low_runs']} / {c['informative_runs']} | {c['placements']} |")
    return rows if len(rows) > 2 else ["No lineage has signal in any run."]


def write_outputs(t2_dir, runs, out, k2_rows, reasons, meta):
    def level_row(arm, x):
        v = {**x, "arm": arm, "fungal_reads": x["fungal"], "informative": yn(x["informative"]),
             "cross_talk": yn(x["cross_talk"]), "N": fmt(x["N"]), "fraction": fmt(x["fraction"]),
             "sylph_ani": fmt(x["sylph_ani"], ".2f"), "spread": fmt(x["spread"], ".2f"),
             "spread_cap": yn(x["spread_cap"])}
        return "\t".join(str(v[c]) for c in LEVEL_COLS) + "\n"

    def call_row(arm, c):
        v = {**c, "arm": arm}
        return "\t".join(yn(v[k]) if k in ("signal", "e1", "e2", "e5") else str(v[k]) for k in CALL_COLS) + "\n"

    tmp = {"levels.tsv": "\t".join(LEVEL_COLS) + "\n" + "".join(level_row(arm, x) for arm in ARMS
                                                               for x in out[arm]["levels"]),
           "calls.tsv": "\t".join(CALL_COLS) + "\n" + "".join(call_row(arm, c) for arm in ARMS
                                                             for c in out[arm]["calls"])}
    units = list(dict.fromkeys(r["unit"] for r in runs))
    md = ["# Test T2: unknown lineages", "",
          f"Made by run_t2.py v{__version__} (run_t1.py v{T.__version__}, call_lineages.py v{C.__version__}) on "
          f"{meta['made']}, with {meta['assign']} and minimap2 {meta['minimap2']}: the T2 index (refs.fna md5 "
          f"{meta['refs_t2']}) for the held-out arm and the full index (refs.fna md5 {meta['refs_full']}) for the "
          f"kept-in arm; VSEARCH {meta['vsearch']} against SILVA `{meta['silva']}`, without the "
          f"{meta['removed']:,} sequences in `t2/silva_removed.tsv` for the held-out arm; {meta['sylph']}.",
          f"Rule: PROTOCOL.md 4.1, section 5 (T2) and section 6 (K2), with section 10, entries 4 to 11. "
          f"t = {meta['t']:.4g} and m = {meta['m']:.4g} (entry 10). Reads: stage 0a's, the first 20,000,000 spots "
          "of each run.", "",
          f"## K2: {'passed' if not reasons else 'FAILED'}", ""]
    md += [f"- {x}" for x in reasons] or ["No held-out isolate is called absent, and each has a lineage where the "
                                          "MRE should be."]
    md += ["", "## The held-out MRE, per isolate", "",
           "| Isolate | MRE reads (full index) | Found by the route | Stay in the route | Lineage that holds most | "
           "Its MRE reads | Call (held out) | K2 |", "|---|---|---|---|---|---|---|---|"]
    for x in k2_rows:
        md.append(f"| {x['unit']} | {x['step1']} | {x['found']} | {x['labelled']} | {x['lineage'] or 'none'} | "
                  f"{x['lineage_mre']} | {x['call']} | {'passes' if x['passed'] else 'FAILS'} |")
    others = [f"{x['unit']}: " + ", ".join(f"{lin} ({n})" for lin, n in x["others"]) for x in k2_rows if x["others"]]
    md += ["", "MRE reads (full index): the reads the full index counts for CaMg/MRE in the isolate's runs (stage 0a, "
               "entry 10). Found by the route: those of them SSUplex or SingleM found, and the T2 index leaves "
               "routed. Stay in the route: those that get a label and are neither organelle DNA nor split reads."]
    if others:
        md += ["", "Other lineages that hold MRE reads: " + "; ".join(others) + "."]
    for arm in ARMS:
        md += ["", f"## Calls, {ARM_NAME[arm]}", ""] + calls_table(out[arm]["calls"], units)
    md += ["", "## SSU placement", "",
           f"SSU pieces placed against all of SILVA: {meta['full_pieces']:,} ({meta['full_before']:,} of them placed "
           f"before, in stage 0a or T1). For the held-out arm, {meta['t2_pieces']:,} pieces: {meta['derived']:,} "
           f"follow from their full-SILVA placement and {meta['searched']:,} were searched again against the T2 "
           f"database, as their top hits were all removed. Of the derived placements, {meta['checked']:,} were "
           "also searched directly: all matched."]
    for arm in ARMS:
        md += ["", f"## The route per run, {ARM_NAME[arm]}", "",
               "| Run | Routed reads found | SSU kept | SSU dropped: outside Bacteria / organelle | "
               "SSU unplaced: several genera / under 60 bp | SingleM kept | Organelle DNA | "
               "Split reads: bacterial / host |", "|---|---|---|---|---|---|---|---|"]
        for r in runs:
            s = out[arm]["data"][r["run"]]["stats"]
            md.append(f"| {r['run']} | {s['found']} | {s['ssu_kept']} | {s['ssu_outside']} / {s['ssu_organelle']} | "
                      f"{s['ssu_several']} / {s['ssu_short']} | {s['sm_kept']} | {s['organelle_dna']} | "
                      f"{s['split_bact']} / {s['split_host']} |")
        capped = [x for x in out[arm]["levels"] if x["spread_cap"]]
        if capped:
            md += ["", "Lowered by the spread check: " + "; ".join(
                f"{x['run']} {x['lineage']} (n = {x['n']}, spread {x['spread']:.2f})" for x in capped) + "."]
    md += ["", "Reads are counted per read (mates separately). `levels.tsv` gives every lineage in every run of both "
               "arms; `calls.tsv` gives every call with its evidence; novel lineages list their most common "
               "placements. The kept-in arm is easy by construction and does not count (section 5). sylph is a "
               "cross-check and does not change calls (4.2)."]
    tmp["T2.md"] = "\n".join(md) + "\n"
    for name, text in tmp.items():
        (t2_dir / (name + ".part")).write_text(text, encoding="utf-8")
    for name in tmp:
        os.replace(t2_dir / (name + ".part"), t2_dir / name)


# ----------------------------------------------------------------------- main ---

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--refdir", default="references")
    ap.add_argument("--t2-dir", default="t2", help="holds runs.tsv and silva_removed.tsv; the outputs go here too")
    ap.add_argument("--stage0a-dir", default="stage0a")
    ap.add_argument("--controls", default="controls.tsv")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--vsearch", default="vsearch")
    ap.add_argument("--sylph", default="sylph")
    args = ap.parse_args(argv)
    work, t2_dir, stage_dir = Path(args.workdir).expanduser(), Path(args.t2_dir), Path(args.stage0a_dir)
    wd = work / "t2"
    try:
        runs = T.read_runs(t2_dir / "runs.tsv")
        if any(r["reuse"] != "stage0a" or r["platform"] != "short" or r["kind"] != "DNA" or r["route"] != "ssu+singlem"
               for r in runs):
            raise RuntimeError("every T2 run must be a short-read DNA run of stage 0a (reuse stage0a, route "
                               "ssu+singlem)")
        controls = {x["control_id"] for x in M2.read_tsv(Path(args.controls))}
        missing = sorted({r["control_id"] for r in runs} - controls)
        if missing:
            raise RuntimeError(f"control_id not in {args.controls}: {', '.join(missing)}")
        t, m = T.stage_values(stage_dir)
        overlaps = Path(args.refdir) / "overlaps.tsv"
        lin_rows = M2.read_tsv(Path(args.refdir) / "lineages.tsv")
        mm2 = subprocess.run([args.minimap2, "--version"], capture_output=True, text=True, check=True).stdout.strip()
        ctx = {arm: arm_ctx(work / "index" / INDEX_SET[arm], INDEX_SET[arm], overlaps, mm2, lin_rows) for arm in ARMS}
        if MRE in ctx["held_out"]["bacterial"] or MRE not in ctx["kept_in"]["bacterial"]:
            raise RuntimeError("the T2 index must leave out every CaMg/MRE genome, and the full index hold them")
        sylph_v = DB.tool_output([args.sylph, "--version"]).strip().splitlines()[0]
        if any(sylph_v != ctx[arm]["info"]["sylph"] for arm in ARMS):
            raise RuntimeError(f"sylph here is '{sylph_v}', but the indexes were built with another")
        silva = json.loads((work / "db" / "silva" / "manifest.json").read_text())
        t2_manifest = work / "db" / "silva_t2" / "manifest.json"
        if not t2_manifest.exists():
            raise RuntimeError(f"{t2_manifest} is missing: run make_t2_silva.py first (and commit its two files)")
        silva_t2 = json.loads(t2_manifest.read_text())
        vsearch = GS.vsearch_version(args.vsearch)
        if vsearch != silva["vsearch"] or vsearch != silva_t2["vsearch"]:
            raise RuntimeError(f"VSEARCH {vsearch} is not the version get_silva.py and make_t2_silva.py used")
        removed_text = (t2_dir / "silva_removed.tsv").read_text(encoding="utf-8")
        if md5_text(removed_text) != silva_t2["removed_list_md5"] or silva_t2["source_md5"] != silva["md5"]:
            raise RuntimeError(f"{work / 'db' / 'silva_t2'} was not made from {t2_dir}/silva_removed.tsv and this "
                               "SILVA; run make_t2_silva.py again")
        if AR.md5_of(silva_t2["fasta"]) != silva_t2["fasta_md5"]:
            raise RuntimeError(f"{silva_t2['fasta']} is not the FASTA make_t2_silva.py wrote; run it again")
        removed = {line.split("\t", 1)[0] for line in removed_text.splitlines()[1:] if line}
        db_path = work / "db" / "manifest.json"
        db = json.loads(db_path.read_text()) if db_path.exists() else {}
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.CalledProcessError) as e:
        print(f"Cannot start: {e}", file=sys.stderr)
        return 1
    args_t2 = copy.copy(args)
    args_t2.index, args.index = ctx["held_out"]["index"], ctx["kept_in"]["index"]
    for sub in ("assign", *(f"{x}_{arm}" for x in ("route", "organelle", "split", "sylph") for arm in ARMS)):
        (wd / sub).mkdir(parents=True, exist_ok=True)

    files, counts, adirs, views, failed = {a: {} for a in ARMS}, {a: {} for a in ARMS}, {a: {} for a in ARMS}, \
        {a: {} for a in ARMS}, []
    for r in runs:
        run = r["run"]
        try:
            src = work / "stage0a" / "route" / run
            if not (src / "route.json").exists():
                raise RuntimeError("stage 0a step 2 did not finish it")
            adirs["kept_in"][run] = work / T.STAGE0A_WORK / "assign" / run
            if not (adirs["kept_in"][run] / "assign.json").exists():
                raise RuntimeError(f"stage 0a step 1 with the entry 10 index ({T.STAGE0A_WORK}) did not finish it")
            files["kept_in"][run], counts["kept_in"][run] = T.assigned(adirs["kept_in"][run], "short", ctx["kept_in"])
            for p, _ in files["kept_in"][run]:      # T2 runs no tool: SSUplex and SingleM outputs must be there
                stem = p.name.split(".")[0]
                if not all((src / f"{stem}.{x}").exists() for x in ("ssuplex.json", "singlem_names.txt")):
                    raise RuntimeError(f"stage 0a step 2 has no SSUplex or SingleM output for {stem}")
            adirs["held_out"][run] = wd / "assign" / run
            if not (adirs["held_out"][run] / "counts.tsv").exists():
                print(f"assigning {run} with the T2 index", flush=True)
            T.assign(run, [p for p, _ in files["kept_in"][run]], adirs["held_out"][run], "short", args_t2)
            files["held_out"][run], counts["held_out"][run] = T.assigned(adirs["held_out"][run], "short",
                                                                         ctx["held_out"])
            for arm in ARMS:
                views[arm][run] = wd / f"route_{arm}" / run
                T.route_view(r, files[arm][run], src, views[arm][run], args_t2 if arm == "held_out" else args, db,
                             ctx[arm])
        except (OSError, KeyError, ValueError, RuntimeError, json.JSONDecodeError) as e:
            failed.append(f"{run}: {e}")
            print(f"FAILED: {run}: {e}", flush=True)
    if failed:
        print(f"\n{len(failed)} run(s) not done, so nothing was written to {t2_dir}/. Run the same command again to "
              "finish them:\n" + "\n".join(failed), file=sys.stderr)
        return 1

    try:
        found, pieces, sm_tax = {}, {}, {}
        for arm in ARMS:
            found[arm], pieces[arm], sm_tax[arm] = {}, {}, {}
            for r in runs:
                run = r["run"]
                found[arm][run] = T.routed_found(r, views[arm][run])
                pieces[arm][run] = R.ssu_sequences(views[arm][run], {k for k, v in found[arm][run].items() if v[1]})
                sm_tax[arm][run] = T.singlem_taxonomies(views[arm][run], {k for k, v in found[arm][run].items()
                                                                          if v[2]})
        queries = {arm: {} for arm in ARMS}
        for arm in ARMS:
            for run in pieces[arm]:
                for s in pieces[arm][run].values():
                    queries[arm].setdefault(R.seq_id(s), s)
        both = {**queries["kept_in"], **queries["held_out"]}
        earlier = [work.joinpath(*p) for p in FULL_CACHES]
        full_before = len(read_batches(earlier)[0] & both.keys())
        print(f"SSU pieces against all of SILVA: {len(both):,}, {full_before:,} of them placed before in stage 0a or "
              "T1", flush=True)
        full = placements(both, wd / "place_full" / "batches", silva["fasta"], args, wd / "place_full" / "vsearch.log",
                          earlier)
        t2_hits, todo = derive({q: full[q] for q in queries["held_out"]}, removed)
        print(f"held-out arm: {len(t2_hits):,} placements follow from all of SILVA; {len(todo):,} to search against "
              "the T2 database", flush=True)
        searched = placements({q: both[q] for q in todo}, wd / "place_t2" / "batches", silva_t2["fasta"], args,
                              wd / "place_t2" / "vsearch.log")
        sample = sorted(t2_hits)
        sample = sorted(random.Random(1).sample(sample, CHECK_SAMPLE)) if len(sample) > CHECK_SAMPLE else sample
        print(f"checking {len(sample):,} derived placements with a direct search", flush=True)
        direct = placements({q: both[q] for q in sample}, wd / "place_check" / "batches", silva_t2["fasta"], args,
                            wd / "place_check" / "vsearch.log")
        differ = [q for q in sample if not same(direct[q], t2_hits[q])]
        if differ:
            raise RuntimeError(f"{len(differ)} of {len(sample)} derived placements differ from a direct search against "
                               f"the T2 database (for example {differ[0]}); the derivation cannot be used")
        t2_hits.update(searched)
        hits = {"kept_in": {q: full[q] for q in queries["kept_in"]}, "held_out": t2_hits}
        tax = {}
        for arm, manifest in (("kept_in", silva), ("held_out", silva_t2)):
            targets = {t for q in queries[arm] for t, _ in hits[arm].get(q, [])}
            if arm == "held_out" and targets & removed:
                raise RuntimeError("a held-out placement hits a removed sequence")
            tax[arm] = {k: [x.strip() for x in v.split(";")] for k, v in
                        R.load_taxonomy(manifest["taxonomy"], targets).items()}

        mre = {r["run"]: mre_reads(views["kept_in"][r["run"]]) for r in runs}
        out = {}
        for arm in ARMS:
            data = arm_data(arm, runs, files[arm], counts[arm], adirs[arm], views[arm], ctx[arm], hits[arm],
                            tax[arm], pieces[arm], found[arm], sm_tax[arm], wd, args_t2 if arm == "held_out" else args,
                            mre if arm == "held_out" else None)
            levels, calls, joiner = T.build(runs, data, t, m, ctx[arm])
            out[arm] = {"data": data, "levels": levels, "calls": calls, "joiner": joiner}
        k2_rows, reasons = k2(runs, out["held_out"]["data"], out["held_out"]["calls"], out["held_out"]["joiner"],
                              counts["kept_in"])
    except (OSError, KeyError, ValueError, RuntimeError) as e:
        print(f"Nothing was written to {t2_dir}/: {e}", file=sys.stderr)
        return 1
    a1 = json.loads((adirs["held_out"][runs[0]["run"]] / "assign.json").read_text())
    meta = {"made": now(), "assign": a1["made_by"], "minimap2": mm2, "refs_t2": ctx["held_out"]["refs_md5"],
            "refs_full": ctx["kept_in"]["refs_md5"], "vsearch": vsearch, "silva": silva["file"],
            "removed": len(removed), "sylph": sylph_v, "t": t, "m": m, "full_pieces": len(both),
            "full_before": full_before, "t2_pieces": len(queries["held_out"]), "derived": len(t2_hits) - len(todo),
            "searched": len(todo), "checked": len(sample)}
    write_outputs(t2_dir, runs, out, k2_rows, reasons, meta)
    print(f"\nK2: {'passed' if not reasons else 'FAILED: ' + '; '.join(reasons)}")
    for x in k2_rows:
        print(f"  {x['unit']}: the held-out MRE's reads are mostly in {x['lineage'] or 'no lineage'}, called "
              f"{x['call']}")
    print(f"-> {t2_dir}/levels.tsv, {t2_dir}/calls.tsv, {t2_dir}/T2.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
