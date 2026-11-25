#!/usr/bin/env python3
"""Test T5, external pair: Linnemannia elongata AD073 wild type against the cured line (PROTOCOL.md section 5;
section 10, entry 15). Reported, not a kill test.

The runs in t5/runs.tsv (six JGI metatranscriptomes, three per line, and the wild type's one DNA run) are
fetched (the first 20,000,000 spots), assigned with the full index, routed, placed and called as T1 does
(run_t1.py v0.2, entries 9 and 10), on their own, with t and m of entry 10. SSU pieces placed before (T1 to
T4) are not placed again.

Expected (controls.tsv): the wild type resident for Mycoavidus, the cured line absent. The metatranscriptomes'
prep is not stated (entry 15), so they are not informative and the cured line cannot be called absent: T5
reads its expectation as "not resident" and reports the call it gets. The wild type also carries an MRE that
has no genome in the references, so every lineage called in the wild type and not in the cured line is
reported too, with no expectation (entry 15).

Writes to t5/ (commit these): levels.tsv, calls.tsv and T5.md. Work files go to WORKDIR/t5/. A rerun skips
what is done.

Usage:
  conda activate phase0-reads
  python scripts/run_t5.py --workdir ~/phase0-data --threads 30
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
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

SYMBIONT = "Mycoavidus"                                     # the pair's known symbiont (section 5, T5)
CALLED = ("resident", "unresolved")                         # "called" in the list of entry 15, choice 4
EARLIER = [("t1", "place", "batches"), ("t2", "place_full", "batches"), ("t3", "place", "batches"),
           ("t4", "place_full", "batches")] + list(RT.FULL_CACHES)   # entry 8's command, all of SILVA


def now():
    return T.now()


def pair(runs, controls):
    """(wild-type unit, cured unit) from the expected calls in controls.tsv: one unit expected resident and one
    expected absent, each with one control_id."""
    units = {}
    for r in runs:
        units.setdefault(r["unit"], set()).add(r["control_id"])
    if len(units) != 2 or any(len(c) != 1 for c in units.values()):
        raise RuntimeError("t5/runs.tsv must hold two lines (units), each with one control_id")
    missing = sorted(c for cs in units.values() for c in cs if c not in controls)
    if missing:
        raise RuntimeError(f"control_id not in controls.tsv: {', '.join(missing)}")
    exp = {u: controls[next(iter(c))] for u, c in units.items()}
    wt = [u for u, e in exp.items() if e == "resident"]
    cured = [u for u, e in exp.items() if e == "absent"]
    if len(wt) != 1 or len(cured) != 1:
        raise RuntimeError(f"controls.tsv must expect one line resident and the other absent (it gives {exp})")
    return wt[0], cured[0], exp


def run_data(r, c, adir, ctx, route=None):
    """One run's data, as run_t1.py v0.2 makes it (its main, entry 10). route: (found, SSU pieces, SingleM
    taxonomies, hits, taxonomy, organelle DNA, split reads), or None for a run without the route."""
    d = {"reads": c[("total", "reads in input")], "fungal": c[("fungal", "fungal reads for E2")],
         "mapped": {lin: c.get(("lineage", lin), 0) for lin in ctx["bacterial"]},
         "spread": AR.read_spread(adir / "spread.tsv") if r["kind"] == "DNA" else None,
         "route": None, "split": Counter(), "placements": defaultdict(Counter), "stats": Counter()}
    if route is None:
        return d
    found, pieces, sm_tax, hits, tax, organelle, split = route
    d["route"], st = Counter(), d["stats"]
    for key, (pos, a, b) in found.items():
        st["found"] += 1
        piece = pieces[key] if a else None
        ssu = C.ssu_label(hits.get(R.seq_id(piece), []), tax, len(piece)) if a else None
        sml = C.singlem_label(sm_tax[key]) if b else None
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
        if (key[0], pos) in organelle:                      # entry 10
            st["organelle_dna"] += 1
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
    return d


def write_outputs(t5_dir, runs, levels, calls, wt, cured, exp, meta):
    def level_row(x):
        v = {**x, "fungal_reads": x["fungal"], "informative": T.yn(x["informative"]),
             "cross_talk": T.yn(x["cross_talk"]), "N": T.fmt(x["N"]), "fraction": T.fmt(x["fraction"]),
             "sylph_ani": T.fmt(x["sylph_ani"], ".2f"), "spread": T.fmt(x["spread"], ".2f"),
             "spread_cap": T.yn(x["spread_cap"])}
        return "\t".join(str(v[c]) for c in T.LEVEL_COLS) + "\n"

    out = {"levels.tsv": "\t".join(T.LEVEL_COLS) + "\n" + "".join(level_row(x) for x in levels),
           "calls.tsv": "\t".join(T.CALL_COLS) + "\n" + "".join(
               "\t".join(T.yn(c[k]) if k in ("signal", "e1", "e2", "e5") else str(c[k]) for k in T.CALL_COLS) + "\n"
               for c in calls)}
    by = {(c["unit"], c["lineage"]): c for c in calls}
    a, b = by.get((wt, SYMBIONT), {}), by.get((cured, SYMBIONT), {})
    ok = a.get("call") == "resident" and b.get("call") != "resident"
    why = []
    if a.get("call") != "resident":
        why.append(f"the wild type is called {a.get('call', 'nothing')} for {SYMBIONT}, not resident")
    if b.get("call") == "resident":
        why.append(f"the cured line is called resident for {SYMBIONT}")
    md = ["# Test T5: the external pair, Linnemannia elongata AD073 wild type against the cured line", "",
          f"Made by run_t5.py v{__version__} (run_t1.py v{T.__version__}, call_lineages.py v{C.__version__}) on "
          f"{meta['made']}, with {meta['assign']}, minimap2 {meta['minimap2']} and the full index (refs.fna md5 "
          f"{meta['refs_md5']}); {meta['tools']}; VSEARCH {meta['vsearch']} against SILVA `{meta['silva']}`; "
          f"{meta['sylph']}.",
          f"Rule: PROTOCOL.md 4.1 and section 5 (T5), with section 10, entries 4 to 10 and 15. t = {meta['t']:.4g} "
          f"and m = {meta['m']:.4g} (entry 10). Reported, not a kill test.", "",
          f"## Result: {'as expected' if ok else 'not as expected'}", "",
          (f"The wild type is called resident for {SYMBIONT} and the cured line is not." if ok else
           "Not as expected: " + "; ".join(why) + "."),
          "",
          "| Line | Expected (controls.tsv) | Read as | Call | Basis | Present / low / informative runs |",
          "|---|---|---|---|---|---|"]
    for u, c, read_as in ((wt, a, "resident"), (cured, b, "not resident")):
        md.append(f"| {u} | {exp[u]} | {read_as} | {c.get('call', '')} | {c.get('basis', '')} | "
                  f"{c.get('present_runs', 0)} / {c.get('low_runs', 0)} / {c.get('informative_runs', 0)} |")
    md += ["", "The six metatranscriptomes are of unstated prep (entry 15), so they are not informative: the cured "
               "line cannot be called absent, and T5 reads its expectation as not resident.", "",
           f"## {SYMBIONT} run by run", "",
           "| Run | Line | Kind | Reads | n | N | Fraction | Level | Informative | Cross-talk | Spread | sylph ANI |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for x in levels:
        if x["lineage"] == SYMBIONT and x["source"] == "mapped":
            md.append(f"| {x['run']} | {x['unit']} | {x['kind']} | {x['reads']:,} | {x['n']:,} | {T.fmt(x['N'])} | "
                      f"{T.fmt(x['fraction'])} | {x['level']} | {T.yn(x['informative'])} | {T.yn(x['cross_talk'])} | "
                      f"{T.fmt(x['spread'], '.2f')} | {T.fmt(x['sylph_ani'], '.2f')} |")
    only = [lin for (u, lin), c in by.items() if u == wt and c["call"] in CALLED and
            by.get((cured, lin), {}).get("call") not in CALLED]
    md += ["", "## Called in the wild type and not in the cured line (entry 15: reported, no expectation)", "",
           ("Called (resident or unresolved) in the wild type, and absent, no call or contaminant in the cured line: "
            + ", ".join(f"{lin} ({by[(wt, lin)]['call']})" for lin in only) + "." if only else
            "No lineage is called in the wild type and not in the cured line.")]
    sig = sorted({lin for (u, lin), c in by.items() if c["signal"]},
                 key=lambda lin: (lin != SYMBIONT, lin))
    md += ["", "## Every lineage with signal in either line", "",
           f"| Lineage | Category | {wt} | {cured} |", "|---|---|---|---|"]
    for lin in sig:
        cw, cc = by.get((wt, lin), {}), by.get((cured, lin), {})
        md.append(f"| {lin} | {(cw or cc).get('category', '')} | {cw.get('call', 'none')} ({cw.get('basis', '')}) | "
                  f"{cc.get('call', 'none')} ({cc.get('basis', '')}) |")
    md += ["", "`levels.tsv` has every run and lineage; `calls.tsv` every call."]
    out["T5.md"] = "\n".join(md) + "\n"
    for name, text in out.items():
        (t5_dir / (name + ".part")).write_text(text, encoding="utf-8")
    for name in out:
        os.replace(t5_dir / (name + ".part"), t5_dir / name)
    return ok


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--refdir", default="references")
    ap.add_argument("--t5-dir", default="t5", help="holds runs.tsv; the outputs to commit go here too")
    ap.add_argument("--stage0a-dir", default="stage0a")
    ap.add_argument("--controls", default="controls.tsv")
    ap.add_argument("--work-name", default="t5", help="work files go to WORKDIR/<this> [t5]")
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
    work, t5_dir, stage_dir = Path(args.workdir).expanduser(), Path(args.t5_dir), Path(args.stage0a_dir)
    refdir, wd = Path(args.refdir), work / args.work_name
    args.index = work / "index" / "full"

    try:
        runs = T.read_runs(t5_dir / "runs.tsv")
        controls = {x["control_id"]: x["expected_call"] for x in M2.read_tsv(Path(args.controls))}
        wt, cured, exp = pair(runs, controls)
        t, m = T.stage_values(stage_dir)
        mm2 = subprocess.run([args.minimap2, "--version"], capture_output=True, text=True, check=True).stdout.strip()
        ctx = RT.arm_ctx(args.index, "full", refdir / "overlaps.tsv", mm2, M2.read_tsv(refdir / "lineages.tsv"))
        if SYMBIONT not in ctx["bacterial"]:
            raise RuntimeError(f"the full index has no {SYMBIONT} lineage")
        sylph_v = DB.tool_output([args.sylph, "--version"]).strip().splitlines()[0]
        if sylph_v != ctx["info"]["sylph"]:
            raise RuntimeError(f"sylph here is '{sylph_v}', but the index was built with '{ctx['info']['sylph']}'")
        silva = json.loads((work / "db" / "silva" / "manifest.json").read_text())
        vsearch = GS.vsearch_version(args.vsearch)
        if vsearch != silva["vsearch"]:
            raise RuntimeError(f"VSEARCH {vsearch} is not the {silva['vsearch']} get_silva.py checked")
        args.silva_fasta = silva["fasta"]
        db = json.loads((work / "db" / "manifest.json").read_text())
        if any(r["route"] != "none" for r in runs):
            versions = DB.check_versions(args)
            if versions != db["versions"]:
                raise RuntimeError(f"the tools here ({versions}) are not those get_marker_databases.py checked "
                                   f"({db['versions']})")
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.CalledProcessError) as e:
        print(f"Cannot start: {e}", file=sys.stderr)
        return 1
    for sub in ("reads", "assign", "route", "organelle", "split", "sylph"):
        (wd / sub).mkdir(parents=True, exist_ok=True)

    # 1. fetch, assign and route, as T1 does
    def one_fetch(r):
        try:
            info, new = T.fetch(r["run"], wd / "reads", T.SPOTS[r["platform"]], args)
            print(f"{'fetched' if new else 'already fetched'}: {r['run']} ({info['spots']:,} spots)", flush=True)
            return None
        except (OSError, RuntimeError) as e:
            return f"{r['run']}: {e}"

    failed, files, counts, adirs, rdirs = [], {}, {}, {}, {}
    print(f"fetching {len(runs)} runs, {args.jobs} at a time", flush=True)
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        failed += [e for e in pool.map(one_fetch, runs) if e]
    for r in runs:
        run = r["run"]
        try:
            fdir = wd / "reads" / run
            if not (fdir / "fetch.json").exists():
                continue
            fi = json.loads((fdir / "fetch.json").read_text())
            adirs[run] = wd / "assign" / run
            if not (adirs[run] / "counts.tsv").exists():
                print(f"assigning {run}", flush=True)
            T.assign(run, [fdir / f for f in fi["files"]], adirs[run], r["platform"], args)
            files[run], counts[run] = T.assigned(adirs[run], r["platform"], ctx)
            if r["route"] == "none":
                continue
            rdirs[run] = wd / "route" / run
            if not (rdirs[run] / "route.json").exists():
                print(f"route: {run}", flush=True)
            T.route_new(r, files[run], rdirs[run], args, db, ctx)
        except (OSError, KeyError, ValueError, RuntimeError, json.JSONDecodeError,
                subprocess.CalledProcessError) as e:
            rdirs.pop(run, None)
            failed.append(f"{run}: {e}")
            print(f"FAILED: {run}: {e}", flush=True)
    if failed:
        print(f"\n{len(failed)} run(s) not done, so nothing was written to {t5_dir}/. Run the same command again "
              "to finish them:\n" + "\n".join(failed), file=sys.stderr)
        return 1

    # 2. placement, data, levels and calls
    try:
        found, pieces, sm_tax = {}, {}, {}
        for run, rdir in rdirs.items():
            r = next(x for x in runs if x["run"] == run)
            found[run] = T.routed_found(r, rdir)
            pieces[run] = R.ssu_sequences(rdir, {k for k, v in found[run].items() if v[1]})
            sm_tax[run] = T.singlem_taxonomies(rdir, {k for k, v in found[run].items() if v[2]})
        queries = {}
        for ps in pieces.values():
            for s in ps.values():
                queries.setdefault(R.seq_id(s), s)
        hits = RT.placements(queries, wd / "place" / "batches", silva["fasta"], args, wd / "place" / "vsearch.log",
                             [work.joinpath(*p) for p in EARLIER])
        tax = {k: [x.strip() for x in v.split(";")] for k, v in
               R.load_taxonomy(silva["taxonomy"], {t_ for q in queries for t_, _ in hits.get(q, [])}).items()}
        data = {}
        for r in runs:
            run, route = r["run"], None
            if run in rdirs:
                print(f"route reads: {run}", flush=True)
                organelle = T.organelle_check(found[run], rdirs[run], wd / "organelle" / f"{run}.json", args, ctx)
                split = T.split_check(found[run], rdirs[run], wd / "split" / f"{run}.json", args, ctx)
                route = (found[run], pieces[run], sm_tax[run], hits, tax, organelle, split)
            data[run] = run_data(r, counts[run], adirs[run], ctx, route)
            print(f"sylph: {run}", flush=True)
            data[run]["sylph"] = T.sylph(files[run], wd / "sylph" / f"{run}.tsv", args, ctx)
        levels, calls, _ = T.build(runs, data, t, m, ctx)
    except (OSError, KeyError, ValueError, RuntimeError, subprocess.CalledProcessError) as e:
        print(f"Nothing was written to {t5_dir}/: {e}", file=sys.stderr)
        return 1

    a1 = json.loads((adirs[runs[0]["run"]] / "assign.json").read_text())
    meta = {"made": now(), "assign": a1["made_by"], "minimap2": mm2, "refs_md5": ctx["refs_md5"],
            "tools": ", ".join(f"{k} {v}" for k, v in sorted(db["versions"].items())), "vsearch": vsearch,
            "silva": silva["file"], "sylph": sylph_v, "t": t, "m": m}
    t5_dir.mkdir(parents=True, exist_ok=True)
    ok = write_outputs(t5_dir, runs, levels, calls, wt, cured, exp, meta)
    by = {(c["unit"], c["lineage"]): c for c in calls}
    print(f"\nT5 (reported): {'as expected' if ok else 'NOT as expected'}")
    for u in (wt, cured):
        c = by.get((u, SYMBIONT), {})
        print(f"  {u}: {SYMBIONT} {c.get('call', '')} ({c.get('basis', '')})")
    print(f"-> {t5_dir}/levels.tsv, {t5_dir}/calls.tsv, {t5_dir}/T5.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
