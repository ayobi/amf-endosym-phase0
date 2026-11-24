#!/usr/bin/env python3
"""Test T3, contaminant challenge (PROTOCOL.md sections 5 and 6; section 10, entries 9, 10 and 13).

The runs are in t3/runs.tsv, in three sets: the Ottawa set (the Illumina runs of PRJNA754762, PRJNA922099
and PRJNA287285, and the B- Hi-C run SRR36103000) and two sets of bacteria from outside the fungus
(PRJNA943439, PRJNA1203428). Every run is processed as T1 processes its runs (run_t1.py v0.2): the first
20,000,000 spots, step 1 with the full index, the route with entry 10's rules, the organelle and split
checks, SSU placement and sylph. A run T1 already fetched is not fetched again, and the outputs SSUplex and
SingleM made for it in T1 are used again. Calls are made per set.

Spikes (entry 13). InSilicoSeq makes each spike's reads from the genome as downloaded (WORKDIR/genomes/,
checked against references/genome_files.tsv):
  iss generate --draft <genome> --model <model> --mode kde --n_reads <n> --seed <seed> --cpus 1
with n = 2 x round(level x the background run's reads / 2) and the model given for the run in t3/runs.tsv
(NovaSeq for the NextSeq 500 runs, as entry 13 sets). Each spike's reads are assigned on their own
(assign_reads.py v0.3, full index). The reads that this routes go, pooled with those of the other spikes
and named by spike, through the route, the organelle and split checks and SSU placement. A spiked run is
its background run with the spike's results added: every rule in 4.1 judges each read on its own, so this
gives what the run would give with the spike's reads in its files. The spread check is worked out on both
read sets together (spread_bits.json.gz). sylph, which only reports, is not run on spikes.

  single-run spikes   each of the 7 kit-contaminant genomes into each run marked background = yes, at
                      0.001%, 0.01%, 0.1% and 1% of its reads, seeds 1, 2 and 3; the Ottawa set is then
                      called again, with only that run spiked
  lab kit spike       Ralstonia pickettii (GCF_014142625.1) at 0.1% into every run marked lab_spike = yes,
                      seeds 1, 2 and 3; the Ottawa set is called again with all of them spiked
  outside the fungus  each set called as it is

K3 (section 6, as applied in entry 13). A single-run spike fails if the spiked genus (its mapped lineage, or
a route lineage in that genus) is called resident in the spiked run's isolate, or if any other lineage
there is resident with the spike but not without it. The lab kit spike fails if any isolate is called
resident for Ralstonia, or if any lineage becomes resident with the spike. An outside set fails if any of
its lineages is called resident.

Writes to t3/ (commit these): levels.tsv and calls.tsv (every set, no spikes), spikes.tsv (one row per
single-run spike), lab_spike.tsv (one row per isolate and seed) and T3.md (K3 and the tables). Work files
go to WORKDIR/t3/. Reruns skip what is done.

Needs t and m from entry 10, the full index (make_index.py v0.2), db/silva/, db/manifest.json, InSilicoSeq
(the version in env/phase0-reads.lock.yml) and the downloaded genomes.

Usage:
  conda activate phase0-reads
  python scripts/run_t3.py --workdir ~/phase0-data --threads 30 --jobs 8
"""

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
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
import run_t2 as RT  # noqa: E402  placements made before are used again

__version__ = "0.1"

LEVELS = [("L1", 0.00001), ("L2", 0.0001), ("L3", 0.001), ("L4", 0.01)]   # 0.001%, 0.01%, 0.1%, 1% (section 5)
SEEDS = (1, 2, 3)                                                          # 4.2
LAB_LEVEL, LAB_GENOME = 0.001, "GCF_014142625.1"                           # 0.1%; Ralstonia pickettii (entry 13)
ISS_ARGS = ["--mode", "kde", "--cpus", "1", "--compress", "--quiet"]       # 4.2; one cpu: entry 13
MODELS = ("HiSeq", "NovaSeq", "MiSeq")                                     # entry 13: NextSeq runs use NovaSeq
LAB_PROJECTS = ("PRJNA754762", "PRJNA922099", "PRJNA287285")               # entry 13 (not the PRJNA277029 amplicons)
SETS = ("ottawa", "hyphosphere", "added_bacterium")
SET_NAME = {"ottawa": "Ottawa set", "hyphosphere": "Hyphosphere bacteria (PRJNA943439)",
            "added_bacterium": "Added bacterium (PRJNA1203428)"}
EXTRA_COLS = ["set", "background", "lab_spike", "iss_model", "instrument"]
SEP = "~"                                                                  # spike reads: <spike>~<pair>~<mate>
EARLIER = [("t1", "place", "batches"), ("t2", "place_full", "batches")] + list(RT.FULL_CACHES)
SPIKE_COLS = ["spike", "background", "unit", "genome", "genus", "level_pct", "seed", "model", "spike_reads",
              "spike_reads_counted", "n_before", "n_after", "N", "level_before", "level_after", "spread_after",
              "route_n_after", "call_before", "call_after", "basis_after", "route_calls_after", "changed_calls", "k3"]
LAB_COLS = ["seed", "unit", "runs_spiked", "spike_reads", "present_runs", "call_before", "call_after", "basis_after",
            "e1", "e2", "e5", "route_calls_after", "changed_calls", "k3"]


def now():
    return T.now()


# ------------------------------------------------------------------- inputs ---

def read_runs(path):
    """t3/runs.tsv: T1's columns (run_t1.read_runs checks them) and T3's own."""
    rows, problems = T.read_runs(path), []
    for r in rows:
        if any(c not in r for c in EXTRA_COLS):
            raise ValueError(f"{path} needs the columns {', '.join(T.RUN_COLS + EXTRA_COLS)}")
        run, spiked = r["run"], "yes" in (r["background"], r["lab_spike"])
        if r["set"] not in SETS:
            problems.append(f"{run}: set must be one of {', '.join(SETS)}")
        if r["background"] not in ("yes", "no") or r["lab_spike"] not in ("yes", "no"):
            problems.append(f"{run}: background and lab_spike must be yes or no")
        if spiked != bool(r["iss_model"]) or (r["iss_model"] and r["iss_model"] not in MODELS):
            problems.append(f"{run}: a spiked run needs iss_model ({', '.join(MODELS)}); other runs none")
        if spiked and (r["set"] != "ottawa" or r["platform"] != "short" or r["layout"] != "PAIRED"):
            problems.append(f"{run}: only paired short-read runs of the Ottawa set are spiked")
        if r["background"] == "yes" and r["kind"] != "DNA":
            problems.append(f"{run}: single-run spikes go into DNA runs only (section 5: informative negatives)")
        if r["lab_spike"] == "yes" and r["bioproject"] not in LAB_PROJECTS:
            problems.append(f"{run}: the lab kit spike goes into {', '.join(LAB_PROJECTS)} only (entry 13)")
        if r["reuse"]:
            problems.append(f"{run}: leave reuse empty; T3 finds T1's files itself")
    for unit in {r["unit"] for r in rows}:
        if len({r["set"] for r in rows if r["unit"] == unit}) != 1:
            problems.append(f"{unit}: its runs are in more than one set")
    if problems:
        raise ValueError("\n".join(problems))
    return rows


def kit_genomes(refdir, work, ctx):
    """The kit-contaminant genomes (references/contaminants.tsv), as downloaded, with their md5 checked."""
    md5 = {x["accession"]: x["downloaded_md5"] for x in M2.read_tsv(refdir / "genome_files.tsv")}
    out = []
    for x in M2.read_tsv(refdir / "contaminants.tsv"):
        acc, genus = x["accession"], x["organism"].split()[0]
        if genus not in ctx["bacterial"]:
            raise RuntimeError(f"{genus} ({acc}) is not a lineage in {refdir}/lineages.tsv")
        path = work / "genomes" / f"{acc}.fna"
        if acc not in md5 or not path.exists() or AR.md5_of(path) != md5[acc]:
            raise RuntimeError(f"{path} is missing or is not the genome clean_references.py downloaded (md5 in "
                               f"{refdir}/genome_files.tsv)")
        out.append({"accession": acc, "genus": genus, "organism": x["organism"], "path": path, "md5": md5[acc]})
    if LAB_GENOME not in {g["accession"] for g in out}:
        raise RuntimeError(f"the lab kit spike genome {LAB_GENOME} is not in {refdir}/contaminants.tsv")
    return out


def spike_plan(runs, genomes, reads_of):
    """Every spike: its key, part (single or lab), background run, genome, level and seed, model and reads."""
    plan, lab = [], next(g for g in genomes if g["accession"] == LAB_GENOME)
    for r in runs:
        todo = []
        if r["background"] == "yes":
            todo += [(f"{r['run']}_{g['genus']}_{code}_s{s}", "single", g, code, level, s)
                     for g in genomes for code, level in LEVELS for s in SEEDS]
        if r["lab_spike"] == "yes":
            todo += [(f"{r['run']}_lab_s{s}", "lab", lab, "lab", LAB_LEVEL, s) for s in SEEDS]
        for key, part, g, code, level, s in todo:
            pairs = round(level * reads_of[r["run"]] / 2)
            if pairs < 1:
                raise RuntimeError(f"{key}: {level:%} of {reads_of[r['run']]:,} reads is not one read pair")
            plan.append({"key": key, "part": part, "run": r["run"], "genome": g, "code": code, "level": level,
                         "seed": s, "model": r["iss_model"], "n_reads": 2 * pairs})
    return plan


# ------------------------------------------------------------------- spikes ---

def rename_pairs(r1, r2, outs, key):
    """Copy InSilicoSeq's pairs, naming them <key>~<pair>~<mate>. Returns the pairs."""
    pairs = 0
    parts = [o.with_name(o.name + ".part") for o in outs]
    with gzip.open(parts[0], "wb", compresslevel=1) as o1, gzip.open(parts[1], "wb", compresslevel=1) as o2:
        it2 = AR.records(r2)
        for _h1, s1, q1 in AR.records(r1):
            h2, s2, q2 = next(it2, (None, None, None))
            if h2 is None:
                raise RuntimeError(f"{key}: InSilicoSeq wrote more mate-1 than mate-2 reads")
            o1.write(b"@%s%s%d%s1\n%s\n+\n%s\n" % (key.encode(), SEP.encode(), pairs, SEP.encode(), s1, q1))
            o2.write(b"@%s%s%d%s2\n%s\n+\n%s\n" % (key.encode(), SEP.encode(), pairs, SEP.encode(), s2, q2))
            pairs += 1
        if next(it2, None) is not None:
            raise RuntimeError(f"{key}: InSilicoSeq wrote more mate-2 than mate-1 reads")
    for p, o in zip(parts, outs, strict=True):
        os.replace(p, o)
    return pairs


def make_spike(p, sdir, args, iss_version):
    """InSilicoSeq for one spike, unless it is done. Returns its spike.json."""
    done = sdir / "spike.json"
    if done.exists():
        return json.loads(done.read_text())
    tmp = sdir / "iss.tmp"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    cmd = [args.iss, "generate", "--draft", str(p["genome"]["path"]), "--model", p["model"], "--n_reads",
           str(p["n_reads"]), "--seed", str(p["seed"]), *ISS_ARGS, "--output", str(tmp / "iss")]
    with open(sdir / "iss.log", "w") as log:
        log.write("$ " + " ".join(cmd) + "\n")
        log.flush()
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        raise RuntimeError(f"{p['key']}: InSilicoSeq failed (exit {rc}); see {sdir / 'iss.log'}")
    outs = [sdir / f"{p['key']}_1.fastq.gz", sdir / f"{p['key']}_2.fastq.gz"]
    pairs = rename_pairs(tmp / "iss_R1.fastq.gz", tmp / "iss_R2.fastq.gz", outs, p["key"])
    shutil.rmtree(tmp)
    if pairs * 2 != p["n_reads"]:
        raise RuntimeError(f"{p['key']}: InSilicoSeq wrote {pairs * 2} reads, not {p['n_reads']}")
    meta = {"made_by": f"run_t3.py v{__version__}", "made": now(), "insilicoseq": iss_version, "spike": p["key"],
            "part": p["part"], "background": p["run"], "genome": p["genome"]["accession"],
            "genome_md5": p["genome"]["md5"], "genus": p["genome"]["genus"], "level": p["level"], "seed": p["seed"],
            "model": p["model"], "reads": pairs * 2, "command": " ".join(cmd[:-1] + ["<tmp>/iss"]),
            "files": [o.name for o in outs]}
    M2.write_atomic(done, json.dumps(meta, indent=1) + "\n")
    return meta


def make_pool(group, keys, route, sdirs, pdir):
    """The routed reads of these spikes, one file per mate, for one run of the route tools. Returns pool.json."""
    done = pdir / "pool.json"
    if done.exists():
        meta = json.loads(done.read_text())
        if meta["spikes"] == keys and meta["route"] == route:
            return meta
    shutil.rmtree(pdir, ignore_errors=True)          # made from other spikes: this script's own work files
    pdir.mkdir(parents=True)
    outs = [pdir / f"{group}_1.fastq.gz", pdir / f"{group}_2.fastq.gz"]
    parts = [o.with_name(o.name + ".part") for o in outs]
    n = [0, 0]
    with gzip.open(parts[0], "wb", compresslevel=1) as o1, gzip.open(parts[1], "wb", compresslevel=1) as o2:
        for key in keys:
            for h, s, q in AR.records(sdirs[key] / "assign" / "routed.fastq.gz"):
                name = M2.first_word(h).decode("ascii")
                if not name.startswith(key + SEP) or name[-2:] not in (SEP + "1", SEP + "2"):
                    raise RuntimeError(f"{key}: unexpected read name {name} in its routed reads")
                k = int(name[-1]) - 1
                (o1, o2)[k].write(h + b"\n" + s + b"\n+\n" + q + b"\n")
                n[k] += 1
    for p_, o in zip(parts, outs, strict=True):
        os.replace(p_, o)
    meta = {"made_by": f"run_t3.py v{__version__}", "spikes": keys, "route": route, "reads": n,
            "files": [o.name for o in outs]}
    M2.write_atomic(done, json.dumps(meta, indent=1) + "\n")
    return meta


# -------------------------------------------------------------- route reads ---

def new_data(c, kind):
    return {"reads": c[("total", "reads in input")], "fungal": c[("fungal", "fungal reads for E2")],
            "mapped": Counter(), "spread": {} if kind == "DNA" else None, "route": Counter(), "split": Counter(),
            "placements": defaultdict(Counter), "stats": Counter(), "sylph": {}}


def tally_route(found, pieces, sm_tax, hits, tax, organelle, split, owner_of, data):
    """T1's route rules (run_t1.py v0.2, as run_t2.arm_data uses them) for each found read, counted for the
    run or spike that owns it."""
    for key, (pos, a, b) in found.items():
        d = data[owner_of(key)]
        st = d["stats"]
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
        if (key[0], pos) in organelle:
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


@lru_cache(maxsize=8)
def read_bits(path):
    """A spread_bits.json.gz, read when a spiked run needs it (bit sets can be large; a few are kept)."""
    return AR.read_bits(path)


def combine(d, s, bits_bg, bits_sp):
    """A spiked run: its background's data with the spike's added (every rule judges each read on its own);
    the spread check on both read sets together (bits_bg, bits_sp: paths of spread_bits.json.gz). sylph stays
    the background's."""
    placements = defaultdict(Counter)
    for src in (d["placements"], s["placements"]):
        for lin, c in src.items():
            placements[lin].update(c)
    return {"reads": d["reads"] + s["reads"], "fungal": d["fungal"] + s["fungal"],
            "mapped": {lin: d["mapped"][lin] + s["mapped"].get(lin, 0) for lin in d["mapped"]},
            "spread": AR.combined_spread([read_bits(bits_bg), read_bits(bits_sp)]) if d["spread"] is not None
            else None,
            "route": d["route"] + s["route"], "split": d["split"] + s["split"], "placements": placements,
            "stats": d["stats"] + s["stats"], "sylph": d["sylph"]}


# -------------------------------------------------------------------- calls ---

def by_unit(calls, unit):
    return {c["lineage"]: c for c in calls if c["unit"] == unit}


def in_genus(genus, lineage):
    """The spiked genus's lineages: its mapped lineage, and route lineages in that genus."""
    return lineage == genus or (lineage.endswith(" (route)") and genus in C.genus_parts(lineage[:-len(" (route)")]))


def compare(before, after, genus):
    """(the genus's lineages, K3 failures, calls that changed) for one unit, without and with a spike."""
    mine = sorted(x for x in after if in_genus(genus, x))
    fails = [f"{x} is called resident" for x in mine if after[x]["call"] == "resident"]
    fails += [f"{x} becomes resident ({before.get(x, {}).get('call', 'no lineage')} without the spike)"
              for x, c in sorted(after.items())
              if x not in mine and c["call"] == "resident" and before.get(x, {}).get("call") != "resident"]
    changed = [f"{x}: {before.get(x, {}).get('call', 'none')} -> {c['call']}" for x, c in sorted(after.items())
               if before.get(x, {}).get("call") != c["call"]]
    return mine, fails, changed


def pct(level):
    return f"{level * 100:g}"


def fmt(v, spec=".4g"):
    return "" if v in ("", None) else format(v, spec)


def yn(v):
    return "" if v is None else ("yes" if v else "no")


def single_rows(plan, runs_by, ott, data, spikes, bits, base, t, m, ctx):
    """Call the Ottawa set again for each single-run spike. Returns (rows, K3 failures)."""
    base_levels = {(x["run"], x["lineage"]): x for x in base[0]}
    base_calls = base[1]
    rows, fails = [], []
    for i, p in enumerate(q for q in plan if q["part"] == "single"):
        run, genus, unit = p["run"], p["genome"]["genus"], runs_by[p["run"]]["unit"]
        scn = dict(data)
        scn[run] = combine(data[run], spikes[p["key"]], bits.get(run), bits[p["key"]])
        levels, calls, _ = T.build(ott, scn, t, m, ctx)
        lv = {(x["run"], x["lineage"]): x for x in levels}
        mine, k3, changed = compare(by_unit(base_calls, unit), by_unit(calls, unit), genus)
        after = by_unit(calls, unit)
        a, b0 = lv.get((run, genus), {}), base_levels.get((run, genus), {})
        route_n = sum(lv[(run, x)]["n"] for x in mine if x != genus and (run, x) in lv)
        rows.append({"spike": p["key"], "background": run, "unit": unit, "genome": p["genome"]["accession"],
                     "genus": genus, "level_pct": pct(p["level"]), "seed": p["seed"], "model": p["model"],
                     "spike_reads": spikes[p["key"]]["reads"],
                     "spike_reads_counted": spikes[p["key"]]["mapped"].get(genus, 0),
                     "n_before": b0.get("n", 0), "n_after": a.get("n", 0), "N": fmt(a.get("N")),
                     "level_before": b0.get("level", ""), "level_after": a.get("level", ""),
                     "spread_after": fmt(a.get("spread"), ".2f"), "route_n_after": route_n,
                     "call_before": by_unit(base_calls, unit).get(genus, {}).get("call", ""),
                     "call_after": after.get(genus, {}).get("call", ""),
                     "basis_after": after.get(genus, {}).get("basis", ""),
                     "route_calls_after": "; ".join(f"{x}: {after[x]['call']}" for x in mine if x != genus),
                     "changed_calls": "; ".join(changed), "k3": "FAILS" if k3 else "passes"})
        fails += [f"{p['key']} ({unit}): {x}" for x in k3]
        if (i + 1) % 50 == 0:
            print(f"single-run spikes: {i + 1} called", flush=True)
    return rows, fails


def lab_rows(plan, runs_by, ott, data, spikes, bits, base, t, m, ctx):
    """Call the Ottawa set again with the lab kit spike in every lab_spike run, per seed."""
    genus = next(p["genome"]["genus"] for p in plan if p["part"] == "lab")
    rows, fails = [], []
    for s in SEEDS:
        mine_p = [p for p in plan if p["part"] == "lab" and p["seed"] == s]
        scn = dict(data)
        for p in mine_p:
            scn[p["run"]] = combine(data[p["run"]], spikes[p["key"]], bits.get(p["run"]), bits[p["key"]])
        _, calls, _ = T.build(ott, scn, t, m, ctx)
        for unit in dict.fromkeys(runs_by[p["run"]]["unit"] for p in mine_p):
            before, after = by_unit(base[1], unit), by_unit(calls, unit)
            mine, k3, changed = compare(before, after, genus)
            c = after.get(genus, {})
            spk = [p for p in mine_p if runs_by[p["run"]]["unit"] == unit]
            rows.append({"seed": s, "unit": unit, "runs_spiked": len(spk),
                         "spike_reads": sum(spikes[p["key"]]["reads"] for p in spk),
                         "present_runs": c.get("present_runs", ""),
                         "call_before": before.get(genus, {}).get("call", ""),
                         "call_after": c.get("call", ""), "basis_after": c.get("basis", ""),
                         "e1": yn(c.get("e1")), "e2": yn(c.get("e2")), "e5": yn(c.get("e5")),
                         "route_calls_after": "; ".join(f"{x}: {after[x]['call']}" for x in mine if x != genus),
                         "changed_calls": "; ".join(changed), "k3": "FAILS" if k3 else "passes"})
            fails += [f"lab kit spike, seed {s} ({unit}): {x}" for x in k3]
    return rows, fails


# ------------------------------------------------------------------ outputs ---

def tsv(cols, rows):
    return "\t".join(cols) + "\n" + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows)


def calls_table(calls, units, only_signal=True):
    out = ["| Unit | Lineage | Source | Call | Basis | Present / low / informative runs |", "|---|---|---|---|---|---|"]
    for unit in units:
        for c in calls:
            if c["unit"] == unit and (c["signal"] or not only_signal):
                out.append(f"| {unit} | {c['lineage']} | {c['source']} | {c['call']} | {c['basis']} | "
                           f"{c['present_runs']} / {c['low_runs']} / {c['informative_runs']} |")
    return out if len(out) > 2 else ["No lineage has signal."]


def write_outputs(t3_dir, runs, base, srows, lrows, fails, meta):
    def level_row(x):
        v = {**x, "fungal_reads": x["fungal"], "informative": yn(x["informative"]),
             "cross_talk": yn(x["cross_talk"]), "N": fmt(x["N"]), "fraction": fmt(x["fraction"]),
             "sylph_ani": fmt(x["sylph_ani"], ".2f"), "spread": fmt(x["spread"], ".2f"),
             "spread_cap": yn(x["spread_cap"])}
        return "\t".join(str(v[c]) for c in T.LEVEL_COLS) + "\n"

    def call_row(c):
        return "\t".join(yn(c[k]) if k in ("signal", "e1", "e2", "e5") else str(c[k]) for k in T.CALL_COLS) + "\n"

    levels = [x for s in SETS if s in base for x in base[s][0]]
    calls = [c for s in SETS if s in base for c in base[s][1]]
    out = {"levels.tsv": "\t".join(T.LEVEL_COLS) + "\n" + "".join(level_row(x) for x in levels),
           "calls.tsv": "\t".join(T.CALL_COLS) + "\n" + "".join(call_row(c) for c in calls),
           "spikes.tsv": tsv(SPIKE_COLS, srows), "lab_spike.tsv": tsv(LAB_COLS, lrows)}
    units = {s: list(dict.fromkeys(r["unit"] for r in runs if r["set"] == s)) for s in SETS}
    md = ["# Test T3: contaminant challenge", "",
          f"Made by run_t3.py v{__version__} (run_t1.py v{T.__version__}, call_lineages.py v{C.__version__}) on "
          f"{meta['made']}, with {meta['assign']}, minimap2 {meta['minimap2']} and the full index (refs.fna md5 "
          f"{meta['refs_md5']}); {meta['tools']}; VSEARCH {meta['vsearch']} against SILVA `{meta['silva']}`; "
          f"{meta['sylph']}; InSilicoSeq {meta['iss']} (`--mode kde --cpus 1`).",
          f"Rule: PROTOCOL.md 4.1, section 5 (T3) and section 6 (K3), with section 10, entries 4 to 13. "
          f"t = {meta['t']:.4g} and m = {meta['m']:.4g} (entry 10). Reads: the first 20,000,000 spots of each run.",
          "", f"## K3: {'passed' if not fails else 'FAILED'}", ""]
    md += [f"- {x}" for x in fails[:200]] or ["No spike and no outside set is called resident, and no lineage "
                                              "becomes resident with a spike."]
    if len(fails) > 200:
        md.append(f"- ... and {len(fails) - 200} more (spikes.tsv, lab_spike.tsv)")

    levels_pct = [pct(lv) for _, lv in LEVELS]
    md += ["", "## Single-run spikes", "",
           f"{len(srows):,} spikes: each kit-contaminant genome into each background run at "
           + ", ".join(f"{x}%" for x in levels_pct) + f" of its reads, seeds {', '.join(map(str, SEEDS))}. Each cell "
           "gives the spiked genus's call in the isolate over the seeds and, in brackets, its level in the spiked "
           "run. Expected: never resident.", "",
           "| Background | Isolate | Genus | " + " | ".join(f"{x}%" for x in levels_pct) + " | K3 |",
           "|---|---|---|" + "---|" * len(levels_pct) + "---|"]
    cells = defaultdict(lambda: defaultdict(list))
    for x in srows:
        cells[(x["background"], x["unit"], x["genus"])][x["level_pct"]].append(x)
    for (bg, unit, genus), by in cells.items():
        def cell(xs):
            c = Counter(x["call_after"] for x in xs)
            lv = Counter(x["level_after"] for x in xs)
            return ", ".join(f"{k} {v}" for k, v in sorted(c.items())) + \
                " (" + ", ".join(f"{k} {v}" for k, v in sorted(lv.items())) + ")"
        k3 = "FAILS" if any(x["k3"] == "FAILS" for xs in by.values() for x in xs) else "passes"
        md.append(f"| {bg} | {unit} | {genus} | " + " | ".join(cell(by[x]) for x in levels_pct) + f" | {k3} |")
    md += ["", "Calls of other lineages that changed with a spike are in `spikes.tsv` (column changed_calls)."]

    md += ["", "## Lab kit spike", "",
           f"*Ralstonia pickettii* ({LAB_GENOME}) at {pct(LAB_LEVEL)}% of the reads of every run marked lab_spike, "
           f"seeds {', '.join(map(str, SEEDS))}. Expected: contaminant, through E5.", "",
           "| Isolate | Runs spiked | Seed | Call without | Call with the spike | Basis | E1 / E2 / E5 | K3 |",
           "|---|---|---|---|---|---|---|---|"]
    for x in lrows:
        md.append(f"| {x['unit']} | {x['runs_spiked']} | {x['seed']} | {x['call_before'] or 'none'} | "
                  f"{x['call_after']} | {x['basis_after']} | {x['e1']} / {x['e2']} / {x['e5']} | {x['k3']} |")

    for s in ("hyphosphere", "added_bacterium"):
        md += ["", f"## Outside the fungus: {SET_NAME[s]}", "", "Expected: nothing resident.", ""]
        md += calls_table(base[s][1], units[s]) if s in base else ["No runs."]
    md += ["", "## The Ottawa set without spikes", "",
           "Reported, not part of K3: the calls the backgrounds give on their own.", ""]
    md += calls_table(base["ottawa"][1], units["ottawa"]) if "ottawa" in base else ["No runs."]
    md += ["", "## The route per run (no spikes)", "",
           "| Run | Set | Routed reads found | SSU kept | SSU unplaced: several genera / under 60 bp | SingleM kept | "
           "Organelle DNA | Split reads: bacterial / host |", "|---|---|---|---|---|---|---|---|"]
    for r in runs:
        st = meta["stats"][r["run"]]
        md.append(f"| {r['run']} | {r['set']} | {st['found']} | {st['ssu_kept']} | {st['ssu_several']} / "
                  f"{st['ssu_short']} | {st['sm_kept'] if r['route'] == 'ssu+singlem' else 'n/a'} | "
                  f"{st['organelle_dna']} | {st['split_bact']} / {st['split_host']} |")
    md += ["", "Spikes are added to their background run's results (entry 13): every rule in 4.1 judges each read on "
               "its own, and the spread check is worked out on both read sets together. sylph is not run on spikes. "
               "`levels.tsv` and `calls.tsv` give every set without spikes; `spikes.tsv` and `lab_spike.tsv` give "
               "each spike's calls. Reads are counted per read (mates separately)."]
    out["T3.md"] = "\n".join(md) + "\n"
    for name, text in out.items():
        (t3_dir / (name + ".part")).write_text(text, encoding="utf-8")
    for name in out:
        os.replace(t3_dir / (name + ".part"), t3_dir / name)


# --------------------------------------------------------------------- main ---

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--refdir", default="references")
    ap.add_argument("--t3-dir", default="t3", help="holds runs.tsv; the outputs to commit go here too")
    ap.add_argument("--stage0a-dir", default="stage0a")
    ap.add_argument("--controls", default="controls.tsv")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--jobs", type=int, default=3, help="runs fetched, and spikes simulated, at the same time [3]")
    ap.add_argument("--fastq-dump", default="fastq-dump")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--ssuplex", default="ssuplex")
    ap.add_argument("--conda", default=os.environ.get("CONDA_EXE") or shutil.which("conda"))
    ap.add_argument("--singlem-env", default="phase0-markers")
    ap.add_argument("--vsearch", default="vsearch")
    ap.add_argument("--sylph", default="sylph")
    ap.add_argument("--iss", default="iss")
    args = ap.parse_args(argv)
    work, t3_dir, stage_dir = Path(args.workdir).expanduser(), Path(args.t3_dir), Path(args.stage0a_dir)
    refdir, wd = Path(args.refdir), work / "t3"
    args.index = work / "index" / "full"

    try:
        runs = read_runs(t3_dir / "runs.tsv")
        controls = {x["control_id"] for x in M2.read_tsv(Path(args.controls))}
        missing = sorted({r["control_id"] for r in runs} - controls)
        if missing:
            raise RuntimeError(f"control_id not in {args.controls}: {', '.join(missing)}")
        t, m = T.stage_values(stage_dir)
        info, lineage_of, _, hosts = AR.load_index(args.index)
        if info.get("set") != "full" or not (args.index / "organelles.fna").exists():
            raise RuntimeError(f"{args.index} is not the full index of make_index.py v0.2")
        overlaps = refdir / "overlaps.tsv"
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
        lin_rows = M2.read_tsv(refdir / "lineages.tsv")
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.CalledProcessError) as e:
        print(f"Cannot start: {e}", file=sys.stderr)
        return 1
    bact_kinds = {x["lineage"]: x["kind"] for x in lin_rows if x["kind"] in ("endosymbiont", "contaminant")}
    order = [x for x in ("CaGg", "CaMg/MRE", "Mycoavidus") if x in bact_kinds] + \
        sorted(x for x in bact_kinds if x not in ("CaGg", "CaMg/MRE", "Mycoavidus"))
    ctx = {"mmi": args.index / AR.INDEX_FILE["short"], "lineage_of": lineage_of, "hosts": hosts,
           "mask": AR.load_mask(overlaps), "minimap2": mm2, "refs_md5": info["refs_md5"],
           "overlaps_md5": AR.md5_of(overlaps), "bacterial": set(bact_kinds), "bacterial_order": order,
           "syldb": args.index / "bacteria.syldb", "acc_lineage": {x["accession"]: x["lineage"] for x in lin_rows},
           "organelles": args.index / "organelles.fna"}
    try:
        genomes = kit_genomes(refdir, work, ctx)
    except (OSError, KeyError, RuntimeError) as e:
        print(f"Cannot start: {e}", file=sys.stderr)
        return 1
    runs_by = {r["run"]: r for r in runs}
    for sub in ("reads", "assign", "route", "organelle", "split", "sylph", "spikes", "spike_route"):
        (wd / sub).mkdir(parents=True, exist_ok=True)
    tools_checked = []

    def check_tools():
        if not tools_checked:
            versions = DB.check_versions(args)
            if versions != db["versions"]:
                raise RuntimeError(f"the tools here ({versions}) are not those get_marker_databases.py checked "
                                   f"({db['versions']})")
            tools_checked.append(True)

    # 1. the runs, as T1 processes them
    def reads_dir(run):
        return work / "t1" / "reads" if (work / "t1" / "reads" / run / "fetch.json").exists() else wd / "reads"

    def one_fetch(r):
        try:
            info_, now_ = T.fetch(r["run"], reads_dir(r["run"]), T.SPOTS["short"], args)
            print(f"{'fetched' if now_ else 'already fetched'}: {r['run']} ({info_['spots']:,} spots)", flush=True)
            return None
        except (OSError, RuntimeError) as e:
            return f"{r['run']}: {e}"

    failed, files, counts = [], {}, {}
    print(f"fetching {len(runs)} runs, {args.jobs} at a time", flush=True)
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        failed += [e for e in pool.map(one_fetch, runs) if e]
    for r in runs:
        run = r["run"]
        try:
            fdir = reads_dir(run) / run
            if not (fdir / "fetch.json").exists():
                continue
            fi = json.loads((fdir / "fetch.json").read_text())
            if not (wd / "assign" / run / "counts.tsv").exists():
                print(f"assigning {run}", flush=True)
            T.assign(run, [fdir / f for f in fi["files"]], wd / "assign" / run, "short", args)
            files[run], counts[run] = T.assigned(wd / "assign" / run, "short", ctx)
        except (OSError, KeyError, RuntimeError, json.JSONDecodeError) as e:
            failed.append(f"{run}: {e}")
            print(f"FAILED: {run}: {e}", flush=True)
    for r in runs:
        run = r["run"]
        if run not in files:
            continue
        src, view = work / T.WORK_NAME / "route" / run, wd / "route" / run
        try:
            if not (view / "route.json").exists():
                if not (src / "route.json").exists():
                    check_tools()
                print(f"route: {run}" + (" (T1's tool outputs used again)" if (src / "route.json").exists() else ""),
                      flush=True)
            T.route_view(r, files[run], src, view, args, db, ctx)
        except (OSError, KeyError, ValueError, RuntimeError, subprocess.CalledProcessError) as e:
            failed.append(f"{run}: {e}")
            print(f"FAILED: {run}: {e}", flush=True)
    if failed:
        print(f"\n{len(failed)} run(s) not done, so nothing was written to {t3_dir}/. Run the same command again to "
              "finish them:\n" + "\n".join(failed), file=sys.stderr)
        return 1

    # 2. the spikes
    try:
        plan = spike_plan(runs, genomes, {run: c[("total", "reads in input")] for run, c in counts.items()})
        sdirs = {p["key"]: wd / "spikes" / p["key"] for p in plan}
        todo = [p for p in plan if not (sdirs[p["key"]] / "spike.json").exists()]
        locked = DB.locked_version("phase0-reads", "insilicoseq")
        if todo:
            words = DB.tool_output([args.iss, "--version"]).split()
            iss_version = words[-1] if words else "missing"
            if iss_version != locked:
                raise RuntimeError(f"InSilicoSeq here is {iss_version}; env/phase0-reads.lock.yml has {locked}")
        else:
            iss_version = locked
        made = {json.loads((sdirs[p["key"]] / "spike.json").read_text())["insilicoseq"] for p in plan
                if p not in todo}
        if made - {locked}:
            raise RuntimeError(f"spikes in {wd / 'spikes'} were made with InSilicoSeq {', '.join(sorted(made))}, not "
                               f"{locked}")
    except (OSError, RuntimeError, KeyError, json.JSONDecodeError) as e:
        print(f"Cannot make the spikes: {e}", file=sys.stderr)
        return 1
    print(f"{len(plan):,} spikes ({len(todo):,} to simulate, {args.jobs} at a time)", flush=True)

    def one_spike(p):
        try:
            make_spike(p, sdirs[p["key"]], args, iss_version)
            return None
        except (OSError, RuntimeError, ValueError) as e:
            return f"{p['key']}: {e}"

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        for i, e in enumerate(pool.map(one_spike, todo)):
            if e:
                failed.append(e)
                print(f"FAILED: {e}", flush=True)
            if (i + 1) % 50 == 0:
                print(f"spikes simulated: {i + 1} of {len(todo)}", flush=True)
    smeta, scounts, sfiles = {}, {}, {}
    for p in plan:
        key, sdir = p["key"], sdirs[p["key"]]
        if not (sdir / "spike.json").exists():
            continue
        try:
            smeta[key] = json.loads((sdir / "spike.json").read_text())
            T.assign(key, [sdir / f for f in smeta[key]["files"]], sdir / "assign", "short", args)
            sfiles[key], scounts[key] = T.assigned(sdir / "assign", "short", ctx)
        except (OSError, KeyError, RuntimeError, json.JSONDecodeError) as e:
            failed.append(f"{key}: {e}")
            print(f"FAILED: {key}: {e}", flush=True)
    groups = defaultdict(list)
    for p in plan:
        r = runs_by[p["run"]]
        g = f"single_{p['run']}" if p["part"] == "single" else f"lab_{'dna' if r['route'] == 'ssu+singlem' else 'rna'}"
        groups[(g, r["route"])].append(p["key"])
    pools = {}
    if not failed:
        for (g, route), keys in groups.items():
            pdir = wd / "spike_route" / g
            try:
                pm = make_pool(g, keys, route, sdirs, pdir)
                pf = [(pdir / f, n) for f, n in zip(pm["files"], pm["reads"], strict=True) if n > 0]
                if pf and not (pdir / "route" / "route.json").exists():
                    check_tools()
                    print(f"route: spikes of {g} ({sum(pm['reads']):,} routed reads)", flush=True)
                if pf:
                    T.route_new({"run": g, "route": route}, pf, pdir / "route", args, db, ctx)
                pools[g] = ({"run": g, "route": route}, pdir / "route" if pf else None)
            except (OSError, KeyError, ValueError, RuntimeError, subprocess.CalledProcessError) as e:
                failed.append(f"spikes of {g}: {e}")
                print(f"FAILED: spikes of {g}: {e}", flush=True)
    if failed:
        print(f"\n{len(failed)} run(s) or spike(s) not done, so nothing was written to {t3_dir}/. Run the same "
              "command again to finish them:\n" + "\n".join(failed[:100]), file=sys.stderr)
        return 1

    # 3. route reads, placement, data
    try:
        found, pieces, sm_tax = {}, {}, {}
        sources = {run: (runs_by[run], wd / "route" / run) for run in files} | \
            {g: v for g, v in pools.items() if v[1] is not None}
        for name, (r, rdir) in sources.items():
            found[name] = T.routed_found(r, rdir)
            pieces[name] = R.ssu_sequences(rdir, {k for k, v in found[name].items() if v[1]})
            sm_tax[name] = T.singlem_taxonomies(rdir, {k for k, v in found[name].items() if v[2]})
        queries = {}
        for ps in pieces.values():
            for s in ps.values():
                queries.setdefault(R.seq_id(s), s)
        earlier = [work.joinpath(*p) for p in EARLIER]
        before = len(RT.read_batches(earlier)[0] & queries.keys())
        print(f"SSU pieces: {len(queries):,}, {before:,} of them placed before", flush=True)
        hits = RT.placements(queries, wd / "place" / "batches", silva["fasta"], args, wd / "place" / "vsearch.log",
                             earlier)
        targets = {t_ for q in queries for t_, _ in hits.get(q, [])}
        tax = {k: [x.strip() for x in v.split(";")] for k, v in R.load_taxonomy(silva["taxonomy"], targets).items()}

        data, bits, stats = {}, {}, {}
        for r in runs:
            run, c = r["run"], counts[r["run"]]
            d = data[run] = new_data(c, r["kind"])
            d["mapped"] = {lin: c.get(("lineage", lin), 0) for lin in ctx["bacterial"]}
            if r["kind"] == "DNA":
                d["spread"] = AR.read_spread(wd / "assign" / run / "spread.tsv")
                bits[run] = str(wd / "assign" / run / "spread_bits.json.gz")
                if not Path(bits[run]).exists():
                    raise RuntimeError(f"{wd / 'assign' / run} has no spread_bits.json.gz: assigned before v0.3?")
            print(f"route reads: {run}", flush=True)
            organelle = T.organelle_check(found[run], wd / "route" / run, wd / "organelle" / f"{run}.json", args, ctx)
            split = T.split_check(found[run], wd / "route" / run, wd / "split" / f"{run}.json", args, ctx)
            tally_route(found[run], pieces[run], sm_tax[run], hits, tax, organelle, split, lambda k, x=run: x, data)
            print(f"sylph: {run}", flush=True)
            d["sylph"] = T.sylph(files[run], wd / "sylph" / f"{run}.tsv", args, ctx)
            stats[run] = d["stats"]
        spikes = {}
        for p in plan:
            key, c = p["key"], scounts[p["key"]]
            spikes[key] = new_data(c, runs_by[p["run"]]["kind"])
            spikes[key]["mapped"] = {lin: c.get(("lineage", lin), 0) for lin in ctx["bacterial"]}
            bits[key] = str(sdirs[key] / "assign" / "spread_bits.json.gz")
            if not Path(bits[key]).exists():
                raise RuntimeError(f"{sdirs[key] / 'assign'} has no spread_bits.json.gz: assigned before v0.3?")
        for g, (_, rdir) in pools.items():
            if rdir is None:
                continue
            print(f"route reads: spikes of {g}", flush=True)
            organelle = T.organelle_check(found[g], rdir, rdir.parent / "organelle.json", args, ctx)   # pool's
            split = T.split_check(found[g], rdir, rdir.parent / "split.json", args, ctx)             # own folder
            tally_route(found[g], pieces[g], sm_tax[g], hits, tax, organelle, split,
                        lambda k: k[1].split(SEP)[0], spikes)

        # 4. calls
        base = {}
        for s in SETS:
            mine = [r for r in runs if r["set"] == s]
            if mine:
                base[s] = T.build(mine, {r["run"]: data[r["run"]] for r in mine}, t, m, ctx)
        ott = [r for r in runs if r["set"] == "ottawa"]
        ott_data = {r["run"]: data[r["run"]] for r in ott}
        print("calling the single-run spikes", flush=True)
        srows, fails = single_rows(plan, runs_by, ott, ott_data, spikes, bits, base["ottawa"], t, m, ctx)
        print("calling the lab kit spike", flush=True)
        lrows, lfails = lab_rows(plan, runs_by, ott, ott_data, spikes, bits, base["ottawa"], t, m, ctx)
        ofails = [f"{SET_NAME[s]}: {c['lineage']} is called resident in {c['unit']}"
                  for s in ("hyphosphere", "added_bacterium") if s in base for c in base[s][1]
                  if c["call"] == "resident"]
    except (OSError, KeyError, ValueError, RuntimeError, subprocess.CalledProcessError) as e:
        print(f"Nothing was written to {t3_dir}/: {e}", file=sys.stderr)
        return 1
    a1 = json.loads((wd / "assign" / runs[0]["run"] / "assign.json").read_text())
    meta = {"made": now(), "assign": a1["made_by"], "minimap2": mm2, "refs_md5": info["refs_md5"],
            "tools": ", ".join(f"{k} {v}" for k, v in sorted(db["versions"].items())), "vsearch": vsearch,
            "silva": silva["file"], "sylph": sylph_v, "iss": iss_version, "t": t, "m": m, "stats": stats}
    write_outputs(t3_dir, runs, base, srows, lrows, fails + lfails + ofails, meta)
    all_fails = fails + lfails + ofails
    print(f"\nK3: {'passed' if not all_fails else 'FAILED'}" + "".join(f"\n  {x}" for x in all_fails[:20]))
    print(f"-> {t3_dir}/levels.tsv, calls.tsv, spikes.tsv, lab_spike.tsv, T3.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
