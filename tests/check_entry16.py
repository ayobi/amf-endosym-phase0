"""Check entry 16: the two lists of known clades in call_lineages.py v0.3, and rescore_entry16.py v0.1
(PROTOCOL.md section 10, entry 16).

Part 1, the caller. Phase 0's list is as locked and is the default; the survey's list adds only CGIF01, as
CaMg/MRE; category_of_route uses it only when given it; JAYUTB01 and the other Mycoplasmatales genera stay
other known genera; SingleM's taxonomy for CGIF01 gives the genus CGIF01; no Phase 0 script uses the
survey's list.

Part 2, the call with the survey's list as rescore_entry16.py works it out from a row of a calls table,
against call_lineages.call() with the new category: hand-made cases and 5,000 random ones. Rows it must not
take (a category Phase 0's list does not give, an unknown call) are refused.

Part 3, rescore_entry16.py on synthetic T1 to T5 tables, written by run_t1.py to run_t5.py's own writers from
run_t1.build: which calls change (call, or basis only), K1 both ways (plant-only roots resident with the
survey's list), K2, K3 for the outside sets both ways, the Ottawa lines at risk from spikes and the spikes
that changed a CGIF01 call, T4's K1 on Logan and its comparison with raw reads both ways, T5. Tests not run
yet are listed. Refusals, with nothing written: a table made with the survey's list, a report that disagrees
with its tables, a compare.tsv that does not match, a call twice, a line in two T3 sets, an unknown arm, a
missing column, T4 without T1, no tables at all. A rerun gives the same files.

Needs only Python (no tools).
"""
import contextlib
import csv
import io
import random
import shutil
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

here = Path(__file__).parent
sys.path.insert(0, str(here.parent / "scripts"))
import call_lineages as C  # noqa: E402
import rescore_entry16 as R  # noqa: E402
import run_t1 as T  # noqa: E402
import run_t2 as T2  # noqa: E402
import run_t3 as T3  # noqa: E402
import run_t4 as T4  # noqa: E402
import run_t5 as T5  # noqa: E402

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)


# ===================================================================== part 1 ===

check(C.__version__ == "0.3", "call_lineages.py is v0.3")
check(C.KNOWN_CLADE_GENERA == {"Glomeribacter": "CaGg", "Moeniiplasma": "CaMg/MRE", "Mycoavidus": "Mycoavidus"},
      "Phase 0's list is as locked")
check(C.SURVEY_CLADE_GENERA == {**C.KNOWN_CLADE_GENERA, "CGIF01": "CaMg/MRE"}, "the survey's list adds only CGIF01")
check(C.KNOWN_CLADES == {"CaGg", "CaMg/MRE", "Mycoavidus"}, "the mapped known clades are as locked")
check(C.category_of_route(["CGIF01"]) == "other known genus", "CGIF01: other known genus with Phase 0's list")
check(C.category_of_route(["CGIF01"], C.KNOWN_CLADE_GENERA) == "other known genus", "Phase 0's list given")
check(C.category_of_route(["CGIF01"], C.SURVEY_CLADE_GENERA) == "known endosymbiont clade", "the survey's list")
for g in ("JAYUTB01", "Spiroplasma", "Hepatoplasma", "Ureaplasma", "Mycoplasmoides"):
    check(C.category_of_route([g]) == C.category_of_route([g], C.SURVEY_CLADE_GENERA) == "other known genus",
          f"{g} stays an other known genus")
KNOWN = "known endosymbiont clade"
for gs, cat in ((["Moeniiplasma"], KNOWN), (["Glomeribacter", "Mycoavidus"], KNOWN),
                (["Burkholderia-Caballeronia-Paraburkholderia"], "kit-contaminant genus"), ([], "novel")):
    check(C.category_of_route(gs) == C.category_of_route(gs, C.SURVEY_CLADE_GENERA) == cat, f"{gs}: {cat} both ways")
check(C.clean_genus("g__CGIF01") == "CGIF01" and C.clean_genus("g__CGIF01_A") == "CGIF01", "CGIF01 as a genus name")
check(C.singlem_label(["Root; d__Bacteria; p__Bacillota_I; c__Bacilli_A; o__Mycoplasmatales; f__CGIF01; g__CGIF01"])
      == ("kept", "CGIF01", "genus CGIF01"), "SingleM's CGIF01 taxonomy gives the genus CGIF01")
for p in sorted((here.parent / "scripts").glob("*.py")):
    if p.name not in ("call_lineages.py", "rescore_entry16.py"):
        check("SURVEY_CLADE_GENERA" not in p.read_text(encoding="utf-8"), f"{p.name} does not use the survey's list")


# ===================================================================== part 2 ===

def row_of(res, category="other known genus", lineage="CGIF01 (route)", members="CGIF01", source="route"):
    return {"unit": "u", "lineage": lineage, "source": source, "members": members, "category": category,
            "call": res["call"], "basis": res["basis"], "present_runs": str(res["present_runs"])}


def rec(level, kind="DNA", bioproject="P1", year="", n=100, fungal=10_000, informative=True, fraction=0.001):
    return {"level": level, "kind": kind, "bioproject": bioproject, "year": year, "n": n, "fungal": fungal,
            "informative": informative, "fraction": fraction}


def both(recs, bg=None):
    r0 = C.call(recs, "other known genus", bg or {})
    r1 = C.call(recs, "known endosymbiont clade", bg or {})
    got = R.survey_call(row_of(r0))
    check(got == ("known endosymbiont clade", r1["call"], r1["basis"]),
          f"survey call {got} is not call()'s {r1['call']} ({r1['basis']}) for {recs} {bg}")
    return r0["call"], r1["call"]


check(both([rec("present", "RNA"), rec("present", "RNA")]) == ("unresolved", "resident"), "unresolved, 2 present")
check(both([rec("present"), rec("low")]) == ("unresolved", "unresolved"), "unresolved, 1 present")
check(both([rec("present", bioproject="P1"), rec("present", bioproject="P2")]) == ("resident", "resident"),
      "E1 and E2: resident both ways")
check(both([rec("present"), rec("present")], {"a": [0.001], "b": [0.002]}) == ("contaminant", "contaminant"),
      "E5 fails: contaminant both ways")
check(both([rec("absent"), rec("absent")]) == ("absent", "absent"), "absent both ways")
check(both([rec("absent", informative=False, n=0)]) == ("no call", "no call"), "no call both ways")
r0 = C.call([rec("present", bioproject="P1"), rec("present", bioproject="P2")], "other known genus", {})
check(R.survey_call(row_of(r0))[2] == "known endosymbiont clade, present in 2 runs", "E1 and E2: the new basis")

rng = random.Random(16)
seen = Counter()
for _ in range(5000):
    recs = [rec(rng.choice(["present", "present", "low", "absent"]), rng.choice(["DNA", "DNA", "RNA"]),
                rng.choice(["P1", "P2"]), rng.choice(["", "", "2019", "2020"]), rng.randint(0, 400),
                rng.choice([0, 5_000, rng.randint(1, 100_000)]), rng.random() < 0.7, rng.random() * 0.01)
            for _ in range(rng.randint(1, 5))]
    bg = {f"i{j}": [rng.random() * 0.01 for _ in range(rng.randint(1, 2))] for j in range(rng.choice([0, 0, 1, 2, 3]))}
    seen[both(recs, bg)] += 1
for pair in (("unresolved", "resident"), ("unresolved", "unresolved"), ("resident", "resident"),
             ("contaminant", "contaminant"), ("absent", "absent"), ("no call", "no call")):
    check(seen[pair] > 0, f"random cases cover {pair}")
check(set(seen) <= {("unresolved", "resident"), ("unresolved", "unresolved"), ("resident", "resident"),
                    ("contaminant", "contaminant"), ("absent", "absent"), ("no call", "no call")},
      f"no other change: {sorted(seen)}")

unchanged = [row_of({"call": "unresolved", "basis": "x", "present_runs": 3}, "other known genus", "JAYUTB01 (route)",
                    "JAYUTB01"),
             row_of({"call": "resident", "basis": "y", "present_runs": 2}, "known endosymbiont clade",
                    "Moeniiplasma (route)", "Moeniiplasma"),
             row_of({"call": "resident", "basis": "z", "present_runs": 4}, "known endosymbiont clade", "CaMg/MRE", "",
                    "mapped"),
             row_of({"call": "contaminant", "basis": "w", "present_runs": 2}, "kit-contaminant genus", "Ralstonia", "",
                    "mapped"),
             row_of({"call": "unresolved", "basis": "v", "present_runs": 2}, "novel", "novel (SingleM)", "")]
for x in unchanged:
    check(R.survey_call(x) == (x["category"], x["call"], x["basis"]), f"{x['lineage']} does not change")
for x, why in ((row_of({"call": "resident", "basis": "b", "present_runs": 2}, "known endosymbiont clade"),
                "a CGIF01 row filed as a known clade"),
               (row_of({"call": "maybe", "basis": "b", "present_runs": 2}), "an unknown call"),
               (row_of({"call": "unresolved", "basis": "b", "present_runs": 2}, KNOWN, "Pseudomonas", "", "mapped"),
                "a mapped lineage with the wrong category")):
    try:
        R.survey_call(x)
        fails.append(f"survey_call took {why}")
    except ValueError:
        pass


# ===================================================================== part 3 ===

t, m = 0.001, 0.25                     # N = 100 mapped and 25 route reads in a run of 100,000 reads
STATS = dict.fromkeys(["found", "ssu_kept", "ssu_outside", "ssu_organelle", "ssu_several", "ssu_short", "sm_kept",
                       "organelle_dna", "split_bact", "split_host", "labelled"], 0)
CTX = {"bacterial_order": ["CaGg", "CaMg/MRE", "Mycoavidus", "Burkholderia", "Ralstonia"]}
CTX_HELD = {"bacterial_order": ["CaGg", "Mycoavidus", "Burkholderia", "Ralstonia"]}
CG, JA, MO = ("genus", "CGIF01"), ("genus", "JAYUTB01"), ("genus", "Moeniiplasma")
META = {"assign": "assign_reads.py v0.3", "minimap2": "2.31", "refs_md5": "x", "tools": "stand-ins",
        "vsearch": "2.32.0", "silva": "s", "sylph": "sylph not run", "t": t, "m": m, "made": "2026-11-25 00:00 UTC",
        "refs_full": "x", "refs_t2": "y", "removed": 0, "full_pieces": 0, "full_before": 0, "t2_pieces": 0,
        "derived": 0, "searched": 0, "checked": 0, "iss": "2.0", "logan_url": "u", "zstd": "zstd"}


def run(run_id, unit, kind, bioproject, lab, isolate, control_id="", prep="", scored="yes", **extra):
    prep = prep or ("dna" if kind == "DNA" else "total RNA")
    return {"run": run_id, "unit": unit, "control_id": control_id or unit.replace(" ", "_"), "scored": scored,
            "kind": kind, "platform": "short", "prep": prep, "layout": "PAIRED", "bioproject": bioproject, "year": "",
            "lab": lab, "isolate": isolate, "route": "ssu+singlem" if kind == "DNA" else "ssu", "reuse": "",
            "note": "", **extra}


def data(kind, mapped=None, route=None, reads=100_000):
    mapped = Counter(mapped or {})
    return {"reads": reads, "fungal": 50_000, "mapped": mapped, "route": Counter(route or {}),
            "spread": {x: 0.9 for x in mapped} if kind == "DNA" else None, "split": {}, "sylph": {},
            "placements": {}, "stats": dict(STATS)}


def write_runs(path, runs, cols):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(runs)


def read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def make_world(root):
    """Synthetic T1 to T5 tables in root, made by the tests' own writers. Returns the controls.tsv path."""
    dirs = {x: root / x for x in ("t1", "t2", "t3", "t4", "t5")}
    for x in dirs.values():
        x.mkdir(parents=True, exist_ok=True)
    # T1: the B+ line is resident for CGIF01 by E1 and E2; the plant-only roots have CGIF01 in two RNA runs
    t1 = [run("P1", T.B_PLUS, "DNA", "PRJ1", "Ottawa", "BEG34", T.K1_PRESENT[0]),
          run("P2", T.B_PLUS, "DNA", "PRJ2", "Ottawa", "BEG34", T.K1_PRESENT[1]),
          run("P3", T.B_PLUS, "RNA", "PRJ3", "Ottawa", "BEG34", "GMAR_BEG34_RNA_Bplus"),
          run("M1", T.B_MINUS, "DNA", "PRJ1", "Ottawa", "BEG34", "GMAR_BEG34_HIC_Bminus"),
          run("M2", T.B_MINUS, "RNA", "PRJ3", "Ottawa", "BEG34", "GMAR_BEG34_RNA_Bminus"),
          run("L1", T.PLANT, "RNA", "PRJ4", "Kyoto", "none", "LJAP_PLANT"),
          run("L2", T.PLANT, "RNA", "PRJ4", "Kyoto", "none", "LJAP_PLANT")]
    d1 = {"P1": data("DNA", {"CaGg": 5000}, {CG: 60}), "P2": data("DNA", {"CaGg": 4000}, {CG: 55}),
          "P3": data("RNA", {"CaGg": 300}), "M1": data("DNA", {}, {CG: 3}), "M2": data("RNA"),
          "L1": data("RNA", {}, {CG: 30}), "L2": data("RNA", {}, {CG: 30})}
    write_runs(dirs["t1"] / "runs.tsv", t1, T.RUN_COLS)
    runs1 = T.read_runs(dirs["t1"] / "runs.tsv")
    lv1, cl1, _ = T.build(runs1, d1, t, m, CTX)
    T.write_outputs(dirs["t1"], runs1, d1, lv1, cl1, T.k1(cl1, lv1, runs1), {}, META)
    # T2: CGIF01 in the held-out MRE positives; JAYUTB01 too; Moeniiplasma (route) in the kept-in arm
    t2 = [run(f"D{i}", "D. epigaea S1", "DNA", "PRJD", "Cornell", "S1") for i in (1, 2, 3)] + \
         [run(f"A{i}", "R. irregularis A2", "DNA", "PRJM", "Bonn", "A2") for i in (1, 2)] + \
         [run(f"Q{i}", "Rhizophagus sp. LPA8", "DNA", "PRJM", "Bonn", "LPA8") for i in (1, 2)]
    held = {"D1": {CG: 100, JA: 30}, "D2": {CG: 100, JA: 30}, "D3": {CG: 100},
            "A1": {CG: 80, ("novel", "SingleM"): 200}, "A2": {CG: 80}, "Q1": {CG: 60}, "Q2": {CG: 4}}
    kept = {"D1": ({"CaMg/MRE": 3000}, {CG: 40, MO: 30}), "D2": ({"CaMg/MRE": 3000}, {CG: 40, MO: 30}),
            "D3": ({"CaMg/MRE": 3000}, {}), "A1": ({"CaMg/MRE": 2000}, {}), "A2": ({"CaMg/MRE": 2000}, {}),
            "Q1": ({"CaMg/MRE": 1000}, {}), "Q2": ({"CaMg/MRE": 1000}, {})}
    write_runs(dirs["t2"] / "runs.tsv", t2, T.RUN_COLS)
    runs2 = T.read_runs(dirs["t2"] / "runs.tsv")
    d2 = {"held_out": {r: data("DNA", {}, x) for r, x in held.items()},
          "kept_in": {r: data("DNA", *x) for r, x in kept.items()}}
    out2 = {}
    for arm, ctx in (("held_out", CTX_HELD), ("kept_in", CTX)):
        lv, cl, _ = T.build(runs2, d2[arm], t, m, ctx)
        out2[arm] = {"levels": lv, "calls": cl, "data": d2[arm]}
    k2_rows = [{"unit": u, "step1": 0, "found": 0, "labelled": 0, "organelle_dna": 0, "split": 0,
                "lineage": "CGIF01 (route)", "lineage_mre": 0, "call": "unresolved", "basis": "", "passed": True,
                "others": []}
               for u in dict.fromkeys(r["unit"] for r in runs2)]
    T2.write_outputs(dirs["t2"], runs2, out2, k2_rows, [], META)
    # T3: CGIF01 in the Ottawa set (contaminant through E5, and low in C2), and resident in the hyphosphere set
    t3 = [run(f"O{u}{i}", f"R. irregularis {u}", "DNA", p, "Ottawa", u, set="ottawa", background="yes", lab_spike="yes",
              iss_model="HiSeq", instrument="HiSeq")
          for u, p in (("4401", "PRJO1"), ("A1", "PRJO1"), ("B3", "PRJO2"), ("C2", "PRJO2")) for i in (1, 2)] + \
         [run(f"H{i}", "MUCL 43194 hyphosphere", "DNA", "PRJH", "Louvain", "MUCL43194", set="hyphosphere",
              background="no", lab_spike="no", iss_model="", instrument="") for i in (1, 2)] + \
         [run("X1", "co-cultured bacterium", "DNA", "PRJX", "Lab2", "X1", set="added_bacterium", background="no",
              lab_spike="no", iss_model="", instrument="")]
    d3 = {r["run"]: data("DNA", {}, {CG: 3 if "C2" in r["run"] else 100}) for r in t3 if r["set"] == "ottawa"}
    d3.update({"H1": data("DNA", {}, {CG: 50}), "H2": data("DNA", {}, {CG: 50}),
               "X1": data("DNA", {"Burkholderia": 2000})})
    write_runs(dirs["t3"] / "runs.tsv", t3, T.RUN_COLS + T3.EXTRA_COLS)
    base = {}
    for s in T3.SETS:
        mine = [r for r in t3 if r["set"] == s]
        base[s] = T.build(mine, {r["run"]: d3[r["run"]] for r in mine}, t, m, CTX)
    srow = dict.fromkeys(T3.SPIKE_COLS, "")
    srow.update({"spike": "s1", "background": "O44011", "unit": "R. irregularis 4401", "genus": "Burkholderia",
                 "level_pct": "0.1", "seed": 1, "call_after": "contaminant", "level_after": "present", "k3": "passes",
                 "changed_calls": "CGIF01 (route): contaminant -> unresolved; Burkholderia: absent -> contaminant"})
    lrow = dict.fromkeys(T3.LAB_COLS, "")
    lrow.update({"seed": 1, "unit": "R. irregularis A1", "k3": "passes",
                 "changed_calls": "Ralstonia: absent -> contaminant"})
    T3.write_outputs(dirs["t3"], t3, base, [srow], [lrow], [], {**META, "stats": {r["run"]: STATS for r in t3}})
    # T4: Logan lacks M2; on Logan the plant-only roots and D. epigaea have CGIF01 present in one run only
    have = {r["run"] for r in t1 + t2} - {"M2"}
    lg1 = {**d1, "L2": data("RNA", {}, {CG: 10})}
    lg_held = {**d2["held_out"], "D2": data("DNA", {}, {CG: 10, JA: 30}), "D3": data("DNA", {}, {CG: 10})}
    built = {}
    for (test, arm), runs, dd, ctx in ((("T1", "full"), runs1, lg1, CTX),
                                       (("T2", "kept_in"), runs2, d2["kept_in"], CTX),
                                       (("T2", "held_out"), runs2, lg_held, CTX_HELD)):
        mine = [r for r in runs if r["run"] in have]
        built[(test, arm)] = (mine, *T.build(mine, {r["run"]: dd[r["run"]] for r in mine}, t, m, ctx))
    raw = {("T1", "full"): (read(dirs["t1"] / T.OUT_LEVELS), read(dirs["t1"] / T.OUT_CALLS))}
    for arm in ("kept_in", "held_out"):
        raw[("T2", arm)] = ([x for x in read(dirs["t2"] / "levels.tsv") if x["arm"] == arm],
                            [x for x in read(dirs["t2"] / "calls.tsv") if x["arm"] == arm])
    left, left_projects = defaultdict(list), set()
    for r in runs1 + runs2:
        if r["run"] not in have:
            left[r["unit"]].append((r["run"], "not in Logan"))
            left_projects.add(r["bioproject"])
    cmp_rows = []
    for (test, arm), (_, lv, cl, _) in built.items():
        cmp_rows += T4.compare(test, arm, raw[(test, arm)], (lv, cl), have, left, left_projects)
    logan_rows = [{**dict.fromkeys(T4.LOGAN_COLS, ""), "run": r["run"], "test": test, "unit": r["unit"],
                   "in_logan": "yes" if r["run"] in have else "no"}
                  for test, runs in (("T1", runs1), ("T2", runs2)) for r in runs]
    mine1, lvl1, cll1, _ = built[("T1", "full")]
    T4.write_outputs(dirs["t4"], logan_rows, {k: v[1:] for k, v in built.items()}, cmp_rows, ([], None, []),
                     T.k1(cll1, lvl1, mine1), None, META)
    # T5: Mycoavidus resident in the wild type, no call in the cured line; CGIF01 in the wild type
    t5 = [run(f"W{i}", "AD073 wild type", "RNA", "PRJJ", "JGI", "AD073", "LELO_AD073_WT", "unstated")
          for i in (1, 2)] + \
         [run("W3", "AD073 wild type", "DNA", "PRJJD", "JGI", "AD073", "LELO_AD073_WT")] + \
         [run(f"K{i}", "AD073 cured", "RNA", "PRJJ", "JGI", "AD073", "LELO_AD073_CURED", "unstated") for i in (1, 2)]
    d5 = {"W1": data("RNA", {"Mycoavidus": 500}, {CG: 30}), "W2": data("RNA", {"Mycoavidus": 500}, {CG: 30}),
          "W3": data("DNA", {"Mycoavidus": 2000}, {CG: 40}), "K1": data("RNA"), "K2": data("RNA")}
    write_runs(dirs["t5"] / "runs.tsv", t5, T.RUN_COLS)
    runs5 = T.read_runs(dirs["t5"] / "runs.tsv")
    lv5, cl5, _ = T.build(runs5, d5, t, m, CTX)
    controls = root / "controls.tsv"
    write_runs(controls, [{"control_id": "LELO_AD073_WT", "expected_call": "resident"},
                          {"control_id": "LELO_AD073_CURED", "expected_call": "absent"}],
               ["control_id", "expected_call"])
    T5.write_outputs(dirs["t5"], runs5, lv5, cl5, "AD073 wild type", "AD073 cured",
                     {"AD073 wild type": "resident", "AD073 cured": "absent"}, META)
    return controls


def rescore(root, controls, out="entry16", **dirs):
    args = []
    for x in ("t1", "t2", "t3", "t4", "t5"):
        args += [f"--{x}-dir", str(dirs.get(x, root / x))]
    err = io.StringIO()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
        rc = R.main(args + ["--controls", str(controls), "--out-dir", str(root / out)])
    return rc, err.getvalue()


def outputs(root, out="entry16"):
    d = root / out
    return {p.name: p.read_text(encoding="utf-8") for p in sorted(d.glob("*"))} if d.exists() else {}


def edit(path, old, new):
    text = path.read_text(encoding="utf-8")
    check(old in text, f"test setup: {old!r} is in {path.name}")
    path.write_text(text.replace(old, new, 1), encoding="utf-8")


tmp = Path(tempfile.mkdtemp(prefix="check_entry16_"))
try:
    ctrl = make_world(tmp)
    rc, err = rescore(tmp, ctrl)
    check(rc == 0, f"rescore_entry16.py runs on the world: {err}")
    rows = read(tmp / "entry16" / R.OUT_TSV)
    md = (tmp / "entry16" / R.OUT_MD).read_text(encoding="utf-8")
    by = {(r["test"], r["part"], r["unit"], r["lineage"]): r for r in rows}

    def row(test, part, unit, lineage):
        return by.get((test, part, unit, lineage), {})

    def expect(test, part, unit, lineage, call, survey, changed):
        r = row(test, part, unit, lineage)
        check((r.get("call"), r.get("call_survey"), r.get("changed")) == (call, survey, changed),
              f"{test} {part} {unit} {lineage}: {r.get('call')} -> {r.get('call_survey')} ({r.get('changed')}), "
              f"not {call} -> {survey} ({changed})")

    expect("T1", "all runs", T.B_PLUS, "CGIF01 (route)", "resident", "resident", "basis only")
    check(row("T1", "all runs", T.B_PLUS, "CGIF01 (route)").get("basis_survey") ==
          "known endosymbiont clade, present in 2 runs", "the B+ line's new basis")
    expect("T1", "all runs", T.PLANT, "CGIF01 (route)", "unresolved", "resident", "call")
    expect("T1", "all runs", T.B_MINUS, "CGIF01 (route)", "unresolved", "unresolved", "no")
    expect("T1", "all runs", T.B_PLUS, "CaGg", "resident", "resident", "no")
    expect("T2", "held_out", "D. epigaea S1", "CGIF01 (route)", "unresolved", "resident", "call")
    expect("T2", "held_out", "R. irregularis A2", "CGIF01 (route)", "unresolved", "resident", "call")
    expect("T2", "held_out", "Rhizophagus sp. LPA8", "CGIF01 (route)", "unresolved", "unresolved", "no")
    expect("T2", "held_out", "D. epigaea S1", "JAYUTB01 (route)", "unresolved", "unresolved", "no")
    expect("T2", "held_out", "R. irregularis A2", "novel (SingleM)", "unresolved", "unresolved", "no")
    expect("T2", "kept_in", "D. epigaea S1", "CGIF01 (route)", "unresolved", "resident", "call")
    expect("T2", "kept_in", "D. epigaea S1", "Moeniiplasma (route)", "resident", "resident", "no")
    expect("T2", "kept_in", "D. epigaea S1", "CaMg/MRE", "resident", "resident", "no")
    for u in ("4401", "A1", "B3"):
        expect("T3", "ottawa", f"R. irregularis {u}", "CGIF01 (route)", "contaminant", "contaminant", "no")
    expect("T3", "ottawa", "R. irregularis C2", "CGIF01 (route)", "unresolved", "unresolved", "no")
    expect("T3", "hyphosphere", "MUCL 43194 hyphosphere", "CGIF01 (route)", "unresolved", "resident", "call")
    expect("T3", "added_bacterium", "co-cultured bacterium", "Burkholderia", "contaminant", "contaminant", "no")
    expect("T4", "T1 full", T.PLANT, "CGIF01 (route)", "unresolved", "unresolved", "no")
    expect("T4", "T2 held_out", "D. epigaea S1", "CGIF01 (route)", "unresolved", "unresolved", "no")
    expect("T4", "T2 held_out", "R. irregularis A2", "CGIF01 (route)", "unresolved", "resident", "call")
    expect("T5", "all runs", "AD073 wild type", "CGIF01 (route)", "unresolved", "resident", "call")
    expect("T5", "all runs", "AD073 wild type", "Mycoavidus", "resident", "resident", "no")
    tables = [tmp / "t1" / T.OUT_CALLS] + [tmp / x / "calls.tsv" for x in ("t2", "t3", "t4", "t5")]
    want = sum(r["signal"] == "yes" for f in tables for r in read(f))
    check(len(rows) == want, f"both_ways.tsv has every call with signal ({len(rows)}, not {want})")
    check(all(r["changed"] == "no" for r in rows if "CGIF01" not in r["members"]), "only CGIF01 lineages change")

    verdict = {}
    for line in md.splitlines():
        cells = [x.strip() for x in line.strip("|").split("|")]
        if len(cells) == 3 and line.startswith("| "):
            verdict[cells[0]] = (cells[1], cells[2])
    check(verdict.get("K1 (T1)") == ("passed", "FAILED: the plant-only roots are called resident for CGIF01 (route)"),
          f"K1 both ways: {verdict.get('K1 (T1)')}")
    check(verdict.get("K2 (T2)") == ("passed", R.SAME_K2), f"K2: {verdict.get('K2 (T2)')}")
    check(verdict.get("K3 (T3), outside sets") ==
          ("passed", "FAILED: resident: CGIF01 (route) in MUCL 43194 hyphosphere"),
          f"K3: {verdict.get('K3 (T3), outside sets')}")
    k3s = verdict.get("K3 (T3), spikes and the lab kit spike", ("", ""))
    check(k3s[0].startswith("passed") and "R. irregularis 4401, R. irregularis A1, R. irregularis B3, R. irregularis C2"
          in k3s[1] and k3s[1].endswith("spikes that changed a CGIF01 call: 1"), f"K3 spikes: {k3s}")
    check(verdict.get("T4: K1 on Logan (reported)") == ("would pass", "would pass"), "K1 on Logan both ways")
    n_cmp = sum(x["what"] == "call" and x["same"] != "n/a" for x in read(tmp / "t4" / "compare.tsv"))
    n_same = sum(x["what"] == "call" and x["same"] == "yes" for x in read(tmp / "t4" / "compare.tsv"))
    check(verdict.get("T4: calls the same as on raw reads") == (f"{n_same} of {n_cmp}", f"{n_same - 2} of {n_cmp}"),
          f"T4's comparison both ways: {verdict.get('T4: calls the same as on raw reads')} ({n_same} of {n_cmp})")
    check(n_cmp > 5 and any(x["same"] == "n/a" for x in read(tmp / "t4" / "compare.tsv")),
          "the T4 world compares calls and leaves some out")
    t5v = ("as expected (wild type resident, cured line no call)",) * 2
    check(verdict.get("T5 (reported): Mycoavidus") == t5v, f"T5: {verdict.get('T5 (reported): Mycoavidus')}")
    check("## Not run yet" not in md and "None." not in md.split("## Calls that change")[1].split("##")[0],
          "every test has run, and calls change")
    check(md.count("| T2 | held_out |") >= 2 and "R. irregularis A2: CGIF01 (route)" in md, "the summary table")

    check(verdict.get("T4: the switch") == ("see t4/T4.md", "the same: it uses levels only")
          and verdict.get("T4: K2 on Logan (reported)") == ("see t4/T4.md", R.SAME_K2), "T4 without a switch or K2")
    w4 = tmp / "w4"
    shutil.copytree(tmp, w4, ignore=shutil.ignore_patterns("entry16", "w4"))
    edit(w4 / "t4" / "T4.md", "## K1 and K2 on Logan", "Logan misses 1 of the 4 positive runs that are present in raw "
         "reads and that Logan holds (25%): more than 10%, so **the survey runs on raw reads**.\n\nK2: would pass.\n\n"
         "## K1 and K2 on Logan")
    rc, err = rescore(w4, w4 / "controls.tsv")
    v4 = {x.split(" | ")[0].strip("| "): x for x in (w4 / "entry16" / R.OUT_MD).read_text(encoding="utf-8").splitlines()
          if x.startswith("| T4: ")}
    check(rc == 0 and "| the survey runs on raw reads | the same" in v4.get("T4: the switch", "")
          and "| would pass | the same: K2" in v4.get("T4: K2 on Logan (reported)", ""),
          f"T4's switch and K2 read: {v4}")
    shutil.rmtree(w4)

    # a rerun gives the same files (but for the date)
    first = outputs(tmp)
    rc, err = rescore(tmp, ctrl)
    again = outputs(tmp)
    check(rc == 0 and first[R.OUT_TSV] == again[R.OUT_TSV], "a rerun gives the same both_ways.tsv")
    check([x for x in first[R.OUT_MD].splitlines() if not x.startswith("Made by")] ==
          [x for x in again[R.OUT_MD].splitlines() if not x.startswith("Made by")], "a rerun gives the same report")

    # tests not run yet
    part = tmp / "partial"
    for x in ("t1", "t2"):
        shutil.copytree(tmp / x, part / x)
    rc, err = rescore(part, ctrl)
    pmd = (part / "entry16" / R.OUT_MD).read_text(encoding="utf-8") if rc == 0 else ""
    check(rc == 0 and "## Not run yet" in pmd and all(f"- {x} (no" in pmd for x in ("T3", "T4", "T5"))
          and "| K1 (T1) |" in pmd and "K3 (T3)" not in pmd, f"T3 to T5 listed as not run: {err}")

    # refusals: nothing is written
    def refused(why, said, setup):
        w = tmp / ("r_" + str(len(fails)) + "_" + str(abs(hash(why))))
        shutil.copytree(tmp, w, ignore=shutil.ignore_patterns("entry16", "partial", "r_*"))
        setup(w)
        rc, err = rescore(w, w / "controls.tsv")
        check(rc == 1 and not (w / "entry16").exists() and "Nothing was written" in err and said in err,
              f"refused: {why} ({err.strip()})")

    refused("a table made with the survey's list", "but Phase 0's list gives",
            lambda w: edit(w / "t1" / T.OUT_CALLS, "CGIF01 (route)\troute\tother known genus\tunresolved",
                           "CGIF01 (route)\troute\tknown endosymbiont clade\tunresolved"))
    refused("a report that disagrees with its tables", "T1_ENTRY10.md says K1 FAILED",
            lambda w: edit(w / "t1" / T.OUT_MD, "## K1: passed", "## K1: FAILED"))
    refused("T4.md's K1 against its tables", "T4.md says K1 would fail",
            lambda w: edit(w / "t4" / "T4.md", "K1: would pass.", "K1: would fail: x."))
    refused("a compare.tsv that does not match", "does not match the calls tables",
            lambda w: edit(w / "t4" / "compare.tsv", "\tunresolved\tunresolved\t\t\tyes",
                           "\tresident\tunresolved\t\t\tno"))

    def twice(w):
        p = w / "t5" / "calls.tsv"
        lines = p.read_text(encoding="utf-8").splitlines(keepends=True)
        p.write_text("".join(lines + [lines[1]]), encoding="utf-8")
    refused("a call twice", "is in the table twice", twice)
    refused("a line in two T3 sets", "in two sets",
            lambda w: edit(w / "t3" / "runs.tsv", "\tottawa\t", "\thyphosphere\t"))
    refused("an unknown arm", "an arm that T2 does not make",
            lambda w: edit(w / "t2" / "calls.tsv", "\nheld_out\t", "\nheld_in\t"))
    refused("a missing column", "has no column category",
            lambda w: edit(w / "t4" / "calls.tsv", "\tcategory\t", "\tkind\t"))
    refused("T4 without T1", "T4 is compared with them", lambda w: shutil.rmtree(w / "t1"))
    refused("no tables at all", "no test has a calls table yet",
            lambda w: [shutil.rmtree(w / x) for x in ("t1", "t2", "t3", "t4", "t5")])
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if fails:
    print("FAIL")
    for f in fails:
        print("  -", f)
    sys.exit(1)
print("PASS: check_entry16 (the two lists, the call with the survey's list against call(), and rescore_entry16.py "
      "on synthetic T1 to T5 tables)")
