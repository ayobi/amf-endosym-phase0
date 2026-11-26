#!/usr/bin/env python3
"""Phase 0's calls with both lists of known endosymbiont clades (PROTOCOL.md section 10, entry 16).

Phase 0's list is call_lineages.py's KNOWN_CLADE_GENERA, as locked. The survey's list, SURVEY_CLADE_GENERA,
adds GTDB's genus CGIF01, the MRE group, as CaMg/MRE. It changes only the category of a route lineage in
CGIF01, from "other known genus" to "known endosymbiont clade". In the order of PROTOCOL.md 4.1 that changes
one step only: such a lineage present in 2 or more runs of its line is resident. So an unresolved call with 2
or more present runs becomes resident, a resident call (by E1 and E2) stays resident on the new basis, and
every other call stays as it is. tests/check_entry16.py checks this against call_lineages.call() itself.
Phase 0's verdicts are the ones made with Phase 0's list; this report changes none of them.

Reads only the committed tables of the tests that have run:
  T1  t1/calls_entry10.tsv. K1 both ways (t1/levels_entry10.tsv, t1/runs.tsv), and checked against
      t1/T1_ENTRY10.md
  T2  t2/calls.tsv, both arms. K2 reads the same both ways: it fails only on absent, no call or no lineage,
      which the list never changes
  T3  t3/calls.tsv, the three sets without spikes (sets from t3/runs.tsv). K3 for the outside sets both ways.
      Spikes are not called again; listed are the Ottawa lines where a spike could make a CGIF01 call
      resident, and how many spikes changed a CGIF01 call (t3/spikes.tsv, t3/lab_spike.tsv)
  T4  t4/calls.tsv. K1 on Logan both ways (t4/levels.tsv), checked against t4/T4.md; the calls compared with
      raw reads (t4/compare.tsv), both ways. The switch uses levels only, so it reads the same
  T5  t5/calls.tsv. T5's reading for Mycoavidus both ways (t5/runs.tsv, controls.tsv)

Writes entry16/both_ways.tsv (every call with signal, with both lists) and entry16/ENTRY16.md; a test that has
not run is listed as not run. Commit both once T5 is done. Nothing else is read or written, so this can run
at any time, also while a test runs.

Usage:
  python scripts/rescore_entry16.py
"""

import argparse
import csv
import hashlib
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import call_lineages as C  # noqa: E402  both lists
import run_stage0a_m as M2  # noqa: E402  tables
import run_t1 as T  # noqa: E402  K1, runs.tsv, the calls table
import run_t3 as T3  # noqa: E402  the sets of T3
import run_t4 as T4  # noqa: E402  the arms of T4, compare.tsv
import run_t5 as T5  # noqa: E402  the pair of T5

__version__ = "0.1"

CALLS = ("no call", "absent", "contaminant", "resident", "unresolved")
OUT_TSV, OUT_MD = "both_ways.tsv", "ENTRY16.md"
KEPT = ["unit", "scored", "lineage", "source", "members", "present_runs", "low_runs", "informative_runs", "category",
        "call", "basis"]                                     # from the calls tables, as they are
COLS = ["test", "part"] + KEPT + ["category_survey", "call_survey", "basis_survey", "changed"]
OUTSIDE = ("hyphosphere", "added_bacterium")
SAME_K2 = "the same: K2 fails only on absent, no call or no lineage, which the list never changes"
SWITCH = r".*\*\*(the survey (?:runs on raw reads|can run on Logan unitigs))\*\*"   # in t4/T4.md


# ------------------------------------------------------------------- calls ---

def survey_call(row):
    """(category, call, basis) with the survey's list, for one row of a calls table made with Phase 0's list."""
    cat, members = row["category"], [g for g in row["members"].split(", ") if g]
    if row["call"] not in CALLS:
        raise ValueError(f"{row['unit']}, {row['lineage']}: unknown call {row['call']!r}")
    if row["source"] == "mapped":
        old, new = C.category_of_mapped(row["lineage"]), None
    elif members:
        old, new = C.category_of_route(members), C.category_of_route(members, C.SURVEY_CLADE_GENERA)
    else:
        old, new = "novel", None
    if cat != old:
        raise ValueError(f"{row['unit']}, {row['lineage']}: the table gives {cat}, but Phase 0's list gives {old}")
    if new is None or new == cat:
        return cat, row["call"], row["basis"]
    if (cat, new) != ("other known genus", "known endosymbiont clade"):
        raise ValueError(f"{row['lineage']}: {cat} would become {new}, which this script does not handle")
    present = int(row["present_runs"])
    if row["call"] in ("unresolved", "resident") and present >= 2:
        return new, "resident", f"known endosymbiont clade, present in {present} runs"
    return new, row["call"], row["basis"]


def rescored(test, part, rows):
    """The rows of one table part with both lists."""
    seen, out = set(), []
    for r in rows:
        key = (r["unit"], r["lineage"])
        if key in seen:
            raise ValueError(f"{test}, {part}: {r['unit']}, {r['lineage']} is in the table twice")
        seen.add(key)
        cat, call, basis = survey_call(r)
        changed = "call" if call != r["call"] else ("basis only" if basis != r["basis"] else "no")
        out.append({**{c: r[c] for c in KEPT}, "test": test, "part": part, "signal": r["signal"],
                    "category_survey": cat, "call_survey": call, "basis_survey": basis, "changed": changed})
    return out


def as_calls(rows, survey):
    return [{"unit": r["unit"], "lineage": r["lineage"], "source": r["source"], "members": r["members"],
             "call": r["call_survey"] if survey else r["call"]} for r in rows]


def survey_only(members):
    """True if a route lineage of these genus names is a known clade with the survey's list only."""
    return C.category_of_route(members) != C.category_of_route(members, C.SURVEY_CLADE_GENERA)


# ------------------------------------------------------------------ tables ---

def note(path, inputs):
    """Records a file read and its md5; returns the path."""
    inputs.append((str(path), hashlib.md5(Path(path).read_bytes()).hexdigest()))
    return path


def table(path, cols, inputs):
    """The rows of a committed table, which must have these columns."""
    with open(note(path, inputs), newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in cols if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path} has no column {', '.join(missing)}")
        return list(reader)


def heading(path, pattern, inputs):
    """The first group of the first line of a committed report that matches, or None (no report or no line)."""
    if not path.exists():
        return None
    for line in note(path, inputs).read_text(encoding="utf-8").splitlines():
        hit = re.match(pattern, line)
        if hit:
            return hit.group(1)
    return None


def k1_text(reasons):
    return "passed" if not reasons else "FAILED: " + "; ".join(reasons)


def words(xs):
    return ", ".join(xs[:-1]) + " or " + xs[-1] if len(xs) > 1 else "".join(xs)


def by_key(rows):
    return {(r["unit"], r["lineage"]): r for r in rows}


# ------------------------------------------------------------------- tests ---

def rescore(d, controls_path):
    """Every test that has run, with both lists. Returns a dict for write_outputs."""
    inputs, parts, missing, verdicts, raw = [], [], [], [], {}

    def add(test, part, rows):
        out = rescored(test, part, rows)
        parts.append((test, part, out))
        return out

    # T1 (entry 10)
    if (d["t1"] / "calls_entry10.tsv").exists():
        rows = add("T1", "all runs", table(d["t1"] / "calls_entry10.tsv", T.CALL_COLS, inputs))
        raw[("T1", "full")] = rows
        levels = table(d["t1"] / "levels_entry10.tsv", ["run", "lineage", "level"], inputs)
        runs = T.read_runs(note(d["t1"] / "runs.tsv", inputs))
        old, new = T.k1(as_calls(rows, False), levels, runs), T.k1(as_calls(rows, True), levels, runs)
        said = heading(d["t1"] / "T1_ENTRY10.md", r"^## K1: (passed|FAILED)$", inputs)
        if said is not None and said != ("passed" if not old else "FAILED"):
            raise ValueError(f"t1/T1_ENTRY10.md says K1 {said}, but its tables give {k1_text(old)}")
        verdicts.append(("K1 (T1)", k1_text(old), k1_text(new)))
    else:
        missing.append(("T1", f"{d['t1']}/calls_entry10.tsv"))

    # T2
    if (d["t2"] / "calls.tsv").exists():
        rows = table(d["t2"] / "calls.tsv", ["arm"] + T.CALL_COLS, inputs)
        other = sorted({r["arm"] for r in rows} - {"held_out", "kept_in"})
        if other:
            raise ValueError(f"t2/calls.tsv has an arm that T2 does not make: {', '.join(other)}")
        for arm in ("held_out", "kept_in"):
            raw[("T2", arm)] = add("T2", arm, [r for r in rows if r["arm"] == arm])
        said = heading(d["t2"] / "T2.md", r"^## K2: (passed|FAILED)$", inputs)
        verdicts.append(("K2 (T2)", said or "see t2/T2.md", SAME_K2))
    else:
        missing.append(("T2", f"{d['t2']}/calls.tsv"))

    # T3, without spikes
    if (d["t3"] / "calls.tsv").exists():
        rows = table(d["t3"] / "calls.tsv", T.CALL_COLS, inputs)
        set_of = {}
        for r in table(d["t3"] / "runs.tsv", ["run", "unit", "set"], inputs):
            if set_of.setdefault(r["unit"], r["set"]) != r["set"]:
                raise ValueError(f"t3/runs.tsv puts {r['unit']} in two sets")
        lost = sorted({r["unit"] for r in rows if r["unit"] not in set_of})
        if lost:
            raise ValueError(f"t3/calls.tsv has lines that t3/runs.tsv does not: {', '.join(lost)}")
        out3 = {s: add("T3", s, [r for r in rows if set_of[r["unit"]] == s]) for s in T3.SETS}
        res = {k: [f"{r['lineage']} in {r['unit']}" for s in OUTSIDE for r in out3[s] if r[k] == "resident"]
               for k in ("call", "call_survey")}
        verdicts.append(("K3 (T3), outside sets", *("passed" if not res[k] else "FAILED: resident: " + "; ".join(res[k])
                                                    for k in ("call", "call_survey"))))
        members = {r["lineage"]: [g for g in r["members"].split(", ") if g] for r in out3["ottawa"]}
        at_risk = sorted({r["unit"] for r in out3["ottawa"] if r["category_survey"] != r["category"]
                          and r["signal"] == "yes" and r["call_survey"] != "resident"})
        touched = 0
        for name in ("spikes.tsv", "lab_spike.tsv"):
            if (d["t3"] / name).exists():
                for x in table(d["t3"] / name, ["changed_calls"], inputs):
                    names = [y.rsplit(": ", 1)[0] for y in x["changed_calls"].split("; ") if y]
                    touched += any(n.endswith(" (route)") and survey_only(members.get(n, [n[:-len(" (route)")]]))
                                   for n in names)
        said = heading(d["t3"] / "T3.md", r"^## K3: (passed|FAILED)$", inputs)
        verdicts.append(("K3 (T3), spikes and the lab kit spike", f"{said or 'see t3/T3.md'} (t3/T3.md, all of K3)",
                         "not called again. Ottawa lines where a spike could make a CGIF01 call resident: "
                         + (", ".join(at_risk) or "none") + f"; spikes that changed a CGIF01 call: {touched}"))
    else:
        missing.append(("T3", f"{d['t3']}/calls.tsv"))

    # T4, on Logan
    if (d["t4"] / "calls.tsv").exists():
        if ("T1", "full") not in raw or ("T2", "held_out") not in raw:
            raise ValueError("t4/ is there, but T1 (entry 10) or T2 is not: T4 is compared with them")
        rows = table(d["t4"] / "calls.tsv", ["test", "arm"] + T.CALL_COLS, inputs)
        other = sorted({(r["test"], r["arm"]) for r in rows} - set(T4.ARMS))
        if other:
            raise ValueError(f"t4/calls.tsv has a part that T4 does not make: {other}")
        logan = {a: add("T4", f"{a[0]} {a[1]}", [r for r in rows if (r["test"], r["arm"]) == a])
                 for a in T4.ARMS}
        levels = [x for x in table(d["t4"] / "levels.tsv", ["test", "arm", "run", "lineage", "level"], inputs)
                  if (x["test"], x["arm"]) == ("T1", "full")]
        if levels:
            held = {x["run"] for x in levels}
            mine = [r for r in T.read_runs(d["t1"] / "runs.tsv") if r["run"] in held]
            k1 = [k1_text(T.k1(as_calls(logan[("T1", "full")], s), levels, mine)) for s in (False, True)]
        else:
            k1 = ["FAILED: no T1 run is in Logan"] * 2
        said = heading(d["t4"] / "T4.md", r"^K1: (would pass|would fail)", inputs)
        if said is not None and (said == "would pass") != (k1[0] == "passed"):
            raise ValueError(f"t4/T4.md says K1 {said}, but its tables give {k1[0]}")
        verdicts.append(("T4: K1 on Logan (reported)", *[x.replace("passed", "would pass", 1) if x == "passed"
                                                         else x.replace("FAILED", "would fail", 1) for x in k1]))
        said2 = heading(d["t4"] / "T4.md", r"^K2: (would pass|would fail)", inputs)
        verdicts.append(("T4: K2 on Logan (reported)", said2 or "see t4/T4.md", SAME_K2))
        same, total = [0, 0], 0
        raw_by = {a: by_key(rows) for a, rows in raw.items()}
        logan_by = {a: by_key(rows) for a, rows in logan.items()}
        for x in table(d["t4"] / "compare.tsv", T4.COMPARE_COLS, inputs):
            if x["what"] != "call" or x["same"] == "n/a":
                continue
            a, key = (x["test"], x["arm"]), (x["run_or_unit"], x["lineage"])
            ra, rb = raw_by[a].get(key), logan_by[a].get(key)
            calls = [(ra["call"] if ra else "none", rb["call"] if rb else "none"),
                     (ra["call_survey"] if ra else "none", rb["call_survey"] if rb else "none")]
            if calls[0] != (x["raw"], x["logan"]) or (x["same"] == "yes") != (x["raw"] == x["logan"]):
                raise ValueError(f"t4/compare.tsv does not match the calls tables at {a}, {key}")
            total += 1
            same = [n + (c[0] == c[1]) for n, c in zip(same, calls)]
        verdicts.append(("T4: calls the same as on raw reads", f"{same[0]} of {total}", f"{same[1]} of {total}"))
        said3 = heading(d["t4"] / "T4.md", SWITCH, inputs)
        verdicts.append(("T4: the switch", said3 or "see t4/T4.md", "the same: it uses levels only"))
    else:
        missing.append(("T4", f"{d['t4']}/calls.tsv"))

    # T5
    if (d["t5"] / "calls.tsv").exists():
        out5 = by_key(add("T5", "all runs", table(d["t5"] / "calls.tsv", T.CALL_COLS, inputs)))
        controls = {x["control_id"]: x["expected_call"] for x in M2.read_tsv(note(controls_path, inputs))}
        wt, cured, _ = T5.pair(T.read_runs(note(d["t5"] / "runs.tsv", inputs)), controls)
        reading = []
        for k in ("call", "call_survey"):
            a, b = out5.get((wt, T5.SYMBIONT), {}).get(k), out5.get((cured, T5.SYMBIONT), {}).get(k)
            reading.append(("as expected" if a == "resident" and b != "resident" else "not as expected")
                           + f" (wild type {a or 'no call made'}, cured line {b or 'no call made'})")
        verdicts.append((f"T5 (reported): {T5.SYMBIONT}", *reading))
    else:
        missing.append(("T5", f"{d['t5']}/calls.tsv"))

    if not parts:
        raise ValueError("no test has a calls table yet")
    return {"inputs": inputs, "parts": parts, "missing": missing, "verdicts": verdicts}


# ----------------------------------------------------------------- outputs ---

def write_outputs(out_dir, res):
    rows = [r for _, _, p in res["parts"] for r in p if r["signal"] == "yes"]
    tsv = "\t".join(COLS) + "\n" + "".join("\t".join(r[c] for c in COLS) + "\n" for r in rows)
    md = ["# Entry 16: Phase 0's calls with both lists of known clades", "",
          f"Made by rescore_entry16.py v{__version__} (call_lineages.py v{C.__version__}) on {T.now()}, from the "
          "committed tables listed at the end.",
          "Phase 0's list (`KNOWN_CLADE_GENERA`, as locked): a route lineage in " + words(sorted(C.KNOWN_CLADE_GENERA))
          + " is a known endosymbiont clade. The survey's list "
          "(`SURVEY_CLADE_GENERA`, entry 16) adds " +
          ", ".join(f"{g} as {c}" for g, c in sorted(C.SURVEY_CLADE_GENERA.items())
                    if g not in C.KNOWN_CLADE_GENERA) +
          ". Phase 0's verdicts are the ones made with Phase 0's list; this report changes none of them.", "",
          "## What changes", "",
          "| Test | Part | Calls with signal | Call changes | Basis only | Now resident |", "|---|---|---|---|---|---|"]
    for test, part, p in res["parts"]:
        sig = [r for r in p if r["signal"] == "yes"]
        now_res = [f"{r['unit']}: {r['lineage']}" for r in sig
                   if r["changed"] == "call" and r["call_survey"] == "resident"]
        md.append(f"| {test} | {part} | {len(sig)} | {sum(r['changed'] == 'call' for r in sig)} | "
                  f"{sum(r['changed'] == 'basis only' for r in sig)} | {'; '.join(now_res) or '-'} |")
    md += ["", "## Verdicts", "", "| Check | Phase 0's list | The survey's list |", "|---|---|---|"]
    md += [f"| {a} | {b} | {c} |" for a, b, c in res["verdicts"]]
    changed = [r for r in rows if r["changed"] != "no"]
    md += ["", "## Calls that change", ""]
    if changed:
        md += ["| Test | Part | Line | Lineage | Present / low / informative runs | Phase 0's list | "
               "The survey's list | Basis with the survey's list |", "|---|---|---|---|---|---|---|---|"]
        md += [f"| {r['test']} | {r['part']} | {r['unit']} | {r['lineage']} | {r['present_runs']} / {r['low_runs']} / "
               f"{r['informative_runs']} | {r['call']} | {r['call_survey']} | {r['basis_survey']} |" for r in changed]
    else:
        md.append("None.")
    if res["missing"]:
        md += ["", "## Not run yet", ""] + [f"- {t} (no `{p}`)" for t, p in res["missing"]]
    md += ["", "## Inputs", "", "| File | md5 |", "|---|---|"]
    md += [f"| `{p}` | {h} |" for p, h in dict.fromkeys(res["inputs"])]
    md += ["", f"`{OUT_TSV}` gives every call with signal with both lists (columns ending in _survey: the survey's "
               "list). A route lineage in CGIF01 changes category; in the order of PROTOCOL.md 4.1 only the step "
               "\"known endosymbiont clade, present in 2 or more runs\" can then change its call."]
    out_dir.mkdir(parents=True, exist_ok=True)
    texts = {OUT_TSV: tsv, OUT_MD: "\n".join(md) + "\n"}
    for name, text in texts.items():
        (out_dir / (name + ".part")).write_text(text, encoding="utf-8")
    for name in texts:
        os.replace(out_dir / (name + ".part"), out_dir / name)
    return rows, changed


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for t in ("t1", "t2", "t3", "t4", "t5"):
        ap.add_argument(f"--{t}-dir", default=t)
    ap.add_argument("--controls", default="controls.tsv")
    ap.add_argument("--out-dir", default="entry16", help="where both_ways.tsv and ENTRY16.md go [entry16]")
    args = ap.parse_args(argv)
    d = {t: Path(getattr(args, f"{t}_dir")) for t in ("t1", "t2", "t3", "t4", "t5")}
    try:
        res = rescore(d, Path(args.controls))
        rows, changed = write_outputs(Path(args.out_dir), res)
    except (OSError, KeyError, ValueError, RuntimeError) as e:
        print(f"Nothing was written to {args.out_dir}/: {e}", file=sys.stderr)
        return 1
    print(f"{len(rows)} calls with signal, {sum(r['changed'] == 'call' for r in changed)} change with the survey's "
          f"list" + "".join(f"\n  {a}: {b} | {c}" for a, b, c in res["verdicts"])
          + "".join(f"\n  {t}: not run yet" for t, _ in res["missing"]))
    print(f"-> {args.out_dir}/{OUT_TSV}, {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
