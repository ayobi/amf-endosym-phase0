"""Check tag_runs.py on mock tables in all three input formats, with and without overrides."""
import subprocess
import sys
from pathlib import Path

import pandas as pd

here = Path(__file__).parent
script = here.parent / "scripts" / "tag_runs.py"
inputs = [str(here / "mock_SraRunTable.csv"), str(here / "mock_runinfo.csv"),
          str(here / "mock_ena.tsv")]
fails = []


def run(outdir, *extra, ok=True):
    cmd = [sys.executable, str(script), *inputs, "-o", str(outdir), *extra]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if ok and p.returncode != 0:
        sys.exit(f"tag_runs.py failed:\n{p.stderr}")
    return p


def load(outdir):
    return pd.read_csv(outdir / "runs_tagged.tsv", sep="\t", dtype=str,
                       keep_default_na=False).set_index("run")


def check(t, expect, contains, label):
    for r, fields in expect.items():
        for k, v in fields.items():
            if t.loc[r, k] != v:
                fails.append(f"[{label}] {r}.{k}: got {t.loc[r, k]!r}, expected {v!r}")
    for r, fields in contains.items():
        for k, words in fields.items():
            got = t.loc[r, k]
            for w in words:
                if w not in got.split(";") and w not in got:
                    fails.append(f"[{label}] {r}.{k}: {w!r} missing from {got!r}")


# ------------------------------------------------------------ plain run ---
out1 = here / "out_plain"
print(run(out1).stdout)
t = load(out1)
assert len(t) == 34, f"expected 34 unique runs after de-duplication, got {len(t)}"
assert "Run" not in t.index, "repeated EDirect header line was kept as a run"
if t.loc["SRR0000001", "found_in"] != "mock_SraRunTable.csv;mock_runinfo.csv":
    fails.append(f"SRR0000001.found_in: got {t.loc['SRR0000001', 'found_in']!r}")
summary = (out1 / "inventory_summary.md").read_text()
if "not in mock_runinfo.csv: ERR0000001" not in summary:
    fails.append("summary does not list runs missing from mock_runinfo.csv")
if "## Review flags" not in summary or "PRJNA999007" not in summary.split("## Review flags")[1]:
    fails.append("summary does not list the flagged BioProject PRJNA999007")

expect = {
    "SRR0000001": dict(data_class="wgs_short", screen_value="high", isolate_hint="BEG34",
                       pair_hint="B+", role="screen"),
    "SRR0000002": dict(data_class="rna", rna_prep="polyA", screen_value="low", isolate_hint="BEG34"),
    "SRR0000003": dict(data_class="single_nucleus_wga", screen_value="uninformative",
                       isolate_hint="CL356"),
    "SRR0000004": dict(data_class="amplicon", amplicon_target="fungal_marker",
                       screen_value="not_applicable"),
    "SRR0000005": dict(data_class="amplicon", amplicon_target="bacterial_16S",
                       screen_value="relevant_16S"),
    "SRR0000006": dict(data_class="rad_reduced", screen_value="low", isolate_hint="A2"),
    "SRR0000007": dict(data_class="wgs_long", screen_value="high", isolate_hint="DAOM197198"),
    "SRR0000008": dict(data_class="hic", screen_value="medium"),
    "SRR0000009": dict(data_class="rna", rna_prep="rRNA_depleted", screen_value="medium"),
    "SRR0000010": dict(data_class="rna", rna_prep="small_rna", screen_value="very_low"),
    # metatranscriptome with no prep stated: medium, not check_prep
    "SRR0000011": dict(data_class="rna", rna_prep="unknown", screen_value="medium",
                       organism="fungus metagenome"),
    # Hi-C cut with a restriction enzyme must not become RAD
    "SRR0000012": dict(data_class="hic", screen_value="medium", pair_hint="Bplus"),
    "SRR0000013": dict(data_class="hic", screen_value="medium", pair_hint="Bminus"),
    # RNA-seq submitted as "Reduced Representation" stays RNA
    "SRR0000014": dict(data_class="rna", rna_prep="unknown", screen_value="check_prep"),
    # OTHER strategy with cDNA selection is RNA
    "SRR0000015": dict(data_class="rna", screen_value="check_prep"),
    "SRR0000016": dict(data_class="rna", rna_prep="polyA", screen_value="low"),
    # amplicon with no marker named: check, not not_applicable
    "SRR0000017": dict(data_class="amplicon", amplicon_target="unknown", screen_value="check"),
    "SRR0000018": dict(data_class="amplicon", amplicon_target="fungal_marker",
                       screen_value="not_applicable"),
    # marker not named, but filed under an AMF taxon: most likely fungal
    "SRR0000031": dict(data_class="amplicon", amplicon_target="unknown",
                       screen_value="not_applicable"),
    "SRR0000019": dict(data_class="wgs_short", screen_value="high"),
    "SRR0000029": dict(data_class="wgs_short", screen_value="medium", role="screen"),
    # filed under the endobacterium's name
    "SRR0000030": dict(data_class="dna_metagenome", screen_value="high", isolate_hint="LPA8"),
    "ERR0000001": dict(data_class="dna_metagenome", screen_value="high", isolate_hint="IT104",
                       platform_group="long"),
    "ERR0000002": dict(data_class="rna", screen_value="check_prep"),
    "ERR0000003": dict(data_class="rna", screen_value="check_prep"),
}
contains = {
    "SRR0000002": dict(pair_hint=["B-", "cured"], culture_hint=["spores", "germinating"]),
    "SRR0000004": dict(culture_hint=["in_planta", "field"]),
    "SRR0000006": dict(culture_hint=["in_vitro_root_culture"]),
    "SRR0000007": dict(culture_hint=["in_vitro_root_culture"]),
    "SRR0000009": dict(culture_hint=["in_planta"]),
    "SRR0000011": dict(culture_hint=["in_planta"]),
    "SRR0000019": dict(review_flag=["10_libraries_one_biosample"]),
    "SRR0000028": dict(review_flag=["10_libraries_one_biosample"]),
    "SRR0000029": dict(review_flag=["human"]),
    "SRR0000030": dict(culture_hint=["in_vitro_root_culture"]),
    "ERR0000002": dict(pair_hint=["WT", "wild-type"]),
    "ERR0000003": dict(pair_hint=["cured"]),
}
check(t, expect, contains, "plain")
for r in ["SRR0000001", "SRR0000012", "SRR0000030"]:
    if t.loc[r, "review_flag"]:
        fails.append(f"[plain] {r} should have no review flag, got {t.loc[r, 'review_flag']!r}")

# -------------------------------------------------------- with overrides ---
out2 = here / "out_overrides"
p = run(out2, "--overrides", str(here / "mock_overrides.tsv"))
print(p.stdout)
if "PRJNA000000" not in p.stderr:
    fails.append("no warning for the override that matches no run")
t2 = load(out2)
expect2 = {
    "SRR0000017": dict(amplicon_target="bacterial_16S", screen_value="relevant_16S",
                       override_note="spore 16S amplicons (test)"),
    "SRR0000019": dict(data_class="single_nucleus_wga", screen_value="uninformative",
                       override_note="one BioSample, one library per nucleus (test)"),
    "SRR0000028": dict(data_class="single_nucleus_wga", screen_value="uninformative"),
    "SRR0000029": dict(role="excluded", data_class="wgs_short"),
    "SRR0000016": dict(screen_value="check", reason="set by override", rna_prep="polyA"),
    "SRR0000018": dict(amplicon_target="fungal_marker", override_note=""),  # only_class=amplicon
                                                                            # but other BioProject
    "SRR0000001": dict(role="screen", override_note=""),
}
check(t2, expect2, {"SRR0000019": dict(review_flag=["10_libraries_one_biosample"])}, "overrides")
s2 = (out2 / "inventory_summary.md").read_text()
for needle in ["## Roles", "## Overrides applied", "matched no run", "PRJNA000000"]:
    if needle not in s2:
        fails.append(f"[overrides] summary is missing {needle!r}")

# ------------------------------------------------------ bad overrides file ---
bad = here / "out_bad_overrides.tsv"
bad.write_text("accession\tonly_class\tfield\tvalue\tnote\n"
               "PRJNA999008\t-\tcolour\tblue\tnot a field\n"
               "PRJNA999007\t-\tdata_class\tsingle_nucleus\ttypo in value\n"
               "PRJNA999009\t-\trole\texcluded\t\n")
p = run(here / "out_bad", "--overrides", str(bad), ok=False)
if p.returncode == 0:
    fails.append("a bad overrides file was accepted")
else:
    for needle in ["colour", "single_nucleus", "note is empty"]:
        if needle not in p.stderr:
            fails.append(f"[bad overrides] error message does not mention {needle!r}")

if fails:
    print("FAIL\n" + "\n".join(fails))
    sys.exit(1)
print(f"PASS: {len(expect)} + {len(expect2)} label checks, flags, overrides and bad-input checks")
