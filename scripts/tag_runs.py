#!/usr/bin/env python3
"""Label public AMF sequencing runs for endosymbiont screening.

Reads one or more run tables, in any of these formats (they can be mixed):
  - NCBI SRA Run Selector "Metadata" download (SraRunTable.csv)
  - EDirect runinfo CSV (efetch -db sra -format runinfo)
  - ENA Portal API TSV (result=read_run)

Merges them (one row per run) and writes to --outdir:
  runs_tagged.tsv        key metadata plus the labels below
  inventory_summary.md   counts, review flags, and the BioProjects most worth screening

Labels added to each run:
  data_class       wgs_short, wgs_long, dna_metagenome, rna, rad_reduced, hic,
                   amplicon, single_nucleus_wga, other
  platform_group   short, long, other
  rna_prep         polyA, rRNA_depleted, small_rna, unknown  (blank if not RNA)
  amplicon_target  bacterial_16S, fungal_marker, unknown     (blank if not amplicon)
  screen_value     high, medium, low, very_low, relevant_16S, check_prep,
                   check, not_applicable, uninformative
  reason           one line explaining screen_value
  role             screen (default), control_only, excluded  (set in overrides)
  review_flag      reasons to look at a run by hand (many libraries on one
                   BioSample, odd organism field)
  override_note    notes from the overrides file, if any row touched this run
  culture_hint     in_vitro_root_culture, in_planta, spores, germinating, field
  isolate_hint     isolate or strain name, if one is stated or recognisable
  pair_hint        words that suggest a with/without-endosymbiont design
                   (B+, B-, cured, WT ...)

These labels come from metadata only. Metadata is often incomplete, so treat
'check_prep' and 'check' as "read the paper", not as "bad data".

Overrides (optional, --overrides FILE.tsv) fix labels by hand. Columns:
  accession   a BioProject, BioSample, Experiment or Run accession
  only_class  apply only to runs that currently have this data_class
              ("-" or blank: any class)
  field       data_class, rna_prep, amplicon_target, screen_value or role
  value       the new value
  note        why (required; it is copied into override_note)
Rows apply in file order. Changing data_class, rna_prep or amplicon_target
re-derives screen_value, unless screen_value is also set for that run.

Usage:
  python tag_runs.py SraRunTable_A.csv SraRunTable_B.csv SraRunTable_C.csv \\
      --overrides overrides.tsv -o inventory
"""

import argparse
import csv
import datetime as dt
import math
import re
import sys
from pathlib import Path

import pandas as pd

__version__ = "0.3"

# ---------------------------------------------------------------- columns ---

# Canonical field -> column names used by Run Selector, EDirect runinfo, ENA.
# Matching ignores case, spaces and underscores.
CANDIDATES = {
    "run": ["Run", "run_accession"],
    "bioproject": ["BioProject", "study_accession"],
    "biosample": ["BioSample", "sample_accession"],
    "experiment": ["Experiment", "experiment_accession"],
    "organism": ["Organism", "ScientificName", "scientific_name"],
    "strategy": ["Assay Type", "LibraryStrategy", "library_strategy"],
    "source": ["LibrarySource", "library_source"],
    "selection": ["LibrarySelection", "library_selection"],
    "layout": ["LibraryLayout", "library_layout"],
    "platform": ["Platform", "instrument_platform"],
    "instrument": ["Instrument", "Model", "instrument_model"],
    "bases": ["Bases", "bases", "base_count"],
    "library_name": ["Library Name", "LibraryName", "library_name"],
    "sample_name": ["Sample Name", "SampleName", "sample_title", "sample_alias"],
    "center": ["Center Name", "CenterName", "center_name"],
    "release": ["ReleaseDate", "first_public"],
}
ISOLATE_COLS = ["isolate", "strain", "fungal_strain", "culture_collection"]
EMPTY_VALUES = {"", "missing", "not applicable", "not collected", "na", "n/a",
                "none", "unknown", "not provided", "nan", "unspecified"}

LONG_PLATFORMS = {"PACBIO_SMRT", "OXFORD_NANOPORE"}
SHORT_PLATFORMS = {"ILLUMINA", "BGISEQ", "DNBSEQ", "ION_TORRENT", "LS454",
                   "ELEMENT", "ULTIMA", "COMPLETE_GENOMICS"}
RNA_STRATEGIES = {"RNA-SEQ", "MIRNA-SEQ", "NCRNA-SEQ", "FL-CDNA", "EST",
                  "SSRNA-SEQ", "RIBO-SEQ"}
WGS_STRATEGIES = {"WGS", "WCS", "SYNTHETIC-LONG-READ", "POOLCLONE", "CLONE"}
# LibrarySelection values that only make sense for RNA. ("size fractionation"
# is left out on purpose: long-read DNA libraries use it too.)
RNA_SELECTIONS = {"cdna", "cdna_randompriming", "cdna_oligo_dt", "polya",
                  "oligo-dt", "inverse rrna", "inverse rrna selection", "rt-pcr"}
RR_SELECTIONS = {"restriction digest", "reduced representation"}

DATA_CLASSES = ["wgs_short", "wgs_long", "dna_metagenome", "rna", "rad_reduced", "hic",
                "amplicon", "single_nucleus_wga", "other"]
SCREEN_VALUES = ["high", "medium", "low", "very_low", "relevant_16S", "check_prep",
                 "check", "not_applicable", "uninformative"]
OVERRIDE_VALUES = {
    "data_class": set(DATA_CLASSES),
    "rna_prep": {"polyA", "rRNA_depleted", "small_rna", "unknown"},
    "amplicon_target": {"bacterial_16S", "fungal_marker", "unknown"},
    "screen_value": set(SCREEN_VALUES),
    "role": {"screen", "control_only", "excluded"},
}
MANY_LIBRARIES = 10   # short-read DNA libraries on one BioSample before we flag it

# ---------------------------------------------------------- text patterns ---

I = re.IGNORECASE
# Word edges that also treat "_" as a separator, so library names such as
# HiC_Bplus or ddRAD_A2 still match.
L = r"(?<![A-Za-z0-9])"
R = r"(?![A-Za-z0-9])"

SINGLE_NUC = re.compile(L + r"(single[ _-]?nucle(us|i)|single[ _-]?cell|snuc|"
                        r"whole[ _-]genome[ _-]amplif)", I)
RAD = re.compile(L + r"(dd[ _-]?rad|rad[ _-]?seq|restriction[ _-]site[ _-]associated)" + R, I)
HIC = re.compile(L + r"(hi-?c|omni-?c|chromatin[ _-]conformation|proximity[ _-]ligation)" + R, I)
SMALL_RNA = re.compile(L + r"(small[ _-]?rna|srna|mirna|micro[ _-]?rna)" + R, I)
POLYA = re.compile(r"poly-?\(?a\)?[ +_-]*(select|enrich|captur|pull)|"
                   r"oligo[ _-]?d\(?t\)?|truseq[ _-]stranded[ _-]mrna|mrna[ _-]?seq|"
                   r"mrna[ _-]librar", I)
DEPLETED = re.compile(r"rrna[ _-]?deplet|ribo-?zero|ribominus|ribo[ _-]?depl|"
                      r"total[ _-]rna[ _-]librar", I)
# Amplicon targets. Words like "AMF" or "arbuscular" are left out on purpose:
# every run here comes from an AMF study, so they say nothing about the marker.
AMP_16S = re.compile(L + r"(16s|v3[ _-]?v4|v4|515f|806r|341f|785r|799f|1193r|27f|1492r|"
                     r"bacteri\w*|endobacteri\w*|prokaryot\w*)" + R, I)
AMP_FUNGAL = re.compile(L + r"((?-i:ITS)[12]?|18s|ssu|lsu|28s|nu-ssu|ns31|aml[12]|amv4\.5nf|"
                        r"amdgr|wanda|ssum[ac]f|lsum[ab]r|its1f|its3|its4|gits7|lr3|lr5)" + R, I)
ODD_ORGANISM = re.compile(r"homo sapiens|\bhuman|patient|\bgut\b|fec(?:al|es)|faec|stool|"
                          r"mus musculus|\bmouse\b", I)

CULTURE = {
    "in_vitro_root_culture": re.compile(
        r"daucus|carrot|hairy[ _-]?root|ri[ _-]?t-?dna|root[ _-]organ[ _-]cultur|"
        r"\broc\b|in[ _-]?vitro|monoxenic|chicory|cichorium", I),
    "in_planta": re.compile(
        r"medicago|lotus japonicus|\boryza\b|\brice\b|zea mays|\bmaize\b|sorghum|"
        r"tomato|solanum|brachypodium|\ballium\b|\bleek\b|plantago|trifolium|"
        r"clover|petunia|nicotiana|populus|poplar|colonized root|"
        r"mycorrhizal root|inoculated root|co-?culture", I),
    "spores": re.compile(r"\bspores?\b|sporocarp", I),
    "germinating": re.compile(r"germinat", I),
    "field": re.compile(r"\bfield\b|\bsoil\b|rhizosphere", I),
}

ISOLATE_CODE = re.compile(r"\b(DAOM[ _-]?\d+|BEG[ _-]?\d+|MUCL[ _-]?\d+|"
                          r"INVAM[ _-]?[A-Z0-9]+)\b", I)
ISOLATE_WORD = re.compile(r"\bisolate[:=\s]+([A-Za-z0-9][A-Za-z0-9._-]{0,20})", I)
# B+/B-/WT stay case-sensitive ((?-i:...)); the word forms do not.
PAIR = re.compile(L + r"(?-i:B)\+|" + L + r"(?-i:B)-" + R + "|"
                  + L + r"(?-i:B)[ _-]?(?:plus|minus)" + R +
                  r"|\bcured\b|endobacteri\w*[ _-]?free|"
                  r"without (?:the |its )?endobacteri\w*|"
                  r"with (?:the |its )?endobacteri\w*|" + L + r"(?-i:WT)" + R + r"|wild[ _-]?type", I)


# ---------------------------------------------------------------- reading ---

def _key(name):
    return re.sub(r"[\s_]", "", str(name)).lower()


def read_table(path):
    """Read a CSV or TSV table as strings, sniffing the separator."""
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        header = fh.readline()
    sep = "\t" if header.count("\t") > header.count(",") else ","
    df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False,
                     quoting=csv.QUOTE_MINIMAL, encoding_errors="replace")
    df = df.fillna("").astype(str)
    df.columns = [c.strip() for c in df.columns]
    return df


def normalise(df, source_file):
    """Map a raw table onto canonical columns; keep a text blob for patterns."""
    by_key = {}
    for col in df.columns:
        by_key.setdefault(_key(col), col)
    out = pd.DataFrame(index=df.index)
    for field, names in CANDIDATES.items():
        col = next((by_key[_key(n)] for n in names if _key(n) in by_key), None)
        out[field] = df[col].str.strip() if col is not None else ""

    iso_cols = [by_key[_key(c)] for c in ISOLATE_COLS if _key(c) in by_key]
    if iso_cols:
        def first_isolate(row):
            for c in iso_cols:
                v = str(row[c]).strip()
                if v.lower() not in EMPTY_VALUES:
                    return v
            return ""
        out["isolate_field"] = df.apply(first_isolate, axis=1)
    else:
        out["isolate_field"] = ""

    out["text"] = df.apply(lambda r: " | ".join(v for v in r.values if v), axis=1)
    out["input_file"] = Path(source_file).name
    if (out["run"] == "").all():
        sys.exit(f"{source_file}: no run accession column found "
                 f"(expected one of {CANDIDATES['run']}).")
    return out


def parse_bases(value):
    s = str(value).strip().replace(",", "")
    m = re.match(r"^([0-9]*\.?[0-9]+)\s*([KMGT]?)", s, I)
    if not m:
        return math.nan
    mult = {"": 1, "K": 1e3, "M": 1e6, "G": 1e9, "T": 1e12}[m.group(2).upper()]
    return float(m.group(1)) * mult


# --------------------------------------------------------------- labelling ---

def classify(r):
    """Data class from strategy, source, selection and text. Order matters."""
    strat = r["strategy"].upper()
    src = r["source"].upper()
    sel = r["selection"].lower().strip()
    text = r["text"]

    if strat == "WGA" or "SINGLE CELL" in src or sel == "mda" or SINGLE_NUC.search(text):
        return "single_nucleus_wga"
    if strat == "AMPLICON":
        return "amplicon"
    # Hi-C before RAD: Hi-C libraries are cut with restriction enzymes, so their
    # LibrarySelection is often "Restriction Digest".
    if (strat == "HI-C" or (strat in ("OTHER", "") and HIC.search(text))
            or (strat in ("WGS", "OTHER", "") and HIC.search(r["library_name"]))):
        return "hic"
    # RNA before RAD: some RNA-seq is submitted with selection "Reduced Representation".
    if ("TRANSCRIPTOMIC" in src or strat in RNA_STRATEGIES
            or (strat in ("OTHER", "") and sel in RNA_SELECTIONS)):
        return "rna"
    if (strat == "RAD-SEQ" or sel in RR_SELECTIONS
            or (strat in ("OTHER", "") and RAD.search(text))):
        return "rad_reduced"
    if "METAGENOMIC" in src:
        return "dna_metagenome"
    if strat in WGS_STRATEGIES or ("GENOMIC" in src and strat in ("OTHER", "")):
        return "wgs_long" if r["platform"].upper() in LONG_PLATFORMS else "wgs_short"
    return "other"


def rna_prep(r):
    sel = r["selection"].lower().strip()
    strat = r["strategy"].upper()
    text = r["text"]
    if strat in ("MIRNA-SEQ", "NCRNA-SEQ") or sel == "size fractionation" or SMALL_RNA.search(text):
        return "small_rna"
    if sel in ("polya", "oligo-dt", "cdna_oligo_dt") or POLYA.search(text):
        return "polyA"
    if sel.startswith("inverse rrna") or DEPLETED.search(text):
        return "rRNA_depleted"
    return "unknown"


def amplicon_target(r):
    if AMP_16S.search(r["text"]):
        return "bacterial_16S"
    if AMP_FUNGAL.search(r["text"]):
        return "fungal_marker"
    return "unknown"


def screen(r):
    c = r["data_class"]
    if c == "single_nucleus_wga":
        return "uninformative", "single nuclei / WGA: sorting should strip cytoplasmic bacteria"
    if c == "amplicon":
        t = r["amplicon_target"]
        if t == "bacterial_16S":
            return "relevant_16S", "bacterial marker amplicon: direct but PCR-biased"
        if t == "fungal_marker":
            return "not_applicable", "fungal marker amplicon: no bacterial signal expected"
        # Marker not named. Filed under an AMF taxon, it is almost always a fungal
        # marker; filed as a metagenome, it may well be 16S (spore microbiomes).
        if "metagenome" in r["organism"].lower():
            return "check", "amplicon filed as a metagenome, marker not stated: read the paper"
        return "not_applicable", "amplicon filed under an AMF name, marker not stated: most likely fungal"
    if c == "rad_reduced":
        return "low", "reduced representation: detection possible but sparse"
    if c == "hic":
        return "medium", "Hi-C: shotgun-like, and contacts may link symbiont to host"
    if c == "rna":
        prep = r["rna_prep"]
        if prep == "small_rna":
            return "very_low", "small RNA only"
        if prep == "polyA":
            return "low", "poly-A selection drops most bacterial mRNA"
        if prep == "rRNA_depleted":
            return "medium", "total RNA: bacterial transcripts kept unless removed"
        if "METATRANSCRIPTOMIC" in r["source"].upper():
            return "medium", "metatranscriptome, prep not stated: confirm before relying on it"
        return "check_prep", "RNA prep not stated in metadata: read the paper"
    if c in ("wgs_short", "wgs_long", "dna_metagenome"):
        b = r["bases_num"]
        if isinstance(b, float) and not math.isnan(b) and b >= 1e9:
            return "high", "shotgun DNA, >= 1 Gbp"
        if isinstance(b, float) and not math.isnan(b) and b > 0:
            return "medium", "shotgun DNA, < 1 Gbp"
        return "medium", "shotgun DNA, size not stated"
    return "check", "strategy not recognised"


def culture_hint(text):
    return ";".join(k for k, rx in CULTURE.items() if rx.search(text))


def isolate_hint(r):
    if r["isolate_field"]:
        return r["isolate_field"]
    m = ISOLATE_CODE.search(r["text"]) or ISOLATE_WORD.search(r["text"])
    return m.group(1) if m else ""


def pair_hint(text):
    found, seen = [], set()
    for m in PAIR.finditer(text):
        w = m.group(0).strip()
        key = re.sub(r"[\s_-]", "", w.lower())
        if key not in seen:
            seen.add(key)
            found.append(w)
    return ";".join(found)


def platform_group(p):
    p = p.upper()
    if p in LONG_PLATFORMS:
        return "long"
    if p in SHORT_PLATFORMS:
        return "short"
    return "other"


def review_flags(df):
    """Flags for things metadata cannot settle. They never change a label."""
    flags = {idx: [] for idx in df.index}

    # Single-nucleus projects often look like ordinary WGS: one BioSample with
    # many separate short-read libraries (one per nucleus).
    dna = df[(df["data_class"].isin(["wgs_short", "dna_metagenome"]))
             & (df["platform_group"] == "short") & (df["biosample"] != "")]
    for bs, g in dna.groupby("biosample"):
        names = set(g["library_name"]) - {""}
        n = len(names) if names else len(g)
        if n >= MANY_LIBRARIES:
            for idx in g.index:
                flags[idx].append(f"{n}_libraries_one_biosample: single nuclei or pooled design?")

    odd = df["text"].str.contains(ODD_ORGANISM)
    for idx in df.index[odd]:
        flags[idx].append("metadata mentions human, patient, gut or mouse: is this AMF data?")
    return pd.Series({idx: "; ".join(v) for idx, v in flags.items()}).reindex(df.index)


def tag(df):
    df = df.copy()
    df["bases_num"] = df["bases"].map(parse_bases)
    df["data_class"] = df.apply(classify, axis=1)
    df["platform_group"] = df["platform"].map(platform_group)
    df["rna_prep"] = df.apply(lambda r: rna_prep(r) if r["data_class"] == "rna" else "", axis=1)
    df["amplicon_target"] = df.apply(
        lambda r: amplicon_target(r) if r["data_class"] == "amplicon" else "", axis=1)
    sv = df.apply(screen, axis=1)
    df["screen_value"] = [s for s, _ in sv]
    df["reason"] = [why for _, why in sv]
    df["role"] = "screen"
    df["review_flag"] = review_flags(df)
    df["override_note"] = ""
    df["culture_hint"] = df["text"].map(culture_hint)
    df["isolate_hint"] = df.apply(isolate_hint, axis=1)
    df["pair_hint"] = df["text"].map(pair_hint)
    return df


# --------------------------------------------------------------- overrides ---

def load_overrides(path):
    ov = read_table(path)
    ov.columns = [c.strip().lower() for c in ov.columns]
    if "only_class" not in ov.columns:
        ov["only_class"] = ""
    missing = [c for c in ("accession", "field", "value", "note") if c not in ov.columns]
    if missing:
        sys.exit(f"{path}: missing column(s) {missing}; need accession, field, value, note "
                 f"(only_class optional)")
    ov = ov[ov["accession"].str.strip() != ""].copy()
    ov = ov[~ov["accession"].str.strip().str.startswith("#")]
    # "-" (or blank) in only_class means "any class".
    ov["only_class"] = ov["only_class"].str.strip().replace("-", "")
    problems = []
    for i, row in ov.iterrows():
        f, v = row["field"].strip(), row["value"].strip()
        line = i + 2  # header is line 1
        if f not in OVERRIDE_VALUES:
            problems.append(f"line {line}: field {f!r} is not one of {sorted(OVERRIDE_VALUES)}")
        elif v not in OVERRIDE_VALUES[f]:
            problems.append(f"line {line}: {v!r} is not a valid {f} ({sorted(OVERRIDE_VALUES[f])})")
        oc = row["only_class"].strip()
        if oc and oc not in DATA_CLASSES:
            problems.append(f"line {line}: only_class {oc!r} is not a data class")
        if not row["note"].strip():
            problems.append(f"line {line}: note is empty (say why the label is changed)")
    if problems:
        sys.exit(f"{path}:\n  " + "\n  ".join(problems))
    return ov


def apply_overrides(df, ov):
    """Apply override rows in order. Returns the table and a per-row report."""
    df = df.copy()
    before = df[["data_class", "rna_prep", "amplicon_target"]].copy()
    screen_set = pd.Series(False, index=df.index)
    report = []
    for _, row in ov.iterrows():
        acc, field, value = row["accession"].strip(), row["field"].strip(), row["value"].strip()
        only, note = row["only_class"].strip(), row["note"].strip()
        mask = ((df["run"] == acc) | (df["experiment"] == acc)
                | (df["biosample"] == acc) | (df["bioproject"] == acc))
        if only:
            mask &= df["data_class"] == only
        n = int(mask.sum())
        report.append(dict(accession=acc, only_class=only or "-", field=field, value=value,
                           runs=n, note=note))
        if not n:
            continue
        df.loc[mask, field] = value
        df.loc[mask, "override_note"] = df.loc[mask, "override_note"].map(
            lambda s: f"{s}; {note}" if s else note)
        if field == "screen_value":
            screen_set |= mask
            df.loc[mask, "reason"] = "set by override"

    # Keep the dependent labels consistent with a changed data_class.
    changed_class = df["data_class"] != before["data_class"]
    to_rna = changed_class & (df["data_class"] == "rna") & (df["rna_prep"] == "")
    df.loc[to_rna, "rna_prep"] = "unknown"
    df.loc[changed_class & (df["data_class"] != "rna"), "rna_prep"] = ""
    to_amp = changed_class & (df["data_class"] == "amplicon") & (df["amplicon_target"] == "")
    df.loc[to_amp, "amplicon_target"] = "unknown"
    df.loc[changed_class & (df["data_class"] != "amplicon"), "amplicon_target"] = ""

    relabelled = (df[["data_class", "rna_prep", "amplicon_target"]] != before).any(axis=1)
    redo = relabelled & ~screen_set
    if redo.any():
        sv = df[redo].apply(screen, axis=1)
        df.loc[redo, "screen_value"] = [s for s, _ in sv]
        df.loc[redo, "reason"] = [why for _, why in sv]
    return df, pd.DataFrame(report)


# ----------------------------------------------------------------- summary ---

def md_table(frame):
    cols = [str(c) for c in frame.columns]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for row in frame.itertuples(index=False):
        cells = [str(v).replace("|", "/").replace("\n", " ") for v in row]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def gbp(series):
    return round(series.fillna(0).sum() / 1e9, 2)


def top_values(series, n=3):
    vals = [v for v in series if v]
    counts = pd.Series(vals).value_counts() if vals else pd.Series(dtype=int)
    return "; ".join(f"{k} ({v})" for k, v in counts.head(n).items())


def top_tokens(series, n=4, sep=";"):
    """Like top_values, but splits joined labels and counts each one."""
    vals = [tok.strip() for v in series if v for tok in v.split(sep) if tok.strip()]
    counts = pd.Series(vals).value_counts() if vals else pd.Series(dtype=int)
    return "; ".join(f"{k} ({v})" for k, v in counts.head(n).items())


def summarise(df, inputs, ov_report=None):
    out = ["# AMF run inventory", "",
           f"Made by tag_runs.py v{__version__} on {dt.date.today().isoformat()}.", ""]
    out.append("## Inputs")
    files_per_run = df["found_in"].str.split(";")
    missing = {}
    for name, n in inputs:
        has = files_per_run.map(lambda fs: name in fs)
        missing[name] = sorted(df.loc[~has, "run"])
        out.append(f"- {name}: {n} rows, {int(has.sum())} unique runs; "
                   f"{len(missing[name])} runs found in other inputs but not in this one")
    out.append(f"- Unique runs after merging: {len(df)}; BioProjects: "
               f"{df['bioproject'].replace('', pd.NA).nunique()}; "
               f"organisms: {df['organism'].replace('', pd.NA).nunique()}")
    out.append("")
    if len(inputs) > 1 and any(missing.values()):
        out.append("## Runs missing from some inputs (first 25 per file)")
        out.append("Expected when the inputs come from different queries; worth a look "
                   "when they are two exports of the same query.")
        for name, runs in missing.items():
            if runs:
                shown = ", ".join(runs[:25]) + (" ..." if len(runs) > 25 else "")
                out.append(f"- not in {name}: {shown}")
        out.append("")

    roles = df["role"].value_counts()
    if set(roles.index) != {"screen"}:
        out.append("## Roles (set in the overrides file)")
        out += [md_table(roles.rename_axis("role").reset_index(name="runs")), ""]
        out.append("Everything below counts only runs with role 'screen', unless it says otherwise.")
        out.append("")
    s = df[df["role"] == "screen"]

    out.append("## Runs by data class and platform")
    ct = pd.crosstab(s["data_class"], s["platform_group"], margins=True, margins_name="total")
    ct = ct.reset_index().rename(columns={"data_class": "data class"})
    out += [md_table(ct), ""]

    out.append("## Screening value")
    sv = (s.groupby("screen_value")
            .agg(runs=("run", "size"), gbp=("bases_num", gbp))
            .sort_values("runs", ascending=False).reset_index())
    out += [md_table(sv), ""]

    rna = s[s["data_class"] == "rna"]
    if len(rna):
        out.append("## RNA library prep (from metadata)")
        rp = rna["rna_prep"].value_counts().rename_axis("rna_prep").reset_index(name="runs")
        out += [md_table(rp), ""]

    amp = s[s["data_class"] == "amplicon"]
    if len(amp):
        out.append("## Amplicon targets (guessed from text)")
        at = amp["amplicon_target"].value_counts().rename_axis("target").reset_index(name="runs")
        out += [md_table(at), ""]

    out.append("## Top organisms")
    org = s["organism"].replace("", "(blank)").value_counts().head(15)
    out += [md_table(org.rename_axis("organism").reset_index(name="runs")), ""]

    useful = s[s["screen_value"].isin(["high", "medium", "relevant_16S"])]
    out.append("## BioProjects most worth screening (high / medium / 16S)")
    if len(useful):
        bp = (useful.groupby("bioproject")
                .agg(runs=("run", "size"),
                     gbp=("bases_num", gbp),
                     organisms=("organism", top_values),
                     classes=("data_class", top_values),
                     culture=("culture_hint", top_tokens),
                     isolates=("isolate_hint", top_values),
                     pair_hints=("pair_hint", top_tokens),
                     flags=("review_flag", lambda x: "yes" if (x != "").any() else ""))
                .sort_values(["runs", "gbp"], ascending=False)
                .head(40).reset_index())
        out += [md_table(bp), ""]
    else:
        out += ["None found.", ""]

    pairs = s[s["pair_hint"] != ""]
    out.append("## Possible with/without-endosymbiont designs (text mentions B+/B-, cured, WT)")
    if len(pairs):
        pp = (pairs.groupby("bioproject")
                .agg(runs=("run", "size"),
                     organisms=("organism", top_values),
                     classes=("data_class", top_values),
                     hints=("pair_hint", top_tokens))
                .sort_values("runs", ascending=False).reset_index())
        out += [md_table(pp), ""]
    else:
        out += ["None found.", ""]

    flagged = df[df["review_flag"] != ""]
    out.append("## Review flags (all roles; labels are not changed by flags)")
    if len(flagged):
        def flag_kinds(x):
            kinds = [re.sub(r"^\d+_", "n_", tok.split(":")[0].strip())
                     for v in x for tok in v.split(";") if tok.strip()]
            return "; ".join(sorted(set(kinds)))
        fl = (flagged.groupby("bioproject")
                .agg(runs=("run", "size"),
                     gbp=("bases_num", gbp),
                     organisms=("organism", top_values),
                     classes=("data_class", top_values),
                     flags=("review_flag", flag_kinds),
                     role=("role", top_values),
                     overridden=("override_note", lambda x: "yes" if (x != "").any() else ""))
                .sort_values(["runs", "gbp"], ascending=False).reset_index())
        out += [md_table(fl), ""]
        out.append("A flagged BioProject with no override still needs a decision: "
                   "add a row to the overrides file, or leave it and say why.")
        out.append("")
    else:
        out += ["None.", ""]

    if ov_report is not None and len(ov_report):
        out.append("## Overrides applied")
        out += [md_table(ov_report[["accession", "only_class", "field", "value", "runs", "note"]]), ""]
        unused = ov_report[ov_report["runs"] == 0]
        if len(unused):
            out.append("Rows that matched no run (check the accession, or the run is not in "
                       "these inputs): " + ", ".join(unused["accession"]))
            out.append("")

    out += ["## Caveats",
            "- Labels come from metadata and free text only. 'check_prep' and 'check' "
            "mean read the paper, not that the data are unusable.",
            "- 'uninformative' runs cannot show absence: no endosymbiont call should be "
            "made from them.",
            "- AMF data can be filed under 'fungus metagenome', plant taxa, or the "
            "endobacterium's own name; an organism=Glomeromycota query misses all three. "
            "See PROTOCOL.md section 3.1 for the extra queries.", ""]
    return "\n".join(out)


# -------------------------------------------------------------------- main ---

KEEP = ["run", "bioproject", "biosample", "experiment", "organism", "strategy", "source",
        "selection", "layout", "platform", "instrument", "bases", "library_name",
        "sample_name", "center", "release", "input_file", "found_in",
        "data_class", "platform_group", "rna_prep", "amplicon_target", "screen_value",
        "reason", "role", "review_flag", "override_note",
        "culture_hint", "isolate_hint", "pair_hint"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tables", nargs="+", help="SraRunTable.csv, runinfo CSV or ENA TSV files")
    ap.add_argument("-o", "--outdir", default="inventory", help="output folder (default: inventory)")
    ap.add_argument("--overrides", help="TSV of manual label fixes (see above)")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args(argv)

    ov = load_overrides(args.overrides) if args.overrides else None

    frames, inputs = [], []
    for t in args.tables:
        raw = read_table(t)
        inputs.append((Path(t).name, len(raw)))
        frames.append(normalise(raw, t))
    merged = pd.concat(frames, ignore_index=True)
    merged = merged[merged["run"] != ""]
    # EDirect runinfo output can repeat its header line between chunks.
    merged = merged[~merged["run"].str.lower().isin({"run", "run_accession"})]
    found_in = merged.groupby("run")["input_file"].agg(lambda s: ";".join(sorted(set(s))))
    merged = merged.drop_duplicates(subset="run", keep="first").copy()
    merged["found_in"] = merged["run"].map(found_in)
    tagged = tag(merged)

    ov_report = None
    if ov is not None:
        tagged, ov_report = apply_overrides(tagged, ov)
        for _, r in ov_report[ov_report["runs"] == 0].iterrows():
            print(f"warning: override for {r['accession']} ({r['field']}={r['value']}) "
                  f"matched no run", file=sys.stderr)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tagged[KEEP].to_csv(outdir / "runs_tagged.tsv", sep="\t", index=False)
    (outdir / "inventory_summary.md").write_text(summarise(tagged, inputs, ov_report),
                                                 encoding="utf-8")

    s = tagged[tagged["role"] == "screen"]
    counts = s["screen_value"].value_counts()
    print(f"{len(tagged)} unique runs -> {outdir}/runs_tagged.tsv, {outdir}/inventory_summary.md")
    if len(s) != len(tagged):
        print("role: " + ", ".join(f"{k}={v}" for k, v in tagged["role"].value_counts().items()))
    print("screen_value (role=screen): " + ", ".join(f"{k}={v}" for k, v in counts.items()))
    n_flag = int((tagged["review_flag"] != "").sum())
    if n_flag:
        print(f"{n_flag} runs have review flags: see 'Review flags' in the summary")


if __name__ == "__main__":
    main()
