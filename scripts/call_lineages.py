#!/usr/bin/env python3
"""The lineage caller (PROTOCOL.md 4.1; section 10, entries 7 to 10 and 16). Functions only; run_t1.py uses them.

Route reads get a label:
  ssu_label       from a read's top SILVA hits (entry 8's VSEARCH command): dropped if a top hit lies
                  outside Bacteria or in chloroplasts or mitochondria (entry 7); a genus when all top hits
                  share one at 94.5% identity or more (4.1); unplaced when the top hits are at 94.5% or
                  more but in several genera (entry 10); otherwise novel, placed at the lowest rank all top
                  hits share, if the SSU piece is at least 60 bp, and unplaced if it is shorter (entry 10)
  singlem_label   from SingleM's taxonomy: dropped unless Bacteria (entry 7); a genus if it reaches one;
                  otherwise novel
Route lineages (entry 9): one per genus, matched across routes (no "Candidatus", no GTDB suffix such as
"_A", no SILVA suffix "--other"; a SILVA name joining genera with hyphens matches each of them), and one
novel lineage per route, "novel (SSU)" and "novel (SingleM)". A route lineage in Glomeribacter, Moeniiplasma
or Mycoavidus is a known endosymbiont clade (E6). That is Phase 0's list (KNOWN_CLADE_GENERA, the default).
The survey's list (SURVEY_CLADE_GENERA, entry 16) adds CGIF01, GTDB's name for the MRE group, as CaMg/MRE;
category_of_route uses it only when given it.

Levels per run (4.1): N = t x reads for a mapped lineage, t x m x reads for a route lineage; present if
n >= N, low if n >= N/10, otherwise absent. A run is informative when N >= 10 and its library can show
absence (DNA or total RNA; entry 9). Cross-talk only lowers "present" to "low": within one BioProject,
when the run's fraction is below 2% of the highest fraction in another run of that BioProject. The
spread check also only lowers "present" to "low" (entry 10): in DNA runs, a mapped lineage none of whose
genomes has a spread ratio of 0.5 or more (assign_reads.py's spread.tsv) holds its reads in a few places.

Calls per isolate or line, in order: no call (no informative run and no signal); absent (every informative
run absent); contaminant (present somewhere and E5 fails, or a kit-contaminant genus without both E1 and
E2); resident (a known endosymbiont clade present in 2 or more runs, or E1 and E2); unresolved.
E1 and E2 use DNA runs only (entry 9). E5 counts unrelated fungal isolates from the same lab where the
lineage is present at a read fraction of at least a tenth of the line's own highest.
"""

import re
from itertools import combinations

__version__ = "0.3"

KNOWN_CLADE_GENERA = {"Glomeribacter": "CaGg", "Moeniiplasma": "CaMg/MRE", "Mycoavidus": "Mycoavidus"}  # Phase 0
SURVEY_CLADE_GENERA = {**KNOWN_CLADE_GENERA, "CGIF01": "CaMg/MRE"}  # entry 16: the survey only (GTDB R232's name)
KNOWN_CLADES = set(KNOWN_CLADE_GENERA.values())       # mapped lineages that are known endosymbiont clades
KIT_GENERA = {"Bradyrhizobium", "Burkholderia", "Herbaspirillum", "Pseudomonas", "Ralstonia", "Sphingomonas",
              "Stenotrophomonas"}                      # 4.1, E6; one reference genome each (section 3.3)
GENUS_ID = 94.5                                        # 4.1 (Yarza et al. 2014)
CROSS_TALK = 0.02                                      # 4.1
NOT_A_NAME = re.compile(r"^(incertae sedis|uncultured|unidentified|unknown|unclassified|metagenome)", re.I)
ORGANELLE = re.compile(r"^(chloroplast|mitochondri)", re.I)                      # entry 7
SILVA_RANKS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus"]  # SILVA 144 (entry 9)
GTDB_RANKS = {"d": "domain", "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus", "s": "species"}
CATEGORIES = ("known endosymbiont clade", "kit-contaminant genus", "other known genus", "novel")
NOVEL_SSU, NOVEL_SINGLEM = "novel (SSU)", "novel (SingleM)"
NOVEL_MIN_BP = 60                                      # entry 10: a novel SSU piece is at least 60 bp
UNPLACED = "unplaced"                                  # entry 10
SEVERAL_GENERA, TOO_SHORT = "best hits in several genera", "piece under 60 bp"
SPREAD_MIN = 0.5                                       # entry 10


# ------------------------------------------------------------------ names ---

def clean_genus(name):
    """A genus name as matched across routes, or None if it is not a genus name."""
    n = (name or "").strip()
    if n.startswith("g__"):
        n = n[3:]
    n = re.sub(r"^Candidatus\s+", "", n)
    n = re.sub(r"--other$", "", n)
    n = re.sub(r"_[A-Z]+$", "", n)
    if not n or NOT_A_NAME.match(n):
        return None
    return n


def genus_parts(name):
    """The genera a name stands for: a SILVA name joining genera with hyphens stands for each of them."""
    parts = name.split("-")
    if len(parts) > 1 and all(p[:1].isupper() for p in parts):
        return parts
    return [name]


def category_of_route(genera, clades=None):
    """E6 for a route lineage made of these genus names (None: a novel lineage). clades: the genus names of
    the known endosymbiont clades; None means Phase 0's, KNOWN_CLADE_GENERA (the survey's: SURVEY_CLADE_GENERA)."""
    if not genera:
        return "novel"
    parts = {p for g in genera for p in genus_parts(g)}
    if parts & (KNOWN_CLADE_GENERA if clades is None else clades).keys():
        return "known endosymbiont clade"
    if parts & KIT_GENERA:
        return "kit-contaminant genus"
    return "other known genus"


def category_of_mapped(lineage):
    return "known endosymbiont clade" if lineage in KNOWN_CLADES else "kit-contaminant genus"


# ----------------------------------------------------------------- labels ---

def ssu_label(hits, tax, piece_len=None):
    """A read's label from its top hits. hits: [(target, identity)], all equally good; tax: target -> the
    SILVA path as a list (Domain to Genus, then the organism); piece_len: the SSU piece's length (entry 10;
    None leaves the length out). Returns (status, genus, placement) with status 'kept', 'outside Bacteria',
    'organelle' or 'unplaced' (then placement says why); genus None means novel."""
    short = piece_len is not None and piece_len < NOVEL_MIN_BP
    if not hits:
        return (UNPLACED, None, TOO_SHORT) if short else ("kept", None, "no hit")
    paths = [tax[t] for t, _ in hits]
    if any(p[0] != "Bacteria" for p in paths):
        return "outside Bacteria", None, None
    if any(ORGANELLE.match(x) for p in paths for x in p):
        return "organelle", None, None
    ranks = [p[:len(SILVA_RANKS)] for p in paths]
    shared = []
    for i in range(min(len(r) for r in ranks)):
        names = {r[i] for r in ranks}
        if len(names) != 1:
            break
        shared.append(names.pop())
    top = max(i for _, i in hits)
    if len(shared) == len(SILVA_RANKS):
        g = clean_genus(shared[-1])
        if g and top >= GENUS_ID:
            return "kept", g, f"genus {g}"
        shared = shared[:-1]               # below 94.5% the genus is not a placement (4.1)
    if top >= GENUS_ID and len({tuple(r) for r in ranks}) > 1:
        return UNPLACED, None, SEVERAL_GENERA                          # entry 10
    if short:
        return UNPLACED, None, TOO_SHORT                               # entry 10
    named = [(SILVA_RANKS[i], x) for i, x in enumerate(shared) if not NOT_A_NAME.match(x)]
    rank, name = named[-1] if named else ("domain", "Bacteria")
    return "kept", None, f"{rank} {re.sub(r'--other$', '', name)}"


def singlem_label(taxonomies):
    """A read's label from the SingleM taxonomies its OTUs got (usually one). Returns (status, genus,
    placement): status 'kept' only for Bacteria; a genus label is preferred, then the deepest placement."""
    best = None
    for tx in taxonomies:
        parts = [p.strip() for p in tx.split(";")]
        if len(parts) < 2 or parts[1] != "d__Bacteria":
            continue
        g = next((clean_genus(p) for p in parts if p.startswith("g__")), None)
        named = [p for p in parts[1:] if len(p) > 3 and p[1:3] == "__" and p[3:]]
        last = named[-1]
        placement = f"genus {g}" if g else f"{GTDB_RANKS.get(last[0], last[0])} {last[3:]}"
        key = (g is not None, len(named))
        if best is None or key > best[0]:
            best = (key, g, placement)
    if best is None:
        return "not Bacteria", None, None
    return "kept", best[1], best[2]


def read_label(ssu, singlem):
    """One routed read's lineage key: ('genus', name), ('novel', route) or None. A read both routes find
    counts once, under its SSU label when the SSU route keeps it (entry 9)."""
    if ssu and ssu[0] == "kept":
        return ("genus", ssu[1]) if ssu[1] else ("novel", "SSU")
    if singlem and singlem[0] == "kept":
        return ("genus", singlem[1]) if singlem[1] else ("novel", "SingleM")
    return None


class Joiner:
    """Genus names into route lineages: a composite SILVA name joins every one of its genera seen."""

    def __init__(self, names):
        self.parent = {n: n for n in names}
        for n in names:
            for p in genus_parts(n):
                if p != n and p in self.parent:
                    self._union(n, p)
        groups = {}
        for n in names:
            groups.setdefault(self._find(n), []).append(n)
        self.members, self.label = {}, {}
        for members in groups.values():
            name = sorted(members, key=lambda x: (-len(genus_parts(x)), x))[0] + " (route)"
            for n in members:
                self.label[n] = name
            self.members[name] = sorted(members)

    def _find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def _union(self, a, b):
        ra, rb = self._find(a), self._find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def lineage_name(key, joiner):
    if key[0] == "novel":
        return NOVEL_SSU if key[1] == "SSU" else NOVEL_SINGLEM
    return joiner.label[key[1]]


# ----------------------------------------------------------------- levels ---

def level(n, N):
    if n >= N:
        return "present"
    if n * 10 >= N:
        return "low"
    return "absent"


def run_levels(records):
    """records: dicts with run, bioproject, n, N, reads, absent_ok, and optionally spread (the lineage's highest
    spread ratio in the run, or None where the check does not apply). Adds fraction, level, informative,
    spread_cap and cross_talk; returns them. All records must be for one lineage."""
    for x in records:
        x["fraction"] = x["n"] / x["reads"] if x["reads"] else 0.0
        x["level"] = level(x["n"], x["N"])
        x["informative"] = x["absent_ok"] and x["N"] >= 10
        x["cross_talk"] = False
        x["spread_cap"] = x["level"] == "present" and x.get("spread") is not None and x["spread"] < SPREAD_MIN
        if x["spread_cap"]:
            x["level"] = "low"
    for x in records:
        if x["level"] != "present":
            continue
        others = [y["fraction"] for y in records if y["bioproject"] == x["bioproject"] and y["run"] != x["run"]]
        if others and x["fraction"] < CROSS_TALK * max(others):
            x["level"], x["cross_talk"] = "low", True
    return records


# ------------------------------------------------------------------ calls ---

def independent(a, b):
    """E1: runs in different BioProjects, or in one BioProject with collection years that differ."""
    return a["bioproject"] != b["bioproject"] or bool(a["year"] and b["year"] and a["year"] != b["year"])


def call(records, category, background):
    """The call for one lineage in one isolate or line. records: the line's runs for the lineage (with
    level, informative, kind, bioproject, year, n, fungal, fraction); background: the read fractions of
    the lineage in present runs of unrelated isolates from the same labs, as {isolate: [fraction, ...]}.
    Returns a dict with the call, its basis and the evidence."""
    present = [x for x in records if x["level"] == "present"]
    signal = any(x["level"] in ("present", "low") for x in records)
    informative = [x for x in records if x["informative"]]
    dna = [x for x in present if x["kind"] == "DNA"]
    e1 = any(independent(a, b) for a, b in combinations(dna, 2))
    ratios = [x["n"] / x["fungal"] for x in dna if x["fungal"] > 0]
    e2 = bool(dna) and len(ratios) == len(dna) and max(ratios) <= 10 * min(ratios)
    top = max((x["fraction"] for x in present), default=0.0)
    e5_hits = sorted(i for i, fs in background.items() if present and any(f * 10 >= top for f in fs))
    e5 = len(e5_hits) < 2
    out = {"signal": signal, "informative_runs": len(informative), "present_runs": len(present),
           "low_runs": sum(x["level"] == "low" for x in records), "e1": e1, "e2": e2, "e5": e5,
           "e5_isolates": e5_hits, "category": category}
    if not signal and not informative:
        return {**out, "call": "no call", "basis": "no informative run and no signal"}
    if informative and all(x["level"] == "absent" for x in informative):
        return {**out, "call": "absent", "basis": f"all {len(informative)} informative runs absent"}
    if present and not e5:
        return {**out, "call": "contaminant", "basis": "E5 fails: present at a similar level in "
                                                       + ", ".join(e5_hits)}
    if category == "kit-contaminant genus" and not (e1 and e2):
        return {**out, "call": "contaminant", "basis": "kit-contaminant genus without both E1 and E2"}
    if category == "known endosymbiont clade" and len(present) >= 2:
        return {**out, "call": "resident", "basis": f"known endosymbiont clade, present in {len(present)} runs"}
    if e1 and e2:
        return {**out, "call": "resident", "basis": "E1 and E2 (DNA runs)"}
    return {**out, "call": "unresolved", "basis": "signal, but no rule above applies"}
