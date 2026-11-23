"""Check call_lineages.py and run_t1.py (PROTOCOL.md section 10, entry 9).

Part 1, the caller on hand-made inputs: genus names (Candidatus, GTDB and SILVA suffixes, names that
are not genera, hyphenated SILVA names), SSU labels (no hit, outside Bacteria, organelle, genus at and
below 94.5%, placement), SingleM labels, one count per read with SSU first, joining route lineages,
levels, cross-talk (lowers present, never raises absent), and every step of the call: no call, absent,
contaminant by E5 (one-sided) and by kit genus, resident by known clade (RNA runs count) and by E1 and
E2, unresolved; E1 and E2 on DNA runs only; collection years; then K1.

Part 2, run_t1.py end to end on synthetic genomes. minimap2, assign_reads.py, make_index.py and VSEARCH
are real; stand-ins play fastq-dump, SSUplex, nhmmer, SingleM (through conda) and sylph. Expected numbers
are worked out here from the planted reads: mapped counts, route counts per lineage, levels and the
cross-talk cap, the split reads of a Hi-C-like chimera, SSU reads dropped as organelle or outside
Bacteria, SingleM reads dropped as not Bacteria, every call, sylph, K1. A run marked reuse=stage0a calls
no tool but sylph; a rerun calls none; refusals (no entry 8 m, a bad runs.tsv, another sylph).

Needs minimap2 and VSEARCH (activate phase0-reads, or set MINIMAP2= and VSEARCH=).
"""
import contextlib
import csv
import gzip
import io
import json
import os
import random
import shutil
import stat
import sys
import tempfile
from pathlib import Path

here = Path(__file__).parent
sys.path.insert(0, str(here.parent / "scripts"))
import call_lineages as C  # noqa: E402
import make_index as MI  # noqa: E402
import run_t1 as T  # noqa: E402

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)


# ===================================================================== part 1 ===

check(C.clean_genus("g__Pseudomonas_E") == "Pseudomonas", "GTDB suffix")
check(C.clean_genus("Candidatus Glomeribacter") == "Glomeribacter", "Candidatus")
check(C.clean_genus("Pseudomonas--other") == "Pseudomonas", "SILVA --other")
check(C.clean_genus("Incertae Sedis") is None and C.clean_genus("uncultured") is None and C.clean_genus("") is None,
      "names that are not genera")
check(C.genus_parts("Burkholderia-Caballeronia-Paraburkholderia") == ["Burkholderia", "Caballeronia",
                                                                      "Paraburkholderia"], "hyphenated SILVA name")
check(C.genus_parts("Pseudomonas") == ["Pseudomonas"] and C.genus_parts("abc-def") == ["abc-def"], "plain names")
check(C.category_of_route(["Glomeribacter"]) == "known endosymbiont clade" and
      C.category_of_route(["Burkholderia-Caballeronia-Paraburkholderia"]) == "kit-contaminant genus" and
      C.category_of_route(["Duganella"]) == "other known genus" and C.category_of_route([]) == "novel",
      "E6 for route lineages")
check(C.category_of_mapped("CaGg") == "known endosymbiont clade" and
      C.category_of_mapped("Ralstonia") == "kit-contaminant genus", "E6 for mapped lineages")


def path(*ranks):
    return list(ranks)


BURK = ["Bacteria", "Pseudomonadati", "Pseudomonadota", "Gammaproteobacteria", "Burkholderiales"]
TAX = {"G1": path(*BURK, "Burkholderiaceae", "Candidatus Glomeribacter", "org"),
       "G2": path(*BURK, "Burkholderiaceae", "Candidatus Glomeribacter", "org2"),
       "M1": path(*BURK, "Burkholderiaceae", "Mycoavidus", "org"),
       "I1": path(*BURK, "Burkholderiaceae", "Incertae Sedis", "org"),
       "O1": path(*BURK[:4], "Burkholderiales--other", "Incertae Sedis", "Incertae Sedis", "org"),
       "O2": path(*BURK[:4], "Burkholderiales--other", "Oxalobacteraceae", "Massilia", "org"),
       "BCP": path(*BURK, "Burkholderiaceae", "Burkholderia-Caballeronia-Paraburkholderia", "org"),
       "CH": path("Bacteria", "Bacillati", "Cyanobacteriota", "Cyanophyceae", "Chloroplast", "Incertae Sedis",
                  "Incertae Sedis", "Lotus"),
       "MI": path("Bacteria", "Pseudomonadati", "Pseudomonadota", "Alphaproteobacteria", "Rickettsiales",
                  "Mitochondria", "Incertae Sedis", "flax"),
       "EU": path("Eukaryota", "Amorphea", "Obazoa", "Opisthokonta", "Holomycota", "Fungi", "Glomeromycota", "x")}
check(C.ssu_label([], TAX) == ("kept", None, "no hit"), "SSU: no hit is kept, novel")
check(C.ssu_label([("G1", 100.0), ("EU", 100.0)], TAX)[0] == "outside Bacteria", "SSU: a eukaryotic top hit drops")
check(C.ssu_label([("G1", 99.0), ("CH", 99.0)], TAX)[0] == "organelle" and
      C.ssu_label([("MI", 99.0)], TAX)[0] == "organelle", "SSU: chloroplast or mitochondrion drops")
check(C.ssu_label([("G1", 99.0), ("G2", 99.0)], TAX) == ("kept", "Glomeribacter", "genus Glomeribacter"),
      "SSU: one genus at 99%")
check(C.ssu_label([("G1", 94.5)], TAX)[1] == "Glomeribacter" and C.ssu_label([("G1", 94.4)], TAX) ==
      ("kept", None, "family Burkholderiaceae"), "SSU: 94.5% is the genus threshold; below it, the family")
check(C.ssu_label([("G1", 99.0), ("M1", 99.0)], TAX) == ("kept", None, "family Burkholderiaceae"),
      "SSU: two genera give the family")
check(C.ssu_label([("I1", 99.0)], TAX) == ("kept", None, "family Burkholderiaceae"), "SSU: Incertae Sedis genus")
check(C.ssu_label([("O1", 99.0), ("O2", 99.0)], TAX) == ("kept", None, "order Burkholderiales"),
      "SSU: --other dropped from the placement name")
check(C.ssu_label([("BCP", 100.0)], TAX)[1] == "Burkholderia-Caballeronia-Paraburkholderia", "SSU: hyphenated genus")

GT = "Root; d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria"
check(C.singlem_label([GT + "; o__Burkholderiales; f__Burkholderiaceae; g__Pseudomonas_E"]) ==
      ("kept", "Pseudomonas", "genus Pseudomonas"), "SingleM: genus with suffix")
check(C.singlem_label([GT]) == ("kept", None, "class Gammaproteobacteria"), "SingleM: novel, class")
check(C.singlem_label(["Root; d__Eukaryota"])[0] == "not Bacteria" and C.singlem_label(["Root"])[0] == "not Bacteria",
      "SingleM: not Bacteria")
check(C.singlem_label([GT, GT + "; o__X; f__Y; g__Massilia", "Root; d__Archaea"])[1] == "Massilia",
      "SingleM: a genus label is preferred")
gen, nov = ("kept", "Massilia", "genus Massilia"), ("kept", None, "family X")
check(C.read_label(gen, None) == ("genus", "Massilia") and C.read_label(nov, ("kept", "Duganella", "g")) ==
      ("novel", "SSU") and C.read_label(("organelle", None, None), ("kept", None, "c")) == ("novel", "SingleM") and
      C.read_label(("outside Bacteria", None, None), ("not Bacteria", None, None)) is None,
      "one count per read, SSU first")
j = C.Joiner(["Burkholderia", "Burkholderia-Caballeronia-Paraburkholderia", "Paraburkholderia", "Pseudomonas"])
check(j.label["Burkholderia"] == j.label["Paraburkholderia"] == "Burkholderia-Caballeronia-Paraburkholderia (route)"
      and j.label["Pseudomonas"] == "Pseudomonas (route)" and len(j.members) == 2, "joining route lineages")

check([C.level(n, 40) for n in (40, 39, 4, 3)] == ["present", "low", "low", "absent"] and
      C.level(0, 0.5) == "absent" and C.level(1, 0.5) == "present", "levels")


def rec(run, bp, n, N=40, reads=4000, kind="DNA", absent_ok=True, fungal=1000, year=""):
    return {"run": run, "bioproject": bp, "n": n, "N": N, "reads": reads, "kind": kind, "absent_ok": absent_ok,
            "fungal": fungal, "year": year}


xs = C.run_levels([rec("A", "P1", 2400), rec("B", "P1", 44), rec("D", "P1", 2), rec("E", "P2", 44)])
lv = {x["run"]: (x["level"], x["cross_talk"]) for x in xs}
check(lv == {"A": ("present", False), "B": ("low", True), "D": ("absent", False), "E": ("present", False)},
      f"cross-talk lowers present only, within one BioProject: {lv}")
xs = C.run_levels([rec("A", "P1", 50, N=5), rec("B", "P1", 50, N=5, absent_ok=False)])
check([x["informative"] for x in xs] == [False, False], "informative needs N >= 10 and a library that can show absence")


def call(recs, category="other known genus", background=None):
    return C.call(C.run_levels(recs), category, background or {})


c = call([rec("A", "P1", 0, absent_ok=False)])
check(c["call"] == "no call", "no informative run, no signal: no call")
c = call([rec("A", "P1", 0), rec("B", "P2", 39, kind="RNA", absent_ok=False)])
check(c["call"] == "absent", "every informative run absent: absent, whatever the other runs show")
c = call([rec("A", "P1", 39)], category="kit-contaminant genus")
check(c["call"] == "contaminant" and "kit" in c["basis"], "kit genus with signal, no E1: contaminant")
c = call([rec("A", "P1", 400), rec("B", "P2", 800)], category="kit-contaminant genus")
check(c["call"] == "resident" and c["e1"] and c["e2"], "kit genus with E1 and E2: resident")
c = call([rec("A", "P1", 400, kind="RNA", absent_ok=False), rec("B", "P1", 400, kind="RNA", absent_ok=False)],
         category="known endosymbiont clade")
check(c["call"] == "resident" and "known" in c["basis"], "known clade present in 2 RNA runs: resident")
c = call([rec("A", "P1", 400), rec("B", "P2", 10)], category="known endosymbiont clade")
check(c["call"] == "unresolved", "known clade present in 1 run: unresolved")
c = call([rec("A", "P1", 400), rec("B", "P2", 400, fungal=11000)])
check(c["call"] == "unresolved" and c["e1"] and not c["e2"], "E2 fails beyond 10-fold: unresolved")
c = call([rec("A", "P1", 400, kind="RNA", absent_ok=False), rec("B", "P2", 400, kind="RNA", absent_ok=False)])
check(c["call"] == "unresolved" and not c["e1"], "E1 and E2 use DNA runs only")
c = call([rec("A", "P1", 400, year="2019"), rec("B", "P1", 400, year="2021")])
check(c["e1"] and c["call"] == "resident", "one BioProject, different years: independent")
c = call([rec("A", "P1", 400, year="2019"), rec("B", "P1", 400)])
check(not c["e1"], "a missing year counts as the same year")
c = call([rec("A", "P1", 400), rec("B", "P2", 400)], background={"iso1": [0.011], "iso2": [0.5]})
check(c["call"] == "contaminant" and not c["e5"] and c["e5_isolates"] == ["iso1", "iso2"],
      "E5 fails with 2 isolates at a tenth or more of the line's highest, also when far higher")
c = call([rec("A", "P1", 400), rec("B", "P2", 400)], background={"iso1": [0.5], "iso2": [0.009]})
check(c["call"] == "resident" and c["e5"], "E5 holds with 1 isolate at a similar level")

# routed reads from step 2's found.tsv.gz: routed only, and SingleM only in DNA runs (entry 7)
with tempfile.TemporaryDirectory() as td:
    fz = Path(td) / "found.tsv.gz"
    with gzip.open(fz, "wt") as fh:
        fh.write("mate\tposition\tread\tssuplex\tsinglem\toutcome\tdetail\n"
                 "m_1\t0\tr0\t1\t0\trouted\tno alignment\nm_1\t1\tr1\t0\t1\trouted\tno alignment\n"
                 "m_1\t2\tr2\t1\t1\tlineage\tCaGg\nm_2\t0\tr0\t1\t1\trouted\tno alignment passes the filter\n")
    dna, rna = T.routed_found({"route": "ssu+singlem"}, Path(td)), T.routed_found({"route": "ssu"}, Path(td))
    check(dna == {("m_1", "r0"): (0, True, False), ("m_1", "r1"): (1, False, True), ("m_2", "r0"): (0, True, True)}
          and rna == {("m_1", "r0"): (0, True, False), ("m_2", "r0"): (0, True, False)},
          f"routed reads: DNA {dna}, RNA {rna}")

# K1 on hand-made calls
base = [{"unit": T.B_PLUS, "lineage": "CaGg", "source": "mapped", "members": "", "call": "resident"},
        {"unit": T.B_MINUS, "lineage": "CaGg", "source": "mapped", "members": "", "call": "unresolved"},
        {"unit": T.PLANT, "lineage": "Ralstonia", "source": "mapped", "members": "", "call": "contaminant"}]
lvl = [{"lineage": "CaGg", "run": "R1", "level": "present"}, {"lineage": "CaGg", "run": "R2", "level": "present"}]
k1runs = [{"run": "R1", "control_id": "GMAR_BEG34_HIC_Bplus"}, {"run": "R2", "control_id": "GMAR_BEG34_GENOME"}]
check(T.k1(base, lvl, k1runs) == [], "K1 passes")
bad = base + [{"unit": T.B_MINUS, "lineage": "Mycoavidus (route)", "source": "route", "members": "Mycoavidus",
               "call": "resident"},
              {"unit": T.PLANT, "lineage": "novel (SSU)", "source": "route", "members": "", "call": "resident"}]
bad[0] = {**bad[0], "call": "unresolved"}
r = T.k1(bad, [lvl[0], {**lvl[1], "level": "low"}], k1runs)
check(len(r) == 4 and any("B- line" in x and "Mycoavidus" in x for x in r) and any("plant-only" in x for x in r) and
      any("B+ line" in x for x in r) and any("R2" in x for x in r), f"K1 reasons: {r}")

# ===================================================================== part 2 ===

MM2 = os.environ.get("MINIMAP2") or shutil.which("minimap2")
VS = os.environ.get("VSEARCH") or shutil.which("vsearch")
if not (MM2 and VS):
    print(("FAIL\n" + "\n".join(fails)) if fails else "PASS (part 1); SKIP part 2: minimap2 or VSEARCH not found")
    sys.exit(1 if fails else 0)

rng = random.Random(11)
COMP = str.maketrans("ACGT", "TGCA")


def rand(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def rc(s):
    return s.translate(COMP)[::-1]


def mutate(s, rate):
    return "".join(rng.choice("ACGT".replace(b, "")) if rng.random() < rate else b for b in s)


MOT = {k: rand(30) for k in ("K1", "K4", "K5", "KE")}


def motifs(g, places):
    s = list(g)
    for pos, k in places:
        s[pos:pos + 30] = MOT[k]
    return "".join(s)


S_C, S_N, S_P, S_O2, S_X = (rand(450) for _ in range(5))
G = {"host": rand(60000),
     "cagg": motifs(rand(5000) + S_C + rand(14550), [(12000, "K1")]),       # SSU at 5000-5450
     "cont": rand(10000),
     "org": rand(5000),
     "nov": motifs(rand(2000) + S_N + rand(9550), [(6000, "K4")]),          # no reference; SSU at 2000-2450
     "dug": motifs(rand(8000), [(3000, "K5")]),                             # no reference
     "pse": rand(1000) + S_P + rand(4550),                                  # no reference; SSU at 1000-1450
     "chl": rand(500) + S_O2 + rand(3050),                                  # no reference; SSU at 500-950
     "eux": motifs(rand(500) + S_X + rand(3050), [(2050, "KE")]),           # no reference
     "plant": rand(30000)}                                                  # no reference, nothing to find
ACC = {"host": "GCA_900000301.1", "cagg": "GCA_800000301.1", "cont": "GCF_700000301.1", "org": "NC_900301.1"}
SSU = [[S_C, "bacteria"], [S_N, "bacteria"], [S_P, "bacteria"], [S_O2, "bacteria"], [S_X, "bacteria"]]
SM_TAX = {"K1": "Root; d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria; o__Burkholderiales; "
                "f__Burkholderiaceae; g__Glomeribacter",
          "K4": "Root; d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria",
          "K5": "Root; d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria; o__Burkholderiales; "
                "f__Burkholderiaceae; g__Duganella",
          "KE": "Root; d__Eukaryota"}
SILVA = [("CAGG16S.1.450", ";".join(TAX["G1"][:7]) + ";Candidatus Glomeribacter gigasporarum", S_C),
         ("NOVREL.1.450", ";".join(TAX["O2"][:7]) + ";Massilia sp.", mutate(S_N, 0.15)),
         ("PSEA.1.450", "Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Pseudomonadales;"
                        "Pseudomonadaceae;Pseudomonas;Pseudomonas fluorescens", S_P),
         ("PSEB.1.450", "Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Pseudomonadales;"
                        "Pseudomonadaceae;Pseudomonas;Pseudomonas putida", mutate(S_P, 0.02)),
         ("CHL1.1.450", ";".join(TAX["CH"]), S_O2),
         ("EUKX.1.450", ";".join(TAX["EU"]), S_X)] + \
        [(f"DEC{i}.1.1200", ";".join(TAX["M1"][:7]) + f";decoy {i}", rand(1200)) for i in range(5)]


def spots(*groups):
    out = [src for src, k in groups for _ in range(k)]
    random.Random(len(out)).shuffle(out)
    return out


CHIM = "cagg@5000+host@30000"
RUNS = {"SRRP1": ("PAIRED", spots(("cagg", 1200), ("host", 700), ("nov@2000", 20), ("nov@5950", 15), ("dug@2950", 30),
                                  ("pse@1000", 10), ("chl@500", 5), ("eux@500", 5), ("eux@2000", 5), (CHIM, 10))),
        "SRRP2": ("PAIRED", spots(("cagg", 1000), ("host", 960), ("dug@2950", 30), ("nov@2000", 10))),
        "SRRP3": ("SINGLE", spots(("cagg", 2400), ("host", 1580), ("pse@1000", 20))),
        "SRRM1": ("PAIRED", spots(("cagg", 22), ("host", 1978))),
        "SRRM2": ("SINGLE", spots(("cagg", 45), ("host", 3955))),
        "SRRN1": ("PAIRED", spots(("org", 800), ("nov@2000", 20), ("pse@1000", 15), ("plant", 1165))),
        "SRRH1": ("SINGLE", spots(("cagg", 300), ("host", 200)), 3000)}     # HiFi: 3 kb reads
COLS = T.RUN_COLS
ROWS = [["SRRP1", T.B_PLUS, "GMAR_BEG34_HIC_Bplus", "yes", "DNA", "short", "dna", "PAIRED", "PRJ1", "", "L1",
         "BEG34", "ssu+singlem", "stage0a", "Hi-C"],
        ["SRRP2", T.B_PLUS, "GMAR_BEG34_GENOME", "yes", "DNA", "short", "dna", "PAIRED", "PRJ2", "", "L2", "BEG34",
         "ssu+singlem", "", "genome"],
        ["SRRP3", T.B_PLUS, "GMAR_BEG34_RNA_Bplus", "yes", "RNA", "short", "poly-A", "SINGLE", "PRJ3", "", "L2",
         "BEG34", "ssu", "", "spores"],
        ["SRRM1", T.B_MINUS, "GMAR_BEG34_HIC_Bminus", "yes", "DNA", "short", "dna", "PAIRED", "PRJ1", "", "L1",
         "BEG34", "ssu+singlem", "", "Hi-C"],
        ["SRRM2", T.B_MINUS, "GMAR_BEG34_RNA_Bminus", "yes", "RNA", "short", "poly-A", "SINGLE", "PRJ3", "", "L2",
         "BEG34", "ssu", "", "spores"],
        ["SRRN1", T.PLANT, "LOTUS_NM_ROOTS", "yes", "RNA", "short", "poly-A", "PAIRED", "PRJ3", "", "L2", "none",
         "ssu", "", "roots"],
        ["SRRH1", "BEG34 HiFi", "GMAR_BEG34_HIFI", "no", "DNA", "hifi", "dna", "SINGLE", "PRJ1", "", "L1", "BEG34",
         "none", "", "HiFi"]]
CONTROLS = {"GMAR_BEG34_HIC_Bplus": "resident", "GMAR_BEG34_GENOME": "resident", "GMAR_BEG34_RNA_Bplus": "resident",
            "GMAR_BEG34_HIC_Bminus": "absent", "GMAR_BEG34_RNA_Bminus": "absent",
            "LOTUS_NM_ROOTS": "absent (any bacterial lineage)", "GMAR_BEG34_HIFI": "reported, not scored"}
T_, M_ = 0.01, 0.25


def stand_in(p, body):
    p.write_text(f"#!{sys.executable}\n" + body)
    p.chmod(p.stat().st_mode | stat.S_IXUSR)
    return str(p)


FAKE_FASTQ_DUMP = '''import gzip, json, os, random, sys
a = sys.argv[1:]
if a == ["--version"]:
    print("\\nfastq-dump : 3.4.1\\n"); sys.exit(0)
spec = json.load(open(os.environ["FAKE_SPEC"]))
run, x, out = a[-1], int(a[a.index("-X") + 1]), a[a.index("-O") + 1]
open(os.environ["FAKE_LOG"], "a").write(f"fastq-dump {run}\\n")
layout, spots, *length = spec["runs"][run]
length = length[0] if length else 150
rng = random.Random(run)
comp = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(comp)[::-1]
os.makedirs(out, exist_ok=True)
names = [f"{run}_1.fastq.gz"] + ([f"{run}_2.fastq.gz"] if layout == "PAIRED" else [])
fhs = [gzip.open(os.path.join(out, n), "wt") for n in names]
g = spec["genomes"]
for i, spot in enumerate(spots[:x], 1):
    if "+" in spot:
        (sa, a_), (sb, b_) = (p.split("@") for p in spot.split("+"))
        a_, b_ = int(a_), int(b_)
        reads = [g[sa][a_:a_ + 75] + g[sb][b_:b_ + 75], rc(g[sb][b_ + 200:b_ + 350])]
    else:
        src, _, at = spot.partition("@")
        s = int(at) if at else rng.randint(0, len(g[src]) - length - 250)
        reads = [g[src][s:s + length], rc(g[src][s + 200:s + 350])]
    for fh, read in zip(fhs, reads):
        fh.write(f"@{run}.{i} {i} length={len(read)}\\n{read}\\n+{run}.{i} {i} length={len(read)}\\n{'I' * len(read)}\\n")
for fh in fhs:
    fh.close()
n = len(spots[:x])
sys.stderr.write(f"Read {n} spots for {run}\\nWritten {n} spots for {run}\\n")
'''

FAKE_SSUPLEX = '''import json, os, sys
a = sys.argv[1:]
if a == ["--version"]:
    print("ssuplex 0.2.0"); sys.exit(0)
opt = lambda k: a[a.index(k) + 1]
inp, prefix = opt("-i"), opt("-o")
open(os.environ["FAKE_LOG"], "a").write(f"ssuplex {os.path.basename(inp)}\\n")
spec = json.load(open(os.environ["FAKE_SPEC"]))
comp = str.maketrans("ACGT", "TGCA")
kmers = {}
for seq, origin in spec["ssu"]:
    for s in (seq, seq.translate(comp)[::-1]):
        for i in range(len(s) - 24):
            kmers[s[i:i + 25]] = origin
origins = ["bacteria", "archaea", "eukaryota", "mitochondria", "chloroplast", "unclassified"]
counts = dict.fromkeys(origins, 0)
fhs, name = {}, None
reads = []
for line in open(inp):
    line = line.rstrip("\\n")
    if line.startswith(">"):
        reads.append([line[1:].split()[0], ""])
    else:
        reads[-1][1] += line
for rid, seq in reads:
    origin = next((kmers[seq[i:i + 25]] for i in range(len(seq) - 24) if seq[i:i + 25] in kmers), "unclassified")
    counts[origin] += 1
    fh = fhs.get(origin) or fhs.setdefault(origin, open(f"{prefix}.{origin}.fasta", "w"))
    fh.write(f">{rid} origin={origin} score=50.0 coords=1-{len(seq)}\\n")
    for i in range(0, len(seq), 80):
        fh.write(seq[i:i + 80] + "\\n")
for fh in fhs.values():
    fh.close()
open(f"{prefix}.extraction.tsv", "w").write("read_id\\n")
with open(f"{prefix}.summary.txt", "w") as s:
    s.write(f"SSUplex v0.2.0 run summary\\n\\ntotal input reads:  {len(reads)}\\n\\nPer-origin counts:\\n")
    for o in origins:
        s.write(f"  {o:<14} {counts[o]:>10}\\n")
'''

FAKE_CONDA = '''import json, os, sys
a = sys.argv[1:]
rest = a[3:]
if rest[:1] == ["--no-capture-output"]:
    rest = rest[1:]
cmd = rest[1:]
if cmd == ["--version"]:
    print("0.21.4"); sys.exit(0)
opt = lambda k: cmd[cmd.index(k) + 1]
fwd, otu = opt("--forward"), opt("--otu-table")
open(os.environ["FAKE_LOG"], "a").write(f"singlem {os.path.basename(fwd)}\\n")
spec = json.load(open(os.environ["FAKE_SPEC"]))
comp = str.maketrans("ACGT", "TGCA")
hits, name = {}, None
for line in open(fwd):
    line = line.rstrip("\\n")
    if line.startswith(">"):
        name = line[1:].split()[0]
        continue
    for k, m in spec["motifs"].items():
        if m in line or m.translate(comp)[::-1] in line:
            hits.setdefault(k, []).append(name)
cols = ["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "read_names", "nucleotides_aligned",
        "taxonomy_by_known?", "read_unaligned_sequences", "equal_best_hit_taxonomies"]
with open(otu, "w") as out:
    out.write("\\t".join(cols) + "\\n")
    for k, names in sorted(hits.items()):
        out.write("\\t".join([k, "s", "ACG", str(len(names)), "1.00", spec["sm_tax"][k], " ".join(sorted(names)),
                              "60", "False", "", ""]) + "\\n")
'''

FAKE_SYLPH = '''import json, os, sys
a = sys.argv[1:]
if a == ["--version"]:
    print(os.environ.get("FAKE_SYLPH_VERSION", "sylph 1.0.0")); sys.exit(0)
if a[0] == "sketch":
    open(a[a.index("-o") + 1] + ".syldb", "w").write("x\\n"); sys.exit(0)
out = a[a.index("-o") + 1]
reads = a[a.index("-1") + 1] if "-1" in a else a[-1]
run = os.path.basename(os.path.dirname(reads))
open(os.environ["FAKE_LOG"], "a").write(f"sylph {run}\\n")
spec = json.load(open(os.environ["FAKE_SPEC"]))
with open(out, "w") as fh:
    fh.write("Sample_file\\tGenome_file\\tTaxonomic_abundance\\tSequence_abundance\\tAdjusted_ANI\\tEff_cov\\tContig_name\\n")
    for acc, ani in spec["sylph"].get(run, []):
        fh.write(f"{reads}\\t/x/clean/bacteria/{acc}.fna\\t50\\t50\\t{ani}\\t3.0\\tc1\\n")
'''


def quiet(fn, argv):
    err = io.StringIO()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
        code = fn(argv)
    return code, err.getvalue()


def write_tsv(p, cols, rows):
    p.write_text("\t".join(cols) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))


def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build(tmp):
    refdir, work = tmp / "references", tmp / "work"
    for sub in ("clean/hosts", "clean/bacteria", "db/silva", "db/hmm", "db/mp"):
        (work / sub).mkdir(parents=True)
    refdir.mkdir()
    (work / "clean" / "hosts" / f"{ACC['host']}.fna").write_text(f">{ACC['host']}|c1\n{G['host']}\n")
    for k in ("cagg", "cont"):
        (work / "clean" / "bacteria" / f"{ACC[k]}.fna").write_text(f">{ACC[k]}|c1\n{G[k]}\n")
    (work / "clean" / "organelles.fna").write_text(f">{ACC['org']}|{ACC['org']}\n{G['org']}\n")
    write_tsv(refdir / "hosts.tsv", ["accession", "organism", "note"], [[ACC["host"], "Gigaspora testii", "host"]])
    write_tsv(refdir / "endosymbionts.tsv", ["accession", "organism", "note"],
              [[ACC["cagg"], "Candidatus Glomeribacter testii", ""]])
    write_tsv(refdir / "contaminants.tsv", ["accession", "organism", "note"], [[ACC["cont"], "Ralstonia testii", ""]])
    write_tsv(refdir / "organelles.tsv", ["plant", "organelle", "accession", "title", "length_bp"],
              [["Lotus testii", "chloroplast", ACC["org"], "x", 5000]])
    (refdir / "overlaps.tsv").write_text("\t".join(["kind", "host_assembly", "host_organism", "host_contig",
                                                    "host_length", "host_start", "host_end", "host_covered",
                                                    "bacterial_genome", "bacterial_organism", "bacterial_contig",
                                                    "bacterial_length", "bacterial_start", "bacterial_end",
                                                    "bacterial_covered", "identity", "aligned_bp", "matched"]) + "\n")
    bin_ = tmp / "bin"
    bin_.mkdir()
    tools = {"fastq_dump": stand_in(bin_ / "fastq-dump", FAKE_FASTQ_DUMP),
             "ssuplex": stand_in(bin_ / "ssuplex", FAKE_SSUPLEX), "conda": stand_in(bin_ / "conda", FAKE_CONDA),
             "sylph": stand_in(bin_ / "sylph", FAKE_SYLPH)}
    stand_in(bin_ / "nhmmer", 'print("# HMMER 3.4 (Aug 2023); http://hmmer.org/")\n')
    os.environ["PATH"] = f"{bin_}{os.pathsep}{os.environ['PATH']}"
    os.environ["FAKE_SPEC"], os.environ["FAKE_LOG"] = str(tmp / "spec.json"), str(tmp / "tools.log")
    (tmp / "spec.json").write_text(json.dumps({
        "genomes": G, "runs": RUNS, "ssu": SSU, "motifs": {k: MOT[k] for k in SM_TAX}, "sm_tax": SM_TAX,
        "sylph": {r: [[ACC["cagg"], 99.12]] for r in ("SRRP1", "SRRP2", "SRRP3")}}))
    code, err = quiet(MI.main, ["--workdir", str(work), "--refdir", str(refdir), "--minimap2", MM2,
                                "--sylph", tools["sylph"], "--threads", "2"])
    if code != 0:
        print("FAIL\nmake_index: " + err)
        sys.exit(1)
    (work / "db" / "silva" / "silva.fasta").write_text("".join(f">{a} {p}\n{s}\n" for a, p, s in SILVA))
    write_tsv(work / "db" / "silva" / "taxonomy.tsv", ["accession", "taxonomy"], [[a, p] for a, p, _ in SILVA])
    (work / "db" / "silva" / "manifest.json").write_text(json.dumps({
        "file": "silva.fasta", "fasta": str(work / "db" / "silva" / "silva.fasta"), "vsearch": "2.32.0",
        "taxonomy": str(work / "db" / "silva" / "taxonomy.tsv")}))
    (work / "db" / "manifest.json").write_text(json.dumps({
        "versions": {"ssuplex": "0.2.0", "hmmer": "3.4", "singlem": "0.21.4"},
        "ssuplex_hmm_dir": str(work / "db" / "hmm"), "singlem_metapackage": str(work / "db" / "mp")}))
    stage = tmp / "stage0a"
    stage.mkdir()
    write_tsv(stage / "fractions.tsv", ["run", "sets_t", "fraction"], [["X1", "yes", T_ * 10], ["X2", "no", 0.0001]])
    write_tsv(stage / "marker_share_entry8.tsv", ["run", "sets_m", "share"], [["X1", "yes", M_], ["X2", "no", 0.01]])
    write_tsv(tmp / "controls.tsv", ["control_id", "expected_call"], list(CONTROLS.items()))
    return refdir, work, stage, tools


def args_for(tmp, refdir, work, stage, t1, tools, extra=()):
    return ["--workdir", str(work), "--refdir", str(refdir), "--t1-dir", str(t1), "--stage0a-dir", str(stage),
            "--controls", str(tmp / "controls.tsv"), "--threads", "2", "--jobs", "2",
            "--fastq-dump", tools["fastq_dump"], "--minimap2", MM2, "--ssuplex", tools["ssuplex"], "--conda",
            tools["conda"], "--vsearch", VS, "--sylph", tools["sylph"], *extra]


def log_lines(tmp):
    p = tmp / "tools.log"
    return p.read_text().splitlines() if p.exists() else []


with tempfile.TemporaryDirectory() as td:
    tmp = Path(td)
    refdir, work, stage, tools = build(tmp)
    t1 = tmp / "t1"
    t1.mkdir()

    # refusals
    write_tsv(t1 / "runs.tsv", COLS, ROWS)
    (stage / "marker_share_entry8.tsv").rename(stage / "hidden.tsv")
    code, err = quiet(T.main, args_for(tmp, refdir, work, stage, t1, tools))
    check(code == 1 and "entry 8" in err, f"no entry 8 m: refused ({err.strip()[:120]})")
    (stage / "hidden.tsv").rename(stage / "marker_share_entry8.tsv")
    write_tsv(t1 / "runs.tsv", COLS, [ROWS[2][:12] + ["ssu+singlem"] + ROWS[2][13:]])
    code, err = quiet(T.main, args_for(tmp, refdir, work, stage, t1, tools))
    check(code == 1 and "route must be" in err, f"RNA run with SingleM: refused ({err.strip()[:120]})")
    os.environ["FAKE_SYLPH_VERSION"] = "sylph 0.9.0"
    write_tsv(t1 / "runs.tsv", COLS, ROWS)
    code, err = quiet(T.main, args_for(tmp, refdir, work, stage, t1, tools))
    check(code == 1 and "sylph" in err, f"another sylph: refused ({err.strip()[:120]})")
    del os.environ["FAKE_SYLPH_VERSION"]
    check(not log_lines(tmp), "refusals call no tool")

    # stage 0a's files for SRRP1, made by the same functions, then moved where stage 0a keeps them
    pre = tmp / "pre_t1"
    pre.mkdir()
    write_tsv(pre / "runs.tsv", COLS, [ROWS[0][:13] + [""] + ROWS[0][14:]])
    code, err = quiet(T.main, args_for(tmp, refdir, work, stage, pre, tools))
    check(code == 0, f"making SRRP1's stage 0a files: {err[-300:]}")
    for sub in ("reads", "assign", "route"):
        (work / "stage0a" / sub).mkdir(parents=True)
        shutil.move(str(work / "t1" / sub / "SRRP1"), str(work / "stage0a" / sub / "SRRP1"))
    aj = work / "stage0a" / "assign" / "SRRP1" / "assign.json"
    aj.write_text(aj.read_text().replace(f"{work}/t1/reads/", f"{work}/stage0a/reads/"))
    shutil.rmtree(work / "t1")
    (tmp / "tools.log").unlink()

    # the T1 run
    write_tsv(t1 / "runs.tsv", COLS, ROWS)
    code, err = quiet(T.main, args_for(tmp, refdir, work, stage, t1, tools))
    check(code == 0, f"run_t1 failed: {err[-600:]}")
    first_log = log_lines(tmp)
    check(not [x for x in first_log if "SRRP1" in x and not x.startswith("sylph")],
          "the stage 0a run calls no tool but sylph: " + str([x for x in first_log if "SRRP1" in x]))
    check(not [x for x in first_log if x.startswith("singlem") and ("SRRP3" in x or "SRRM2" in x or "SRRN1" in x)],
          "SingleM is not run on RNA runs")
    check(not [x for x in first_log if ("ssuplex" in x or "singlem" in x) and "SRRH1" in x],
          "the HiFi run (route none) is mapped only")
    if code == 0:
        lv = {(x["run"], x["lineage"]): x for x in read_tsv(t1 / "levels.tsv")}
        calls = {(x["unit"], x["lineage"]): x for x in read_tsv(t1 / "calls.tsv")}
        md = (t1 / "T1.md").read_text()

        def reads_of(run):
            layout, sp = RUNS[run][:2]
            return len(sp) * (2 if layout == "PAIRED" else 1)

        def mates(run):
            return 2 if RUNS[run][0] == "PAIRED" else 1

        def k(run, src):
            return RUNS[run][1].count(src)

        # mapped CaGg, N, levels and the cross-talk cap
        for run in RUNS:
            x = lv[(run, "CaGg")]
            n = k(run, "cagg") * mates(run)
            N = T_ * reads_of(run)
            check(int(x["n"]) == n and abs(float(x["N"]) - N) < 1e-6 * N, f"{run} CaGg n/N: {x['n']}/{x['N']}, "
                                                                          f"want {n}/{N}")
        want_lv = {"SRRP1": "present", "SRRP2": "present", "SRRP3": "present", "SRRM1": "low", "SRRM2": "low",
                   "SRRN1": "absent", "SRRH1": "present"}
        got_lv = {r: lv[(r, "CaGg")]["level"] for r in RUNS}
        check(got_lv == want_lv, f"CaGg levels {got_lv}")
        check(lv[("SRRM1", "CaGg")]["cross_talk"] == "yes" and lv[("SRRM2", "CaGg")]["cross_talk"] == "yes" and
              lv[("SRRP1", "CaGg")]["cross_talk"] == "no", "cross-talk flags")
        check(lv[("SRRM1", "CaGg")]["informative"] == "yes" and lv[("SRRM2", "CaGg")]["informative"] == "no" and
              lv[("SRRN1", "CaGg")]["informative"] == "no", "informative: DNA yes, poly-A RNA no")

        # route lineages: counts worked out from the planted reads
        want_route = {("SRRP1", "novel (SSU)"): 2 * k("SRRP1", "nov@2000"),
                      ("SRRP2", "novel (SSU)"): 2 * k("SRRP2", "nov@2000"),
                      ("SRRN1", "novel (SSU)"): 2 * k("SRRN1", "nov@2000"),
                      ("SRRP1", "novel (SingleM)"): k("SRRP1", "nov@5950"),
                      ("SRRP1", "Duganella (route)"): k("SRRP1", "dug@2950"),
                      ("SRRP2", "Duganella (route)"): k("SRRP2", "dug@2950"),
                      ("SRRP1", "Pseudomonas (route)"): 2 * k("SRRP1", "pse@1000"),
                      ("SRRP3", "Pseudomonas (route)"): k("SRRP3", "pse@1000"),
                      ("SRRN1", "Pseudomonas (route)"): 2 * k("SRRN1", "pse@1000")}
        route_lins = {"novel (SSU)", "novel (SingleM)", "Duganella (route)", "Pseudomonas (route)"}
        got_lins = {x["lineage"] for x in read_tsv(t1 / "levels.tsv") if x["source"] == "route"}
        check(got_lins == route_lins, f"route lineages {sorted(got_lins)}")
        for (run, lin), x in lv.items():
            if lin in route_lins:
                check(int(x["n"]) == want_route.get((run, lin), 0), f"{run} {lin}: n {x['n']}, "
                                                                    f"want {want_route.get((run, lin), 0)}")
                check(abs(float(x["N"]) - T_ * M_ * reads_of(run)) < 1e-9, f"{run} {lin}: route N")
        check(("SRRN1", "novel (SingleM)") not in lv and ("SRRH1", "novel (SSU)") not in lv,
              "SingleM lineages are not measured in RNA runs, route lineages not in mapped-only runs")
        check(int(lv[("SRRP1", "CaGg")]["split_reads"]) == k("SRRP1", CHIM) and
              int(lv[("SRRP2", "CaGg")]["split_reads"]) == 0, "split reads of the Hi-C-like chimera")
        check(lv[("SRRP1", "CaGg")]["sylph_ani"] == "99.12" and lv[("SRRM1", "CaGg")]["sylph_ani"] == "",
              "sylph ANI")

        # calls
        want_calls = {(T.B_PLUS, "CaGg"): "resident", (T.B_MINUS, "CaGg"): "unresolved", (T.PLANT, "CaGg"): "no call",
                      (T.B_PLUS, "Ralstonia"): "absent", (T.PLANT, "Ralstonia"): "no call",
                      (T.B_PLUS, "Duganella (route)"): "resident", (T.B_PLUS, "novel (SSU)"): "resident",
                      (T.B_PLUS, "novel (SingleM)"): "unresolved", (T.B_PLUS, "Pseudomonas (route)"): "contaminant",
                      (T.PLANT, "Pseudomonas (route)"): "contaminant", (T.PLANT, "novel (SSU)"): "unresolved",
                      (T.B_MINUS, "novel (SSU)"): "absent", (T.B_MINUS, "Pseudomonas (route)"): "absent",
                      ("BEG34 HiFi", "CaGg"): "unresolved"}
        got_calls = {key: calls[key]["call"] if key in calls else "missing" for key in want_calls}
        check(got_calls == want_calls, "calls: " + str({a: b for a, b in got_calls.items() if want_calls[a] != b}))
        check(calls[(T.B_PLUS, "Duganella (route)")]["category"] == "other known genus" and
              calls[(T.B_PLUS, "Pseudomonas (route)")]["category"] == "kit-contaminant genus" and
              calls[(T.B_PLUS, "novel (SSU)")]["category"] == "novel", "categories")
        check("family Oxalobacteraceae" in calls[(T.B_PLUS, "novel (SSU)")]["placements"] and
              "class Gammaproteobacteria" in calls[(T.B_PLUS, "novel (SingleM)")]["placements"],
              "novel placements: " + calls[(T.B_PLUS, "novel (SSU)")]["placements"])
        check((T.PLANT, "novel (SingleM)") not in calls and ("BEG34 HiFi", "novel (SSU)") not in calls,
              "no call row where a lineage is not measured")

        # the route table, K1, the confusion table
        p1 = next(x for x in md.splitlines() if x.startswith("| SRRP1 |") and x.count("|") == 7)
        ssu_found = 2 * sum(k("SRRP1", s) for s in ("nov@2000", "pse@1000", "chl@500", "eux@500")) + k("SRRP1", CHIM)
        sm_found = sum(k("SRRP1", s) for s in ("nov@5950", "dug@2950", "eux@2000"))
        want_p1 = (f"| SRRP1 | {ssu_found + sm_found} | {ssu_found - 4 * 5} | 10 / 10 | "
                   f"{k('SRRP1', 'nov@5950') + k('SRRP1', 'dug@2950')} | {k('SRRP1', CHIM)} / 0 |")
        check(p1 == want_p1, f"route table for SRRP1: {p1}, want {want_p1}")
        check("## K1: passed" in md, "K1 passes")
        check("| BEG34 B- | yes | absent | unresolved |" in md, "confusion table row for the B- line")

        # a rerun calls no tool and gives the same tables
        before = {n: (t1 / n).read_text() for n in ("levels.tsv", "calls.tsv", "T1.md")}
        code, err = quiet(T.main, args_for(tmp, refdir, work, stage, t1, tools))
        check(code == 0 and log_lines(tmp) == first_log, "a rerun calls no tool")
        after = {n: (t1 / n).read_text() for n in before}
        check(after["levels.tsv"] == before["levels.tsv"] and after["calls.tsv"] == before["calls.tsv"] and
              after["T1.md"].split("\n", 3)[3] == before["T1.md"].split("\n", 3)[3], "a rerun gives the same tables")

if fails:
    print("FAIL\n" + "\n".join(fails))
    sys.exit(1)
print("PASS: caller rules (names, SSU and SingleM labels, one count per read, joining, levels, cross-talk, no call, "
      "absent, contaminant by E5 and kit genus, resident by known clade and by E1 and E2, DNA-only E1 and E2, "
      "years, K1); run_t1 end to end (mapped and route counts, N, cross-talk cap, split reads, SSU and SingleM "
      "drops, every call, placements, sylph, K1, the stage 0a run reused, no SingleM on RNA, HiFi mapped only, "
      "rerun calls no tool, refusals)")
