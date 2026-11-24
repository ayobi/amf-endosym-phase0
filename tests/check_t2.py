"""Check make_t2_silva.py and run_t2.py (PROTOCOL.md section 5, T2; section 6, K2; section 10, entry 11).

Needs minimap2 and VSEARCH (activate phase0-reads, or set MINIMAP2= and VSEARCH=).

Part 1, the SILVA removal rules on hand-made headers: each taxon name in a rank or in the name, in any
case; Moeniiplasma and Mycoplasma-related (followed by any word) in the name; MRE in capitals, including
inside a name such as GvMRE, but not AMRE1 or MREB1.

Part 2, end to end on synthetic genomes. Stage 0a is made by run_t1.py with stand-ins for fastq-dump,
SSUplex and SingleM, and real minimap2, assign_reads.py, make_index.py and VSEARCH; its folders are then
moved where stage 0a keeps them. Then make_t2_silva.py and run_t2.py. Planted:
  Iso A (2 runs)  an MRE whose SSU is in SILVA under Mollicutes, plus a relative at 85% under Bacillus:
                  held out, its reads go to novel (SSU), called unresolved (K2 passes); kept in, CaMg/MRE
                  is resident. A bacterium without a reference whose SSU is in SILVA twice, at 100%, once
                  in Mollicutes and once in Bacillus: unplaced when kept in, "Bacillus (route)" when held
                  out (the placement that follows from all of SILVA, with the removed hit left out).
  Iso B (1 run)   an MRE filed in SILVA under another class but named as one: its SSU reads go to novel
                  (SSU), its SingleM reads to novel (SingleM); unresolved (K2 passes).
  Iso C (1 run)   the MRE of Iso A at a low level: novel (SSU) is absent (K2 fails).
  Iso D (1 run)   an MRE whose reads no tool finds: no lineage holds them (K2 fails).
Checked: the removal list and the T2 database; the held-out arm assigned with the T2 index and holding no
CaMg/MRE lineage; placements that follow from all of SILVA, searched again only when every top hit was
removed, and checked against a direct search (a wrong derivation stops the run); K2 per isolate, with the
lineage that holds most of the MRE's reads and the others; both arms in levels.tsv and calls.tsv; a rerun
aligns and places nothing again; with Iso A and B only, K2 passes; refusals (a removal list that is not
the one the T2 database was made from, a T2 database whose FASTA was changed).
"""
import contextlib
import csv
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
import make_index as MI  # noqa: E402
import make_t2_silva as MS  # noqa: E402
import run_t1 as T  # noqa: E402
import run_t2 as RT  # noqa: E402

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)


# ----------------------------------------------------------- part 1: the rules ---

B = "Bacteria;Bacillati;"
cases = {f"{B}Mycoplasmatota;Mollicutes;Mycoplasmatales;Mycoplasmataceae;Mycoplasma;Mycoplasma sp.":
         ["Mollicutes", "Mycoplasmatota", "Mycoplasmatales"],
         f"{B}Bacillota;Bacilli;Bacillales;Bacillaceae;uncultured;uncultured mollicutes bacterium": ["Mollicutes"],
         f"{B}Tenericutes;x;y;z;Acholeplasma;Acholeplasma sp.": ["Tenericutes"],
         f"{B}Bacillota;Bacilli;Entomoplasmatales--other;x;y;z": ["Entomoplasmatales"],
         f"{B}a;b;Anaeroplasmatales;c;d;e": ["Anaeroplasmatales"],
         f"{B}a;b;Haloplasmatales;c;d;e": ["Haloplasmatales"],
         f"{B}a;b;Acholeplasmatales;c;d;e": ["Acholeplasmatales"],
         f"{B}a;b;c;d;uncultured;Candidatus Moeniiplasma glomeromycotorum": ["Moeniiplasma"],
         f"{B}a;b;c;d;uncultured;Mycoplasma-related endobacterium of Dentiscutata": ["Mycoplasma-related"],
         f"{B}a;b;c;d;uncultured;mycoplasma related endobacteria clone 3": ["Mycoplasma-related"],
         f"{B}a;b;c;d;uncultured;Mycoplasma-related endosymbiont of Gigaspora": ["Mycoplasma-related"],
         f"{B}a;b;c;d;Mycoplasma;Mycoplasma sp. related to M. hominis": [],
         f"{B}a;b;c;d;uncultured;endosymbiont GvMRE of Glomus versiforme": ["MRE"],
         f"{B}a;b;c;d;uncultured;MRE-2 of Gigaspora": ["MRE"],
         f"{B}a;b;c;d;Ralstonia;Ralstonia sp. AMRE1": [],
         f"{B}a;b;c;d;Bacillus;Bacillus sp. MREB1": [],
         f"{B}a;b;c;d;Bacillus;Bacillus sp. mre1": [],
         "Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Burkholderiales;Burkholderiaceae;Ralstonia;R.": []}
for path, want in cases.items():
    got = MS.matches(path)
    check(got == want, f"matches({path.split(';')[-1]!r}): {got}, want {want}")

MM2 = os.environ.get("MINIMAP2") or shutil.which("minimap2")
VS = os.environ.get("VSEARCH") or shutil.which("vsearch")
if not (MM2 and VS):
    print(("FAIL\n" + "\n".join(fails)) if fails else "PASS (part 1); SKIP part 2: minimap2 or VSEARCH not found")
    sys.exit(1 if fails else 0)

# --------------------------------------------------------------- part 2: world ---

rng = random.Random(29)
COMP = str.maketrans("ACGT", "TGCA")


def rand(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(s, rate):
    return "".join(rng.choice("ACGT".replace(b, "")) if rng.random() < rate else b for b in s)


S_A, S_B, S_N, S_X, K_B = rand(450), rand(450), rand(450), rand(150), rand(30)
mre_b = rand(3000) + S_B + rand(16550)
mre_b = mre_b[:9000] + K_B + mre_b[9030:]                                 # a SingleM marker at 9000
G = {"host": rand(60000),
     "mreA": rand(3000) + S_A + rand(16550),                              # SSU at 3000-3450
     "mreB": mre_b,                                                       # SSU at 3000-3450
     "cont": rand(10000), "org": rand(5000), "mito": rand(6000),
     "unk": rand(1000) + S_X + rand(1000),                                # no reference; SSU at 1000-1150
     "nov": rand(2000) + S_N + rand(1000)}                                # no reference; SSU at 2000-2450
ACC = {"host": "GCA_900000401.1", "mreA": "GCA_800000401.1", "mreB": "GCA_800000402.1", "cont": "GCF_700000401.1",
       "org": "NC_900401.1", "mito": "MN900401.1"}
SSU = [[S_A, "bacteria"], [S_B, "bacteria"], [S_N, "bacteria"], [S_X, "bacteria"]]
SM_TAX = {"KB": "Root; d__Bacteria; p__Bacillota; c__Bacilli"}
MOL = "Bacteria;Bacillati;Mycoplasmatota;Mollicutes;Mycoplasmatales;Mycoplasmataceae"
BAC = "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae"
SILVA = [("MREA.1.450", f"{MOL};Candidatus Moeniiplasma;Candidatus Moeniiplasma glomeromycotorum", S_A),
         ("MREC.1.450", "Bacteria;Bacillati;Bacillota;Bacilli;Lactobacillales;uncultured;uncultured;"
                        "endosymbiont GvMRE of Glomus versiforme", mutate(S_A, 0.01)),
         ("MREB.1.450", "Bacteria;Bacillati;Bacillota;Bacilli;Lactobacillales;Streptococcaceae;uncultured;"
                        "Mycoplasma-related endobacterium of Gigaspora", S_B),
         ("RELA.1.450", f"{BAC};Bacillus;Bacillus sp. AMRE1", mutate(S_A, 0.15)),
         ("RELB.1.450", "Bacteria;Bacillati;Bacillota;Erysipelotrichia;Erysipelotrichales;Erysipelotrichaceae;"
                        "Erysipelothrix;Erysipelothrix sp. MREB1", mutate(S_B, 0.15)),
         ("MIXR.1.350", f"{MOL};Mycoplasma;Mycoplasma sp.", rand(100) + S_X + rand(100)),
         ("MIXK.1.350", f"{BAC};Bacillus;Bacillus sp.", rand(100) + S_X + rand(100))] + \
        [(f"DEC{i}.1.900", f"{BAC};Bacillus;decoy {i}", rand(900)) for i in range(4)] + \
        [("SPIR.1.900", f"{B}Bacillota;Bacilli;Spiroplasmales;Spiroplasmataceae;Spiroplasma;Spiroplasma sp.",
          rand(900))]                                                        # a gap the report's check must show
REMOVED = {"MREA.1.450", "MREC.1.450", "MREB.1.450", "MIXR.1.350"}


def mre_spots(src, lo, hi, step):
    return [f"{src}@{p}" for p in range(lo, hi, step)]


def spots(*groups):
    out = [x for g in groups for x in g]
    random.Random(len(out)).shuffle(out)
    return out


def rep(x, k):
    return [x] * k


RUNS = {"SRRA1": ("PAIRED", spots(mre_spots("mreA", 4000, 19500, 50), rep("mreA@3000", 15), rep("unk@1000", 20),
                                  rep("nov@2000", 8), rep("host", 1647))),
        "SRRA2": ("PAIRED", spots(mre_spots("mreA", 4025, 19500, 60), rep("mreA@3010", 12), rep("host", 1730))),
        "SRRB1": ("PAIRED", spots(mre_spots("mreB", 9500, 19500, 40), mre_spots("mreB", 3600, 8500, 40),
                                  rep("mreB@3000", 15), rep("mreB@8990", 15), rep("host", 1597))),
        "SRRC1": ("PAIRED", spots(rep("mreA@10000", 30), rep("mreA@3100", 1), rep("host", 7969))),
        "SRRD1": ("PAIRED", spots(mre_spots("mreB", 12000, 19000, 70), rep("host", 1900)))}
UNIT = {"SRRA1": ("Iso A", "CTRL_A", "PRJA", "L1", "A"), "SRRA2": ("Iso A", "CTRL_A", "PRJA", "L1", "A"),
        "SRRB1": ("Iso B", "CTRL_B", "PRJB", "L2", "B"), "SRRC1": ("Iso C", "CTRL_C", "PRJB", "L2", "C"),
        "SRRD1": ("Iso D", "CTRL_D", "PRJD", "L3", "D")}
T_, M_ = 0.01, 0.25


def rows(reuse, only=None):
    return [[run, u, c, "yes", "DNA", "short", "dna", "PAIRED", p, "2017", lab, iso, "ssu+singlem", reuse, "plate"]
            for run, (u, c, p, lab, iso) in UNIT.items() if only is None or u in only]


def stand_in(p, body):
    p.write_text(f"#!{sys.executable}\n" + body)
    p.chmod(p.stat().st_mode | stat.S_IXUSR)
    return str(p)


def quiet(fn, argv):
    err = io.StringIO()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
        code = fn(argv)
    return code, err.getvalue()


def write_tsv(p, cols, rows_):
    p.write_text("\t".join(cols) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows_))


def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


FAKE_FASTQ_DUMP = '''import gzip, json, os, random, sys
a = sys.argv[1:]
if a == ["--version"]:
    print("\\nfastq-dump : 3.4.1\\n"); sys.exit(0)
spec = json.load(open(os.environ["FAKE_SPEC"]))
run, x, out = a[-1], int(a[a.index("-X") + 1]), a[a.index("-O") + 1]
layout, spots = spec["runs"][run]
rng = random.Random(run)
comp = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(comp)[::-1]
os.makedirs(out, exist_ok=True)
names = [f"{run}_1.fastq.gz"] + ([f"{run}_2.fastq.gz"] if layout == "PAIRED" else [])
fhs = [gzip.open(os.path.join(out, n), "wt") for n in names]
g = spec["genomes"]
for i, spot in enumerate(spots[:x], 1):
    src, _, at = spot.partition("@")
    s = int(at) if at else rng.randint(0, len(g[src]) - 400)
    for fh, read in zip(fhs, [g[src][s:s + 150], rc(g[src][s + 200:s + 350])]):
        fh.write(f"@{run}.{i} {i} length=150\\n{read}\\n+{run}.{i} {i} length=150\\n{'I' * 150}\\n")
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
spec = json.load(open(os.environ["FAKE_SPEC"]))
comp = str.maketrans("ACGT", "TGCA")
kmers = {}
for seq, origin in spec["ssu"]:
    for s in (seq, seq.translate(comp)[::-1]):
        for i in range(len(s) - 24):
            kmers[s[i:i + 25]] = origin
origins = ["bacteria", "archaea", "eukaryota", "mitochondria", "chloroplast", "unclassified"]
counts = dict.fromkeys(origins, 0)
fhs, reads = {}, []
for line in open(inp):
    line = line.rstrip("\\n")
    if line.startswith(">"):
        reads.append([line[1:].split()[0], ""])
    else:
        reads[-1][1] += line
for rid, seq in reads:
    at = [i for i in range(len(seq) - 24) if seq[i:i + 25] in kmers]
    origin = kmers[seq[at[0]:at[0] + 25]] if at else "unclassified"
    lo, hi = (at[0], at[-1] + 25) if at else (0, len(seq))
    counts[origin] += 1
    fh = fhs.get(origin) or fhs.setdefault(origin, open(f"{prefix}.{origin}.fasta", "w"))
    fh.write(f">{rid} origin={origin} score=50.0 coords={lo + 1}-{hi}\\n")
    for i in range(lo, hi, 80):
        fh.write(seq[i:min(i + 80, hi)] + "\\n")
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

FAKE_SYLPH = '''import os, sys
a = sys.argv[1:]
if a == ["--version"]:
    print("sylph 1.0.0"); sys.exit(0)
if a[0] == "sketch":
    open(a[a.index("-o") + 1] + ".syldb", "w").write("x\\n"); sys.exit(0)
cols = ["Sample_file", "Genome_file", "Taxonomic_abundance", "Sequence_abundance", "Adjusted_ANI", "Eff_cov",
        "Contig_name"]
with open(a[a.index("-o") + 1], "w") as fh:
    fh.write("\\t".join(cols) + "\\n")
'''

WRAP = '''import os, sys
real, a = {real!r}, sys.argv[1:]
if a != ["--version"]:
    open(os.environ["FAKE_CALLS"], "a").write({name!r} + "\\n")
os.execv(real, [real] + a)
'''


def build(tmp):
    refdir, work, stage, pre = tmp / "references", tmp / "work", tmp / "stage0a", tmp / "pre"
    for sub in ("clean/hosts", "clean/bacteria", "db/silva", "db/hmm", "db/mp"):
        (work / sub).mkdir(parents=True)
    for d in (refdir, stage, pre):
        d.mkdir()
    (work / "clean" / "hosts" / f"{ACC['host']}.fna").write_text(f">{ACC['host']}|c1\n{G['host']}\n")
    for k in ("mreA", "mreB", "cont"):
        (work / "clean" / "bacteria" / f"{ACC[k]}.fna").write_text(f">{ACC[k]}|c1\n{G[k]}\n")
    (work / "clean" / "organelles.fna").write_text(f">{ACC['org']}|{ACC['org']}\n{G['org']}\n"
                                                   f">{ACC['mito']}|{ACC['mito']}\n{G['mito']}\n")
    write_tsv(refdir / "hosts.tsv", ["accession", "organism", "note"], [[ACC["host"], "Rhizophagus testii", "host"]])
    write_tsv(refdir / "endosymbionts.tsv", ["accession", "organism", "note"],
              [[ACC["mreA"], "Candidatus Moeniiplasma testii A", "held_out_in_T2"],
               [ACC["mreB"], "Candidatus Moeniiplasma testii B", "held_out_in_T2"]])
    write_tsv(refdir / "contaminants.tsv", ["accession", "organism", "note"], [[ACC["cont"], "Ralstonia testii", ""]])
    write_tsv(refdir / "organelles.tsv", ["plant", "organelle", "accession", "title", "length_bp"],
              [["Lotus testii", "chloroplast", ACC["org"], "x", 5000]])
    write_tsv(refdir / "fungal_mitochondria.tsv", ["fungus", "organelle", "accession", "title", "length_bp",
                                                   "source", "query"],
              [["Rhizophagus testii", "mitochondrion", ACC["mito"], "x", 6000, "GenBank", "q"]])
    (refdir / "overlaps.tsv").write_text("\t".join(["kind", "host_assembly", "host_organism", "host_contig",
                                                    "host_length", "host_start", "host_end", "host_covered",
                                                    "bacterial_genome", "bacterial_organism", "bacterial_contig",
                                                    "bacterial_length", "bacterial_start", "bacterial_end",
                                                    "bacterial_covered", "identity", "aligned_bp", "matched"]) + "\n")
    bin_ = tmp / "bin"
    bin_.mkdir()
    tools = {"fastq_dump": stand_in(bin_ / "fastq-dump", FAKE_FASTQ_DUMP),
             "ssuplex": stand_in(bin_ / "ssuplex", FAKE_SSUPLEX), "conda": stand_in(bin_ / "conda", FAKE_CONDA),
             "sylph": stand_in(bin_ / "sylph", FAKE_SYLPH),
             "mm2": stand_in(bin_ / "mm2wrap", WRAP.format(real=MM2, name="minimap2")),
             "vsearch": stand_in(bin_ / "vswrap", WRAP.format(real=VS, name="vsearch"))}
    stand_in(bin_ / "nhmmer", 'print("# HMMER 3.4 (Aug 2023); http://hmmer.org/")\n')
    os.environ["PATH"] = f"{bin_}{os.pathsep}{os.environ['PATH']}"
    os.environ["FAKE_SPEC"], os.environ["FAKE_LOG"] = str(tmp / "spec.json"), str(tmp / "tools.log")
    os.environ["FAKE_CALLS"] = str(tmp / "calls.log")
    (tmp / "spec.json").write_text(json.dumps({"genomes": G, "runs": RUNS, "ssu": SSU, "motifs": {"KB": K_B},
                                               "sm_tax": SM_TAX}))
    code, err = quiet(MI.main, ["--workdir", str(work), "--refdir", str(refdir), "--minimap2", MM2,
                                "--sylph", tools["sylph"], "--threads", "2"])
    if code != 0:
        print("FAIL\nmake_index: " + err)
        sys.exit(1)
    fasta = work / "db" / "silva" / "silva.fasta"
    fasta.write_text("".join(f">{a} {p}\n{s}\n" for a, p, s in SILVA))
    write_tsv(work / "db" / "silva" / "taxonomy.tsv", ["accession", "taxonomy"], [[a, p] for a, p, _ in SILVA])
    (work / "db" / "silva" / "manifest.json").write_text(json.dumps({
        "file": "silva.fasta", "md5": "0123", "sequences": len(SILVA), "fasta": str(fasta), "vsearch": "2.32.0",
        "taxonomy": str(work / "db" / "silva" / "taxonomy.tsv")}))
    (work / "db" / "manifest.json").write_text(json.dumps({
        "versions": {"ssuplex": "0.2.0", "hmmer": "3.4", "singlem": "0.21.4"},
        "ssuplex_hmm_dir": str(work / "db" / "hmm"), "singlem_metapackage": str(work / "db" / "mp")}))
    for name in ("fractions_entry10.tsv", "marker_share_entry10.tsv"):
        col, val = ("fraction", T_ * 10) if name.startswith("fr") else ("share", M_)
        write_tsv(stage / name, ["run", "sets_t" if col == "fraction" else "sets_m", col], [["X1", "yes", val]])
    write_tsv(tmp / "controls.tsv", ["control_id", "expected_call"],
              [[c, "resident"] for c in ("CTRL_A", "CTRL_B", "CTRL_C", "CTRL_D")])

    # Stage 0a, made by run_t1.py and moved where stage 0a keeps its files.
    write_tsv(pre / "runs.tsv", T.RUN_COLS, rows(""))
    code, err = quiet(T.main, ["--workdir", str(work), "--refdir", str(refdir), "--t1-dir", str(pre),
                               "--stage0a-dir", str(stage), "--controls", str(tmp / "controls.tsv"), "--threads", "2",
                               "--jobs", "2", "--fastq-dump", tools["fastq_dump"], "--minimap2", MM2,
                               "--ssuplex", tools["ssuplex"], "--conda", tools["conda"], "--vsearch", VS,
                               "--sylph", tools["sylph"]])
    if code != 0:
        print("FAIL\nmaking stage 0a with run_t1.py: " + err)
        sys.exit(1)
    for src, dest in (("t1/reads", "stage0a/reads"), ("t1_entry10/assign", "stage0a_entry10/assign"),
                      ("t1_entry10/route", "stage0a/route")):
        (work / dest).mkdir(parents=True, exist_ok=True)
        for run in RUNS:
            shutil.move(str(work / src / run), str(work / dest / run))
    w = work.resolve()
    for aj in (work / "stage0a_entry10" / "assign").glob("*/assign.json"):
        aj.write_text(aj.read_text().replace(f"{w}/t1/reads/", f"{w}/stage0a/reads/"))
    shutil.rmtree(work / "t1_entry10")       # work/t1/place stays: T1's placements, used again by run_t2.py
    return refdir, work, stage, tools


def counts_of(p):
    return {(r["section"], r["name"]): int(r["reads"]) for r in read_tsv(p)}


with tempfile.TemporaryDirectory() as td:
    tmp = Path(td)
    refdir, work, stage, tools = build(tmp)
    t2 = tmp / "t2"
    t2.mkdir()

    # make_t2_silva.py
    code, err = quiet(MS.main, ["--workdir", str(work), "--t2-dir", str(t2), "--vsearch", tools["vsearch"]])
    check(code == 0, f"make_t2_silva failed: {err}")
    removed = {r["sequence"]: r["matched"] for r in read_tsv(t2 / "silva_removed.tsv")}
    check(set(removed) == REMOVED, f"removed: {sorted(removed)}, want {sorted(REMOVED)}")
    check(removed.get("MREB.1.450") == "Mycoplasma-related" and removed.get("MREC.1.450") == "MRE",
          f"matched terms: {removed}")
    man = json.loads((work / "db" / "silva_t2" / "manifest.json").read_text())
    kept = [line[1:].split()[0] for line in open(man["fasta"]) if line.startswith(">")]
    check(man["removed"] == len(REMOVED) and man["sequences"] == len(SILVA) - len(REMOVED) and
          set(kept) == {a for a, _, _ in SILVA} - REMOVED, f"T2 database: {man}, {kept}")
    md = (t2 / "SILVA_T2.md").read_text()
    check("**Removed: 4 sequences**" in md and "| Mollicutes | 2 |" in md and "| MRE | 1 |" in md and
          "Mycoplasmopsis, Phytoplasma, Spiroplasma, Ureaplasma: **1**." in md and "| Spiroplasma | 1 |" in md and
          "| SPIR.1.900 |" in md, f"SILVA_T2.md: {md}")
    bad = json.loads((work / "db" / "silva" / "manifest.json").read_text())
    (work / "db" / "silva" / "manifest.json").write_text(json.dumps({**bad, "sequences": 99}))
    t2x = tmp / "t2x"
    code, err = quiet(MS.main, ["--workdir", str(work), "--t2-dir", str(t2x), "--vsearch", tools["vsearch"]])
    check(code == 1 and not (work / "db" / "silva_t2" / "manifest.json").exists() and
          not (t2x / "silva_removed.tsv").exists(), f"a wrong sequence count: exit {code}, {err!r}")
    (work / "db" / "silva" / "manifest.json").write_text(json.dumps(bad))
    code, err = quiet(MS.main, ["--workdir", str(work), "--t2-dir", str(t2), "--vsearch", tools["vsearch"]])
    check(code == 0, f"make_t2_silva again: {err}")

    # run_t2.py
    write_tsv(t2 / "runs.tsv", T.RUN_COLS, rows("stage0a"))

    def run_t2(t2_dir=t2):
        return quiet(RT.main, ["--workdir", str(work), "--refdir", str(refdir), "--t2-dir", str(t2_dir),
                               "--stage0a-dir", str(stage), "--controls", str(tmp / "controls.tsv"), "--threads", "2",
                               "--minimap2", tools["mm2"], "--vsearch", tools["vsearch"], "--sylph", tools["sylph"]])

    calls_log = tmp / "calls.log"
    calls_log.unlink(missing_ok=True)
    code, err = run_t2()
    check(code == 0, f"run_t2 failed: {err[-2000:]}")
    if code == 0:
        first_calls = calls_log.read_text().split()
        aj = json.loads((work / "t2" / "assign" / "SRRA1" / "assign.json").read_text())
        check(aj["index"]["set"] == "t2", f"held-out assignment index: {aj['index']}")
        lv = {(x["arm"], x["run"], x["lineage"]): x for x in read_tsv(t2 / "levels.tsv")}
        cl = {(x["arm"], x["unit"], x["lineage"]): x for x in read_tsv(t2 / "calls.tsv")}
        check(not any(k[0] == "held_out" and k[2] == "CaMg/MRE" for k in lv) and
              any(k[0] == "kept_in" and k[2] == "CaMg/MRE" for k in lv), "CaMg/MRE: kept in only")
        check(cl[("kept_in", "Iso A", "CaMg/MRE")]["call"] == "resident", "kept in: Iso A CaMg/MRE resident")
        got = {u: cl[("held_out", u, "novel (SSU)")]["call"] for u in ("Iso A", "Iso B", "Iso C")}
        check(got == {"Iso A": "unresolved", "Iso B": "unresolved", "Iso C": "absent"},
              f"held out: novel (SSU) calls {got}")
        check("family Bacillaceae" in cl[("held_out", "Iso A", "novel (SSU)")]["placements"] and
              "family Erysipelotrichaceae" in cl[("held_out", "Iso B", "novel (SSU)")]["placements"],
              f"held out: novel (SSU) placements {cl[('held_out', 'Iso A', 'novel (SSU)')]['placements']}, "
              f"{cl[('held_out', 'Iso B', 'novel (SSU)')]['placements']}")
        check(("held_out", "Iso A", "Bacillus (route)") in cl and ("kept_in", "Iso A", "Bacillus (route)") not in cl,
              "the tied SSU is unplaced when kept in and Bacillus when held out")
        check(cl[("held_out", "Iso B", "novel (SingleM)")]["call"] == "unresolved", "held out: Iso B novel (SingleM)")
        # MRE read tracking, from the stage 0a counts and the planted reads
        kc = {r: counts_of(work / "stage0a_entry10" / "assign" / r / "counts.tsv") for r in RUNS}
        md = (t2 / "T2.md").read_text()
        a_mre = kc["SRRA1"][("lineage", "CaMg/MRE")] + kc["SRRA2"][("lineage", "CaMg/MRE")]
        row_a = next(x for x in md.splitlines() if x.startswith("| Iso A |"))
        check(row_a == f"| Iso A | {a_mre} | {2 * (15 + 12)} | {2 * (15 + 12)} | novel (SSU) | {2 * (15 + 12)} | "
                       "unresolved | passes |", f"T2.md Iso A: {row_a}")
        b_mre = kc["SRRB1"][("lineage", "CaMg/MRE")]
        row_b = next(x for x in md.splitlines() if x.startswith("| Iso B |"))
        check(row_b == f"| Iso B | {b_mre} | {2 * 15 + 15} | {2 * 15 + 15} | novel (SSU) | {2 * 15} | unresolved | "
                       "passes |" and "Iso B: novel (SingleM) (15)" in md, f"T2.md Iso B: {row_b}")
        row_c = next(x for x in md.splitlines() if x.startswith("| Iso C |"))
        check(row_c == f"| Iso C | {kc['SRRC1'][('lineage', 'CaMg/MRE')]} | 2 | 2 | novel (SSU) | 2 | absent | FAILS |",
              f"T2.md Iso C: {row_c}")
        row_d = next(x for x in md.splitlines() if x.startswith("| Iso D |"))
        check(row_d == f"| Iso D | {kc['SRRD1'][('lineage', 'CaMg/MRE')]} | 0 | 0 | none | 0 | none | FAILS |" and
              kc["SRRD1"][("lineage", "CaMg/MRE")] > 0, f"T2.md Iso D: {row_d}")
        check("## K2: FAILED" in md and "- Iso C: novel (SSU), where the held-out MRE is, is called absent" in md and
              "- Iso D: no route lineage holds the held-out MRE's reads" in md, "T2.md K2 reasons")
        # pieces: S_X and two of S_N are routed in both arms and were placed by T1; the MRE's 8 only when held out
        sp = next(x for x in md.splitlines() if x.startswith("SSU pieces placed against all of SILVA"))
        check(sp == "SSU pieces placed against all of SILVA: 11 (3 of them placed before, in stage 0a or T1). For the "
                    "held-out arm, 11 pieces: 3 follow from their full-SILVA placement and 8 were searched again "
                    "against the T2 database, as their top hits were all removed. Of the derived placements, 3 were "
                    "also searched directly: all matched.", f"T2.md placement line: {sp}")
        check(first_calls.count("vsearch") == 3, f"VSEARCH runs: {first_calls.count('vsearch')} (want 3: all of "
                                                 "SILVA, the T2 database and the check)")
        searched = {d: sum(line.startswith(">") for fa in (work / "t2" / d / "batches").glob("batch_*.fa")
                           for line in open(fa)) for d in ("place_full", "place_t2", "place_check")}
        check(searched == {"place_full": 8, "place_t2": 8, "place_check": 3},
              f"pieces searched: {searched} (T1's 3 placements should be used again)")

        # a rerun aligns and places nothing again, and gives the same tables
        before = {n: (t2 / n).read_text() for n in ("levels.tsv", "calls.tsv", "T2.md")}
        calls_log.unlink()
        code, err = run_t2()
        again = calls_log.read_text().split() if calls_log.exists() else []
        check(code == 0 and again == [], f"rerun: exit {code}, calls {again}")
        check(all((t2 / n).read_text() == before[n] for n in ("levels.tsv", "calls.tsv")) and
              (t2 / "T2.md").read_text().split("\n", 3)[3] == before["T2.md"].split("\n", 3)[3],
              "a rerun gives the same tables")

        # with Iso A and B only, K2 passes
        t2b = tmp / "t2b"
        t2b.mkdir()
        shutil.copy(t2 / "silva_removed.tsv", t2b / "silva_removed.tsv")
        write_tsv(t2b / "runs.tsv", T.RUN_COLS, rows("stage0a", {"Iso A", "Iso B"}))
        code, err = run_t2(t2b)
        check(code == 0 and "## K2: passed" in (t2b / "T2.md").read_text(), f"Iso A and B only: {err[-500:]}")

        # a derivation that does not match the direct search stops the run
        real = RT.derive

        def wrong(full, removed_):
            t2_hits, todo = real(full, removed_)
            q = next(q for q, h in t2_hits.items() if h)
            t2_hits[q] = t2_hits[q][:-1] + [("DEC0.1.900", 99.0)]
            return t2_hits, todo

        RT.derive = wrong
        code, err = run_t2()
        RT.derive = real
        check(code == 1 and "differ from a direct search" in err, f"a wrong derivation: exit {code}, {err[-300:]!r}")

    # refusals
    lst = t2 / "silva_removed.tsv"
    keep = lst.read_text()
    lst.write_text(keep.replace("MRE\t", "MRE \t", 1))
    code, err = run_t2()
    check(code == 1 and "was not made from" in err, f"an edited removal list: exit {code}, {err!r}")
    lst.write_text(keep)
    fa = Path(json.loads((work / "db" / "silva_t2" / "manifest.json").read_text())["fasta"])
    seqs = fa.read_text()
    fa.write_text(seqs.replace(">MIXK.1.350", ">MIXQ.1.350"))
    code, err = run_t2()
    check(code == 1 and "FASTA" in err, f"a changed T2 FASTA: exit {code}, {err!r}")
    fa.write_text(seqs)

if fails:
    print("FAIL\n" + "\n".join(fails))
    sys.exit(1)
print(f"PASS: {len(cases)} removal-rule cases; make_t2_silva (removal list, T2 database, counts, a wrong sequence "
      "count writes nothing); run_t2 end to end (the held-out arm assigned with the T2 index and holding no "
      "CaMg/MRE lineage, the kept-in arm calling CaMg/MRE resident, MRE reads routed to novel (SSU) and novel "
      "(SingleM), a tie unplaced when kept in and placed in Bacillus when held out, placements derived from all of "
      "SILVA and checked, K2 per isolate: two pass, one absent and one with no lineage fail), a rerun aligns and "
      "places nothing, K2 passes with Iso A and B only, a wrong derivation stops the run, refusals")
