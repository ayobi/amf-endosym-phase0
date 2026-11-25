"""Check run_t4.py (PROTOCOL.md section 5, T4; section 6, the Logan switch; section 10, entry 14).

Needs minimap2, VSEARCH and zstd (activate phase0-reads, or set MINIMAP2=, VSEARCH= and ZSTD=). SSUplex,
SingleM, fastq-dump and sylph are stand-ins, and Logan is a local web server.

Part 1, cutting unitigs into pieces, on random and edge lengths: every piece as long as the reads (or the
whole unitig), every 31-mer in a piece and counted once, and weights that add up to ka x the unitig's 31-mers
/ the 31-mers in a read.

Part 2, end to end on synthetic genomes. Stage 0a, T1 and T2 are made with run_t1.py, make_t2_silva.py and
run_t2.py, as in check_t2.py. Logan's unitigs are made from the same reads: each read is one unitig with
ka:f:1, so each piece is a read with weight 1. T4 must then give what raw reads gave: the same levels (n, N,
fraction, level) and calls in T1 and in both arms of T2, K2 failing for the same isolate, and the switch at
0%. Then a second Logan, with changes:
  SRRP2 (a positive) is missing: listed and left out; the B+ line's calls are not compared, and the runs of
                its BioProject are marked (cross-talk)
  SRRA2         its unitigs lack the MRE's reads, so CaMg/MRE is not present there: a miss. With SRRP2 left
                out, 1 of 5 counted positives is missed (20%), so the survey runs on raw reads
  SRRP1, SRRA1  unitigs at ka:f:3 and 2: those weights, in mapped and route counts, and the same levels
  SRRN1         a long host unitig is added: pieces and weights as cut
  SRRH1         a HiFi run added to T1's list: left out
Checked as well: a rerun downloads, unpacks, aligns and places nothing; a cut-off download is tried three
times, reported, and fetched on the next run; refusals (no zstd, T2 not done, a positive not in the run lists).
"""
import contextlib
import csv
import functools
import gzip
import hashlib
import http.server
import io
import json
import os
import random
import shutil
import stat
import subprocess
import sys
import tempfile
import threading
from pathlib import Path

here = Path(__file__).parent
sys.path.insert(0, str(here.parent / "scripts"))
import make_index as MI  # noqa: E402
import make_t2_silva as MS  # noqa: E402
import run_t1 as T  # noqa: E402
import run_t2 as RT  # noqa: E402
import run_t4 as T4  # noqa: E402

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)


# -------------------------------------------------------------- part 1: pieces ---

prng = random.Random(7)
n_cases = 0
for L, r in [(30, 150), (31, 150), (149, 150), (150, 150), (151, 150), (270, 150), (271, 150), (1000, 150)] + \
        [(prng.randint(31, 1500), prng.choice([75, 100, 150, 151, 250])) for _ in range(400)]:
    seq = bytes(prng.choice(b"ACGT") for _ in range(L))
    ka, step, n_k = prng.choice([1.0, 2.5, 17.0]), r - 30, L - 30
    ps = T4.unitig_pieces(seq, ka, r)
    n_cases += 1
    if n_k <= 0:
        check(ps == [], f"L={L}: pieces of a unitig shorter than 31 bases")
        continue
    full = [j * step for j in range(len(ps))]
    if L > r and len(ps) > 1 and full[-1] + r > L:
        full[-1] = L - r
    counted = [w * step / ka for _, w in ps]
    covered = set()
    for (p, _), s in zip(ps, full, strict=True):
        covered.update(range(s, s + len(p) - 30))
    check(all(len(p) == min(r, L) and p == seq[s:s + len(p)] for (p, _), s in zip(ps, full, strict=True)),
          f"L={L}, r={r}: piece lengths or places")
    check(covered == set(range(n_k)), f"L={L}, r={r}: not every 31-mer is in a piece")
    check(all(abs(c - round(c)) < 1e-9 and 1 <= round(c) <= step for c in counted) and
          abs(sum(counted) - n_k) < 1e-6, f"L={L}, r={r}: k-mers counted {counted}, want {n_k} in all")
    check(abs(sum(w for _, w in ps) - ka * n_k / step) < 1e-9, f"L={L}, r={r}: weights")

MM2 = os.environ.get("MINIMAP2") or shutil.which("minimap2")
VS = os.environ.get("VSEARCH") or shutil.which("vsearch")
ZSTD = os.environ.get("ZSTD") or shutil.which("zstd")
if not (MM2 and VS and ZSTD):
    print(("FAIL\n" + "\n".join(fails)) if fails else
          f"PASS (part 1, {n_cases} unitigs); SKIP part 2: minimap2, VSEARCH or zstd not found")
    sys.exit(1 if fails else 0)

# --------------------------------------------------------------- part 2: world ---

rng = random.Random(31)
COMP = str.maketrans("ACGT", "TGCA")
T4.RETRY_WAIT = 0


def rand(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def rc(s):
    return s.translate(COMP)[::-1]


def mutate(s, rate):
    return "".join(rng.choice("ACGT".replace(b, "")) if rng.random() < rate else b for b in s)


S_A, S_B, S_N, S_X, K_B, S_G = rand(450), rand(450), rand(450), rand(150), rand(30), rand(450)
mre_b = rand(3000) + S_B + rand(16550)
mre_b = mre_b[:9000] + K_B + mre_b[9030:]                                 # a SingleM marker at 9000
G = {"host": rand(60000), "mreA": rand(3000) + S_A + rand(16550), "mreB": mre_b,
     "gg": rand(3000) + S_G + rand(16550),                                # CaGg, SSU at 3000
     "cont": rand(10000), "org": rand(5000), "mito": rand(6000),
     "unk": rand(1000) + S_X + rand(1000), "nov": rand(2000) + S_N + rand(1000)}
ACC = {"host": "GCA_900000401.1", "mreA": "GCA_800000401.1", "mreB": "GCA_800000402.1", "gg": "GCA_800000403.1",
       "cont": "GCF_700000401.1", "org": "NC_900401.1", "mito": "MN900401.1"}
SSU = [[S_A, "bacteria"], [S_B, "bacteria"], [S_N, "bacteria"], [S_X, "bacteria"], [S_G, "bacteria"]]
SM_TAX = {"KB": "Root; d__Bacteria; p__Bacillota; c__Bacilli"}
MOL = "Bacteria;Bacillati;Mycoplasmatota;Mollicutes;Mycoplasmatales;Mycoplasmataceae"
BAC = "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae"
BURK = "Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Burkholderiales;Burkholderiaceae"
SILVA = [("MREA.1.450", f"{MOL};Candidatus Moeniiplasma;Candidatus Moeniiplasma glomeromycotorum", S_A),
         ("MREC.1.450", "Bacteria;Bacillati;Bacillota;Bacilli;Lactobacillales;uncultured;uncultured;"
                        "endosymbiont GvMRE of Glomus versiforme", mutate(S_A, 0.01)),
         ("MREB.1.450", "Bacteria;Bacillati;Bacillota;Bacilli;Lactobacillales;Streptococcaceae;uncultured;"
                        "Mycoplasma-related endobacterium of Gigaspora", S_B),
         ("RELA.1.450", f"{BAC};Bacillus;Bacillus sp. AMRE1", mutate(S_A, 0.15)),
         ("RELB.1.450", "Bacteria;Bacillati;Bacillota;Erysipelotrichia;Erysipelotrichales;Erysipelotrichaceae;"
                        "Erysipelothrix;Erysipelothrix sp. MREB1", mutate(S_B, 0.15)),
         ("MIXR.1.350", f"{MOL};Mycoplasma;Mycoplasma sp.", rand(100) + S_X + rand(100)),
         ("MIXK.1.350", f"{BAC};Bacillus;Bacillus sp.", rand(100) + S_X + rand(100)),
         ("GLOM.1.450", f"{BURK};Candidatus Glomeribacter;Candidatus Glomeribacter gigasporarum", S_G)] + \
        [(f"DEC{i}.1.900", f"{BAC};Bacillus;decoy {i}", rand(900)) for i in range(4)]


def at(src, lo, hi, step):
    return [f"{src}@{p}" for p in range(lo, hi, step)]


def spots(*groups):
    out = [x for g in groups for x in g]
    random.Random(len(out)).shuffle(out)
    return out


def rep(x, k):
    return [x] * k


RUNS = {"SRRA1": ("PAIRED", spots(at("mreA", 4000, 19500, 50), rep("mreA@3000", 15), rep("unk@1000", 20),
                                  rep("nov@2000", 8), rep("host", 1647))),
        "SRRA2": ("PAIRED", spots(at("mreA", 4025, 19500, 60), rep("mreA@3010", 12), rep("host", 1730))),
        "SRRB1": ("PAIRED", spots(at("mreB", 9500, 19500, 40), at("mreB", 3600, 8500, 40), rep("mreB@3000", 15),
                                  rep("mreB@8990", 15), rep("host", 1597))),
        "SRRC1": ("PAIRED", spots(rep("mreA@10000", 30), rep("mreA@3100", 1), rep("host", 7969))),
        "SRRP1": ("PAIRED", spots(at("gg", 200, 19500, 100), rep("nov@2000", 10), rep("host", 1797))),
        "SRRP2": ("PAIRED", spots(at("gg", 230, 19500, 110), rep("host", 1825))),
        "SRRP3": ("PAIRED", spots(at("gg", 300, 19000, 300), rep("host", 1937))),
        "SRRN1": ("PAIRED", spots(rep("host", 2000))),
        "SRRN2": ("SINGLE", spots(rep("host", 1500))),
        "SRRL1": ("PAIRED", spots(rep("host", 1000))),
        "SRRJ1": ("PAIRED", spots(at("gg", 500, 19000, 500), rep("host", 963)))}
T2_ROWS = [["SRRA1", "Iso A", "CTRL_A", "PRJA", "L1", "A"], ["SRRA2", "Iso A", "CTRL_A", "PRJA", "L1", "A"],
           ["SRRB1", "Iso B", "CTRL_B", "PRJB", "L2", "B"], ["SRRC1", "Iso C", "CTRL_C", "PRJB", "L2", "C"]]
#         run      unit                      control_id               scored kind   prep      layout    project year
T1_ROWS = [["SRRP1", "BEG34 B+", "GMAR_BEG34_HIC_Bplus", "yes", "DNA", "dna", "PAIRED", "PRJ1", "2019"],
           ["SRRP2", "BEG34 B+", "GMAR_BEG34_GENOME", "yes", "DNA", "dna", "PAIRED", "PRJ1", "2020"],
           ["SRRP3", "BEG34 B+", "GMAR_BEG34_RNA_Bplus", "yes", "RNA", "poly-A", "PAIRED", "PRJ3", "2018"],
           ["SRRN1", "BEG34 B-", "GMAR_BEG34_HIC_Bminus", "yes", "DNA", "dna", "PAIRED", "PRJ1", "2019"],
           ["SRRN2", "BEG34 B-", "GMAR_BEG34_RNA_Bminus", "yes", "RNA", "poly-A", "SINGLE", "PRJ3", "2018"],
           ["SRRL1", "Lotus plant-only roots", "LOTUS_NM_ROOTS", "yes", "RNA", "poly-A", "PAIRED", "PRJ4", "2021"],
           ["SRRJ1", "BEG34 JGI MT", "GMAR_BEG34_JGI_MT", "no", "RNA", "unstated", "PAIRED", "PRJ5", "2017"]]
HIFI = ["SRRH1", "BEG34 HiFi", "GMAR_BEG34_HIFI", "no", "DNA", "hifi", "dna", "SINGLE", "PRJ1", "2022", "Torino",
        "BEG34", "none", "", "HiFi"]
POSITIVES = [("SRRP1", "T1", "CaGg"), ("SRRP2", "T1", "CaGg"), ("SRRP3", "T1", "CaGg"), ("SRRA1", "T2", "CaMg/MRE"),
             ("SRRA2", "T2", "CaMg/MRE"), ("SRRB1", "T2", "CaMg/MRE"), ("SRRC1", "T2", "CaMg/MRE")]
T_, M_ = 0.01, 0.25


def t1_rows():
    out = []
    for run, unit, ctrl, scored, kind, prep, layout, bp, year in T1_ROWS:
        lab, iso = {"Lotus plant-only roots": ("Kyoto", "none"), "BEG34 JGI MT": ("JGI", "BEG34")}.get(
            unit, ("Torino", "BEG34"))
        route = "none" if unit == "BEG34 JGI MT" else ("ssu+singlem" if kind == "DNA" else "ssu")
        out.append([run, unit, ctrl, scored, kind, "short", prep, layout, bp, year, lab, iso, route, "", "test"])
    return out


def t2_rows(reuse):
    return [[run, u, c, "yes", "DNA", "short", "dna", "PAIRED", p, "2017", lab, iso, "ssu+singlem", reuse, "plate"]
            for run, u, c, p, lab, iso in T2_ROWS]


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


class Logan(http.server.SimpleHTTPRequestHandler):
    """A local Logan: HEAD and GET as S3 answers them; paths in `cut` send half of the file."""
    requests, cut = [], set()

    def log_message(self, *a):
        pass

    def do_HEAD(self):
        Logan.requests.append(("HEAD", self.path))
        super().do_HEAD()

    def do_GET(self):
        Logan.requests.append(("GET", self.path))
        if self.path in Logan.cut:
            data = Path(self.translate_path(self.path)).read_bytes()
            self.send_response(200)
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data[:len(data) // 2])
            self.close_connection = True
            return
        super().do_GET()


def build(tmp):
    refdir, work, stage, pre = tmp / "references", tmp / "work", tmp / "stage0a", tmp / "pre"
    for sub in ("clean/hosts", "clean/bacteria", "db/silva", "db/hmm", "db/mp"):
        (work / sub).mkdir(parents=True)
    for d in (refdir, stage, pre):
        d.mkdir()
    (work / "clean" / "hosts" / f"{ACC['host']}.fna").write_text(f">{ACC['host']}|c1\n{G['host']}\n")
    for k in ("mreA", "mreB", "gg", "cont"):
        (work / "clean" / "bacteria" / f"{ACC[k]}.fna").write_text(f">{ACC[k]}|c1\n{G[k]}\n")
    (work / "clean" / "organelles.fna").write_text(f">{ACC['org']}|{ACC['org']}\n{G['org']}\n"
                                                   f">{ACC['mito']}|{ACC['mito']}\n{G['mito']}\n")
    write_tsv(refdir / "hosts.tsv", ["accession", "organism", "note"], [[ACC["host"], "Gigaspora testii", "host"]])
    write_tsv(refdir / "endosymbionts.tsv", ["accession", "organism", "note"],
              [[ACC["mreA"], "Candidatus Moeniiplasma testii A", "held_out_in_T2"],
               [ACC["mreB"], "Candidatus Moeniiplasma testii B", "held_out_in_T2"],
               [ACC["gg"], "Candidatus Glomeribacter testii", ""]])
    write_tsv(refdir / "contaminants.tsv", ["accession", "organism", "note"], [[ACC["cont"], "Ralstonia testii", ""]])
    write_tsv(refdir / "organelles.tsv", ["plant", "organelle", "accession", "title", "length_bp"],
              [["Lotus testii", "chloroplast", ACC["org"], "x", 5000]])
    write_tsv(refdir / "fungal_mitochondria.tsv", ["fungus", "organelle", "accession", "title", "length_bp",
                                                   "source", "query"],
              [["Gigaspora testii", "mitochondrion", ACC["mito"], "x", 6000, "GenBank", "q"]])
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
             "vsearch": stand_in(bin_ / "vswrap", WRAP.format(real=VS, name="vsearch")),
             "zstd": stand_in(bin_ / "zstdwrap", WRAP.format(real=ZSTD, name="zstd"))}
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
    write_tsv(stage / "fractions_entry10.tsv", ["run", "sets_t", "fraction"], [["X1", "yes", T_ * 10]])
    write_tsv(stage / "marker_share_entry10.tsv", ["run", "sets_m", "share"], [["X1", "yes", M_]])
    write_tsv(stage / "runs.tsv", ["run", "test", "kind", "line", "symbiont_lineage", "sets_t", "bioproject", "layout",
                                   "note"],
              [[run, test, "DNA", "x", lin, "yes", "P", "PAIRED", ""] for run, test, lin in POSITIVES])
    write_tsv(tmp / "controls.tsv", ["control_id", "expected_call"],
              [[r[2], "resident"] for r in T1_ROWS + [HIFI]] +
              [[c, "resident"] for c in ("CTRL_A", "CTRL_B", "CTRL_C")])

    def t1(d, rows_):
        write_tsv(d / "runs.tsv", T.RUN_COLS, rows_)
        return quiet(T.main, ["--workdir", str(work), "--refdir", str(refdir), "--t1-dir", str(d), "--stage0a-dir",
                              str(stage), "--controls", str(tmp / "controls.tsv"), "--threads", "2", "--jobs", "2",
                              "--fastq-dump", tools["fastq_dump"], "--minimap2", MM2, "--ssuplex", tools["ssuplex"],
                              "--conda", tools["conda"], "--vsearch", VS, "--sylph", tools["sylph"]])

    # Stage 0a for the T2 runs, made by run_t1.py and moved where stage 0a keeps its files (as in check_t2.py)
    code, err = t1(pre, t2_rows(""))
    if code != 0:
        print("FAIL\nmaking stage 0a with run_t1.py: " + err)
        sys.exit(1)
    for src, dest in (("t1/reads", "stage0a/reads"), ("t1_entry10/assign", "stage0a_entry10/assign"),
                      ("t1_entry10/route", "stage0a/route")):
        (work / dest).mkdir(parents=True, exist_ok=True)
        for run, *_ in T2_ROWS:
            shutil.move(str(work / src / run), str(work / dest / run))
    w = work.resolve()
    for aj in (work / "stage0a_entry10" / "assign").glob("*/assign.json"):
        aj.write_text(aj.read_text().replace(f"{w}/t1/reads/", f"{w}/stage0a/reads/"))
    shutil.rmtree(work / "t1_entry10")
    # T1 and T2 on raw reads
    t1_dir, t2_dir = tmp / "t1", tmp / "t2"
    t1_dir.mkdir()
    t2_dir.mkdir()
    code, err = t1(t1_dir, t1_rows())
    if code != 0:
        print("FAIL\nrun_t1.py: " + err)
        sys.exit(1)
    code, err = quiet(MS.main, ["--workdir", str(work), "--t2-dir", str(t2_dir), "--vsearch", tools["vsearch"]])
    if code != 0:
        print("FAIL\nmake_t2_silva.py: " + err)
        sys.exit(1)
    write_tsv(t2_dir / "runs.tsv", T.RUN_COLS, t2_rows("stage0a"))
    code, err = quiet(RT.main, ["--workdir", str(work), "--refdir", str(refdir), "--t2-dir", str(t2_dir),
                                "--stage0a-dir", str(stage), "--controls", str(tmp / "controls.tsv"), "--threads", "2",
                                "--minimap2", MM2, "--vsearch", VS, "--sylph", tools["sylph"]])
    if code != 0:
        print("FAIL\nrun_t2.py: " + err)
        sys.exit(1)
    return refdir, work, stage, tools, t1_dir, t2_dir


def raw_reads(work, run):
    """Every read of a run's raw files, mate 1 then mate 2."""
    for d in (work / "t1" / "reads" / run, work / "stage0a" / "reads" / run):
        if (d / "fetch.json").exists():
            out = []
            for f in json.loads((d / "fetch.json").read_text())["files"]:
                with gzip.open(d / f, "rt") as fh:
                    out += fh.read().split("\n")[1::4]
            return out
    raise RuntimeError(run)


def put(root, run, items):
    """Logan's file for a run: [(sequence, ka)], one unitig each, zstd-compressed."""
    d = root / run
    d.mkdir(parents=True, exist_ok=True)
    fa = d / f"{run}.unitigs.fa"
    fa.write_text("".join(f">{run}_{i} ka:f:{ka:g}   L:+:{i + 1}:-\n{s}\n" for i, (s, ka) in enumerate(items)))
    subprocess.run([ZSTD, "-q", "-f", str(fa), "-o", str(d / f"{run}.unitigs.fa.zst")], check=True)
    fa.unlink()


def fnum(v):
    return float(v) if v != "" else None


LEVEL_SAME = ("level", "N", "fraction", "informative", "cross_talk", "spread_cap")
CALL_SAME = ("category", "call", "basis", "signal", "informative_runs", "present_runs", "low_runs", "e1", "e2", "e5",
             "members", "placements")


def same_as_raw(raw_lv, raw_cl, lv, cl, where):
    """Every raw level and call, the same on Logan."""
    bad = []
    for x in raw_lv:
        y = lv.get((x["run"], x["lineage"]))
        if y is None or fnum(x["n"]) != fnum(y["n"]) or fnum(x["reads"]) != fnum(y["reads"]) or \
                any(x[c] != y[c] for c in LEVEL_SAME):
            bad.append(f"{x['run']} {x['lineage']}: raw {[x[c] for c in ('n',) + LEVEL_SAME]}, Logan "
                       f"{[y[c] for c in ('n',) + LEVEL_SAME] if y else None}")
    for x in raw_cl:
        y = cl.get((x["unit"], x["lineage"]))
        if y is None or any(x[c] != y[c] for c in CALL_SAME):
            bad.append(f"{x['unit']} {x['lineage']}: raw {[x[c] for c in CALL_SAME]}, Logan "
                       f"{[y[c] for c in CALL_SAME] if y else None}")
    check(not bad and len(lv) == len(raw_lv) and len(cl) == len(raw_cl),
          f"{where}: Logan differs from raw reads ({len(lv)} levels and {len(cl)} calls; raw {len(raw_lv)} and "
          f"{len(raw_cl)}): " + "; ".join(bad[:6]))


with tempfile.TemporaryDirectory() as td:
    tmp = Path(td)
    refdir, work, stage, tools, t1_dir, t2_dir = build(tmp)
    os.environ["no_proxy"] = os.environ["NO_PROXY"] = "127.0.0.1,localhost"
    site = tmp / "logan"
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(Logan, directory=str(site)))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{server.server_address[1]}"
    calls_log = tmp / "calls.log"
    short = [r[0] for r in T1_ROWS] + [r[0] for r in T2_ROWS]
    reads = {run: raw_reads(work, run) for run in short}

    def run_t4(name, t1=t1_dir, t2=t2_dir, stage_dir=stage, extra=()):
        return quiet(T4.main, ["--workdir", str(work), "--refdir", str(refdir), "--t1-dir", str(t1), "--t2-dir",
                               str(t2), "--t4-dir", str(tmp / name), "--stage0a-dir", str(stage_dir), "--work-name",
                               name, "--threads", "2", "--jobs", "3", "--minimap2", tools["mm2"],
                               "--ssuplex", tools["ssuplex"], "--conda", tools["conda"], "--vsearch", tools["vsearch"],
                               "--zstd", tools["zstd"], "--logan-url", f"{base}/{name}", *extra])

    raw1 = {"lv": read_tsv(t1_dir / "levels_entry10.tsv"), "cl": read_tsv(t1_dir / "calls_entry10.tsv")}
    raw2 = {"lv": read_tsv(t2_dir / "levels.tsv"), "cl": read_tsv(t2_dir / "calls.tsv")}
    raw_present = {(x["run"], x["lineage"]) for x in raw1["lv"] + [y for y in raw2["lv"] if y["arm"] == "kept_in"]
                   if x["level"] == "present"}
    check({p for p, _, lin in POSITIVES if (p, lin) in raw_present} ==
          {"SRRP1", "SRRP2", "SRRP3", "SRRA1", "SRRA2", "SRRB1"} and "## K1: passed" in (t1_dir / "T1_ENTRY10.md")
          .read_text() and "## K2: FAILED" in (t2_dir / "T2.md").read_text(),
          "the world: raw positives present, K1 passes, K2 fails for Iso C")

    # ---- Logan 1: each read one unitig with ka 1: the same results as raw reads
    for run in short:
        put(site / "t4a", run, [(s, 1.0) for s in reads[run]])
    calls_log.unlink(missing_ok=True)
    code, err = run_t4("t4a")
    check(code == 0, f"run_t4 (Logan = reads) failed: {err[-2500:]}")
    if code == 0:
        out = tmp / "t4a"
        lg = {x["run"]: x for x in read_tsv(out / "logan.tsv")}
        check(all(lg[r]["in_logan"] == "yes" and lg[r]["read_length"] == "150" and
                  int(lg[r]["pieces"]) == len(reads[r]) and float(lg[r]["weight"]) == len(reads[r]) for r in short),
              f"logan.tsv: {[(r, lg[r]['pieces'], lg[r]['weight']) for r in short]}")
        md5 = hashlib.md5((site / "t4a" / "SRRP1" / "SRRP1.unitigs.fa.zst").read_bytes()).hexdigest()
        check(lg["SRRP1"]["unitigs_md5"] == md5, "logan.tsv: the md5 of the download")
        check(not any((work / "t4a" / "unitigs").iterdir()), "the downloads are removed once cut")
        lv4, cl4 = read_tsv(out / "levels.tsv"), read_tsv(out / "calls.tsv")
        for test, arm, rl, rc_ in (("T1", "full", raw1["lv"], raw1["cl"]),
                                   ("T2", "kept_in", [x for x in raw2["lv"] if x["arm"] == "kept_in"],
                                    [x for x in raw2["cl"] if x["arm"] == "kept_in"]),
                                   ("T2", "held_out", [x for x in raw2["lv"] if x["arm"] == "held_out"],
                                    [x for x in raw2["cl"] if x["arm"] == "held_out"])):
            same_as_raw(rl, rc_, {(x["run"], x["lineage"]): x for x in lv4 if (x["test"], x["arm"]) == (test, arm)},
                        {(x["unit"], x["lineage"]): x for x in cl4 if (x["test"], x["arm"]) == (test, arm)},
                        f"{test} {arm}")
        cmp_ = read_tsv(out / "compare.tsv")
        check(cmp_ and all(x["same"] == "yes" for x in cmp_),
              f"compare.tsv: {[x for x in cmp_ if x['same'] != 'yes']}")
        md = (out / "T4.md").read_text()
        check("Logan misses 0 of the 6 positive runs that are present in raw reads and that Logan holds (0%): 10% or "
              "less, so **the survey can run on Logan unitigs**." in md, "T4.md: the switch at 0%")
        check("K1: would pass." in md and "K2: would fail: Iso C: novel (SSU), where the held-out MRE is, is called "
              "absent." in md, "T4.md: K1 and K2 on Logan")
        raw_k2 = [x for x in (t2_dir / "T2.md").read_text().splitlines() if x.startswith("| Iso ") and
                  x.endswith(("passes |", "FAILS |"))]
        k2_rows = [x for x in md.splitlines() if x.startswith("| Iso ") and x.endswith(("passes |", "FAILS |"))]
        conv = ["| " + " | ".join(f"{float(c):.1f}" if c.strip().isdigit() else c.strip() for c in x.strip("|")
                                  .split("|")) + " |" for x in raw_k2]
        check(k2_rows == conv, f"T4.md K2 rows {k2_rows}, want {conv}")

        # a rerun downloads, unpacks, aligns and places nothing
        before = {n: (out / n).read_text() for n in ("logan.tsv", "levels.tsv", "calls.tsv", "compare.tsv")}
        n_req = len(Logan.requests)
        calls_log.unlink(missing_ok=True)
        code, err = run_t4("t4a")
        again = calls_log.read_text().split() if calls_log.exists() else []
        check(code == 0 and again == [] and len(Logan.requests) == n_req,
              f"rerun: exit {code}, tools {again}, requests {Logan.requests[n_req:]}")
        check(all((out / n).read_text() == before[n] for n in before), "a rerun gives the same tables")

    # ---- Logan 2: one run missing, one without the MRE, one at ka 3, a long unitig, a HiFi run, a cut download
    for run in short:
        if run == "SRRP2":
            continue
        items = [(s, 1.0) for s in reads[run]]
        if run == "SRRA2":
            items = [(s, 1.0) for s in reads[run] if s not in G["mreA"] and s not in rc(G["mreA"])]
        elif run in ("SRRP1", "SRRA1"):
            items = [(s, 3.0 if run == "SRRP1" else 2.0) for s in reads[run]]
        elif run == "SRRN1":
            items.append((G["host"][5000:6000], 2.0))
        put(site / "t4b", run, items)
    t1b = tmp / "t1b"
    shutil.copytree(t1_dir, t1b)
    write_tsv(t1b / "runs.tsv", T.RUN_COLS, t1_rows() + [HIFI])
    Logan.cut.add("/t4b/SRRN2/SRRN2.unitigs.fa.zst")
    code, err = run_t4("t4b", t1=t1b)
    gets = [p for m, p in Logan.requests if m == "GET" and p.startswith("/t4b/SRRN2/")]
    check(code == 1 and "SRRN2: " in err and "the download failed 3 times" in err and len(gets) == 3 and
          not (tmp / "t4b" / "T4.md").exists() and not (work / "t4b" / "unitigs" / "SRRN2.unitigs.fa.zst").exists(),
          f"a cut-off download: exit {code}, {len(gets)} tries, {err[-400:]!r}")
    Logan.cut.clear()
    code, err = run_t4("t4b", t1=t1b)
    check(code == 0, f"run_t4 (Logan with changes) failed: {err[-2500:]}")
    heads = [p for m, p in Logan.requests if m == "HEAD" and p.startswith("/t4b/")]
    check(len(heads) == len(set(heads)) == len(short), f"one HEAD request per run: {len(heads)}")
    if code == 0:
        out = tmp / "t4b"
        lg = {x["run"]: x for x in read_tsv(out / "logan.tsv")}
        check(lg["SRRP2"]["in_logan"] == "no" and lg["SRRP2"]["note"] == "not in Logan v1.2 (HTTP 404)" and
              lg["SRRH1"]["note"].startswith("HiFi") and lg["SRRH1"]["in_logan"] == "", "logan.tsv: runs left out")
        pm = json.loads((work / "t4b" / "pieces" / "SRRN1" / "pieces.json").read_text())
        check(pm["unitigs"] == len(reads["SRRN1"]) + 1 and pm["pieces"] == len(reads["SRRN1"]) + 9 and
              abs(pm["weight"] - (len(reads["SRRN1"]) + 2 * 970 / 120)) < 1e-3, f"SRRN1 pieces: {pm}")
        check(float(lg["SRRP1"]["weight"]) == 3 * len(reads["SRRP1"]), "SRRP1: weights of 3")
        lv4 = {(x["test"], x["arm"], x["run"], x["lineage"]): x for x in read_tsv(out / "levels.tsv")}
        r1 = {(x["run"], x["lineage"]): x for x in raw1["lv"]}
        a, b = lv4[("T1", "full", "SRRP1", "CaGg")], r1[("SRRP1", "CaGg")]
        check(a["level"] == b["level"] == "present" and float(a["n"]) == 3 * int(b["n"]) and
              a["fraction"] == b["fraction"], f"SRRP1 at ka 3: {a['n']} {a['fraction']}, raw {b['n']} {b['fraction']}")
        r2 = {(x["arm"], x["run"], x["lineage"]): x for x in raw2["lv"]}
        doubled = [(arm, lin) for arm, lin in (("kept_in", "CaMg/MRE"), ("kept_in", "novel (SSU)"),
                                               ("held_out", "novel (SSU)"), ("held_out", "Bacillus (route)"))
                   if float(lv4.get(("T2", arm, "SRRA1", lin), {}).get("n", -1)) ==
                   2 * int(r2[(arm, "SRRA1", lin)]["n"]) > 0 and
                   lv4[("T2", arm, "SRRA1", lin)]["level"] == r2[(arm, "SRRA1", lin)]["level"]]
        check(len(doubled) == 4, f"SRRA1 at ka 2, mapped and route lineages doubled: {doubled}")
        check(not any(k[2] in ("SRRP2", "SRRH1") for k in lv4), "levels.tsv: no run left out")
        check(lv4[("T2", "kept_in", "SRRA2", "CaMg/MRE")]["level"] == "absent", "SRRA2 without the MRE: absent")
        cmp_ = read_tsv(out / "compare.tsv")
        bplus = [x for x in cmp_ if x["what"] == "call" and x["run_or_unit"] == "BEG34 B+"]
        check(bplus and all(x["same"] == "n/a" and "SRRP2 (not in Logan)" in x["note"] for x in bplus),
              f"B+ calls not compared: {bplus}")
        prj1 = [x for x in cmp_ if x["what"] == "level" and x["run_or_unit"] in ("SRRP1",)]
        check(prj1 and all("cross-talk" in x["note"] for x in prj1), f"PRJ1 levels marked: {prj1}")
        diff = {(x["what"], x["arm"], x["run_or_unit"], x["lineage"]): (x["raw"], x["logan"]) for x in cmp_
                if x["same"] == "no"}
        check(diff.get(("level", "kept_in", "SRRA2", "CaMg/MRE")) == ("present", "absent") and
              diff.get(("call", "kept_in", "Iso A", "CaMg/MRE"), ("", ""))[0] == "resident" and
              all(k[2] in ("SRRA2", "Iso A") for k in diff), f"differences: {diff}")
        md = (out / "T4.md").read_text()
        check("Logan misses 1 of the 5 positive runs that are present in raw reads and that Logan holds (20%): more "
              "than 10%, so **the survey runs on raw reads**." in md and
              "| SRRA2 | T2 | CaMg/MRE | present | absent | yes | yes |" in md and
              "| SRRC1 | T2 | CaMg/MRE | low | low | no | no |" in md and
              "Positive runs left out: SRRP2 (not in Logan; present in raw reads)." in md, "T4.md: the switch at 20%")
        check("Logan holds 10 of the 12 runs" in md and "| SRRP2 | T1 | BEG34 B+ | not in Logan v1.2 (HTTP 404) |" in md
              and "| SRRH1 | T1 | BEG34 HiFi | HiFi: not repeated on Logan (entry 14) |" in md, "T4.md: runs")

    # refusals
    code, err = run_t4("t4c", extra=("--zstd", str(tmp / "no_zstd")))
    check(code == 1 and "Cannot start" in err and "zstd" in err, f"no zstd: exit {code}, {err!r}")
    t2x = tmp / "t2x"
    shutil.copytree(t2_dir, t2x)
    (t2x / "levels.tsv").unlink()
    code, err = run_t4("t4c", t2=t2x)
    check(code == 1 and "levels.tsv is missing: run run_t2.py first" in err, f"T2 not done: exit {code}, {err!r}")
    sx = tmp / "stage_x"
    shutil.copytree(stage, sx)
    with open(sx / "runs.tsv", "a") as fh:
        fh.write("SRRZ9\tT1\tDNA\tx\tCaGg\tyes\tP\tPAIRED\t\n")
    code, err = run_t4("t4c", stage_dir=sx)
    check(code == 1 and "SRRZ9" in err, f"a positive not in the run lists: exit {code}, {err!r}")
    check(not (work / "t4c").exists() or not any((work / "t4c" / "logan").iterdir()), "a refusal fetches nothing")
    server.shutdown()

if fails:
    print("FAIL\n" + "\n".join(fails))
    sys.exit(1)
print(f"PASS: {n_cases} unitigs cut into pieces; run_t4 end to end: with each read as a unitig, the same levels and "
      "calls as raw reads in T1 and both arms of T2, K2 failing for the same isolate and the switch at 0%; a rerun "
      "fetches and runs nothing; a cut-off download tried three times and fetched on the next run; with a positive "
      "missing, one without its symbiont, weights of 3, a long unitig and a HiFi run: runs left out and listed, calls "
      "of the line not compared, the miss found and the switch at 20%; refusals")
