"""Check run_t3.py (PROTOCOL.md section 5, T3; section 6, K3; section 10, entry 13).

Needs minimap2 and VSEARCH (activate phase0-reads, or set MINIMAP2= and VSEARCH=). InSilicoSeq, SSUplex,
SingleM, fastq-dump and sylph are stand-ins. The spike levels are set larger than the protocol's, so the
runs can be small (the code path is the same).

A synthetic lab ("Ottawa") with three isolates, each with runs in two BioProjects of different years, and a
B- line whose run T1 already fetched and routed; and two sets from outside the fungus. Two kit-contaminant
genomes: Ralstonia, whose download holds a contig with an SSU gene that the cleaned reference lacks (its
reads are routed, found by SSUplex and placed in Ralstonia), and Pseudomonas, whose download holds a copy of
part of the CaGg genome (its reads count for CaGg). Planted:
  A1 + Pseudomonas   isolate A has Pseudomonas in its other run (another year): from the third level up, E1
                     and E2 hold, E5 holds, and Pseudomonas is resident (K3 fails)
  B1 + Pseudomonas   isolate B has CaGg in its other run, and some CaGg reads in B1: the spike's CaGg reads
                     make B1 present (the spread check passes on both read sets together), so CaGg, a known
                     clade, becomes resident (K3 fails); Pseudomonas itself is only a contaminant
  Ralstonia spikes   contaminant; their route reads form "Ralstonia (route)"
  lab kit spike      Ralstonia in all six runs of A, B and C: contaminant through E5; with C left out, E5
                     holds for A and B and both are resident (K3 fails)
  outside sets       Mycoavidus, a known clade, in both runs of the hyphosphere set: resident (K3 fails);
                     Pseudomonas in the RNA set: contaminant
Checked as well: a spiked run gives what the same run gives with the spike's reads in its files (run_t1.py
on the concatenated reads); T1's reads and route outputs are used again; a rerun runs no tool; a run list
with only the passing parts passes; refusals (a genome that is not the download, a lab spike outside the
lab projects).
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
import run_t1 as T  # noqa: E402
import run_t3 as T3  # noqa: E402

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)


MM2 = os.environ.get("MINIMAP2") or shutil.which("minimap2")
VS = os.environ.get("VSEARCH") or shutil.which("vsearch")
if not (MM2 and VS):
    print("SKIP: minimap2 or VSEARCH not found")
    sys.exit(0)

T3.LEVELS = [("L1", 0.001), ("L2", 0.01), ("L3", 0.05), ("L4", 0.2)]   # test sizes; the protocol's are 1000x less
T3.LAB_LEVEL, T3.LAB_PROJECTS = 0.05, ("P1", "P2")

rng = random.Random(41)
COMP = str.maketrans("ACGT", "TGCA")


def rand(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


S_G, S_M, S_R = rand(450), rand(450), rand(450)
G = {"host1": rand(60000), "host2": rand(40000), "mre": rand(8000),
     "gg": rand(3000) + S_G + rand(16550),                               # CaGg, 20 kb, SSU at 3000
     "my": rand(2000) + S_M + rand(17550),                               # Mycoavidus, 20 kb
     "ral": rand(15000), "ralY": rand(1000) + S_R + rand(1000),          # Ralstonia; Y only in the download
     "pse": rand(12000), "org": rand(5000), "mito": rand(6000)}
G["pseX"] = G["gg"][8000:11000]                                          # Pseudomonas download: CaGg 8-11 kb
ACC = {"host1": "GCA_900000501.1", "host2": "GCA_900000502.1", "gg": "GCA_800000501.1", "my": "GCA_800000502.1",
       "mre": "GCA_800000503.1", "ral": T3.LAB_GENOME, "pse": "GCF_700000502.1", "org": "NC_900501.1",
       "mito": "MN900501.1"}
RAW = {"ral": [("c1", G["ral"]), ("cY", G["ralY"])], "pse": [("c1", G["pse"]), ("cX", G["pseX"])]}
SSU = [[S_G, "bacteria"], [S_M, "bacteria"], [S_R, "bacteria"]]
BURK = "Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Burkholderiales;Burkholderiaceae"
SILVA = [("RALS.1.450", f"{BURK};Ralstonia;Ralstonia pickettii", S_R),
         ("GLOM.1.450", f"{BURK};Candidatus Glomeribacter;Candidatus Glomeribacter gigasporarum", S_G),
         ("MYCO.1.450", f"{BURK};Mycoavidus;Mycoavidus cysteinexigens", S_M)] + \
        [(f"DEC{i}.1.900", f"{BURK};Burkholderia;decoy {i}", rand(900)) for i in range(3)]
T_, M_ = 0.01, 0.25


def spread(src, n, length, off=0):
    step = (length - 400) // n
    return [f"{src}@{off + i * step}" for i in range(n)]


def rep(x, k):
    return [x] * k


RUNS = {"SRRA1": spread("host1", 2000, 60000),
        "SRRA2": spread("host1", 1940, 60000) + spread("pse", 60, 12000),
        "SRRB1": spread("host1", 1988, 60000) + spread("gg", 12, 20000, 100),
        "SRRB2": spread("host1", 1960, 60000) + spread("gg", 40, 20000, 50),
        "SRRC1": spread("host1", 2000, 60000), "SRRC2": spread("host1", 2000, 60000, 7),
        "SRRM1": spread("host2", 2000, 40000),
        "SRRH1": spread("host1", 1960, 60000) + spread("my", 40, 20000),
        "SRRH2": spread("host1", 1960, 60000, 11) + spread("my", 40, 20000, 77),
        "SRRG1": spread("host1", 1960, 60000) + spread("pse", 40, 12000),
        "SRRG2": spread("host1", 1960, 60000, 5) + spread("pse", 40, 12000, 33)}
#        run      unit                    isolate set              kind   bioproject year  lab     bg     lab   model
LIST = [("SRRA1", "R. testii A", "A", "ottawa", "DNA", "P1", "2019", "Ottawa", "yes", "yes", "NovaSeq"),
        ("SRRA2", "R. testii A", "A", "ottawa", "DNA", "P2", "2021", "Ottawa", "no", "yes", "HiSeq"),
        ("SRRB1", "R. testii B", "B", "ottawa", "DNA", "P1", "2019", "Ottawa", "yes", "yes", "NovaSeq"),
        ("SRRB2", "R. testii B", "B", "ottawa", "DNA", "P2", "2021", "Ottawa", "no", "yes", "HiSeq"),
        ("SRRC1", "R. testii C", "C", "ottawa", "DNA", "P1", "2019", "Ottawa", "no", "yes", "NovaSeq"),
        ("SRRC2", "R. testii C", "C", "ottawa", "RNA", "P2", "2021", "Ottawa", "no", "yes", "NovaSeq"),
        ("SRRM1", "BEG34 B-", "BEG34", "ottawa", "DNA", "P0", "", "Ottawa", "yes", "no", "NovaSeq"),
        ("SRRH1", "MUCL hyphosphere", "MUCL", "hyphosphere", "DNA", "P3", "2020", "CAU", "no", "no", ""),
        ("SRRH2", "MUCL hyphosphere", "MUCL", "hyphosphere", "DNA", "P3", "2020", "CAU", "no", "no", ""),
        ("SRRG1", "MUCL added bacterium", "MUCL", "added_bacterium", "RNA", "P4", "2023", "UCL", "no", "no", ""),
        ("SRRG2", "MUCL added bacterium", "MUCL", "added_bacterium", "RNA", "P4", "2023", "UCL", "no", "no", "")]


def t3_rows(edit=None, keep=None):
    out = []
    for run, unit, iso, set_, kind, bp, year, lab, bg, lab_s, model in LIST:
        if keep and run not in keep:
            continue
        if edit and run in edit:
            bg, lab_s, model = edit[run]
        out.append([run, unit, f"CTRL_{iso}", "yes", kind, "short", "dna" if kind == "DNA" else "unstated", "PAIRED",
                    bp, year, lab, iso, "ssu+singlem" if kind == "DNA" else "ssu", "", "test run", set_, bg, lab_s,
                    model, "test instrument"])
    return out


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


FAKE_FASTQ_DUMP = '''import gzip, json, os, sys
a = sys.argv[1:]
if a == ["--version"]:
    print("\\nfastq-dump : 3.4.1\\n"); sys.exit(0)
spec = json.load(open(os.environ["FAKE_SPEC"]))
run, x, out = a[-1], int(a[a.index("-X") + 1]), a[a.index("-O") + 1]
open(os.environ["FAKE_LOG"], "a").write(f"fastq-dump {run}\\n")
os.makedirs(out, exist_ok=True)
names = [f"{run}_1.fastq.gz", f"{run}_2.fastq.gz"]
if run in spec.get("copy", {}):
    n = 0
    for name, srcs in zip(names, spec["copy"][run]):
        with gzip.open(os.path.join(out, name), "wb") as fh:
            for s in srcs:
                data = gzip.open(s, "rb").read()
                fh.write(data)
                n += data.count(b"\\n") // 4 if name.endswith("_1.fastq.gz") else 0
    sys.stderr.write(f"Read {n} spots for {run}\\nWritten {n} spots for {run}\\n"); sys.exit(0)
comp = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(comp)[::-1]
g, spots = spec["genomes"], spec["runs"][run][:x]
fhs = [gzip.open(os.path.join(out, n), "wt") for n in names]
for i, spot in enumerate(spots, 1):
    src, _, at = spot.partition("@")
    s = int(at)
    for fh, read in zip(fhs, [g[src][s:s + 150], rc(g[src][s + 200:s + 350])]):
        fh.write(f"@{run}.{i} {i} length=150\\n{read}\\n+{run}.{i} {i} length=150\\n{'I' * 150}\\n")
for fh in fhs:
    fh.close()
sys.stderr.write(f"Read {len(spots)} spots for {run}\\nWritten {len(spots)} spots for {run}\\n")
'''

FAKE_ISS = '''import gzip, os, random, sys
a = sys.argv[1:]
if a == ["--version"]:
    print("iss version 2.0.1"); sys.exit(0)
opt = lambda k: a[a.index(k) + 1]
genome, model, n, seed, out = opt("--draft"), opt("--model"), int(opt("--n_reads")), int(opt("--seed")), opt("--output")
assert a[0] == "generate" and opt("--mode") == "kde" and opt("--cpus") == "1" and "--compress" in a
open(os.environ["FAKE_LOG"], "a").write(f"iss {os.path.basename(genome)} {n} {seed} {model}\\n")
L = {"HiSeq": 126, "NovaSeq": 151, "MiSeq": 301}[model]
frag = 2 * L + 100
seqs, name = {}, None
for line in open(genome):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:].split()[0]
        seqs[name] = []
    elif line:
        seqs[name].append(line)
seqs = {k: "".join(v) for k, v in seqs.items()}
total = sum(len(s) for s in seqs.values())
rng = random.Random(seed * 1000003 + n)
comp = str.maketrans("ACGT", "TGCA")
with gzip.open(out + "_R1.fastq.gz", "wt") as f1, gzip.open(out + "_R2.fastq.gz", "wt") as f2:
    for i in range(n // 2):
        x = rng.randrange(total)
        for k, s in seqs.items():
            if x < len(s):
                break
            x -= len(s)
        st = min(x, len(s) - frag)
        r1, r2 = s[st:st + L], s[st + frag - L:st + frag].translate(comp)[::-1]
        f1.write(f"@{k}_{i}_0/1\\n{r1}\\n+\\n{'I' * L}\\n")
        f2.write(f"@{k}_{i}_0/2\\n{r2}\\n+\\n{'I' * L}\\n")
open(out + "_abundance.txt", "w").write(os.path.basename(genome) + "\\t1.0\\n")
open(out + ".iss.tmp.0.vcf", "w").close()
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

FAKE_CONDA = '''import sys
a = sys.argv[1:]
rest = a[3:]
if rest[:1] == ["--no-capture-output"]:
    rest = rest[1:]
cmd = rest[1:]
if cmd == ["--version"]:
    print("0.21.4"); sys.exit(0)
otu = cmd[cmd.index("--otu-table") + 1]
cols = ["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "read_names", "nucleotides_aligned",
        "taxonomy_by_known?", "read_unaligned_sequences", "equal_best_hit_taxonomies"]
open(otu, "w").write("\\t".join(cols) + "\\n")
'''

FAKE_SYLPH = '''import sys
a = sys.argv[1:]
if a == ["--version"]:
    print("sylph 1.0.0"); sys.exit(0)
if a[0] == "sketch":
    open(a[a.index("-o") + 1] + ".syldb", "w").write("x\\n"); sys.exit(0)
cols = ["Sample_file", "Genome_file", "Taxonomic_abundance", "Sequence_abundance", "Adjusted_ANI", "Eff_cov",
        "Contig_name"]
open(a[a.index("-o") + 1], "w").write("\\t".join(cols) + "\\n")
'''

WRAP = '''import os, sys
real, a = {real!r}, sys.argv[1:]
if a != ["--version"]:
    open(os.environ["FAKE_CALLS"], "a").write({name!r} + "\\n")
os.execv(real, [real] + a)
'''


def build(tmp):
    refdir, work, stage = tmp / "references", tmp / "work", tmp / "stage0a"
    for sub in ("clean/hosts", "clean/bacteria", "db/silva", "db/hmm", "db/mp", "genomes"):
        (work / sub).mkdir(parents=True)
    for d in (refdir, stage):
        d.mkdir()
    for k in ("host1", "host2"):
        (work / "clean" / "hosts" / f"{ACC[k]}.fna").write_text(f">{ACC[k]}|c1\n{G[k]}\n")
    for k in ("gg", "my", "mre", "ral", "pse"):
        (work / "clean" / "bacteria" / f"{ACC[k]}.fna").write_text(f">{ACC[k]}|c1\n{G[k]}\n")
    gf = []
    for k in ("ral", "pse"):
        p = work / "genomes" / f"{ACC[k]}.fna"
        p.write_text("".join(f">{ACC[k]}|{c}\n{s}\n" for c, s in RAW[k]))
        gf.append([ACC[k], "contaminant", k, len(RAW[k]), sum(len(s) for _, s in RAW[k]), T3.AR.md5_of(p)])
    write_tsv(refdir / "genome_files.tsv", ["accession", "role", "organism", "downloaded_sequences",
                                            "downloaded_bp", "downloaded_md5"], gf)
    (work / "clean" / "organelles.fna").write_text(f">{ACC['org']}|{ACC['org']}\n{G['org']}\n"
                                                   f">{ACC['mito']}|{ACC['mito']}\n{G['mito']}\n")
    write_tsv(refdir / "hosts.tsv", ["accession", "organism", "note"],
              [[ACC["host1"], "Rhizophagus testii", "host"], [ACC["host2"], "Gigaspora testii", "host"]])
    write_tsv(refdir / "endosymbionts.tsv", ["accession", "organism", "note"],
              [[ACC["gg"], "Candidatus Glomeribacter testii", ""], [ACC["my"], "Mycoavidus testii", ""],
               [ACC["mre"], "Candidatus Moeniiplasma testii", "held_out_in_T2"]])
    write_tsv(refdir / "contaminants.tsv", ["accession", "organism", "note"],
              [[ACC["ral"], "Ralstonia pickettii", ""], [ACC["pse"], "Pseudomonas testii", ""]])
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
             "sylph": stand_in(bin_ / "sylph", FAKE_SYLPH), "iss": stand_in(bin_ / "iss", FAKE_ISS),
             "mm2": stand_in(bin_ / "mm2wrap", WRAP.format(real=MM2, name="minimap2")),
             "vsearch": stand_in(bin_ / "vswrap", WRAP.format(real=VS, name="vsearch"))}
    stand_in(bin_ / "nhmmer", 'print("# HMMER 3.4 (Aug 2023); http://hmmer.org/")\n')
    os.environ["PATH"] = f"{bin_}{os.pathsep}{os.environ['PATH']}"
    os.environ["FAKE_SPEC"], os.environ["FAKE_LOG"] = str(tmp / "spec.json"), str(tmp / "tools.log")
    os.environ["FAKE_CALLS"] = str(tmp / "calls.log")
    (tmp / "spec.json").write_text(json.dumps({"genomes": G, "runs": RUNS, "ssu": SSU}))
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
    write_tsv(tmp / "controls.tsv", ["control_id", "expected_call"],
              [[f"CTRL_{x}", "absent"] for x in ("A", "B", "C", "BEG34", "MUCL", "PHYS")])
    return refdir, work, stage, tools


def t1_run(tmp, refdir, work, stage, tools, rows, name):
    d = tmp / name
    d.mkdir()
    write_tsv(d / "runs.tsv", T.RUN_COLS, rows)
    return quiet(T.main, ["--workdir", str(work), "--refdir", str(refdir), "--t1-dir", str(d), "--stage0a-dir",
                          str(stage), "--controls", str(tmp / "controls.tsv"), "--threads", "2", "--jobs", "2",
                          "--fastq-dump", tools["fastq_dump"], "--minimap2", MM2, "--ssuplex", tools["ssuplex"],
                          "--conda", tools["conda"], "--vsearch", VS, "--sylph", tools["sylph"]]), d


def log_lines(p):
    return p.read_text().splitlines() if p.exists() else []


with tempfile.TemporaryDirectory() as td:
    tmp = Path(td)
    refdir, work, stage, tools = build(tmp)
    # T1 has the B- run already: its reads and route outputs are used again
    (code, err), _ = t1_run(tmp, refdir, work, stage, tools, [r[:15] for r in t3_rows(keep={"SRRM1"})], "t1out")
    check(code == 0, f"T1 on the B- run failed: {err[-1500:]}")

    def run_t3(t3_dir, argv_extra=()):
        return quiet(T3.main, ["--workdir", str(work), "--refdir", str(refdir), "--t3-dir", str(t3_dir),
                               "--stage0a-dir", str(stage), "--controls", str(tmp / "controls.tsv"), "--threads", "2",
                               "--jobs", "3", "--fastq-dump", tools["fastq_dump"], "--minimap2", tools["mm2"],
                               "--ssuplex", tools["ssuplex"], "--conda", tools["conda"], "--vsearch", tools["vsearch"],
                               "--sylph", tools["sylph"], "--iss", tools["iss"], *argv_extra])

    t3 = tmp / "t3"
    t3.mkdir()
    write_tsv(t3 / "runs.tsv", T3.T.RUN_COLS + T3.EXTRA_COLS, t3_rows())
    log0 = len(log_lines(tmp / "tools.log"))
    code, err = run_t3(t3)
    check(code == 0, f"run_t3 failed: {err[-3000:]}")
    if code == 0:
        new_log = log_lines(tmp / "tools.log")[log0:]
        check("fastq-dump SRRM1" not in new_log and not any(x.startswith("ssuplex SRRM1_") for x in new_log) and
              any(x.startswith("ssuplex SRRA1_") for x in new_log),
              f"T1's reads and SSUplex outputs for the B- run are used again: {new_log[:12]}")
        n_iss = sum(x.startswith("iss ") for x in new_log)
        check(n_iss == 3 * 2 * 4 * 3 + 6 * 3, f"InSilicoSeq runs: {n_iss}, want 90")
        sp = {x["spike"]: x for x in read_tsv(t3 / "spikes.tsv")}
        lab = read_tsv(t3 / "lab_spike.tsv")
        md = (t3 / "T3.md").read_text()
        check(len(sp) == 72 and len(lab) == 9, f"rows: {len(sp)} spikes, {len(lab)} lab rows")
        # the reads of a spike: level x the background's reads, from the download (Y and X contigs included)
        s4 = sp["SRRA1_Ralstonia_L4_s1"]
        check(int(s4["spike_reads"]) == 800 and 0 < int(s4["spike_reads_counted"]) < 800 and
              int(s4["route_n_after"]) > 0 and "Ralstonia (route)" in s4["route_calls_after"],
              f"A1 + Ralstonia, 20%: {s4}")
        # planted single-run results
        want = {f"SRRA1_Pseudomonas_L{lv}_s{s}" for lv in (3, 4) for s in (1, 2, 3)} | \
               {f"SRRB1_Pseudomonas_L{lv}_s{s}" for lv in (3, 4) for s in (1, 2, 3)}
        got = {k for k, x in sp.items() if x["k3"] == "FAILS"}
        check(got == want, f"single-run spikes failing K3: {sorted(got ^ want)} differ")
        a3 = sp["SRRA1_Pseudomonas_L3_s1"]
        check(a3["call_before"] == "contaminant" and a3["call_after"] == "resident" and a3["level_after"] == "present"
              and a3["level_before"] == "absent", f"A1 + Pseudomonas, 5%: {a3}")
        b3 = sp["SRRB1_Pseudomonas_L3_s1"]
        check(b3["call_after"] == "contaminant" and "CaGg: unresolved -> resident" in b3["changed_calls"],
              f"B1 + Pseudomonas, 5%: {b3['call_after']}, {b3['changed_calls']}")
        b2 = sp["SRRB1_Pseudomonas_L2_s1"]
        check(b2["k3"] == "passes" and "CaGg:" not in b2["changed_calls"], f"B1 + Pseudomonas, 1%: {b2}")
        check(all(x["call_after"] in ("contaminant", "absent", "unresolved", "") for k, x in sp.items()
                  if "Ralstonia" in k), "Ralstonia single-run spikes: never resident")
        check(all(x["k3"] == "passes" for k, x in sp.items() if k.startswith("SRRM1")), "B- spikes pass")
        # lab kit spike: contaminant through E5 in A, B and C
        check(all(x["call_after"] == "contaminant" and x["e5"] == "no" and x["k3"] == "passes" for x in lab),
              f"lab kit spike: {[(x['unit'], x['seed'], x['call_after'], x['e5']) for x in lab]}")
        check("## K3: FAILED" in md and "Mycoavidus is called resident in MUCL hyphosphere" in md and
              "Pseudomonas is called resident" in md and "CaGg becomes resident (unresolved without the spike)" in md,
              "T3.md K3 reasons")
        calls = {(c["unit"], c["lineage"]): c for c in read_tsv(t3 / "calls.tsv")}
        check(calls[("MUCL added bacterium", "Pseudomonas")]["call"] == "contaminant" and
              calls[("MUCL hyphosphere", "Mycoavidus")]["call"] == "resident" and
              calls[("R. testii B", "CaGg")]["call"] == "unresolved", "baseline calls")

        # a spiked run gives what the run gives with the spike's reads in its files (run_t1.py on both)
        sdir = work / "t3" / "spikes" / "SRRA1_Ralstonia_L4_s1"
        rdir = work / "t3" / "reads" / "SRRA1"
        spec = json.loads((tmp / "spec.json").read_text())
        spec["copy"] = {"SRRPHYS": [[str(rdir / "SRRA1_1.fastq.gz"), str(sdir / "SRRA1_Ralstonia_L4_s1_1.fastq.gz")],
                                    [str(rdir / "SRRA1_2.fastq.gz"), str(sdir / "SRRA1_Ralstonia_L4_s1_2.fastq.gz")]]}
        (tmp / "spec.json").write_text(json.dumps(spec))
        phys = [["SRRPHYS", "PHYS", "CTRL_PHYS", "yes", "DNA", "short", "dna", "PAIRED", "P9", "2019", "Lab9",
                 "PHYS", "ssu+singlem", "", "A1 with a spike in its files"]]
        (code2, err2), d = t1_run(tmp, refdir, work, stage, tools, phys, "t1phys")
        check(code2 == 0, f"T1 on the concatenated run failed: {err2[-1500:]}")
        if code2 == 0:
            lv = {x["lineage"]: x for x in read_tsv(d / "levels_entry10.tsv") if x["run"] == "SRRPHYS"}
            route_n = int(lv.get("Ralstonia (route)", {}).get("n", -1))
            check(int(lv["Ralstonia"]["n"]) == int(s4["n_after"]) and lv["Ralstonia"]["spread"] == s4["spread_after"]
                  and route_n == int(s4["route_n_after"]),
                  f"added vs in the files: n {lv['Ralstonia']['n']} / {s4['n_after']}, spread "
                  f"{lv['Ralstonia']['spread']} / {s4['spread_after']}, route "
                  f"{route_n} / {s4['route_n_after']}")

        # a rerun runs no tool and writes the same tables
        before = {n: (t3 / n).read_text() for n in ("levels.tsv", "calls.tsv", "spikes.tsv", "lab_spike.tsv")}
        (tmp / "calls.log").unlink(missing_ok=True)
        log1 = len(log_lines(tmp / "tools.log"))
        code, err = run_t3(t3)
        again = log_lines(tmp / "calls.log") + log_lines(tmp / "tools.log")[log1:]
        check(code == 0 and again == [], f"rerun: exit {code}, tool calls {again[:5]}, {err[-300:]}")
        check(all((t3 / n).read_text() == before[n] for n in before), "a rerun writes the same tables")

        # the lab kit spike without isolate C: E5 holds for A and B, and both are resident
        t3b = tmp / "t3b"
        t3b.mkdir()
        write_tsv(t3b / "runs.tsv", T3.T.RUN_COLS + T3.EXTRA_COLS,
                  t3_rows(edit={"SRRC1": ("no", "no", ""), "SRRC2": ("no", "no", "")}))
        code, err = run_t3(t3b)
        check(code == 0, f"run_t3 without C in the lab spike failed: {err[-1000:]}")
        if code == 0:
            labb = {(x["unit"], x["seed"]): x for x in read_tsv(t3b / "lab_spike.tsv")}
            check(len(labb) == 6 and all(x["call_after"] == "resident" and x["e5"] == "yes" and x["k3"] == "FAILS"
                                         for x in labb.values()), f"lab spike without C: {labb}")
            check("lab kit spike, seed 1 (R. testii A): Ralstonia is called resident" in (t3b / "T3.md").read_text(),
                  "T3.md: lab kit spike failure")

        # only the passing parts: K3 passes
        t3c = tmp / "t3c"
        t3c.mkdir()
        write_tsv(t3c / "runs.tsv", T3.T.RUN_COLS + T3.EXTRA_COLS,
                  t3_rows(edit={"SRRA1": ("no", "yes", "NovaSeq"), "SRRB1": ("no", "yes", "NovaSeq")},
                          keep={"SRRA1", "SRRA2", "SRRB1", "SRRB2", "SRRC1", "SRRC2", "SRRM1", "SRRG1", "SRRG2"}))
        code, err = run_t3(t3c)
        check(code == 0 and "## K3: passed" in (t3c / "T3.md").read_text(), f"passing run list: {err[-1000:]}")

    # refusals
    t3d = tmp / "t3d"
    t3d.mkdir()
    write_tsv(t3d / "runs.tsv", T3.T.RUN_COLS + T3.EXTRA_COLS, t3_rows(edit={"SRRH1": ("no", "yes", "NovaSeq")}))
    code, err = run_t3(t3d)
    check(code == 1 and "only paired short-read runs of the Ottawa set are spiked" in err and
          not (t3d / "T3.md").exists(), f"a spike outside the Ottawa set: exit {code}, {err[-300:]!r}")
    g = work / "genomes" / f"{ACC['pse']}.fna"
    keep = g.read_text()
    g.write_text(keep.replace("cX", "cZ"))
    code, err = run_t3(t3)
    check(code == 1 and "is not the genome clean_references.py downloaded" in err, f"a changed genome: {err[-300:]!r}")
    g.write_text(keep)

if fails:
    print("FAIL\n" + "\n".join(fails))
    sys.exit(1)
print("PASS: spikes made from the downloaded genomes and added to their background runs (equal to run_t1.py on the "
      "concatenated reads), single-run spikes (a resident spiked genus and a lineage that becomes resident both "
      "fail K3; contaminants pass), the lab kit spike (contaminant through E5; resident when E5 holds), outside sets "
      "(a known clade resident fails K3), T1's reads and route outputs used again, a rerun runs no tool, a passing "
      "run list passes, refusals")
