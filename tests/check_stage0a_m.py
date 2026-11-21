"""Check get_marker_databases.py and run_stage0a_m.py end to end on synthetic genomes. Stand-ins play
SSUplex, nhmmer, SingleM (through a stand-in conda), Zenodo and the Metaxa2 download; minimap2 and
assign_reads.py are real.

Needs minimap2 (activate phase0-reads, or set MINIMAP2=/path/to/minimap2).
The genomes carry planted 30 bp motifs. The SSUplex stand-in sorts a read by the SSU motif it holds
(either strand); the SingleM stand-in reports reads holding a marker motif. The stand-in fastq-dump
records where each read came from, so every expected number is counted here, apart from the scripts.

Databases: versions are checked against env/*.lock.yml; the metapackage comes from its Zenodo record,
its md5 is checked against the value in references/DATABASES.md, and singlem data --verify-only runs;
the five profile files are made as dev/assemble_metaxa2_hmms.sh makes them (mitochondria = M + N, LSU
profiles ignored); the quick start must give the README's counts; a wrong md5, a missing profile file,
a wrong quick start, a wrong HMMER version or a damaged metapackage writes no manifest; a rerun
downloads nothing.
Step 2: per run, the symbiont reads found by each tool, by both and by either, and the share; found
reads that are routed (no reference) or counted for another lineage do not count; m is the lowest
DNA-positive share, and an RNA positive with a lower share does not set it; found reads reach minimap2
under their step 1 names (positions); big SSUplex outputs and the FASTA are removed; a rerun calls no
tool; SingleM not seeing MRE gives "MRE rest on SSU" and still sets m; SSUplex not finding MRE SSU
fails its check and leaves m unset; a tool failure writes nothing and keeps finished work; a changed
read file, the T2 index, a missing manifest, step 1 not finished and an unlocked tool version are refused.
"""
import contextlib
import csv
import gzip
import hashlib
import io
import json
import os
import random
import shutil
import stat
import sys
import tarfile
import tempfile
from collections import Counter
from pathlib import Path

here = Path(__file__).parent
sys.path.insert(0, str(here.parent / "scripts"))
import get_marker_databases as DB  # noqa: E402
import make_index as MI  # noqa: E402
import run_stage0a as S  # noqa: E402
import run_stage0a_m as M  # noqa: E402

MM2 = os.environ.get("MINIMAP2") or shutil.which("minimap2")
if not MM2:
    print("SKIP: minimap2 not found (activate phase0-reads, or set MINIMAP2)")
    sys.exit(0)

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)


rng = random.Random(7)
COMP = str.maketrans("ACGT", "TGCA")


def rand(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def rc(s):
    return s.translate(COMP)[::-1]


MOTIFS = {k: rand(30) for k in ("B1", "B2", "B3", "A", "E", "M", "C", "K1", "K2", "K3")}
SSU_ORIGIN = {"B1": "bacteria", "B2": "bacteria", "B3": "bacteria", "A": "archaea", "E": "eukaryota",
              "M": "mitochondria", "C": "chloroplast"}
MARKERS = ["K1", "K2", "K3"]


def genome(n, pattern):
    """Random sequence with a motif block every 300 bp, cycling through pattern (lists of motif names,
    placed 40 bp apart, so a read can hold both motifs of a block)."""
    s = list(rand(n))
    for j, p in enumerate(range(100, n - 200, 300)):
        for k, name in enumerate(pattern[j % len(pattern)]):
            s[p + 40 * k:p + 40 * k + 30] = MOTIFS[name]
    return "".join(s)


def bact(b, k):
    return [[b], [k], [b, k]]


GENOMES = {"host": ("GCA_900000201.1", "Rhizophagus testii", genome(60000, [["E"], [], []])),
           "cagg": ("GCA_800000201.1", "Candidatus Glomeribacter testii", genome(20000, bact("B1", "K1"))),
           "mre": ("GCA_800000202.1", "Candidatus Moeniiplasma testii", genome(20000, bact("B2", "K2"))),
           "cont": ("GCF_700000201.1", "Ralstonia testii", genome(10000, bact("B3", "K3")))}
NOVEL = genome(10000, bact("B3", "K3"))          # a bacterium with no reference: its reads are routed
ORG = ("NC_900201.1", genome(5000, [["C"], [], []]))
T = "\t"


def stand_in(path, body):
    path.write_text(f"#!{sys.executable}\n" + body)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return str(path)


def quiet(fn, argv):
    err = io.StringIO()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
        code = fn(argv)
    return code, err.getvalue()


def build_world(tmp):
    refdir, work = tmp / "references", tmp / "work"
    (work / "clean" / "hosts").mkdir(parents=True)
    (work / "clean" / "bacteria").mkdir(parents=True)
    refdir.mkdir()
    for key, (acc, _, seq) in GENOMES.items():
        sub = "hosts" if key == "host" else "bacteria"
        (work / "clean" / sub / f"{acc}.fna").write_text(f">{acc}|c1\n{seq}\n")
    (work / "clean" / "organelles.fna").write_text(f">{ORG[0]}|{ORG[0]}\n{ORG[1]}\n")
    h, c, m, k = (GENOMES[x] for x in ("host", "cagg", "mre", "cont"))
    (refdir / "hosts.tsv").write_text(f"accession{T}organism{T}note\n{h[0]}{T}{h[1]}{T}host\n")
    (refdir / "endosymbionts.tsv").write_text(f"accession{T}organism{T}note\n{c[0]}{T}{c[1]}{T}\n"
                                              f"{m[0]}{T}{m[1]}{T}held_out_in_T2\n")
    (refdir / "contaminants.tsv").write_text(f"accession{T}organism{T}note\n{k[0]}{T}{k[1]}{T}\n")
    (refdir / "organelles.tsv").write_text(f"plant{T}organelle{T}accession{T}title{T}length_bp\n"
                                           f"Plantago testii{T}chloroplast{T}{ORG[0]}{T}x{T}5000\n")
    (refdir / "overlaps.tsv").write_text("\t".join(["kind", "host_assembly", "host_organism", "host_contig",
                                                    "host_length", "host_start", "host_end", "host_covered",
                                                    "bacterial_genome", "bacterial_organism", "bacterial_contig",
                                                    "bacterial_length", "bacterial_start", "bacterial_end",
                                                    "bacterial_covered", "identity", "aligned_bp", "matched"]) + "\n")
    sylph = stand_in(tmp / "sylph", 'import sys\na = sys.argv[1:]\nif a == ["--version"]:\n    print("sylph 1.0.0")\n'
                     '    sys.exit(0)\nopen(a[a.index("-o") + 1] + ".syldb", "w").write("x\\n")\n')
    code, err = quiet(MI.main, ["--workdir", str(work), "--refdir", str(refdir), "--minimap2", MM2,
                                "--sylph", sylph, "--threads", "2"])
    if code != 0:
        print("FAIL\nmake_index: " + err)
        sys.exit(1)
    return refdir, work


# fastq-dump: a spot is "source" (random start) or "source@start"; writes sources.tsv (mate, position, source).
FAKE_FASTQ_DUMP = '''import gzip, json, os, random, sys
a = sys.argv[1:]
if a == ["--version"]:
    print("\\nfastq-dump : 3.4.1\\n"); sys.exit(0)
spec = json.load(open(os.environ["FAKE_SPEC"]))
run, x, out = a[-1], int(a[a.index("-X") + 1]), a[a.index("-O") + 1]
r = spec["runs"][run]
rng = random.Random(run)
comp = str.maketrans("ACGT", "TGCA")
os.makedirs(out, exist_ok=True)
names = [f"{run}_1.fastq.gz"] + ([f"{run}_2.fastq.gz"] if r["layout"] == "PAIRED" else [])
fhs = [gzip.open(os.path.join(out, n), "wt") for n in names]
src_out = open(os.path.join(out, "sources.tsv"), "w")
spots = r["spots"][:x]
for i, spot in enumerate(spots, 1):
    src, _, at = spot.partition("@")
    g = spec["genomes"][src]
    s = int(at) if at else rng.randint(0, len(g) - 400)
    reads = [g[s:s + 150], g[s + 200:s + 350][::-1].translate(comp)]
    for mate, (fh, read) in enumerate(zip(fhs, reads), 1):
        fh.write(f"@{run}.{i} {i} length=150\\n{read}\\n+{run}.{i} {i} length=150\\n{'I' * 150}\\n")
        src_out.write(f"{mate}\\t{i - 1}\\t{src}\\n")
for fh in fhs:
    fh.close()
src_out.close()
sys.stderr.write(f"Read {len(spots)} spots for {run}\\nWritten {len(spots)} spots for {run}\\n")
'''

# SSUplex 0.2.0: same options, outputs and summary format; sorts reads by the SSU motif they hold.
FAKE_SSUPLEX = '''import json, os, sys
a = sys.argv[1:]
if a == ["--version"]:
    print("ssuplex " + os.environ.get("FAKE_SSUPLEX_VERSION", "0.2.0")); sys.exit(0)
opt = lambda k: a[a.index(k) + 1]
inp, prefix, hmm = opt("-i"), opt("-o"), opt("--hmm-dir")
base = os.path.basename(inp)
open(os.environ["FAKE_LOG"], "a").write(f"ssuplex {base}\\n")
if os.environ.get("FAKE_SSUPLEX_FAIL") and base.startswith(os.environ["FAKE_SSUPLEX_FAIL"]):
    sys.exit(2)
origins = ["bacteria", "archaea", "eukaryota", "mitochondria", "chloroplast"]
if sorted(os.listdir(hmm)) != sorted(o + ".hmm" for o in origins):
    sys.exit(3)
outs = [f"{prefix}.{o}.fasta" for o in origins + ["unclassified"]] + [f"{prefix}.extraction.tsv", f"{prefix}.summary.txt"]
if any(os.path.exists(p) for p in outs):
    if "--force" not in a:
        sys.exit(4)
    for p in outs:
        if os.path.exists(p):
            os.remove(p)
motifs = json.load(open(os.environ["FAKE_MOTIFS"]))
mito = set(os.environ.get("FAKE_SSU_AS_MITO", "").split(","))
comp = str.maketrans("ACGT", "TGCA")
reads = []
for line in open(inp):
    line = line.rstrip("\\n")
    if line.startswith(">"):
        reads.append([line[1:].split()[0], ""])
    elif reads:
        reads[-1][1] += line
counts = dict.fromkeys(origins + ["unclassified"], 0)
fhs = {}
ext = open(f"{prefix}.extraction.tsv", "w")
ext.write("read_id\\torigin\\tscore\\tn_regions\\tevalue\\tenv_from\\tenv_to\\n")
for rid, seq in reads:
    origin = "unclassified"
    for k, o in motifs["ssu"].items():
        m = motifs["seqs"][k]
        if m in seq or m.translate(comp)[::-1] in seq:
            origin = "mitochondria" if k in mito else o
            break
    counts[origin] += 1
    if origin not in fhs:
        fhs[origin] = open(f"{prefix}.{origin}.fasta", "w")
    fhs[origin].write(f">{rid} origin={origin} score=50.0 coords=1-{len(seq)}\\n{seq}\\n")
    ext.write(f"{rid}\\t{origin}\\t50.0\\t1\\t1e-10\\t1\\t{len(seq)}\\n")
for fh in list(fhs.values()) + [ext]:
    fh.close()
total, unc = len(reads), counts["unclassified"]
with open(f"{prefix}.summary.txt", "w") as s:
    s.write("SSUplex v0.2.0 run summary\\n===========================\\n")
    s.write(f"input:          {inp}\\noutput prefix:  {prefix}\\nthreads:        4\\nE-value cutoff: 1e-5\\n")
    s.write("min bit score:  0\\nmin domains:    1\\n\\n")
    s.write(f"total input reads:  {total}\\nclassified reads:   {total - unc}\\nunclassified reads: {unc}\\n\\n")
    s.write("Per-origin counts:\\n")
    for o in origins + ["unclassified"]:
        s.write(f"  {o:<14} {counts[o]:>10}\\n")
'''

# conda run -n ENV [--no-capture-output] singlem ...: --version, data --verify-only, pipe.
FAKE_CONDA = '''import hashlib, json, os, sys
a = sys.argv[1:]
if a[:2] != ["run", "-n"]:
    sys.exit(9)
rest = a[3:]
if rest[:1] == ["--no-capture-output"]:
    rest = rest[1:]
if rest[:1] != ["singlem"]:
    sys.exit(9)
cmd = rest[1:]
log = os.environ["FAKE_LOG"]
if cmd == ["--version"]:
    print(os.environ.get("FAKE_SINGLEM_VERSION", "0.21.4")); sys.exit(0)
if cmd == ["data", "--verify-only"]:
    open(log, "a").write("singlem verify\\n")
    p = os.environ["SINGLEM_METAPACKAGE_PATH"]
    c = json.load(open(os.path.join(p, "CONTENTS.json")))
    if c["data_version"] != "6.5.0":
        sys.exit(1)
    for f, m in c["md5sums"].items():
        if hashlib.md5(open(os.path.join(p, f), "rb").read()).hexdigest() != m:
            sys.exit(1)
    sys.exit(0)
if cmd[:1] == ["pipe"]:
    opt = lambda k: cmd[cmd.index(k) + 1]
    fwd, otu, mp = opt("--forward"), opt("--otu-table"), opt("--metapackage")
    base = os.path.basename(fwd)
    open(log, "a").write(f"singlem {base}\\n")
    if os.environ.get("FAKE_SINGLEM_FAIL") and base.startswith(os.environ["FAKE_SINGLEM_FAIL"]):
        sys.exit(5)
    if not os.path.exists(os.path.join(mp, "CONTENTS.json")):
        sys.exit(6)
    motifs = json.load(open(os.environ["FAKE_MOTIFS"]))
    off = set(os.environ.get("FAKE_SM_IGNORE", "").split(","))
    comp = str.maketrans("ACGT", "TGCA")
    hits, name = {}, None
    for line in open(fwd):
        line = line.rstrip("\\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            continue
        for k in motifs["markers"]:
            m = motifs["seqs"][k]
            if k not in off and (m in line or m.translate(comp)[::-1] in line):
                hits.setdefault(k, []).append(name)
    extras = "--output-extras" in cmd
    cols = ["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy"]
    cols += ["read_names", "nucleotides_aligned", "taxonomy_by_known?", "read_unaligned_sequences",
             "equal_best_hit_taxonomies"] if extras else []
    with open(otu, "w") as out:
        out.write("\\t".join(cols) + "\\n")
        for k, names in sorted(hits.items()):
            row = [k, base, "ACG-" * 15, str(len(names)), "1.00", "Root; d__Bacteria"]
            row += [" ".join(sorted(names)), "60", "False", "", ""] if extras else []
            out.write("\\t".join(row) + "\\n")
    sys.exit(0)
sys.exit(9)
'''

FAKE_NHMMER = '''import os
print("# nhmmer :: search a DNA model, alignment, or sequence against a DNA database")
print("# HMMER " + os.environ.get("FAKE_HMMER_VERSION", "3.4") + " (Aug 2023); http://hmmer.org/")
'''

# minimap2 that keeps a copy of what it reads from stdin.
MM2_WRAP = f'''import os, subprocess, sys
real, a = {MM2!r}, sys.argv[1:]
if a[-1:] == ["-"]:
    data = sys.stdin.buffer.read()
    with open(os.environ["FAKE_MM2_LOG"], "ab") as fh:
        fh.write(b"### call\\n" + data)
    sys.exit(subprocess.run([real] + a, input=data).returncode)
os.execv(real, [real] + a)
'''

PROFILES = {"A": 2, "B": 3, "C": 1, "E": 2, "M": 2, "N": 1}
QUICK = [("B1", 12), ("A", 10), ("E", 10), ("M", 9), ("C", 9), (None, 50)]


def add_member(tf, name, data):
    info = tarfile.TarInfo(name)
    info.size = len(data)
    tf.addfile(info, io.BytesIO(data))


def metaxa2_tarball(path, skip=()):
    hmm = {c: "".join(f"HMMER3/f [3.1b2 | February 2015]\nNAME  {c}_region{i}\nLENG  60\n//\n" for i in range(n))
           for c, n in PROFILES.items()}
    test = "".join(f">q{i}_{k}\n{rand(80) + (MOTIFS[k] if k else rand(30)) + rand(90)}\n"
                   for k, n in QUICK for i in range(n))
    with tarfile.open(path, "w:gz") as tf:
        for c, text in hmm.items():
            if c not in skip:
                add_member(tf, f"Metaxa2_2.2.3/metaxa2_db/SSU/HMMs/{c}.hmm", text.encode())
        add_member(tf, "Metaxa2_2.2.3/metaxa2_db/LSU/HMMs/B.hmm", b"HMMER3/f\nNAME  LSU_decoy\n//\n")
        add_member(tf, "Metaxa2_2.2.3/test.fasta", test.encode())
    return hmm


def fake_zenodo(root):
    name = DB.METAPACKAGE + ".smpkg.zb"
    payload = b"fake metapackage payload\n"
    contents = {"data_version": "6.5.0", "zenodo_backpack_version": "0.4.0", "payload_directory": "payload_directory",
                "md5sums": {"payload_directory/data.bin": hashlib.md5(payload).hexdigest()}}
    root.mkdir(parents=True, exist_ok=True)
    archive = root / f"{name}.tar.gz"
    with tarfile.open(archive, "w:gz") as tf:
        add_member(tf, f"{name}/CONTENTS.json", json.dumps(contents).encode())
        add_member(tf, f"{name}/payload_directory/data.bin", payload)
    md5 = hashlib.md5(archive.read_bytes()).hexdigest()
    (root / "records").mkdir(exist_ok=True)
    (root / "records" / DB.METAPACKAGE_DOI.rsplit(".", 1)[-1]).write_text(json.dumps(
        {"files": [{"key": archive.name, "size": archive.stat().st_size, "checksum": f"md5:{md5}",
                    "links": {"self": archive.as_uri()}}], "metadata": {"version": "6.5.0"}}))
    return archive, md5


def motif_free(g, n):
    """n spot starts in g whose two reads hold no whole motif."""
    starts = [s for s in range(0, len(g) - 400, 7)
              if not any(m in g[s:s + 150] or m in g[s + 200:s + 350] for m in MOTIFS.values())]
    return starts[::max(1, len(starts) // n)][:n]


cagg = GENOMES["cagg"][2]
SPEC = {
    "SRRT1DNA": {"layout": "PAIRED", "spots": ["cagg"] * 110 + ["host"] * 100 + ["cont"] * 40 + ["novel"] * 30
                 + ["org"] * 20},
    "SRRT2DNA": {"layout": "PAIRED", "spots": ["mre"] * 100 + ["host"] * 120 + ["cont"] * 40 + ["novel"] * 40},
    "SRRT1SE": {"layout": "SINGLE", "spots": ["cagg"] * 90 + ["host"] * 110},
    # RNA: 40 cagg spots without motifs, 4 whose first read holds B1 (block at 100 + 900j; start 50 bp before)
    "SRRT1RNA": {"layout": "PAIRED", "spots": [f"cagg@{s}" for s in motif_free(cagg, 40)]
                 + [f"cagg@{50 + 900 * j}" for j in range(4)] + ["host"] * 256},
}
ROWS = [["SRRT1DNA", "T1", "DNA", "test B+", "CaGg", "yes", "PRJTEST", "PAIRED", "note 1"],
        ["SRRT1RNA", "T1", "RNA", "test B+", "CaGg", "no", "PRJTEST", "PAIRED", "note 2"],
        ["SRRT2DNA", "T2", "DNA", "test MRE", "CaMg/MRE", "yes", "PRJTEST", "PAIRED", "note 3"],
        ["SRRT1SE", "T1", "DNA", "test B+", "CaGg", "yes", "PRJTEST", "SINGLE", "note 4"]]
ORDER = ["SRRT1DNA", "SRRT2DNA", "SRRT1SE", "SRRT1RNA"]          # DNA positives first, then the rest
SYM = {"CaGg": "cagg", "CaMg/MRE": "mre"}


def expected(work, run, lineage, ssu_off=(), marker_off=()):
    """Count, from the reads and their sources, what step 2 should report for one run."""
    ssu = [MOTIFS[k] for k in ("B1", "B2", "B3") if k not in ssu_off]
    mk = [MOTIFS[k] for k in MARKERS if k not in marker_off]
    rdir = work / "stage0a" / "reads" / run
    src = {}
    for line in open(rdir / "sources.tsv"):
        mate, pos, s = line.split()
        src[(int(mate), int(pos))] = s
    e, found = Counter(), {}
    for mate in (1, 2):
        f = rdir / f"{run}_{mate}.fastq.gz"
        if not f.exists():
            continue
        found[f"{run}_{mate}"] = []
        for pos, seq in enumerate(gzip.open(f, "rt").read().split("\n")[1::4]):
            a = any(m in seq or rc(m) in seq for m in ssu)
            b = any(m in seq or rc(m) in seq for m in mk)
            s = src[(mate, pos)]
            if a or b:
                found[f"{run}_{mate}"].append(pos)
                e[f"found_{s}"] += 1
            e["ssu_reads"] += a
            e["singlem_reads"] += b
            if s == SYM[lineage]:
                e["symbiont_reads"] += 1
                e["ssu_symbiont"] += a
                e["singlem_symbiont"] += b
                e["both_symbiont"] += a and b
                e["route_symbiont"] += a or b
    return e, found


def runs_tsv(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(S.RUN_COLS)
        w.writerows(rows)


def log_lines(path):
    return path.read_text().splitlines() if path.exists() else []


S.SPOTS = 300
INT_COLS = ["symbiont_reads", "ssu_reads", "ssu_symbiont", "singlem_reads", "singlem_symbiont", "both_symbiont",
            "route_symbiont"]

with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    bindir = tmp / "bin"
    bindir.mkdir()
    stand_in(bindir / "nhmmer", FAKE_NHMMER)
    os.environ["PATH"] = f"{bindir}{os.pathsep}{os.environ['PATH']}"
    ssuplex, conda = stand_in(tmp / "ssuplex", FAKE_SSUPLEX), stand_in(tmp / "conda", FAKE_CONDA)
    mm2wrap = stand_in(tmp / "minimap2", MM2_WRAP)
    (tmp / "motifs.json").write_text(json.dumps({"seqs": MOTIFS, "ssu": SSU_ORIGIN, "markers": MARKERS}))
    os.environ.update(FAKE_MOTIFS=str(tmp / "motifs.json"), FAKE_LOG=str(tmp / "tools.log"),
                      FAKE_MM2_LOG=str(tmp / "mm2.log"))

    # ---------------------------------------------------------- step 1, for real
    refdir, work = build_world(tmp)
    (tmp / "spec.json").write_text(json.dumps({"genomes": {**{k: v[2] for k, v in GENOMES.items()},
                                                           "novel": NOVEL, "org": ORG[1]}, "runs": SPEC}))
    os.environ["FAKE_SPEC"] = str(tmp / "spec.json")
    S.sra_cache_dir = lambda: tmp / "ncbi_cache"
    sd = tmp / "stage0a"
    sd.mkdir()
    runs_tsv(sd / "runs.tsv", ROWS)
    code, err = quiet(S.main, ["--workdir", str(work), "--refdir", str(refdir), "--stage0a-dir", str(sd),
                               "--threads", "2", "--jobs", "2", "--fastq-dump",
                               stand_in(tmp / "fastq-dump", FAKE_FASTQ_DUMP), "--minimap2", MM2])
    if code != 0:
        print("FAIL\nstep 1: " + err)
        sys.exit(1)

    # ---------------------------------------------------------- databases
    real = {"name": DB.METAPACKAGE, "doi": DB.METAPACKAGE_DOI, "md5": DB.METAPACKAGE_MD5, "url": DB.METAXA2_URL}
    dbmd = (here.parent / "references" / "DATABASES.md").read_text()
    check(all(v in dbmd for v in real.values()), "get_marker_databases.py constants differ from references/DATABASES.md")
    archive, md5 = fake_zenodo(tmp / "zenodo")
    hmm_src = metaxa2_tarball(tmp / "Metaxa2_2.2.3.tar.gz")
    metaxa2_tarball(tmp / "Metaxa2_noN.tar.gz", skip=("N",))
    DB.ZENODO_API = (tmp / "zenodo" / "records").as_uri() + "/"
    DB.METAPACKAGE_MD5 = md5
    DB.METAXA2_URL = (tmp / "Metaxa2_2.2.3.tar.gz").as_uri()

    def get_db(workdir, extra=()):
        return quiet(DB.main, ["--workdir", str(workdir), "--ssuplex", ssuplex, "--conda", conda, *extra])

    code, err = get_db(work)
    manifest_path = work / "db" / "manifest.json"
    check(code == 0 and manifest_path.exists(), f"get_marker_databases failed: {err}")
    if code == 0:
        man = json.loads(manifest_path.read_text())
        hd = work / "db" / "ssuplex_hmms"
        check((hd / "bacteria.hmm").read_text() == hmm_src["B"] and (hd / "archaea.hmm").read_text() == hmm_src["A"]
              and (hd / "eukaryota.hmm").read_text() == hmm_src["E"]
              and (hd / "chloroplast.hmm").read_text() == hmm_src["C"]
              and (hd / "mitochondria.hmm").read_text() == hmm_src["M"] + hmm_src["N"],
              "profile files are not made as dev/assemble_metaxa2_hmms.sh makes them")
        check({o: v["profiles"] for o, v in man["ssuplex_hmms"].items()} ==
              {"bacteria": 3, "archaea": 2, "eukaryota": 2, "mitochondria": 3, "chloroplast": 1}, "profile counts")
        check(man["versions"] == {"ssuplex": "0.2.0", "hmmer": "3.4", "singlem": "0.21.4"}, f"versions {man['versions']}")
        check(Path(man["singlem_metapackage"]).name == DB.METAPACKAGE + ".smpkg.zb"
              and not list((work / "db" / "singlem").glob("*.tar.gz*")), "metapackage folder, or archive left")
        check(man["ssuplex_quickstart"]["per_origin"] == DB.QUICKSTART, "quick start counts in the manifest")
        check(log_lines(tmp / "tools.log").count("singlem verify") == 1, "singlem data --verify-only not run")
        archive.unlink()                                             # a rerun must not download again
        (tmp / "Metaxa2_2.2.3.tar.gz").rename(tmp / "Metaxa2_away.tar.gz")
        code, err = get_db(work)
        check(code == 0, f"rerun downloaded again or failed: {err}")
        (tmp / "Metaxa2_away.tar.gz").rename(tmp / "Metaxa2_2.2.3.tar.gz")
        archive, md5 = fake_zenodo(tmp / "zenodo")
        DB.METAPACKAGE_MD5 = md5

    def db_fails(name, want, extra_env=None, **patch):
        w = tmp / f"dbfail_{name}"
        old = {k: getattr(DB, k) for k in patch}
        for k, v in patch.items():
            setattr(DB, k, v)
        os.environ.update(extra_env or {})
        code, err = get_db(w)
        for k in extra_env or {}:
            os.environ.pop(k)
        for k, v in old.items():
            setattr(DB, k, v)
        check(code == 1 and want in err and not (w / "db" / "manifest.json").exists(), f"{name}: exit {code}, {err!r}")
        return w

    w = db_fails("md5", "references/DATABASES.md gives", METAPACKAGE_MD5="0" * 32)
    check(not list((w / "db" / "singlem").glob("*.part")), "wrong md5: the download was kept")
    db_fails("no_N", "lacks N.hmm", METAXA2_URL=(tmp / "Metaxa2_noN.tar.gz").as_uri())
    db_fails("quickstart", "quick start gave", extra_env={"FAKE_SSU_AS_MITO": "B1"})
    db_fails("hmmer", "hmmer is 3.3.2", extra_env={"FAKE_HMMER_VERSION": "3.3.2"})
    damaged = tmp / "dbfail_damaged"
    shutil.copytree(work / "db", damaged / "db")
    (damaged / "db" / "manifest.json").unlink()
    (damaged / "db" / "singlem" / (DB.METAPACKAGE + ".smpkg.zb") / "payload_directory" / "data.bin").write_text("x")
    code, err = get_db(damaged)
    check(code == 1 and "verify-only failed" in err and not (damaged / "db" / "manifest.json").exists(),
          f"damaged metapackage: exit {code}, {err!r}")

    # ---------------------------------------------------------- step 2
    base_work, base_sd = tmp / "base_work", tmp / "base_stage0a"
    shutil.copytree(work, base_work)
    shutil.copytree(sd, base_sd)

    def scenario(name):
        w, s = tmp / f"work_{name}", tmp / f"stage0a_{name}"
        shutil.copytree(base_work, w)
        shutil.copytree(base_sd, s)
        return w, s

    def step2(w, s, extra=()):
        return quiet(M.main, ["--workdir", str(w), "--refdir", str(refdir), "--stage0a-dir", str(s), "--threads", "2",
                              "--minimap2", mm2wrap, "--ssuplex", ssuplex, "--conda", conda, *extra])

    def share_rows(s):
        return {r["run"]: r for r in csv.DictReader(open(s / "marker_share.tsv"), delimiter="\t")}

    def compare(s, w, ssu_off=(), marker_off=(), label=""):
        rows, exp = share_rows(s), {}
        for run, *_ in ROWS:
            lin = next(r[4] for r in ROWS if r[0] == run)
            e, _ = expected(base_work, run, lin, ssu_off, marker_off)
            exp[run] = e
            for col in INT_COLS:
                check(int(rows[run][col]) == e[col], f"{label}{run} {col}: {rows[run][col]}, want {e[col]}")
            share = e["route_symbiont"] / e["symbiont_reads"]
            check(abs(float(rows[run]["share"]) - share) < 1e-9, f"{label}{run} share {rows[run]['share']}, want {share}")
        return rows, exp

    # A: everything works.
    wA, sA = scenario("A")
    (tmp / "mm2.log").unlink(missing_ok=True)
    calls0 = len(log_lines(tmp / "tools.log"))
    code, err = step2(wA, sA)
    check(code == 0, f"step 2 failed: {err}")
    if code == 0:
        rows, exp = compare(sA, wA)
        fr = {r["run"]: r for r in csv.DictReader(open(base_sd / "fractions.tsv"), delimiter="\t")}
        for run in exp:
            check(int(fr[run]["symbiont_reads"]) == exp[run]["symbiont_reads"], f"{run}: test world and step 1 disagree")
        dna = {r: exp[r]["route_symbiont"] / exp[r]["symbiont_reads"] for r in ("SRRT1DNA", "SRRT2DNA", "SRRT1SE")}
        rna = exp["SRRT1RNA"]["route_symbiont"] / exp["SRRT1RNA"]["symbiont_reads"]
        check(all(exp[r]["ssu_symbiont"] >= 10 for r in dna) and exp["SRRT2DNA"]["singlem_symbiont"] >= 10
              and 0 < rna < min(dna.values()) and all(exp[r]["both_symbiont"] > 0 for r in dna),
              f"test design: the world does not exercise the checks as meant ({dna}, RNA {rna})")
        m, low = min(dna.values()), min(dna, key=dna.get)
        md = (sA / "STAGE0A_M.md").read_text()
        check(f"**m = {m:.4g}**" in md and f"in {low} " in md, f"STAGE0A_M.md does not give m = {m:.4g} from {low}")
        check("SSUplex finds CaGg and MRE SSU: **passed**" in md and "SingleM sees MRE: **yes**" in md, "checks in md")
        check(list(rows) == [r[0] for r in ROWS] and rows["SRRT1RNA"]["sets_m"] == "no", "row order or sets_m")
        found = {}
        for run in exp:
            found.update(expected(base_work, run, next(r[4] for r in ROWS if r[0] == run))[1])
        want_calls = [found[f"{run}_{mate}"] for run in ORDER for mate in (1, 2)
                      if f"{run}_{mate}" in found and found[f"{run}_{mate}"]]
        fed = [[int(line[1:]) for line in chunk.decode().split("\n")[0::4] if line.startswith("@")]
               for chunk in (tmp / "mm2.log").read_bytes().split(b"### call\n")[1:]]
        check(fed == want_calls, "found reads did not reach minimap2 under their step 1 names (positions), in order")
        with gzip.open(wA / "stage0a" / "route" / "SRRT1DNA" / "found.tsv.gz", "rt") as fh:
            outc = Counter((r["outcome"], r["detail"]) for r in csv.DictReader(fh, delimiter="\t"))
        e1 = exp["SRRT1DNA"]
        check(sum(k for (o, _), k in outc.items() if o == "routed") == e1["found_novel"] > 0
              and outc[("lineage", "Ralstonia")] == e1["found_cont"] > 0
              and outc[("lineage", "CaGg")] == e1["route_symbiont"], f"SRRT1DNA found-read outcomes {dict(outc)}")
        leftovers = [p.name for p in (wA / "stage0a" / "route").rglob("*")
                     if p.name.endswith((".extraction.tsv", ".reads.fasta")) or p.name == "tmp"
                     or (".ssuplex." in p.name and p.name.endswith(".fasta") and ".bacteria." not in p.name)]
        check(not leftovers, f"left behind: {leftovers}")
        calls = log_lines(tmp / "tools.log")[calls0:]
        check(sorted(calls) == sorted([f"{tool} {k}.reads.fasta" for k in found for tool in ("ssuplex", "singlem")]),
              f"tool calls {calls}")
        before, n_calls = (sA / "marker_share.tsv").read_bytes(), len(log_lines(tmp / "tools.log"))
        code, err = step2(wA, sA)                                    # rerun: no tool runs again
        check(code == 0 and len(log_lines(tmp / "tools.log")) == n_calls
              and (sA / "marker_share.tsv").read_bytes() == before, f"rerun: exit {code}, {err!r}")

    # B: SingleM does not see MRE.
    wB, sB = scenario("B")
    os.environ["FAKE_SM_IGNORE"] = "K2"
    code, err = step2(wB, sB)
    os.environ.pop("FAKE_SM_IGNORE")
    check(code == 0, f"B: step 2 failed: {err}")
    if code == 0:
        rows, exp = compare(sB, wB, marker_off=("K2",), label="B: ")
        dna = {r: exp[r]["route_symbiont"] / exp[r]["symbiont_reads"] for r in ("SRRT1DNA", "SRRT2DNA", "SRRT1SE")}
        md = (sB / "STAGE0A_M.md").read_text()
        check("SingleM sees MRE: **no, so MRE rest on SSU**" in md and f"**m = {min(dna.values()):.4g}**" in md,
              "B: SingleM not seeing MRE should give 'MRE rest on SSU' and still set m")

    # C: SSUplex calls MRE SSU mitochondrial, so it fails its check.
    wC, sC = scenario("C")
    os.environ["FAKE_SSU_AS_MITO"] = "B2"
    code, err = step2(wC, sC)
    os.environ.pop("FAKE_SSU_AS_MITO")
    check(code == 0, f"C: step 2 failed: {err}")
    if code == 0:
        compare(sC, wC, ssu_off=("B2",), label="C: ")
        md = (sC / "STAGE0A_M.md").read_text()
        check("SSUplex finds CaGg and MRE SSU: **failed**" in md and "**m is not set**" in md and "**m =" not in md,
              "C: SSUplex failing its check should leave m unset")

    # D: SSUplex fails on one mate file; nothing is written, finished work is kept, the rerun finishes.
    wD, sD = scenario("D")
    os.environ["FAKE_SSUPLEX_FAIL"] = "SRRT2DNA_2"
    code, err = step2(wD, sD)
    os.environ.pop("FAKE_SSUPLEX_FAIL")
    route = wD / "stage0a" / "route"
    check(code == 1 and "SRRT2DNA" in err and not (sD / "marker_share.tsv").exists()
          and all((route / r / "route.json").exists() for r in ("SRRT1DNA", "SRRT1SE", "SRRT1RNA"))
          and not (route / "SRRT2DNA" / "route.json").exists() and (route / "SRRT2DNA" / "SRRT2DNA_1.found.tsv").exists(),
          f"D: exit {code}, {err!r}")
    n_calls = len(log_lines(tmp / "tools.log"))
    code, err = step2(wD, sD)
    new = log_lines(tmp / "tools.log")[n_calls:]
    check(code == 0 and (sD / "marker_share.tsv").exists()
          and sorted(new) == ["singlem SRRT2DNA_2.reads.fasta", "ssuplex SRRT2DNA_2.reads.fasta"],
          f"D: rerun after a failure: exit {code}, new tool calls {new}, {err!r}")

    # Refused at the start.
    wE, sE = scenario("E")
    code, err = step2(wE, sE, extra=["--index", str(wE / "index" / "t2")])
    check(code == 1 and "full references" in err, f"T2 index: exit {code}, {err!r}")
    (wE / "db" / "manifest.json").rename(wE / "db" / "manifest.away")
    code, err = step2(wE, sE)
    check(code == 1 and "get_marker_databases.py first" in err, f"missing manifest: exit {code}, {err!r}")
    (wE / "db" / "manifest.away").rename(wE / "db" / "manifest.json")
    (sE / "fractions.tsv").unlink()
    code, err = step2(wE, sE)
    check(code == 1 and "finish step 1 first" in err, f"step 1 not finished: exit {code}, {err!r}")
    shutil.copy(base_sd / "fractions.tsv", sE / "fractions.tsv")
    os.environ["FAKE_SSUPLEX_VERSION"] = "0.2.1"
    code, err = step2(wE, sE)
    os.environ.pop("FAKE_SSUPLEX_VERSION")
    check(code == 1 and "ssuplex is 0.2.1" in err and not (wE / "stage0a" / "route").exists(),
          f"unlocked SSUplex version: exit {code}, {err!r}")

    # A read file that changed after step 1 (last: step 1's assign.json in every copy points at these files).
    with gzip.open(work / "stage0a" / "reads" / "SRRT1DNA" / "SRRT1DNA_1.fastq.gz", "ab") as fh:
        fh.write(b"@extra\nACGT\n+\nIIII\n")
    wF, sF = tmp / "work_F", tmp / "stage0a_F"
    shutil.copytree(base_work, wF)
    shutil.copytree(base_sd, sF)
    code, err = step2(wF, sF)
    check(code == 1 and "not the file step 1 counted" in err and not (sF / "marker_share.tsv").exists(),
          f"changed read file: exit {code}, {err!r}")

if fails:
    print("FAIL\n" + "\n".join(fails))
    sys.exit(1)
print("PASS: databases (locked versions, metapackage from its Zenodo record with the DATABASES.md md5 and SingleM's "
      "own check, profile files as SSUplex's script makes them, quick start as in the README, failures write no "
      "manifest, reruns download nothing); step 2 (per-run counts and shares, routed and other-lineage reads not "
      "counted, m from DNA positives only, found reads reach minimap2 under their step 1 names, big files removed, "
      "rerun calls no tool, SingleM not seeing MRE, SSUplex failing its check, a tool failure writes nothing and "
      "the rerun finishes, refusals)")
