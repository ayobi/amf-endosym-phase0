"""Check run_t5.py (PROTOCOL.md section 5, T5; section 10, entry 15).

Needs minimap2 and VSEARCH (activate phase0-reads, or set MINIMAP2= and VSEARCH=). SSUplex, SingleM, fastq-dump
and sylph are stand-ins.

A synthetic pair: three wild-type metatranscriptomes and one wild-type DNA run with a Mycoavidus (in the
references) and an MRE that has no genome in the references, against three cured metatranscriptomes. Two of
the cured runs carry a few Mycoavidus reads, below 2% of the wild type's level, as hopping between libraries
would leave. Checked:
  T5 gives exactly the tables run_t1.py gives for the same runs (the same code path, entry 15)
  with the six metatranscriptomes as one project, the cross-talk check turns the cured runs' Mycoavidus to
  low, so the cured line is unresolved, not resident: as expected. The wild type is resident, in 4 runs
  the MRE comes out through the route and is listed as called in the wild type and not in the cured line
  with each run its own project (as SRA lists them), the cured line is resident: not as expected
  8 chimeric reads in the DNA run (half Mycoavidus SSU, half host) are split reads, as in T1
  a rerun runs no tool; refusals (a control_id missing from controls.tsv, two lines expected resident, three
  lines)
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
import run_t5 as T5  # noqa: E402

MM2 = os.environ.get("MINIMAP2") or shutil.which("minimap2")
VS = os.environ.get("VSEARCH") or shutil.which("vsearch")
if not (MM2 and VS):
    print("SKIP: minimap2 or VSEARCH not found")
    sys.exit(0)

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)


rng = random.Random(53)


def rand(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


S_M, S_A, S_L, S_N = rand(450), rand(450), rand(450), rand(450)
G = {"host": rand(60000), "myco": rand(3000) + S_M + rand(16550), "mreA": rand(3000) + S_A + rand(16550),
     "mreL": rand(3000) + S_L + rand(16550), "nov": rand(2000) + S_N + rand(1000), "cont": rand(10000),
     "org": rand(5000), "mito": rand(6000)}
ACC = {"host": "GCA_900000501.1", "myco": "GCA_800000501.1", "mreA": "GCA_800000502.1", "cont": "GCF_700000501.1",
       "org": "NC_900501.1", "mito": "MN900501.1"}
SSU = [[S_M, "bacteria"], [S_A, "bacteria"], [S_L, "bacteria"], [S_N, "bacteria"]]
MOL = "Bacteria;Bacillati;Mycoplasmatota;Mollicutes;Mycoplasmatales;Mycoplasmataceae"
BURK = "Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Burkholderiales;Burkholderiaceae"
BAC = "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae"
SILVA = [("MYCO.1.450", f"{BURK};Mycoavidus;Mycoavidus cysteinexigens", S_M),
         ("MREA.1.450", f"{MOL};Candidatus Moeniiplasma;Candidatus Moeniiplasma glomeromycotorum", S_A),
         ("MREL.1.450", f"{MOL};Candidatus Moeniiplasma;Mycoplasma-related endobacterium of Linnemannia", S_L)] + \
        [(f"DEC{i}.1.900", f"{BAC};Bacillus;decoy {i}", rand(900)) for i in range(4)]


CHIM = "myco@3000+host@30000"          # a Hi-C-like chimera: the Mycoavidus SSU and host


def at(src, lo, hi, step):
    return [f"{src}@{p}" for p in range(lo, hi, step)]


def spots(*groups):
    out = [x for g in groups for x in g]
    random.Random(len(out)).shuffle(out)
    return out


def rep(x, k):
    return [x] * k


WT_RNA = lambda: spots(at("myco", 200, 19500, 50), rep("mreL@3000", 12), at("mreL", 4000, 19000, 500),  # noqa: E731
                       rep("nov@2000", 8), rep("host", 1564))
RUNS = {"SRRW1": ("PAIRED", WT_RNA()), "SRRW2": ("PAIRED", WT_RNA()), "SRRW3": ("PAIRED", WT_RNA()),
        "SRRWD": ("PAIRED", spots(at("myco", 150, 19500, 60), rep("mreL@3010", 10), rep(CHIM, 8),
                                  rep("host", 1660))),
        "SRRC1": ("PAIRED", spots(rep("myco@5000", 6), rep("host", 1994))),
        "SRRC2": ("PAIRED", spots(rep("myco@9000", 6), rep("host", 1994))),
        "SRRC3": ("PAIRED", spots(rep("host", 2000)))}
GROUP = "PRJNA368402-PRJNA368407"
#       run      unit               control             kind   prep        bioproject     route          SRA project
ROWS = [["SRRW1", "AD073 wild type", "LELO_AD073_WT", "RNA", "unstated", GROUP, "ssu", "PRJW1"],
        ["SRRW2", "AD073 wild type", "LELO_AD073_WT", "RNA", "unstated", GROUP, "ssu", "PRJW2"],
        ["SRRW3", "AD073 wild type", "LELO_AD073_WT", "RNA", "unstated", GROUP, "ssu", "PRJW3"],
        ["SRRWD", "AD073 wild type", "LELO_AD073_WT", "DNA", "dna", "PRJWD", "ssu+singlem", "PRJWD"],
        ["SRRC1", "AD073 cured", "LELO_AD073_CURED", "RNA", "unstated", GROUP, "ssu", "PRJC1"],
        ["SRRC2", "AD073 cured", "LELO_AD073_CURED", "RNA", "unstated", GROUP, "ssu", "PRJC2"],
        ["SRRC3", "AD073 cured", "LELO_AD073_CURED", "RNA", "unstated", GROUP, "ssu", "PRJC3"]]
T_, M_ = 0.003, 0.25


def rows(split=False):
    return [[run, unit, ctrl, "yes", kind, "short", prep, "PAIRED", sra if split else bp, "", "JGI", "AD073", route,
             "", "test"] for run, unit, ctrl, kind, prep, bp, route, sra in ROWS]


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
    if "+" in spot:                      # a chimera, as Hi-C makes: half bacterium, half host
        (sa, a_), (sb, b_) = (p.split("@") for p in spot.split("+"))
        a_, b_ = int(a_), int(b_)
        reads = [g[sa][a_:a_ + 75] + g[sb][b_:b_ + 75], rc(g[sb][b_ + 200:b_ + 350])]
    else:
        src, _, at = spot.partition("@")
        s = int(at) if at else rng.randint(0, len(g[src]) - 400)
        reads = [g[src][s:s + 150], rc(g[src][s + 200:s + 350])]
    for fh, read in zip(fhs, reads):
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
open(os.environ["FAKE_CALLS"], "a").write("ssuplex\\n")
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

FAKE_CONDA = '''import os, sys
a = sys.argv[1:]
rest = a[3:]
if rest[:1] == ["--no-capture-output"]:
    rest = rest[1:]
cmd = rest[1:]
if cmd == ["--version"]:
    print("0.21.4"); sys.exit(0)
open(os.environ["FAKE_CALLS"], "a").write("singlem\\n")
otu = cmd[cmd.index("--otu-table") + 1]
cols = ["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "read_names", "nucleotides_aligned",
        "taxonomy_by_known?", "read_unaligned_sequences", "equal_best_hit_taxonomies"]
open(otu, "w").write("\\t".join(cols) + "\\n")
'''

FAKE_SYLPH = '''import os, sys
a = sys.argv[1:]
if a == ["--version"]:
    print("sylph 1.0.0"); sys.exit(0)
if a[0] == "sketch":
    open(a[a.index("-o") + 1] + ".syldb", "w").write("x\\n"); sys.exit(0)
if "FAKE_CALLS" in os.environ:
    open(os.environ["FAKE_CALLS"], "a").write("sylph\\n")
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
    refdir, work, stage = tmp / "references", tmp / "work", tmp / "stage0a"
    for sub in ("clean/hosts", "clean/bacteria", "db/silva", "db/hmm", "db/mp"):
        (work / sub).mkdir(parents=True)
    refdir.mkdir()
    stage.mkdir()
    (work / "clean" / "hosts" / f"{ACC['host']}.fna").write_text(f">{ACC['host']}|c1\n{G['host']}\n")
    for k in ("myco", "mreA", "cont"):
        (work / "clean" / "bacteria" / f"{ACC[k]}.fna").write_text(f">{ACC[k]}|c1\n{G[k]}\n")
    (work / "clean" / "organelles.fna").write_text(f">{ACC['org']}|{ACC['org']}\n{G['org']}\n"
                                                   f">{ACC['mito']}|{ACC['mito']}\n{G['mito']}\n")
    write_tsv(refdir / "hosts.tsv", ["accession", "organism", "note"], [[ACC["host"], "Linnemannia testii", "host"]])
    write_tsv(refdir / "endosymbionts.tsv", ["accession", "organism", "note"],
              [[ACC["myco"], "Mycoavidus testii", ""], [ACC["mreA"], "Candidatus Moeniiplasma testii A",
                                                       "held_out_in_T2"]])
    write_tsv(refdir / "contaminants.tsv", ["accession", "organism", "note"], [[ACC["cont"], "Ralstonia testii", ""]])
    write_tsv(refdir / "organelles.tsv", ["plant", "organelle", "accession", "title", "length_bp"],
              [["Lotus testii", "chloroplast", ACC["org"], "x", 5000]])
    write_tsv(refdir / "fungal_mitochondria.tsv", ["fungus", "organelle", "accession", "title", "length_bp",
                                                   "source", "query"],
              [["Linnemannia testii", "mitochondrion", ACC["mito"], "x", 6000, "GenBank", "q"]])
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
    os.environ["FAKE_SPEC"], os.environ["FAKE_CALLS"] = str(tmp / "spec.json"), str(tmp / "calls.log")
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
              [["LELO_AD073_WT", "resident"], ["LELO_AD073_CURED", "absent"]])
    return refdir, work, stage, tools


with tempfile.TemporaryDirectory() as td:
    tmp = Path(td)
    refdir, work, stage, tools = build(tmp)
    calls_log = tmp / "calls.log"

    def common(controls=None):
        return ["--workdir", str(work), "--refdir", str(refdir), "--stage0a-dir", str(stage),
                "--controls", str(controls or tmp / "controls.tsv"), "--threads", "2", "--jobs", "2",
                "--fastq-dump", tools["fastq_dump"], "--minimap2", tools["mm2"], "--ssuplex", tools["ssuplex"],
                "--conda", tools["conda"], "--vsearch", tools["vsearch"], "--sylph", tools["sylph"]]

    def t5(name, rows_=None, controls=None):
        d = tmp / name
        d.mkdir(exist_ok=True)
        write_tsv(d / "runs.tsv", T.RUN_COLS, rows_ if rows_ is not None else rows())
        return quiet(T5.main, common(controls) + ["--t5-dir", str(d)])

    code, err = t5("t5")
    check(code == 0, f"run_t5 failed: {err[-2500:]}")
    if code == 0:
        out = tmp / "t5"
        lv = {(x["run"], x["lineage"]): x for x in read_tsv(out / "levels.tsv")}
        cl = {(x["unit"], x["lineage"]): x for x in read_tsv(out / "calls.tsv")}
        md = (out / "T5.md").read_text()
        myco = {r: lv[(r, "Mycoavidus")] for r in ("SRRW1", "SRRW2", "SRRW3", "SRRWD", "SRRC1", "SRRC2", "SRRC3")}
        check(all(myco[r]["level"] == "present" for r in ("SRRW1", "SRRW2", "SRRW3", "SRRWD")) and
              all(myco[r]["level"] == "low" and myco[r]["cross_talk"] == "yes" for r in ("SRRC1", "SRRC2")) and
              myco["SRRC3"]["level"] == "absent" and myco["SRRWD"]["spread"] != "" and
              myco["SRRWD"]["split_reads"] == "8" and
              all(myco[r]["informative"] == "no" for r in myco if r != "SRRWD"),
              f"Mycoavidus levels: {[(r, x['n'], x['level'], x['cross_talk'], x['spread']) for r, x in myco.items()]}")
        cw, cc = cl[("AD073 wild type", "Mycoavidus")], cl[("AD073 cured", "Mycoavidus")]
        check(cw["call"] == "resident" and cw["basis"] == "known endosymbiont clade, present in 4 runs" and
              cc["call"] == "unresolved", f"calls: wild type {cw['call']} ({cw['basis']}), cured {cc['call']}")
        check("## Result: as expected" in md and
              "| AD073 wild type | resident | resident | resident | known endosymbiont clade, present in 4 runs | "
              "4 / 0 / 1 |" in md and
              "| AD073 cured | absent | not resident | unresolved | signal, but no rule above applies | 0 / 2 / 0 |"
              in md, "T5.md: the result table")
        route = [lin for (u, lin), c in cl.items() if u == "AD073 wild type" and c["source"] == "route" and
                 c["call"] in ("resident", "unresolved")]
        check(any("Moeniiplasma" in lin for lin in route) and
              all(cl.get(("AD073 cured", lin), {}).get("call", "none") not in ("resident", "unresolved")
                  for lin in route), f"the MRE through the route: {route}")
        line = next((x for x in md.splitlines() if x.startswith("Called (resident or unresolved) in the wild type")),
                    "")
        check(any(f"{lin} (" in line for lin in route if "Moeniiplasma" in lin) and "Mycoavidus" not in line,
              f"T5.md: the list of lineages called in the wild type only: {line!r}")

        # the same runs through run_t1.py: the same tables
        t1x = tmp / "t1x"
        t1x.mkdir()
        write_tsv(t1x / "runs.tsv", T.RUN_COLS, rows())
        code1, err1 = quiet(T.main, common() + ["--t1-dir", str(t1x)])
        check(code1 == 0, f"run_t1 on the T5 runs failed: {err1[-800:]}")
        if code1 == 0:
            check((t1x / T.OUT_LEVELS).read_text() == (out / "levels.tsv").read_text() and
                  (t1x / T.OUT_CALLS).read_text() == (out / "calls.tsv").read_text(),
                  "T5's tables differ from run_t1.py's on the same runs")

        # a rerun, and each RNA run its own project (as SRA lists them): no tool runs again
        calls_log.unlink(missing_ok=True)
        code, err = t5("t5")
        code_b, err_b = t5("t5b", rows(split=True))
        again = calls_log.read_text().split() if calls_log.exists() else []
        check(code == 0 and code_b == 0 and again == [], f"reruns: exit {code}, {code_b}, tools {again}")
        if code_b == 0:
            cb = {(x["unit"], x["lineage"]): x for x in read_tsv(tmp / "t5b" / "calls.tsv")}
            mdb = (tmp / "t5b" / "T5.md").read_text()
            check(cb[("AD073 cured", "Mycoavidus")]["call"] == "resident" and "## Result: not as expected" in mdb and
                  "the cured line is called resident for Mycoavidus" in mdb,
                  f"split projects: cured {cb[('AD073 cured', 'Mycoavidus')]['call']}")

    # refusals
    c2 = tmp / "controls_short.tsv"
    write_tsv(c2, ["control_id", "expected_call"], [["LELO_AD073_WT", "resident"]])
    code, err = t5("t5c", controls=c2)
    check(code == 1 and "control_id not in controls.tsv: LELO_AD073_CURED" in err, f"missing control: {err!r}")
    c3 = tmp / "controls_both.tsv"
    write_tsv(c3, ["control_id", "expected_call"], [["LELO_AD073_WT", "resident"], ["LELO_AD073_CURED", "resident"]])
    code, err = t5("t5c", controls=c3)
    check(code == 1 and "one line resident and the other absent" in err, f"both resident: {err!r}")
    three = rows() + [["SRRX1", "AD073 other", "LELO_AD073_WT", "yes", "RNA", "short", "unstated", "PAIRED", GROUP,
                       "", "JGI", "AD073", "ssu", "", "test"]]
    code, err = t5("t5c", three)
    check(code == 1 and "two lines" in err, f"three lines: {err!r}")

if fails:
    print("FAIL\n" + "\n".join(fails))
    sys.exit(1)
print("PASS: run_t5 gives the same tables as run_t1.py on the same runs; with the metatranscriptomes as one project, "
      "the cured runs' few Mycoavidus reads fall to low (cross-talk), so the cured line is unresolved and the wild "
      "type resident in 4 runs: as expected; the MRE comes out through the route and is listed as called in the "
      "wild type only; with each run its own project, the cured line is resident: not as expected; a rerun runs no "
      "tool; refusals")
