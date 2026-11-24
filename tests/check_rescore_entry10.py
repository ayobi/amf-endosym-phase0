"""Check rescore_m_entry10.py (PROTOCOL.md section 10, entry 10) on synthetic genomes and hand-made step 2
work files, with real minimap2, make_index.py, assign_reads.py and VSEARCH, so every expected number is
known in advance.

Needs minimap2 and VSEARCH (activate phase0-reads, or set MINIMAP2= and VSEARCH=).

Step 1 is run again here with an index that holds a plant chloroplast and a fungal mitochondrion. Step 2's
found reads per run, as planted:
  - symbiont reads with SSU pieces that match the symbiont's genus (kept, short ones too), that tie
    across two genera at 100% (unplaced), whose best hits include a chloroplast or a eukaryote (dropped),
    random pieces of 80 bp (kept, no hit) and of 50 bp (unplaced, too short);
  - symbiont reads that SingleM places in Bacteria (kept) or elsewhere (not);
  - symbiont reads holding 100 bp of the chloroplast (organelle DNA: not found);
  - reads that step 2 counted for the symbiont but that the new index counts for the mitochondrion, and
    a read it did not count but the new index does;
  - found reads of the host (ignored).
Checked: per run, every column of marker_share_entry10.tsv; m and the 4.2 checks; an RNA run uses entry
8's sample (seed 1) of the reads still the symbiont's, plus new ones, and applies the kept share; a rerun
aligns and places nothing again; entry 8's VSEARCH batches are used again; refusals (no
fractions_entry10.tsv, an index from before entry 10, step 1 with another index, found reads that do not
agree with step 1).
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
import assign_reads as AR  # noqa: E402
import get_silva as GS  # noqa: E402
import make_index as MI  # noqa: E402
import rescore_m_entry10 as RE  # noqa: E402

MM2 = os.environ.get("MINIMAP2") or shutil.which("minimap2")
VS = os.environ.get("VSEARCH") or shutil.which("vsearch")
if not (MM2 and VS):
    print("SKIP: minimap2 or VSEARCH not found (activate phase0-reads, or set MINIMAP2= and VSEARCH=)")
    sys.exit(0)

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)


rng = random.Random(23)


def rand(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


T = "\t"
S_A, S_B, CHL, MITO = rand(450), rand(450), rand(5000), rand(6000)
CHLP = CHL[1000:1100]
G = {"host": ("GCA_900000401.1", "Rhizophagus testii", rand(40000)),
     "cagg": ("GCA_800000401.1", "Candidatus Glomeribacter testii",
              rand(5000) + S_A + rand(4550) + CHLP + rand(5000)),      # SSU at 5000-5450, chloroplast 10000-10100
     "mre": ("GCA_800000402.1", "Candidatus Moeniiplasma testii", rand(4000) + S_B + rand(10550)),
     "cont": ("GCF_700000401.1", "Ralstonia testii", rand(10000))}
ORG = [("NC_900401.1", CHL), ("MN900401.1", MITO)]
BURK = "Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Burkholderiales;Burkholderiaceae"
SILVA = [("GLO1.1.450", f"{BURK};Candidatus Glomeribacter;Candidatus Glomeribacter gigasporarum", S_A),
         ("TIEX.1.400", f"{BURK};Mycoavidus;Mycoavidus sp.", rand(100) + S_A[200:400] + rand(100)),
         ("CHLS.1.450", "Bacteria;Bacillati;Cyanobacteriota;Cyanophyceae;Chloroplast;Incertae Sedis;Incertae Sedis;"
                        "Lotus japonicus", S_A[0:100] + rand(350)),
         ("MRE1.1.450", "Bacteria;Bacillati;Mycoplasmatota;Mollicutes;Mycoplasmatales;Mycoplasmataceae;"
                        "Candidatus Moeniiplasma;Candidatus Moeniiplasma glomeromycotorum", S_B),
         ("EUK1.1.400", "Eukaryota;Amorphea;Obazoa;Opisthokonta;Holomycota;Fungi;Glomeromycota;Gigaspora sp.",
          rand(150) + S_B[300:400] + rand(150))] + \
        [(f"DEC{i}.1.900", f"{BURK};Ralstonia;decoy {i}", rand(900)) for i in range(4)]


def stand_in(path, body):
    path.write_text(f"#!{sys.executable}\n" + body)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return str(path)


def quiet(fn, argv):
    err = io.StringIO()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
        code = fn(argv)
    return code, err.getvalue()


def from_genome(key, at=None, n=150):
    g = G[key][2]
    s = rng.randint(0, len(g) - n) if at is None else at
    return g[s:s + n]


def cagg_plain():
    """A CaGg read away from its SSU and from the chloroplast stretch."""
    g = G["cagg"][2]
    s = rng.choice([rng.randint(0, 4800), rng.randint(5500, 9800), rng.randint(10200, len(g) - 150)])
    return g[s:s + 150]


# ------------------------------------------------------------------ the runs ---
# Each found read: (mate, read sequence, ssuplex piece or None, singlem domain or None, step 2's outcome, what it is)

def piece_read(a, b):
    """A CaGg read holding S_A[a:b], with that piece."""
    o = max(0, min(a, 300))
    return G["cagg"][2][5000 + o:5150 + o], S_A[a:b]


def run_a():
    """SRRA, a CaGg DNA run (sets m)."""
    found = []
    for j in range(10):
        found.append(("_1", *piece_read(100 + j, 190 + j), None, "CaGg", "genus"))
    for _ in range(4):
        found.append(("_1", *piece_read(110, 150), None, "CaGg", "genus"))           # a 40 bp genus match
    for j in range(5):
        found.append(("_1", *piece_read(220 + j, 370 + j), None, "CaGg", "several"))
    for j in range(3):
        found.append(("_1", *piece_read(5 + j, 95 + j), None, "CaGg", "organelle"))
    for _ in range(3):
        found.append(("_1", cagg_plain(), rand(80), None, "CaGg", "no hit"))
    for _ in range(3):
        found.append(("_1", cagg_plain(), rand(50), None, "CaGg", "short"))
    for j in range(6):
        found.append(("_2", *piece_read(120 + j, 200 + j), None, "CaGg", "genus"))
    for dom in ["d__Bacteria"] * 8 + ["d__Archaea"] * 3:
        found.append(("_1", cagg_plain(), None, dom, "CaGg", "singlem"))
    for j in range(4):                                     # 100 bp of the chloroplast: organelle DNA
        found.append(("_1", G["cagg"][2][9975 + j:10125 + j], S_A[100:190], "d__Bacteria", "CaGg", "chloroplast"))
    for _ in range(5):                                     # the new index counts these for the mitochondrion
        found.append(("_1", from_genome_seq(MITO), None, "d__Bacteria", "CaGg", "mito"))
    for _ in range(5):
        found.append(("_1", from_genome("host"), rand(120), None, "Rhizophagus testii", "host"))
    return {"plain": {"_1": [cagg_plain() for _ in range(60)] + [from_genome("host") for _ in range(200)],
                      "_2": [cagg_plain() for _ in range(40)] + [from_genome("host") for _ in range(150)]},
            "found": found}


def from_genome_seq(seq, n=150):
    s = rng.randint(0, len(seq) - n)
    return seq[s:s + n]


def run_b():
    """SRRB, an MRE DNA run (sets m)."""
    g = G["mre"][2]
    found = []
    for j in range(15):
        found.append(("_1", g[4000 + 10 * j:4150 + 10 * j], S_B[10 * j:10 * j + 120], None, "CaMg/MRE", "genus"))
    for j in range(4):
        found.append(("_1", g[4300:4450], S_B[310 + j:390 + j], None, "CaMg/MRE", "outside"))
    for dom in ["d__Bacteria"] * 10 + ["d__Eukaryota"] * 2:
        found.append(("_2", from_genome("mre", at=rng.randint(4500, 14000)), None, dom, "CaMg/MRE", "singlem"))
    plain = {"_1": [from_genome("mre", at=rng.randint(4500, 14000)) for _ in range(300)] +
                   [from_genome("host") for _ in range(300)],
             "_2": [from_genome("mre", at=rng.randint(4500, 14000)) for _ in range(100)] +
                   [from_genome("host") for _ in range(200)]}
    return {"plain": plain, "found": found}


def run_c():
    """SRRC, a CaGg RNA run (reported only; sampled)."""
    found = []
    for j in range(12):                                    # pieces SRRA does not have, so they are placed anew
        found.append(("_1", *piece_read(101 + j % 9, 189 + j % 9), None, "CaGg", "genus"))
    for j in range(4):
        found.append(("_1", *piece_read(221 + j, 369 + j), None, "CaGg", "several"))
    for _ in range(2):
        found.append(("_1", from_genome_seq(MITO), S_A[130:200], None, "CaGg", "mito"))
    found.append(("_1", *piece_read(106, 194), None, "routed", "new"))          # step 2 did not count it
    return {"plain": {"_1": [cagg_plain() for _ in range(80)] + [from_genome("host") for _ in range(100)]},
            "found": found}


RUNS = [("SRRA", "T1", "DNA", "CaGg", "yes", "PAIRED", run_a()),
        ("SRRB", "T2", "DNA", "CaMg/MRE", "yes", "PAIRED", run_b()),
        ("SRRC", "T1", "RNA", "CaGg", "no", "SINGLE", run_c())]
SHARE8 = {"SRRA": 0.2, "SRRB": 0.05, "SRRC": 0.3}
SAMPLE = 5


def build(tmp):
    refdir, work, stage = tmp / "references", tmp / "work", tmp / "stage0a"
    for sub in ("clean/hosts", "clean/bacteria", "db/silva"):
        (work / sub).mkdir(parents=True)
    refdir.mkdir()
    stage.mkdir()
    for key, (acc, _, seq) in G.items():
        sub = "hosts" if key == "host" else "bacteria"
        (work / "clean" / sub / f"{acc}.fna").write_text(f">{acc}|c1\n{seq}\n")
    (work / "clean" / "organelles.fna").write_text("".join(f">{a}|{a}\n{s}\n" for a, s in ORG))
    h, c, m, k = (G[x] for x in ("host", "cagg", "mre", "cont"))
    (refdir / "hosts.tsv").write_text(f"accession{T}organism{T}note\n{h[0]}{T}{h[1]}{T}host\n")
    (refdir / "endosymbionts.tsv").write_text(f"accession{T}organism{T}note\n{c[0]}{T}{c[1]}{T}\n"
                                              f"{m[0]}{T}{m[1]}{T}held_out_in_T2\n")
    (refdir / "contaminants.tsv").write_text(f"accession{T}organism{T}note\n{k[0]}{T}{k[1]}{T}\n")
    (refdir / "organelles.tsv").write_text(f"plant{T}organelle{T}accession{T}title{T}length_bp\n"
                                           f"Lotus testii{T}chloroplast{T}{ORG[0][0]}{T}x{T}5000\n")
    (refdir / "fungal_mitochondria.tsv").write_text(
        f"fungus{T}organelle{T}accession{T}title{T}length_bp{T}source{T}query\n"
        f"Gigaspora testii{T}mitochondrion{T}{ORG[1][0]}{T}x{T}6000{T}GenBank{T}q\n")
    (refdir / "overlaps.tsv").write_text("\t".join(AR_OVERLAP_COLS) + "\n")
    sylph = stand_in(tmp / "sylph", 'import sys\na = sys.argv[1:]\nif a == ["--version"]:\n    print("sylph 1.0.0")\n'
                     '    sys.exit(0)\nopen(a[a.index("-o") + 1] + ".syldb", "w").write("x\\n")\n')
    code, err = quiet(MI.main, ["--workdir", str(work), "--refdir", str(refdir), "--minimap2", MM2,
                                "--sylph", sylph, "--threads", "2"])
    if code != 0:
        print("FAIL\nmake_index: " + err)
        sys.exit(1)

    (work / "db" / "silva" / "silva.fasta").write_text("".join(f">{a} {p}\n{s}\n" for a, p, s in SILVA))
    (work / "db" / "silva" / "taxonomy.tsv").write_text("accession\ttaxonomy\n" +
                                                        "".join(f"{a}\t{p}\n" for a, p, _ in SILVA))
    (work / "db" / "silva" / "manifest.json").write_text(json.dumps({
        "file": "silva.fasta", "md5": "0123", "fasta": str(work / "db" / "silva" / "silva.fasta"),
        "taxonomy": str(work / "db" / "silva" / "taxonomy.tsv"), "vsearch": GS.vsearch_version(VS)}))

    with open(stage / "runs.tsv", "w") as fh:
        fh.write("run\ttest\tkind\tline\tsymbiont_lineage\tsets_t\tbioproject\tlayout\tnote\n")
        for run, test, kind, lin, sets, layout, _ in RUNS:
            fh.write(f"{run}\t{test}\t{kind}\tline {run}\t{lin}\t{sets}\tPRJTEST\t{layout}\tnote {run}\n")
    with open(stage / "marker_share_entry8.tsv", "w") as fh:
        fh.write("run\tsets_m\tshare\n" + "".join(f"{r[0]}\t{r[4]}\t{SHARE8[r[0]]}\n" for r in RUNS))
    (stage / "fractions_entry10.tsv").write_text("run\tsets_t\tfraction\nSRRA\tyes\t0.1\n")

    expect = {}
    for run, test, kind, lin, sets, layout, spec in RUNS:
        reads_dir, rdir = work / "stage0a" / "reads" / run, work / "stage0a" / "route" / run
        reads_dir.mkdir(parents=True)
        rdir.mkdir(parents=True)
        mates = sorted(spec["plain"])
        per_mate = {m: [(f"{run}{m}.p{i}", s) for i, s in enumerate(spec["plain"][m])] for m in mates}
        rows = []
        for i, (m, seq, piece, dom, old, what) in enumerate(spec["found"]):
            name = f"{run}{m}.f{i}"
            rows.append((f"{run}{m}", len(per_mate[m]), name, seq, piece, dom, old, what))
            per_mate[m].append((name, seq))
        files = []
        for m in mates:
            f = reads_dir / f"{run}{m}.fastq.gz"
            with gzip.open(f, "wt") as fh:
                for name, s in per_mate[m]:
                    fh.write(f"@{name} x\n{s}\n+\n{'I' * len(s)}\n")
            files.append(str(f))
        code, err = quiet(AR.main, ["--index", str(work / "index" / "full"), "--platform", "short", "--reads", *files,
                                    "--out", str(work / "stage0a_entry10" / "assign" / run), "--refdir", str(refdir),
                                    "--threads", "2", "--minimap2", MM2])
        if code != 0:
            print("FAIL\nassign_reads: " + err)
            sys.exit(1)
        # step 2's work files
        by_mate = {}
        for mate, pos, name, seq, piece, dom, old, what in rows:
            by_mate.setdefault(mate, []).append((pos, name, seq, piece, dom))
        with gzip.open(rdir / "found.tsv.gz", "wt") as fh:
            fh.write("mate\tposition\tread\tssuplex\tsinglem\toutcome\tdetail\n")
            for mate, pos, name, seq, piece, dom, old, what in rows:
                oc = ("routed", "no alignment") if old == "routed" else ("lineage", old)
                fh.write(f"{mate}\t{pos}\t{name}\t{int(piece is not None)}\t{int(dom is not None)}\t{oc[0]}\t{oc[1]}\n")
        for mate in sorted(by_mate):
            with gzip.open(rdir / f"{mate}.found.fastq.gz", "wt") as fh:
                for pos, name, seq, _, _ in by_mate[mate]:
                    fh.write(f"@{pos} {name} x\n{seq}\n+\n{'I' * len(seq)}\n")
            with open(rdir / f"{mate}.ssuplex.bacteria.fasta", "w") as fh:
                for pos, name, seq, piece, _ in by_mate[mate]:
                    if piece is not None:
                        fh.write(f">{name} origin=bacteria score=50.0 coords=1-{len(piece)}\n")
                        fh.writelines(piece[i:i + 80] + "\n" for i in range(0, len(piece), 80))
            with gzip.open(rdir / f"{mate}.singlem.otu.tsv.gz", "wt") as fh:
                fh.write("gene\tsample\tsequence\tnum_hits\tcoverage\ttaxonomy\tread_names\n")
                for dom in ("d__Bacteria", "d__Archaea", "d__Eukaryota"):
                    names = [n for _, n, _, _, d in by_mate[mate] if d == dom]
                    if names:
                        fh.write(f"S3.1\t{mate}\tACG\t{len(names)}\t1\tRoot; {dom}; p__X\t{' '.join(names)}\n")
        expect[run] = expected(run, kind, lin, rows, work)
    return refdir, work, stage, expect


AR_OVERLAP_COLS = ["kind", "host_assembly", "host_organism", "host_contig", "host_length", "host_start", "host_end",
                   "host_covered", "bacterial_genome", "bacterial_organism", "bacterial_contig", "bacterial_length",
                   "bacterial_start", "bacterial_end", "bacterial_covered", "identity", "aligned_bp", "matched"]
VERDICT = {"genus": "kept", "no hit": "kept", "new": "kept", "several": "several", "short": "short",
           "organelle": "organelle", "outside": "outside"}


def expected(run, kind, lin, rows, work):
    """What the run should come to, worked out from the design (step 1's n is read from its counts)."""
    counts = {(r["section"], r["name"]): int(r["reads"]) for r in
              csv.DictReader(open(work / "stage0a_entry10" / "assign" / run / "counts.tsv"), delimiter="\t")}
    n = counts[("lineage", lin)]
    sym = [x for x in rows if x[7] not in ("mito", "host")]
    org = [x for x in sym if x[7] == "chloroplast"]
    ssu = [x for x in sym if x[4] is not None and x[7] != "chloroplast"]
    sm = [x for x in sym if x[5] is not None and x[7] != "chloroplast"] if kind == "DNA" else []
    placed = ssu
    if kind != "DNA":
        old = sorted((x[0], x[2]) for x in rows if x[6] == lin and x[4] is not None)
        sample = set(random.Random(1).sample(old, SAMPLE)) if len(old) > SAMPLE else set(old)
        placed = [x for x in ssu if (x[0], x[2]) in sample or (x[0], x[2]) not in set(old)]
    v = [VERDICT[x[7]] for x in placed]
    kept = v.count("kept")
    if kind == "DNA":
        sm_kept = sum(x[5] == "d__Bacteria" for x in sm)
        ssu_kept, route, sm_found = kept, kept + sm_kept, len(sm)
    else:
        ssu_kept = round(len(ssu) * kept / len(placed))
        sm_kept, route, sm_found = "", ssu_kept, ""
    return {"symbiont_reads": n, "found_symbiont": len(sym), "organelle_dna": len(org), "ssu_found": len(ssu),
            "ssu_placed": len(placed), "ssu_outside_bacteria": v.count("outside"), "ssu_organelle": v.count("organelle"),
            "ssu_several_genera": v.count("several"), "ssu_too_short": v.count("short"), "ssu_kept": ssu_kept,
            "singlem_found": sm_found, "singlem_kept": sm_kept, "route_kept": route, "share": route / n}


# minimap2 and VSEARCH that log their calls
WRAP = '''import os, sys
real, a = {real!r}, sys.argv[1:]
if a != ["--version"]:
    open(os.environ["FAKE_LOG"], "a").write({name!r} + "\\n")
if {name!r} == "vsearch" and os.environ.get("FAKE_VSEARCH_FAIL") and "--usearch_global" in a:
    sys.exit(3)
os.execv(real, [real] + a)
'''

with tempfile.TemporaryDirectory() as td:
    tmp = Path(td)
    refdir, work, stage, expect = build(tmp)
    mm2 = stand_in(tmp / "minimap2", WRAP.format(real=MM2, name="minimap2"))
    vs = stand_in(tmp / "vsearch", WRAP.format(real=VS, name="vsearch"))
    os.environ["FAKE_LOG"] = str(tmp / "calls.log")

    def calls():
        p = tmp / "calls.log"
        return p.read_text().split() if p.exists() else []

    def rescore(w=work, s=stage, extra=()):
        return quiet(RE.main, ["--workdir", str(w), "--refdir", str(refdir), "--stage0a-dir", str(s), "--threads", "2",
                               "--minimap2", mm2, "--vsearch", vs, "--rna-sample", str(SAMPLE), *extra])

    def results(s=stage):
        return {x["run"]: x for x in csv.DictReader(open(s / RE.OUT_TSV), delimiter="\t")}

    code, err = rescore()
    check(code == 0, f"rescore failed: {err}")
    if code == 0:
        got = results()
        for run, e in expect.items():
            for col, want in e.items():
                have = got[run][col]
                ok = abs(float(have) - want) < 1e-9 if col == "share" else have == str(want)
                check(ok, f"{run} {col}: {have}, want {want}")
            check(abs(float(got[run]["share_entry8"]) - SHARE8[run]) < 1e-12, f"{run}: entry 8's share")
        check(expect["SRRA"]["organelle_dna"] == 4 and expect["SRRA"]["ssu_several_genera"] == 5 and
              expect["SRRA"]["ssu_too_short"] == 3 and expect["SRRC"]["ssu_placed"] < expect["SRRC"]["ssu_found"],
              "the design plants every case")
        m = min(expect[r]["share"] for r in ("SRRA", "SRRB"))
        low = min(("SRRA", "SRRB"), key=lambda r: expect[r]["share"])
        md = (stage / RE.OUT_MD).read_text()
        check(f"**m = {m:.4g}** (entry 8: 0.05)" in md and f"in {low} " in md and
              "SSUplex finds CaGg and MRE SSU: **passed**" in md and "SingleM sees MRE: **yes**" in md and
              "SRRC" in md.split("--rna-sample")[1][:60] and "section 10, entries 6, 7, 8 and 10" in md,
              f"{RE.OUT_MD}: m, entry 8's m, the checks or the RNA note")
        first = calls()
        check(first.count("vsearch") == 2 and first.count("minimap2") >= 6,
              f"calls: {first.count('vsearch')} VSEARCH (want 2: DNA, then RNA), {first.count('minimap2')} minimap2")

        # a rerun aligns and places nothing again, and writes the same table
        before = (stage / RE.OUT_TSV).read_bytes()
        n0 = len(calls())
        code, err = rescore()
        check(code == 0 and len(calls()) == n0 and (stage / RE.OUT_TSV).read_bytes() == before,
              f"rerun: exit {code}, new calls {calls()[n0:]}, {err!r}")

        # entry 8's batches are used again: here this run's, put where entry 8 kept them
        w2, s2 = tmp / "work2", tmp / "stage2"
        shutil.copytree(work, w2, symlinks=True)
        shutil.copytree(stage, s2)
        shutil.rmtree(w2 / "stage0a" / RE.WORK_NAME)
        shutil.copytree(work / "stage0a" / RE.WORK_NAME / "batches", w2 / "stage0a" / RE.ENTRY8_WORK / "batches")
        for p in (w2 / "stage0a_entry10" / "assign").glob("*/assign.json"):    # the reads now live in work2
            p.write_text(p.read_text().replace(str(work), str(w2)))
        (s2 / RE.OUT_TSV).unlink()
        n1 = len(calls())
        code, err = rescore(w2, s2)
        check(code == 0 and "vsearch" not in calls()[n1:] and (s2 / RE.OUT_TSV).read_bytes() == before,
              f"entry 8's batches: exit {code}, calls {calls()[n1:]}, {err[-300:]!r}")

    # refusals
    (stage / "fractions_entry10.tsv").rename(stage / "hidden.tsv")
    code, err = rescore()
    check(code == 1 and "run_stage0a.py" in err, f"no fractions_entry10.tsv: {err!r}")
    (stage / "hidden.tsv").rename(stage / "fractions_entry10.tsv")
    org = work / "index" / "full" / "organelles.fna"
    org.rename(org.with_suffix(".hidden"))
    code, err = rescore()
    check(code == 1 and "before entry 10" in err, f"an index from before entry 10: {err!r}")
    org.with_suffix(".hidden").rename(org)

    w3, s3 = tmp / "work3", tmp / "stage3"
    shutil.copytree(work, w3, symlinks=True)
    shutil.copytree(stage, s3)
    (s3 / RE.OUT_TSV).unlink(missing_ok=True)
    shutil.rmtree(w3 / "stage0a" / RE.WORK_NAME)
    aj = w3 / "stage0a_entry10" / "assign" / "SRRB" / "assign.json"
    a = json.loads(aj.read_text())
    a["index"]["refs_md5"] = "0" * 32
    aj.write_text(json.dumps(a))
    code, err = rescore(w3, s3)
    check(code == 1 and "another platform, minimap2, index" in err and not (s3 / RE.OUT_TSV).exists(),
          f"step 1 with another index: {err!r}")

    w4, s4 = tmp / "work4", tmp / "stage4"
    shutil.copytree(work, w4, symlinks=True)
    shutil.copytree(stage, s4)
    shutil.rmtree(w4 / "stage0a" / RE.WORK_NAME)
    ct = w4 / "stage0a_entry10" / "assign" / "SRRA" / "counts.tsv"       # 5 found reads are the mitochondrion's
    rows = list(csv.DictReader(open(ct), delimiter="\t"))
    for r in rows:
        if (r["section"], r["name"]) == ("lineage", "organelles"):
            r["reads"] = "1"
    with open(ct, "w") as fh:
        fh.write("section\tname\treads\n" + "".join(f"{r['section']}\t{r['name']}\t{r['reads']}\n" for r in rows))
    code, err = rescore(w4, s4)
    check(code == 1 and "do not agree with step 1" in err, f"found reads that disagree with step 1: {err!r}")

if fails:
    print("FAIL\n" + "\n".join(fails))
    sys.exit(1)
print("PASS: rescore_m_entry10 (found reads aligned again with the entry 10 index, organelle DNA left out, SSU "
      "pieces kept, dropped and unplaced as T1 labels them, SingleM in Bacteria only, every column per run, m and "
      "the 4.2 checks, RNA sampling as in entry 8, a rerun aligns and places nothing, entry 8's VSEARCH batches "
      "used again, refusals)")
