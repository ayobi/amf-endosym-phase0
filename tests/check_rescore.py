"""Check get_silva.py and rescore_stage0a_m.py with real VSEARCH on a small SILVA-like database and
synthetic step 2 work files, so every expected number is known in advance.

Needs VSEARCH (activate phase0-reads, or set VSEARCH=/path/to/vsearch).

The database holds two bacterial sequences (stand-ins for CaGg and an MRE relative), and records that
share a stretch with one of them: a chloroplast, a mitochondrion, a fungus and an archaeon. Symbiont
reads are cut from the bacterial sequences: unique stretches (kept, on the plus or the minus strand),
stretches shared with an organelle (dropped), stretches shared with a fungus or an archaeon (dropped),
and random sequence with no hit (kept). Two more records, a foraminifer and an archaeon, end or start
with the first or last 25-30 bases of some reads: under VSEARCH's default identity these short overlaps
score 100% and the reads were dropped; with --iddef 0 (entry 8) they are kept.

get_silva.py: the release folder is tried before the "current" one; the md5 and the sequence count are
checked; a wrong md5, a wrong count or an unlocked VSEARCH writes no manifest; a rerun downloads nothing.
rescore_stage0a_m.py: per run, the SSU reads kept and dropped (by reason), the SingleM reads kept
(Bacteria only, and none in RNA runs), the share and m; reads on the minus strand are found (--strand
both); an RNA run is sampled and the kept share applied; the two 4.2 checks, passing and failing; a rerun
runs no VSEARCH; a VSEARCH failure writes nothing, keeps finished batches, and the rerun finishes; an
unlocked VSEARCH, a missing manifest and work files that disagree with step 2 are refused.
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
import tempfile
from pathlib import Path

here = Path(__file__).parent
sys.path.insert(0, str(here.parent / "scripts"))
import get_silva as GS  # noqa: E402
import rescore_stage0a_m as R  # noqa: E402

VSEARCH = os.environ.get("VSEARCH") or shutil.which("vsearch")
if not VSEARCH:
    print("SKIP: vsearch not found (activate phase0-reads, or set VSEARCH)")
    sys.exit(0)

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)


rng = random.Random(5)
COMP = str.maketrans("ACGT", "TGCA")


def rand(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def rc(s):
    return s.translate(COMP)[::-1]


A, B = rand(1500), rand(1500)                      # stand-ins for a CaGg 16S and an MRE-relative 16S
SILVA = [("BAC1.1.1500", "Bacteria;Pseudomonadota;Gammaproteobacteria;Burkholderiales;Burkholderiaceae;Glomeribacter;"
          "Candidatus Glomeribacter gigasporarum", A),
         ("BAC2.1.1500", "Bacteria;Bacillota;Bacilli;Mycoplasmatales;Mycoplasmataceae;Mycoplasma;Mycoplasma sp.", B),
         ("CHL1.1.1500", "Bacteria;Cyanobacteria;Cyanobacteriia;Chloroplast;Lotus japonicus", A[:600] + rand(900)),
         ("MIT1.1.1500", "Bacteria;Pseudomonadota;Alphaproteobacteria;Rickettsiales;Mitochondria;Zea mays",
          rand(800) + B[800:]),
         ("EUK1.1.1500", "Eukaryota;Amorphea;Obazoa;Opisthokonta;Nucletmycea;Fungi;Gigaspora margarita",
          rand(1200) + A[1200:]),
         ("ARC1.1.1500", "Archaea;Thermoproteota;Nitrososphaeria;Nitrososphaerales;Nitrososphaeraceae;Nitrososphaera;"
          "Nitrososphaera sp.", rand(300) + B[300:450] + rand(1050)),
         ("EUK2.1.1800", "Eukaryota;Amorphea;Obazoa;Opisthokonta;Holozoa;Metazoa;Homo sapiens", rand(1800))]
SILVA += [("EUK3.1.1325", "Eukaryota;SAR;Rhizaria;Retaria;Foraminifera;Globothalamea;Rotaliida;Cibicidoides sp.",
           rand(1300) + A[900:925]),                # ends with the first 25 bases of A's overlap reads
          ("ARC2.1.1230", "Archaea;Methanobacteriati--other;Halobacteriota;Methanosarcinia;Methanosarcinales;"
           "Methanosarcinaceae;Methanosarcina;metagenome", B[200:230] + rand(1200))]   # starts with B's reads' last 30

# Where each kind of segment starts (the stretch a segment comes from is shared with nothing, or with exactly
# the record named), and the verdict it should get. A: unique 600-1200; B: unique 0-300 and 450-800.
# "overlap" reads start at A[900] or end at B[230], where EUK3 ends or ARC2 starts (see segment()).
STARTS = {"keep+": {"A": 700, "B": 480}, "keep-": {"A": 1000, "B": 500},
          "organelle": {"A": 100, "B": 900}, "outside": {"A": 1250, "B": 300}}
VERDICT = {"keep+": "Bacteria", "keep-": "Bacteria", "organelle": "organelle", "organelle-": "organelle",
           "outside": "outside Bacteria", "outside-": "outside Bacteria", "nohit": "no hit", "overlap": "Bacteria"}
KEEP = ("Bacteria", "no hit")      # entry 7; deliberately not taken from the script


def segment(src, kind, length):
    """An SSU segment of the given kind, cut from A (CaGg) or B (MRE). Runs use different lengths, so no
    two runs share a sequence."""
    if kind == "nohit":
        return rand(length)
    if kind == "overlap":            # the read's first 25 bases end EUK3, or its last 30 start ARC2
        return A[900:900 + length] if src is A else B[230 - length:230]
    key, base = ("A" if src is A else "B"), kind.rstrip("-") if kind != "keep-" else kind
    start = STARTS[base][key] + (0 if (base, key) == ("outside", "B") else rng.randint(0, 40))
    seq = src[start:start + length]
    return rc(seq) if kind.endswith("-") else seq


def stand_in(path, body):
    path.write_text(f"#!{sys.executable}\n" + body)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return str(path)


def quiet(fn, argv):
    err = io.StringIO()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
        code = fn(argv)
    return code, err.getvalue()


# VSEARCH that logs its calls, can fail on a named batch, and can claim another version.
WRAP = f'''import os, sys
real, a = {VSEARCH!r}, sys.argv[1:]
if a == ["--version"] and os.environ.get("FAKE_VSEARCH_VERSION"):
    sys.stderr.write("vsearch v" + os.environ["FAKE_VSEARCH_VERSION"] + "_linux_x86_64\\n"); sys.exit(0)
if "--usearch_global" in a:
    q = a[a.index("--usearch_global") + 1]
    open(os.environ["FAKE_LOG"], "a").write(os.path.basename(q) + "\\n")
    if os.environ.get("FAKE_VSEARCH_FAIL") and os.environ["FAKE_VSEARCH_FAIL"] in q:
        sys.exit(3)
os.execv(real, [real] + a)
'''

RUNS = [  # run, test, kind, lineage, sets, n (step 1), source, SSU kinds, SingleM domains
    ("SRRA", "T1", "DNA", "CaGg", "yes", 10000, A,
     {"keep+": 12, "keep-": 6, "organelle": 4, "organelle-": 2, "outside": 3, "outside-": 2, "nohit": 2, "overlap": 3},
     {"d__Bacteria": 40, "d__Archaea": 5, "d__Eukaryota": 3, "": 2}),
    ("SRRB", "T2", "DNA", "CaMg/MRE", "yes", 20000, B,
     {"keep+": 15, "organelle": 5, "organelle-": 3, "outside": 4, "overlap": 3},
     {"d__Bacteria": 30, "d__Archaea": 4, "": 1}),
    ("SRRC", "T1", "RNA", "CaGg", "no", 500, A,
     {"keep+": 20, "keep-": 10, "organelle": 10, "outside": 5, "nohit": 5, "overlap": 4},
     {"d__Bacteria": 3}),
]
RUN_COLS = ["run", "test", "kind", "line", "symbiont_lineage", "sets_t", "bioproject", "layout", "note"]


def build_route(work, stage_dir, runs):
    """Step 2's work files and outputs for the given runs; returns what each run should come to."""
    stage_dir.mkdir(parents=True, exist_ok=True)
    expect, share_rows = {}, []
    with open(stage_dir / "runs.tsv", "w") as fh:
        fh.write("\t".join(RUN_COLS) + "\n")
        for run, test, kind, lin, sets, *_ in runs:
            fh.write(f"{run}\t{test}\t{kind}\tline {run}\t{lin}\t{sets}\tPRJTEST\tPAIRED\tnote {run}\n")
    for i, (run, test, kind, lin, sets, n, src, ssu_kinds, sm_domains) in enumerate(runs):
        rdir = work / "stage0a" / "route" / run
        rdir.mkdir(parents=True, exist_ok=True)
        found, fasta, otu = [], {f"{run}_1": [], f"{run}_2": []}, {f"{run}_1": {}, f"{run}_2": {}}
        pos, ssu_reads, length = 0, [], 150 - 5 * i
        for kind_, count in ssu_kinds.items():
            for _ in range(count):
                mate, name = f"{run}_{1 + pos % 2}", f"{run}.{pos}"
                seq = segment(src, kind_, length)
                fasta[mate].append((name, seq))
                found.append((mate, pos, name, 1, 0, "lineage", lin))
                ssu_reads.append(((mate, name), VERDICT[kind_]))
                pos += 1
        sm_reads = []
        for dom, count in sm_domains.items():
            for _ in range(count):
                mate, name = f"{run}_{1 + pos % 2}", f"{run}.{pos}"
                otu[mate].setdefault(dom, []).append(name)
                found.append((mate, pos, name, 0, 1, "lineage", lin))
                sm_reads.append(dom)
                pos += 1
        for j in range(7):                       # found reads that are not the symbiont's: ignored
            mate, name = f"{run}_1", f"{run}.{pos}"
            fasta[mate].append((name, segment(src, "keep+", length)))
            found.append((mate, pos, name, 1, 0, "routed" if j % 2 else "lineage", "no alignment" if j % 2 else "Ralstonia"))
            pos += 1
        with gzip.open(rdir / "found.tsv.gz", "wt") as fh:
            fh.write("mate\tposition\tread\tssuplex\tsinglem\toutcome\tdetail\n")
            fh.writelines(f"{m}\t{p}\t{x}\t{a}\t{b}\t{o}\t{d}\n" for m, p, x, a, b, o, d in found)
        for mate, recs in fasta.items():         # SSUplex writes 80 bases a line
            with open(rdir / f"{mate}.ssuplex.bacteria.fasta", "w") as fh:
                for name, seq in recs:
                    fh.write(f">{name} origin=bacteria score=50.0 coords=1-{len(seq)}\n")
                    fh.writelines(seq[i:i + 80] + "\n" for i in range(0, len(seq), 80))
        for mate, groups in otu.items():
            with gzip.open(rdir / f"{mate}.singlem.otu.tsv.gz", "wt") as fh:
                fh.write("gene\tsample\tsequence\tnum_hits\tcoverage\ttaxonomy\tread_names\n")
                for dom, names in groups.items():
                    tax = "Root" + (f"; {dom}; p__X" if dom else "")
                    fh.write(f"S3.1\t{mate}\tACG\t{len(names)}\t1\t{tax}\t{' '.join(names)}\n")
                fh.write("S3.2\tx\tACG\t1\t1\tRoot; d__Eukaryota; Fungi\tnot.a.symbiont.read\n")
        ssu_n, sm_n = sum(ssu_kinds.values()), sum(sm_domains.values())
        row = {"run": run, "test": test, "kind": kind, "line": f"line {run}", "symbiont_lineage": lin,
               "note": f"note {run}", "sets_m": sets, "symbiont_reads": n, "ssu_reads": ssu_n + 7, "ssu_symbiont": ssu_n,
               "singlem_reads": sm_n, "singlem_symbiont": sm_n, "both_symbiont": 0, "route_symbiont": ssu_n + sm_n,
               "share": (ssu_n + sm_n) / n}
        (rdir / "route.json").write_text(json.dumps({"row": row}))
        share_rows.append(row)
        # expected, counted here from the design
        placed = sorted(k for k, _ in ssu_reads)
        if kind == "RNA" and SAMPLE and len(placed) > SAMPLE:
            placed = sorted(random.Random(1).sample(placed, SAMPLE))
        verdict_of = dict(ssu_reads)
        kept = sum(verdict_of[k] in KEEP for k in placed)
        if kind == "DNA":
            ssu_kept, sm_kept = kept, sm_reads.count("d__Bacteria")
            route = ssu_kept + sm_kept
        else:
            ssu_kept, sm_kept = round(ssu_n * kept / len(placed)), ""
            route = ssu_kept
        expect[run] = {"ssu_step2": ssu_n, "ssu_placed": len(placed), "ssu_kept": ssu_kept, "singlem_kept": sm_kept,
                       "ssu_outside_bacteria": sum(verdict_of[k] == "outside Bacteria" for k in placed),
                       "ssu_organelle": sum(verdict_of[k] == "organelle" for k in placed),
                       "ssu_no_hit": sum(verdict_of[k] == "no hit" for k in placed),
                       "route_kept": route, "share": route / n, "share_step2": (ssu_n + sm_n) / n}
    with open(stage_dir / "marker_share.tsv", "w", newline="") as fh:
        cols = ["run", "test", "kind", "line", "symbiont_lineage", "sets_m", "symbiont_reads", "ssu_reads",
                "ssu_symbiont", "singlem_reads", "singlem_symbiont", "both_symbiont", "route_symbiont", "share", "note"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(share_rows)
    return expect


SAMPLE = 20

with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    vs = stand_in(tmp / "vsearch", WRAP)
    os.environ["FAKE_LOG"] = str(tmp / "vsearch_calls.log")
    calls = lambda: (tmp / "vsearch_calls.log").read_text().split() if (tmp / "vsearch_calls.log").exists() else []

    # ------------------------------------------------------------- get_silva.py
    site = tmp / "silva_site" / "current"
    site.mkdir(parents=True)
    fasta_text = "".join(f">{acc} {tax}\n" + "".join(seq.replace("T", "U")[i:i + 60] + "\n" for i in range(0, len(seq), 60))
                         for acc, tax, seq in SILVA)
    with gzip.open(site / GS.FILE, "wt") as fh:
        fh.write(fasta_text)
    md5 = hashlib.md5((site / GS.FILE).read_bytes()).hexdigest()
    (site / (GS.FILE + ".md5")).write_text(f"{md5}  {GS.FILE}\n")
    real = {"file": GS.FILE, "count": GS.EXPECTED_SEQUENCES}
    dbmd = (here.parent / "references" / "DATABASES.md").read_text()
    check("release 144" in dbmd and f"{real['count']:,}" in dbmd and real["file"].startswith("SILVA_144_SSURef_NR99"),
          "get_silva.py's release or sequence count differs from references/DATABASES.md")
    GS.BASES = [(tmp / "silva_site" / "release_144").as_uri() + "/", site.as_uri() + "/"]   # first one is missing
    GS.EXPECTED_SEQUENCES = len(SILVA)
    work = tmp / "work"
    code, err = quiet(GS.main, ["--workdir", str(work), "--vsearch", vs])
    man_path = work / "db" / "silva" / "manifest.json"
    check(code == 0 and man_path.exists(), f"get_silva failed: {err}")
    man = json.loads(man_path.read_text())
    check(man["md5"] == md5 and man["sequences"] == len(SILVA) and man["url"].startswith(site.as_uri()),
          f"manifest {man}")
    tax_rows = (work / "db" / "silva" / "taxonomy.tsv").read_text().splitlines()
    check(tax_rows[1] == f"{SILVA[0][0]}\t{SILVA[0][1]}" and len(tax_rows) == len(SILVA) + 1, "taxonomy.tsv")
    (site / GS.FILE).write_bytes(b"not the file")                 # a rerun must not download again
    code, err = quiet(GS.main, ["--workdir", str(work), "--vsearch", vs])
    check(code == 0, f"get_silva rerun downloaded again or failed: {err}")

    def silva_fails(name, want, env=None, count=None, md5_text=None):
        w = tmp / f"silvafail_{name}"
        with gzip.open(site / GS.FILE, "wt") as fh:
            fh.write(fasta_text)
        good = hashlib.md5((site / GS.FILE).read_bytes()).hexdigest()
        (site / (GS.FILE + ".md5")).write_text((md5_text or good) + f"  {GS.FILE}\n")
        old = GS.EXPECTED_SEQUENCES
        GS.EXPECTED_SEQUENCES = count or old
        os.environ.update(env or {})
        code, err = quiet(GS.main, ["--workdir", str(w), "--vsearch", vs])
        for k in env or {}:
            os.environ.pop(k)
        GS.EXPECTED_SEQUENCES = old
        check(code == 1 and want in err and not (w / "db" / "silva" / "manifest.json").exists(),
              f"get_silva {name}: exit {code}, {err!r}")
        return w

    w = silva_fails("md5", "does not match SILVA's md5", md5_text="0" * 32)
    check(not (w / "db" / "silva" / GS.FILE).exists(), "a download with the wrong md5 was kept")
    silva_fails("count", "references/DATABASES.md says", count=len(SILVA) + 1)
    silva_fails("version", "VSEARCH is 2.30.0", env={"FAKE_VSEARCH_VERSION": "2.30.0"})

    # ------------------------------------------------------- rescore_stage0a_m.py
    def rescore(w, s, extra=()):
        return quiet(R.main, ["--workdir", str(w), "--stage0a-dir", str(s), "--threads", "2", "--vsearch", vs,
                              "--rna-sample", str(SAMPLE), *extra])

    def results(s):
        return {x["run"]: x for x in csv.DictReader(open(s / R.OUT_TSV), delimiter="\t")}

    sd = tmp / "stage0a"
    expect = build_route(work, sd, RUNS)
    (sd / "marker_share_entry7.tsv").write_text("run\tsets_m\tshare\nSRRA\tyes\t0.0123\nSRRB\tyes\t0.0042\n"
                                                "SRRC\tno\t0.0001\n")
    n0 = len(calls())
    code, err = rescore(work, sd)
    check(code == 0, f"rescore failed: {err}")
    if code == 0:
        got = results(sd)
        for run, e in expect.items():
            for col, want in e.items():
                have = got[run][col]
                ok = abs(float(have) - want) < 1e-9 if col.startswith("share") else have == str(want)
                check(ok, f"{run} {col}: {have}, want {want}")
        m = min(expect[r]["share"] for r in ("SRRA", "SRRB"))
        md = (sd / R.OUT_MD).read_text()
        check(f"**m = {m:.4g}**" in md and "in SRRB " in md and "SSUplex finds CaGg and MRE SSU: **passed**" in md
              and "SingleM sees MRE: **yes**" in md and "SRRC" in md.split("random sample")[0][-200:]
              and "entry 7: 0.0042)" in md and "--iddef 0" in md,
              f"{R.OUT_MD}: m, entry 7's m, the rule, the checks or the RNA sampling note")
        check(calls()[n0:] == ["batch_0001.fa", "batch_0002.fa"], f"VSEARCH calls {calls()[n0:]} (DNA, then RNA)")
        before, n1 = (sd / R.OUT_TSV).read_bytes(), len(calls())
        code, err = rescore(work, sd)                                # rerun: no VSEARCH
        check(code == 0 and len(calls()) == n1 and (sd / R.OUT_TSV).read_bytes() == before,
              f"rerun: exit {code}, {err!r}")

    def scenario(name, runs):
        w, s = tmp / f"work_{name}", tmp / f"stage0a_{name}"
        shutil.copytree(work / "db", w / "db")
        return w, s, build_route(w, s, runs)

    # SSUplex fails its check: every MRE SSU read is dropped.
    runs = [list(x) for x in RUNS]
    runs[1][7] = {"organelle": 6, "outside": 6}
    wC, sC, eC = scenario("C", [tuple(x) for x in runs])
    code, err = rescore(wC, sC)
    md = (sC / R.OUT_MD).read_text() if code == 0 else ""
    check(code == 0 and "**failed**" in md and "**m is not set**" in md, f"SSUplex failing its check: {err!r}")

    # SingleM does not see MRE: every MRE SingleM read is placed outside Bacteria.
    runs = [list(x) for x in RUNS]
    runs[1][8] = {"d__Archaea": 20, "": 5}
    wD, sD, eD = scenario("D", [tuple(x) for x in runs])
    code, err = rescore(wD, sD)
    md = (sD / R.OUT_MD).read_text() if code == 0 else ""
    m = min(eD[r]["share"] for r in ("SRRA", "SRRB"))
    check(code == 0 and "no, so MRE rest on SSU" in md and f"**m = {m:.4g}**" in md, f"SingleM not seeing MRE: {err!r}")

    # VSEARCH fails on the RNA batch: nothing written, the DNA batch kept, the rerun places only the RNA runs.
    wE, sE, _ = scenario("E", RUNS)
    os.environ["FAKE_VSEARCH_FAIL"] = "batch_0002"
    code, err = rescore(wE, sE)
    os.environ.pop("FAKE_VSEARCH_FAIL")
    bd = wE / "stage0a" / "rescore_entry8" / "batches"
    check(code == 1 and "VSEARCH failed" in err and not (sE / R.OUT_TSV).exists()
          and (bd / "batch_0001.tsv").exists() and not (bd / "batch_0002.fa").exists(), f"VSEARCH failure: {err!r}")
    n2 = len(calls())
    code, err = rescore(wE, sE)
    check(code == 0 and calls()[n2:] == ["batch_0002.fa"] and (sE / R.OUT_TSV).exists(),
          f"rerun after a failure: exit {code}, calls {calls()[n2:]}, {err!r}")

    # Refused at the start.
    os.environ["FAKE_VSEARCH_VERSION"] = "2.30.0"
    code, err = rescore(work, sd)
    os.environ.pop("FAKE_VSEARCH_VERSION")
    check(code == 1 and "Cannot start" in err and "2.30.0" in err, f"unlocked VSEARCH: {err!r}")
    wF, sF, _ = scenario("F", RUNS)
    (wF / "db" / "silva" / "manifest.json").unlink()
    code, err = rescore(wF, sF)
    check(code == 1 and "get_silva.py" in err, f"missing SILVA manifest: {err!r}")
    wG, sG, _ = scenario("G", RUNS)
    rj = wG / "stage0a" / "route" / "SRRA" / "route.json"
    row = json.loads(rj.read_text())
    row["row"]["ssu_symbiont"] += 1
    rj.write_text(json.dumps(row))
    code, err = rescore(wG, sG)
    check(code == 1 and "does not match step 2" in err, f"work files that disagree with step 2: {err!r}")

if fails:
    print("FAIL\n" + "\n".join(fails))
    sys.exit(1)
print("PASS: get_silva (release folder first, md5 and sequence count checked, failures write no manifest, reruns "
      "download nothing); rescore (SSU reads kept and dropped by reason, short end overlaps no longer tie, SingleM "
      "Bacteria only and none in RNA, minus-strand reads found, RNA sampling, shares and m, both 4.2 checks passing "
      "and failing, rerun runs no VSEARCH, a VSEARCH failure writes nothing and the rerun finishes, refusals)")
