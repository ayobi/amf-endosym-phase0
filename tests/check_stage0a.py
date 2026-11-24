"""Check run_stage0a.py end to end on synthetic genomes: a stand-in fastq-dump writes reads with known
symbiont counts, then real minimap2 and assign_reads.py count them.

Needs minimap2 (activate phase0-reads, or set MINIMAP2=/path/to/minimap2).
Checks: only the first SPOTS spots are used; fractions per run; t is the lowest DNA-positive
fraction / 10, and an RNA positive with a lower fraction does not set it; single-end and paired
runs; reruns fetch and assign nothing again; a failed fetch writes no outputs and keeps the
finished runs; a DNA positive with no symbiont reads stops before t; RNA rows with sets_t = yes
and the T2 index are refused. Measuring t again (section 10, entry 10: --work-name and --suffix)
uses the fetched reads, assigns into the new folder, writes the suffixed outputs with the spread of
the DNA runs and the earlier t, and leaves the earlier outputs alone; --work-name without --suffix
is refused.
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
import run_stage0a as S  # noqa: E402

MM2 = os.environ.get("MINIMAP2") or shutil.which("minimap2")
if not MM2:
    print("SKIP: minimap2 not found (activate phase0-reads, or set MINIMAP2)")
    sys.exit(0)

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)


rng = random.Random(5)


def rand(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


GENOMES = {"host": ("GCA_900000201.1", "Rhizophagus testii", rand(60000)),
           "cagg": ("GCA_800000201.1", "Candidatus Glomeribacter testii", rand(20000)),
           "mre": ("GCA_800000202.1", "Candidatus Moeniiplasma testii", rand(20000)),
           "cont": ("GCF_700000201.1", "Ralstonia testii", rand(10000))}
ORG = ("NC_900201.1", rand(5000))
T = "\t"


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


def stand_in(path, body):
    path.write_text(f"#!{sys.executable}\n" + body)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return str(path)


def quiet(fn, argv):
    err = io.StringIO()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
        code = fn(argv)
    return code, err.getvalue()


# The stand-in fastq-dump: spot sources per run come from a JSON spec; it honours -X and -O,
# writes RUN_1/RUN_2 for paired runs (RUN_1 only for single-end), logs each call, and fails on request.
FAKE = '''import gzip, json, os, random, sys
a = sys.argv[1:]
if a == ["--version"]:
    print("\\nfastq-dump : 3.4.1\\n"); sys.exit(0)
spec = json.load(open(os.environ["FAKE_SPEC"]))
run, x, out = a[-1], int(a[a.index("-X") + 1]), a[a.index("-O") + 1]
open(os.environ["FAKE_LOG"], "a").write(run + "\\n")
if run in os.environ.get("FAKE_FAIL", "").split(","):
    sys.stderr.write("network trouble\\n"); sys.exit(3)
r = spec["runs"][run]
rng = random.Random(run)
comp = str.maketrans("ACGT", "TGCA")
os.makedirs(out, exist_ok=True)
names = [f"{run}_1.fastq.gz"] + ([f"{run}_2.fastq.gz"] if r["layout"] == "PAIRED" else [])
fhs = [gzip.open(os.path.join(out, n), "wt") for n in names]
spots = r["spots"][:x]
for i, src in enumerate(spots, 1):
    g = spec["genomes"][src]
    s = rng.randint(0, len(g) - 400)
    reads = [g[s:s + 150], g[s + 200:s + 350][::-1].translate(comp)]
    for fh, read in zip(fhs, reads):
        fh.write(f"@{run}.{i} {i} length=150\\n{read}\\n+{run}.{i} {i} length=150\\n{'I' * 150}\\n")
for fh in fhs:
    fh.close()
sys.stderr.write(f"Read {len(spots)} spots for {run}\\nWritten {len(spots)} spots for {run}\\n")
'''


def runs_tsv(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(S.RUN_COLS)
        w.writerows(rows)


def spec_for(runs):
    return {"genomes": {k: v[2] for k, v in GENOMES.items()}, "runs": runs}


S.SPOTS = 300
# Spot sources, in order. Only the first 300 count. Expected symbiont reads = symbiont spots in the
# first 300 x reads per spot.
SPEC = {
    "SRRT1DNA": {"layout": "PAIRED", "spots": ["cagg"] * 30 + ["host"] * 270 + ["cagg"] * 50},   # 30 of 300
    "SRRT2DNA": {"layout": "PAIRED", "spots": ["host"] * 290 + ["mre"] * 10 + ["mre"] * 40},     # 10 of 300
    "SRRT1RNA": {"layout": "PAIRED", "spots": ["host"] * 299 + ["cagg"]},                        # 1 of 300
    "SRRT1SE": {"layout": "SINGLE", "spots": ["host"] * 180 + ["cagg"] * 20},                    # 20 of 200
}
ROWS = [["SRRT1DNA", "T1", "DNA", "test B+", "CaGg", "yes", "PRJTEST", "PAIRED", "x"],
        ["SRRT2DNA", "T2", "DNA", "test MRE", "CaMg/MRE", "yes", "PRJTEST", "PAIRED", "x"],
        ["SRRT1RNA", "T1", "RNA", "test B+", "CaGg", "no", "PRJTEST", "PAIRED", "x"],
        ["SRRT1SE", "T1", "DNA", "test B+", "CaGg", "yes", "PRJTEST", "SINGLE", "x"]]
WANT = {"SRRT1DNA": (600, 60), "SRRT2DNA": (600, 20), "SRRT1RNA": (600, 2), "SRRT1SE": (200, 20)}

with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    refdir, work = build_world(tmp)
    fake = stand_in(tmp / "fastq-dump", FAKE)
    (tmp / "spec.json").write_text(json.dumps(spec_for(SPEC)))
    os.environ.update(FAKE_SPEC=str(tmp / "spec.json"), FAKE_LOG=str(tmp / "calls.log"), FAKE_FAIL="")
    S.sra_cache_dir = lambda: tmp / "ncbi_cache"

    def stage(stage_dir, workdir=work, extra=()):
        return quiet(S.main, ["--workdir", str(workdir), "--refdir", str(refdir), "--stage0a-dir", str(stage_dir),
                              "--threads", "2", "--jobs", "2", "--fastq-dump", fake, "--minimap2", MM2, *extra])

    sd = tmp / "stage0a"
    sd.mkdir()
    runs_tsv(sd / "runs.tsv", ROWS)
    code, err = stage(sd)
    check(code == 0, f"stage 0a failed: {err}")
    if code == 0:
        fr = {r["run"]: r for r in csv.DictReader(open(sd / "fractions.tsv"), delimiter="\t")}
        for run, (reads, sym) in WANT.items():
            got = (int(fr[run]["reads"]), int(fr[run]["symbiont_reads"]))
            check(got == (reads, sym), f"{run}: reads and symbiont reads {got}, want {(reads, sym)}")
            check(abs(float(fr[run]["fraction"]) / (sym / reads) - 1) < 1e-6, f"{run}: fraction {fr[run]['fraction']}")
        check(int(fr["SRRT1DNA"]["spots"]) == 300 and int(fr["SRRT1SE"]["spots"]) == 200, "spots column")
        md = (sd / "STAGE0A.md").read_text()
        t = 20 / 600 / 10
        check(f"**t = {t:.4g}**" in md and "in SRRT2DNA" in md, f"STAGE0A.md does not give t = {t:.4g} from SRRT2DNA")
        check("fetched with fastq-dump 3.4.1" in md and "the first 300 spots" in md.replace(",", ""), "STAGE0A.md versions/spots")
        first_md = md
        calls = (tmp / "calls.log").read_text().split()
        check(sorted(calls) == sorted(SPEC), f"fastq-dump calls: {calls}")
        stamp = {r: (work / "stage0a" / "assign" / r / "counts.tsv").stat().st_mtime for r in SPEC}

        code, err = stage(sd)                                       # rerun: nothing fetched or assigned again
        check(code == 0, f"rerun failed: {err}")
        check((tmp / "calls.log").read_text().split() == calls, "rerun fetched again")
        check(all((work / "stage0a" / "assign" / r / "counts.tsv").stat().st_mtime == stamp[r] for r in SPEC),
              "rerun assigned again")

        # Measuring t again (entry 10): the reads are not fetched again; the runs are assigned into the new folder.
        first_tsv = (sd / "fractions.tsv").read_bytes()
        code, err = stage(sd, extra=["--work-name", "stage0a_entry10", "--suffix", "_entry10"])
        check(code == 0, f"entry 10 run failed: {err}")
        if code == 0:
            check((tmp / "calls.log").read_text().split() == calls, "entry 10 run fetched again")
            check(all((work / "stage0a_entry10" / "assign" / r / "counts.tsv").exists() for r in SPEC),
                  "entry 10 run did not assign into stage0a_entry10/assign")
            check(all((work / "stage0a" / "assign" / r / "counts.tsv").stat().st_mtime == stamp[r] for r in SPEC),
                  "entry 10 run touched the earlier assignments")
            check((sd / "STAGE0A.md").read_text() == first_md and (sd / "fractions.tsv").read_bytes() == first_tsv,
                  "entry 10 run changed the earlier outputs")
            fr10 = {r["run"]: r for r in csv.DictReader(open(sd / "fractions_entry10.tsv"), delimiter="\t")}
            check(all((fr10[r]["reads"], fr10[r]["symbiont_reads"]) == (fr[r]["reads"], fr[r]["symbiont_reads"])
                      for r in SPEC), "entry 10 run: counts differ with the same index")
            dna = [r for r in SPEC if r != "SRRT1RNA"]
            check(all(float(fr10[r]["spread"]) >= 0.5 for r in dna) and fr10["SRRT1RNA"]["spread"] == "",
                  f"entry 10 run: spread {[(r, fr10[r]['spread']) for r in SPEC]}")
            md10 = (sd / "STAGE0A_ENTRY10.md").read_text()
            for w in ("(entry10)", f"**t = {t:.4g}**", f"Before this (`fractions.tsv`): t = {t:.4g}.",
                      "`stage0a_entry10/assign/`", "`fractions_entry10.tsv`", "| Spread |"):
                check(w in md10, f"STAGE0A_ENTRY10.md does not have {w!r}")
        code, err = stage(sd, extra=["--work-name", "stage0a_x"])
        check(code == 1 and "together" in err, f"--work-name without --suffix: exit {code}, {err!r}")

    # A failed fetch: the other runs are done and kept, nothing is written, the partial folder is gone.
    sd2, work2 = tmp / "stage0a_fail", tmp / "work_fail"
    sd2.mkdir()
    shutil.copytree(work / "index", work2 / "index")
    runs_tsv(sd2 / "runs.tsv", ROWS)
    os.environ["FAKE_FAIL"] = "SRRT2DNA"
    code, err = stage(sd2, workdir=work2)
    os.environ["FAKE_FAIL"] = ""
    check(code == 1 and "SRRT2DNA" in err and not (sd2 / "fractions.tsv").exists(),
          f"failed fetch: exit {code}, {err!r}")
    check((work2 / "stage0a" / "assign" / "SRRT1DNA" / "counts.tsv").exists(), "failed fetch: finished runs not kept")
    check(not (work2 / "stage0a" / "reads" / "SRRT2DNA.part").exists(), "failed fetch: partial folder left")
    code, err = stage(sd2, workdir=work2)                           # the rerun finishes the job
    check(code == 0 and (sd2 / "fractions.tsv").exists(), f"rerun after a failed fetch: exit {code}, {err!r}")

    # A DNA positive with no symbiont reads stops before t.
    sd3, work3 = tmp / "stage0a_zero", tmp / "work_zero"
    sd3.mkdir()
    shutil.copytree(work / "index", work3 / "index")
    spec0 = dict(SPEC, SRRT2DNA={"layout": "PAIRED", "spots": ["host"] * 300})
    (tmp / "spec0.json").write_text(json.dumps(spec_for(spec0)))
    os.environ["FAKE_SPEC"] = str(tmp / "spec0.json")
    runs_tsv(sd3 / "runs.tsv", ROWS)
    code, err = stage(sd3, workdir=work3)
    os.environ["FAKE_SPEC"] = str(tmp / "spec.json")
    check(code == 1 and "t would be 0" in err and not (sd3 / "fractions.tsv").exists(),
          f"zero symbiont reads: exit {code}, {err!r}")

    # Refused: an RNA run that sets t, and the T2 index.
    sd4 = tmp / "stage0a_bad"
    sd4.mkdir()
    runs_tsv(sd4 / "runs.tsv", [ROWS[0], ["SRRT1RNA", "T1", "RNA", "test B+", "CaGg", "yes", "PRJTEST", "PAIRED", "x"]])
    code, err = stage(sd4)
    check(code == 1 and "RNA runs do not set t" in err, f"RNA sets_t: exit {code}, {err!r}")
    code, err = stage(sd, extra=["--index", str(work / "index" / "t2")])
    check(code == 1 and "full references" in err, f"T2 index: exit {code}, {err!r}")

if fails:
    print("FAIL\n" + "\n".join(fails))
    sys.exit(1)
print("PASS: first SPOTS spots only, fractions per run (paired and single-end), t from DNA positives only "
      "(a lower RNA fraction does not set it), rerun fetches and assigns nothing, measuring t again into a new "
      "folder with suffixed outputs and spread (earlier outputs untouched), failed fetch keeps finished runs "
      "and writes nothing, zero symbiont reads stops before t, RNA sets_t, the T2 index and a lone --work-name "
      "refused")
