"""Check make_index.py and assign_reads.py end to end: synthetic genomes are cleaned
(clean_references.py, with a fake NCBI), indexed, and simulated reads are assigned with real minimap2.

Needs minimap2 (activate phase0-reads, or set MINIMAP2=/path/to/minimap2). sylph is replaced by a
stand-in that records what it was asked to sketch.

Planted cases (150 bp reads; expected outcome under PROTOCOL.md 4.1 and section 10, entry 4):
  unique host, CaGg, MRE, Mycoavidus, contaminant and organelle reads -> their lineage
  reads from a region identical in two Mycoavidus genomes           -> Mycoavidus (tie within a lineage)
  reads from a 10-copy repeat, all in Mycoavidus                     -> Mycoavidus
  reads from a region identical in the two host assemblies          -> none (tied), but fungal for E2
  reads from a region identical in CaGg and a contaminant           -> none (tied across lineages)
  reads from a 60-copy repeat (Mycoavidus and contaminant)          -> none (51 equally best alignments)
  reads from the shared stretch found by the cleaning, either side  -> none (shared stretch)
  random reads                                                       -> routed (no alignment)
  reads with a bad tail, and half-and-half chimeras                  -> routed (no alignment passes)
Also: the T2 index (MRE reads become routed), paired files and interleaved mates with the same
names (counts unchanged), HiFi reads, the routed reads file, and failure cases.
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
import urllib.parse
import zipfile
from pathlib import Path

here = Path(__file__).parent
sys.path.insert(0, str(here.parent / "scripts"))
import clean_references as C  # noqa: E402
import make_index as MI  # noqa: E402
import assign_reads as A  # noqa: E402

MM2 = os.environ.get("MINIMAP2") or shutil.which("minimap2")
if not MM2:
    print("SKIP: minimap2 not found (activate phase0-reads, or set MINIMAP2)")
    sys.exit(0)

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)


# ---------------------------------------------------------------- unit checks of the rule
L = {b"h1": "Host one", b"h2": "Host two", b"m1": "Mycoavidus", b"m2": "Mycoavidus", b"c1": "Ralstonia"}
H = {"Host one", "Host two"}
MASK = {b"h1": [(1000, 2000)]}


def aln(score, target, ident=100, cov=100, ts=5000):
    return (score, ident, 100, 0, cov, 100, target, ts, ts + cov)


unit = [
    ("no alignment", [], ("routed", "no alignment", False)),
    ("nothing passes", [aln(50, b"m1", ident=90)], ("routed", "no alignment passes the filter", False)),
    ("best fails, another passes", [aln(60, b"m1", ident=90), aln(55, b"m2")],
     ("none", "best alignment fails the filter", False)),
    ("51 equal best", [aln(60, b"m1")] * 51, ("none", "51 or more equally best alignments", False)),
    ("50 equal best, one lineage", [aln(60, b"m1")] * 25 + [aln(60, b"m2")] * 25, ("lineage", "Mycoavidus", False)),
    ("shared stretch", [aln(60, b"h1", ts=1500)], ("none", "overlaps a shared stretch", False)),
    ("tie across lineages", [aln(60, b"m1"), aln(60, b"c1")], ("none", "tied across lineages", False)),
    ("tie between hosts", [aln(60, b"h1"), aln(60, b"h2")], ("none", "tied across lineages", True)),
    ("51 equal best in hosts", [aln(60, b"h1")] * 51, ("none", "51 or more equally best alignments", True)),
    ("unique host", [aln(60, b"h1"), aln(40, b"m1")], ("lineage", "Host one", True)),
    ("95% over 80% passes", [aln(60, b"m1", ident=95, cov=80)], ("lineage", "Mycoavidus", False)),
    ("94% fails", [aln(60, b"m1", ident=94)], ("routed", "no alignment passes the filter", False)),
    ("79% of the read fails", [aln(60, b"m1", cov=79)], ("routed", "no alignment passes the filter", False)),
    ("touching the stretch end is outside it", [aln(60, b"h1", ts=2000)], ("lineage", "Host one", True)),
    ("1 bp inside the stretch", [aln(60, b"h1", ts=1999)], ("none", "overlaps a shared stretch", False)),
]
for name, alns, want in unit:
    got = A.evaluate(alns, L, H, MASK)
    check(got == want, f"rule, {name}: got {got}, want {want}")

paf = (b"0\t100\t0\t100\t+\th1\t9000\t5000\t5100\t100\t100\t60\tAS:i:200\n"
       b"0\t100\t0\t100\t+\th2\t9000\t5000\t5100\t100\t100\t0\tAS:i:200\n"
       b"2\t100\t0\t90\t+\tm1\t9000\t10\t100\t88\t90\t60\tAS:i:150\n"
       b"3\t100\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0\n")
got = [(q, len(a)) for q, a in A.parse_paf(io.BytesIO(paf))]
check(got == [(0, 2), (2, 1), (3, 0)], f"parse_paf: {got}")

# ---------------------------------------------------------------- synthetic genomes
rng = random.Random(20261120)


def rand(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq, rate):
    s = list(seq)
    for i in rng.sample(range(len(s)), int(len(s) * rate)):
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def repeats(unit_seq, n):
    return "".join(unit_seq + rand(300) for _ in range(n))


X_hh, X_bc, X_mm = rand(10000), rand(1000), rand(5000)
REP, REP2, E1 = rand(400), rand(400), rand(25000)
H1c1 = rand(20000) + X_hh + rand(20000)
H1c2 = rand(10000) + mutate(E1[5000:8000], 0.01) + rand(27000)       # the shared stretch
H2c1 = rand(15000) + X_hh + rand(15000)
H2c2 = rand(30000)
B1c1 = rand(10000) + X_bc + rand(19000)
M1c1 = rand(10000) + X_mm + rand(10000) + repeats(REP, 30) + repeats(REP2, 6)
M2c1 = rand(8000) + X_mm + rand(8000) + repeats(REP, 25) + repeats(REP2, 4)
C1c1 = rand(12000) + X_bc + rand(12000) + repeats(REP, 5)
O1, O2 = rand(8000), rand(6000)

H1, H2, E1A, B1, M1, M2, C1 = ("GCA_900000101.1", "GCA_900000102.1", "GCA_800000104.1", "GCA_800000101.1",
                                "GCA_800000102.1", "GCA_800000103.1", "GCF_700000101.1")
GENOMES = {H1: [("H1c1.1", H1c1), ("H1c2.1", H1c2)], H2: [("H2c1.1", H2c1), ("H2c2.1", H2c2)],
           E1A: [("E1c1.1", E1)], B1: [("B1c1.1", B1c1)], M1: [("M1c1.1", M1c1)], M2: [("M2c1.1", M2c1)],
           C1: [("C1c1.1", C1c1)]}
ORGANELLES = {"NC_900101.1": O1, "NC_900102.1": O2}


def wrap(seq, width=80):
    return "".join(seq[i:i + width] + "\n" for i in range(0, len(seq), width))


def fake_http_get(url, params=None, dest=None, tries=4):
    params = params or {}
    if url.endswith("/version"):
        return json.dumps({"version": "18.41.1"}).encode()
    if "/genome/accession/" in url:
        acc = urllib.parse.unquote(url.split("/genome/accession/")[1].split("/")[0])
        fasta = "".join(f">{n} synthetic\n" + wrap(s) for n, s in GENOMES[acc])
        member = f"ncbi_dataset/data/{acc}/{acc}_TEST_genomic.fna"
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as z:
            z.writestr(member, fasta)
            z.writestr("md5sum.txt", f"{hashlib.md5(fasta.encode()).hexdigest()}  {member}\n")
        dest.write_bytes(buf.getvalue())
        return dest
    if url.endswith("efetch.fcgi"):
        return "".join(f">{i} synthetic organelle\n" + wrap(ORGANELLES[i]) for i in params["id"].split(",")).encode()
    raise AssertionError(url)


C.http_get = fake_http_get
T = "\t"


def write_lists(refdir):
    refdir.mkdir(parents=True)
    (refdir / "hosts.tsv").write_text(f"accession{T}organism{T}note\n{H1}{T}Rhizophagus testii{T}host\n"
                                      f"{H2}{T}Rhizophagus sp. TEST2{T}host\n")
    (refdir / "endosymbionts.tsv").write_text(
        f"accession{T}organism{T}note\n{E1A}{T}Candidatus Moeniiplasma testii{T}held_out_in_T2\n"
        f"{B1}{T}Candidatus Glomeribacter testii{T}\n{M1}{T}Mycoavidus testii{T}\n{M2}{T}Mycoavidus sp. T2{T}\n")
    (refdir / "contaminants.tsv").write_text(f"accession{T}organism{T}note\n{C1}{T}Ralstonia testii{T}\n")
    (refdir / "organelles.tsv").write_text(f"plant{T}organelle{T}accession{T}title{T}length_bp\n"
                                           f"Plantago testii{T}chloroplast{T}NC_900101.1{T}x{T}8000\n"
                                           f"Plantago testii{T}mitochondrion{T}NC_900102.1{T}x{T}6000\n")


def quiet(fn, argv):
    err = io.StringIO()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
        code = fn(argv)
    return code, err.getvalue()


def stand_in(path, body):
    path.write_text(f"#!{sys.executable}\n" + body)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return str(path)


SYLPH = '''import sys
a = sys.argv[1:]
if a == ["--version"]:
    print("sylph 1.0.0"); sys.exit(0)
assert a[0] == "sketch" and a[1] == "-g", a
o = a[a.index("-o") + 1]
open(o + ".syldb", "w").write("\\n".join(a[2:a.index("-o")]) + "\\n")
'''

# ---------------------------------------------------------------- simulated reads
COMP = str.maketrans("ACGT", "TGCA")


def take(seq, lo, hi, n, length, label, reads, expect):
    for k in range(n):
        s = rng.randint(lo, hi - length)
        r = seq[s:s + length]
        if rng.random() < 0.5:
            r = r[::-1].translate(COMP)
        reads.append((f"{label}_{k}", r))
        expect.append(label)


def short_reads():
    reads, labels = [], []
    take(H1c1, 0, 20000, 20, 150, "H1u", reads, labels)
    take(H1c1, 20000, 30000, 20, 150, "HH", reads, labels)
    take(H2c2, 0, 30000, 20, 150, "H2u", reads, labels)
    take(H1c2, 10200, 12800, 10, 150, "STRh", reads, labels)
    take(E1, 5200, 7800, 10, 150, "STRb", reads, labels)
    take(E1, 10000, 25000, 20, 150, "E1u", reads, labels)
    take(B1c1, 12000, 29000, 20, 150, "B1u", reads, labels)
    take(B1c1, 10000, 11000, 10, 150, "BC", reads, labels)
    take(M1c1, 0, 10000, 20, 150, "M1u", reads, labels)
    take(M1c1, 10000, 15000, 10, 150, "MM", reads, labels)
    take(REP, 0, 400, 10, 150, "REP", reads, labels)
    take(REP2, 0, 400, 10, 150, "REP2", reads, labels)
    take(C1c1, 0, 12000, 20, 150, "C1u", reads, labels)
    take(O1, 0, 8000, 10, 150, "O1", reads, labels)
    take(O2, 0, 6000, 10, 150, "O2", reads, labels)
    for k in range(20):
        reads.append((f"rand_{k}", rand(150)))
    for k in range(10):                      # 100 bp exact, then a mismatch every 3 bp
        s = rng.randint(15000, 18000)
        r = list(B1c1[s:s + 150])
        for i in range(101, 150, 3):
            r[i] = rng.choice([b for b in "ACGT" if b != r[i]])
        reads.append((f"low_{k}", "".join(r)))
    for k in range(10):                      # half Mycoavidus, half random
        s = rng.randint(2000, 9000)
        reads.append((f"chim_{k}", M1c1[s:s + 75] + rand(75)))
    return reads


def write_fastq(path, reads):
    with gzip.open(path, "wt") as fh:
        for name, s in reads:
            fh.write(f"@{name} extra words\n{s}\n+\n{'I' * len(s)}\n")


def counts(outdir):
    rows = list(csv.DictReader(open(outdir / "counts.tsv"), delimiter="\t"))
    return {(r["section"], r["name"]): int(r["reads"]) for r in rows}


def routed_names(outdir):
    with gzip.open(outdir / "routed.fastq.gz", "rt") as fh:
        return [line[1:].split()[0] for i, line in enumerate(fh) if i % 4 == 0]


FULL = {("total", "reads in input"): 260,
        ("lineage", "Rhizophagus testii"): 20, ("lineage", "Rhizophagus sp. TEST2"): 20,
        ("lineage", "CaMg/MRE"): 20, ("lineage", "CaGg"): 20, ("lineage", "Mycoavidus"): 40,
        ("lineage", "Ralstonia"): 20, ("lineage", "plant organelles"): 20,
        ("fungal", "fungal reads for E2"): 60,
        ("none", "best alignment fails the filter"): 0, ("none", "51 or more equally best alignments"): 10,
        ("none", "overlaps a shared stretch"): 20, ("none", "tied across lineages"): 30,
        ("routed", "no alignment"): 20, ("routed", "no alignment passes the filter"): 20}
T2 = {k: v for k, v in FULL.items() if k != ("lineage", "CaMg/MRE")}
T2[("routed", "no alignment")] = 40
ROUTED = {f"{p}_{k}" for p in ("rand",) for k in range(20)} | {f"{p}_{k}" for p in ("low", "chim") for k in range(10)}

with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    refdir, work = tmp / "references", tmp / "work"
    write_lists(refdir)
    code, err = quiet(C.main, ["--workdir", str(work), "--refdir", str(refdir), "--threads", "2", "--minimap2", MM2])
    check(code == 0, f"cleaning failed: {err}")
    ov = list(csv.DictReader(open(refdir / "overlaps.tsv"), delimiter="\t"))
    check([r["kind"] for r in ov] == ["shared_stretch"], f"cleaning should find exactly one shared stretch: {ov}")

    sylph = stand_in(tmp / "sylph", SYLPH)
    code, err = quiet(MI.main, ["--workdir", str(work), "--refdir", str(refdir), "--threads", "2",
                                "--minimap2", MM2, "--sylph", sylph])
    check(code == 0, f"make_index failed: {err}")
    if code == 0:
        lin = {r["accession"]: r for r in csv.DictReader(open(refdir / "lineages.tsv"), delimiter="\t")}
        want = {H1: "Rhizophagus testii", H2: "Rhizophagus sp. TEST2", E1A: "CaMg/MRE", B1: "CaGg", M1: "Mycoavidus",
                M2: "Mycoavidus", C1: "Ralstonia", "NC_900101.1": "plant organelles", "NC_900102.1": "plant organelles"}
        check({a: r["lineage"] for a, r in lin.items()} == want, f"lineages.tsv: {[(a, r['lineage']) for a, r in lin.items()]}")
        check(lin[E1A]["in_t2"] == "no" and lin[B1]["in_t2"] == "yes", "lineages.tsv: in_t2 column")
        for s, n_genomes in (("full", 9), ("t2", 8)):
            info = json.loads((work / "index" / s / "index.json").read_text())
            check(info["genomes"] == n_genomes, f"{s}: {info['genomes']} genomes, want {n_genomes}")
            check(info["refs_md5"] == hashlib.md5((work / "index" / s / "refs.fna").read_bytes()).hexdigest(),
                  f"{s}: refs.fna md5 in index.json")
            sketched = (work / "index" / s / "bacteria.syldb").read_text().split()
            want_b = {B1, M1, M2, C1} | ({E1A} if s == "full" else set())
            check({Path(p).stem for p in sketched} == want_b and all("/clean/bacteria/" in p for p in sketched),
                  f"{s}: sylph was asked to sketch {sketched}")
        index_md = (refdir / "INDEX.md").read_text()
        for w in ("minimap2 2.31", "sylph 1.0.0", "section 10, entry 4", "| t2 | 8 |"):
            check(w in index_md, f"INDEX.md does not mention {w!r}")

        reads = short_reads()
        r_all = tmp / "reads.fq.gz"
        write_fastq(r_all, reads)

        def assign(index_set, platform, files, out, minimap2=MM2):
            return quiet(A.main, ["--index", str(work / "index" / index_set), "--platform", platform,
                                  "--reads", *map(str, files), "--out", str(out), "--refdir", str(refdir),
                                  "--threads", "3", "--minimap2", minimap2])

        out = tmp / "out_full"
        code, err = assign("full", "short", [r_all], out)
        check(code == 0, f"assign (full) failed: {err}")
        if code == 0:
            got = counts(out)
            for k, v in FULL.items():
                check(got.get(k) == v, f"full: {k[0]} / {k[1]}: got {got.get(k)}, want {v}")
            check(set(routed_names(out)) == ROUTED and len(routed_names(out)) == 40,
                  f"full: routed reads differ: {sorted(set(routed_names(out)) ^ ROUTED)[:10]}")
            meta = json.loads((out / "assign.json").read_text())
            check(meta["index"]["set"] == "full" and meta["inputs"][0]["reads"] == 260, "full: assign.json")
            check("--frag=no" in meta["command"] and "-N 50" in meta["command"], f"full: command {meta['command']}")

        out = tmp / "out_t2"
        code, err = assign("t2", "short", [r_all], out)
        check(code == 0, f"assign (t2) failed: {err}")
        if code == 0:
            got = counts(out)
            check(("lineage", "CaMg/MRE") not in got, "t2: CaMg/MRE should not be a lineage")
            for k, v in T2.items():
                check(got.get(k) == v, f"t2: {k[0]} / {k[1]}: got {got.get(k)}, want {v}")

        # Paired files whose mates share names, and one interleaved file with the same names.
        r1, r2, inter = tmp / "r_1.fq.gz", tmp / "r_2.fq.gz", tmp / "inter.fq.gz"
        half = len(reads) // 2
        mates1 = [(f"pair{k}", s) for k, (_, s) in enumerate(reads[:half])]
        mates2 = [(f"pair{k}", s) for k, (_, s) in enumerate(reads[half:])]
        write_fastq(r1, mates1)
        write_fastq(r2, mates2)
        write_fastq(inter, [x for pair in zip(mates1, mates2) for x in pair])
        for label, files in (("paired", [r1, r2]), ("interleaved", [inter])):
            out = tmp / f"out_{label}"
            code, err = assign("full", "short", files, out)
            check(code == 0, f"assign ({label}) failed: {err}")
            if code == 0:
                got = counts(out)
                check(all(got.get(k) == v for k, v in FULL.items()),
                      f"{label}: counts differ from the single-file run: "
                      f"{[(k, got.get(k), v) for k, v in FULL.items() if got.get(k) != v]}")
                check(len(routed_names(out)) == 40, f"{label}: {len(routed_names(out))} routed reads, want 40")

        # HiFi reads.
        hreads, hlabels = [], []
        take(H1c1, 0, 20000, 5, 3000, "hH1u", hreads, hlabels)
        take(E1, 9000, 25000, 5, 3000, "hE1u", hreads, hlabels)
        take(B1c1, 12000, 29000, 5, 3000, "hB1u", hreads, hlabels)
        take(H1c1, 20000, 30000, 3, 3000, "hHH", hreads, hlabels)
        for k in range(3):
            hreads.append((f"hrand_{k}", rand(3000)))
        hf = tmp / "hifi.fq.gz"
        write_fastq(hf, hreads)
        out = tmp / "out_hifi"
        code, err = assign("full", "hifi", [hf], out)
        check(code == 0, f"assign (hifi) failed: {err}")
        if code == 0:
            got = counts(out)
            want_h = {("total", "reads in input"): 21, ("lineage", "Rhizophagus testii"): 5,
                      ("lineage", "CaMg/MRE"): 5, ("lineage", "CaGg"): 5, ("none", "tied across lineages"): 3,
                      ("fungal", "fungal reads for E2"): 8, ("routed", "no alignment"): 3}
            for k, v in want_h.items():
                check(got.get(k) == v, f"hifi: {k[0]} / {k[1]}: got {got.get(k)}, want {v}")
            check(sorted(routed_names(out)) == [f"hrand_{k}" for k in range(3)], "hifi: routed reads")

        # A minimap2 that drops the output line of every second read: those reads count as unaligned.
        real_version = json.loads((work / "index" / "full" / "index.json").read_text())["minimap2"]
        dropper = stand_in(tmp / "mm2_dropper", f'''import subprocess, sys
if sys.argv[1:] == ["--version"]:
    print("{real_version}"); sys.exit(0)
p = subprocess.run([{MM2!r}, *sys.argv[1:]], stdin=sys.stdin, capture_output=True, check=True)
for line in p.stdout.decode().splitlines():
    if int(line.split("\\t")[0]) % 2 == 0:
        print(line)
''')
        out = tmp / "out_dropper"
        code, err = assign("full", "short", [r_all], out, minimap2=dropper)
        got = counts(out) if code == 0 else {}
        check(code == 0 and got.get(("total", "reads in input")) == 260 and got.get(("routed", "no alignment")) >= 130,
              f"missing output lines: exit {code}, {err}, routed {got.get(('routed', 'no alignment'))}")

        # A different minimap2 version is refused; a cut-short FASTQ writes nothing.
        other = stand_in(tmp / "mm2_other", 'import sys\nprint("2.99-r1"); sys.exit(0)\n')
        code, err = assign("full", "short", [r_all], tmp / "out_other", minimap2=other)
        check(code == 1 and "2.99" in err and not (tmp / "out_other" / "counts.tsv").exists(),
              f"version mismatch: exit {code}, {err!r}")
        bad = tmp / "bad.fq"
        bad.write_text("@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\n")
        code, err = assign("full", "short", [bad], tmp / "out_bad")
        check(code == 1 and "cut short" in err and not (tmp / "out_bad" / "counts.tsv").exists(),
              f"cut-short FASTQ: exit {code}, {err!r}")

# An endosymbiont that fits none of the lineages in 4.1 stops the index build before anything is written.
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    refdir = tmp / "references"
    write_lists(refdir)
    (refdir / "endosymbionts.tsv").write_text(f"accession{T}organism{T}note\n{B1}{T}Burkholderia testii{T}\n")
    code, err = quiet(MI.main, ["--workdir", str(tmp / "work"), "--refdir", str(refdir), "--minimap2", MM2,
                                "--sylph", stand_in(tmp / "sylph", SYLPH)])
    check(code == 1 and "fits none of the endosymbiont lineages" in err and not (refdir / "lineages.tsv").exists(),
          f"unknown endosymbiont: exit {code}, {err!r}")

if fails:
    print("FAIL\n" + "\n".join(fails))
    sys.exit(1)
print(f"PASS: {len(unit)} rule checks, PAF parsing, cleaning -> index (lineages, full and T2, sylph stand-in), "
      "short reads (every planted case), T2 index, paired and interleaved mates, HiFi, routed reads, "
      "missing output lines, version check, cut-short FASTQ, unknown endosymbiont stops the index build")
