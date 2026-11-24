"""Check clean_references.py on synthetic genomes, with real minimap2 and a fake NCBI.

Needs minimap2 (activate phase0-reads, or set MINIMAP2=/path/to/minimap2).
Planted cases, and what the rule in PROTOCOL.md 3.4 must do with them:
  - a bacterial contig copied into a host assembly         -> host contig removed
  - a second copy of it in the same host at 3% divergence  -> removed too (needs -p 0: minimap2
                                                              hides this secondary alignment by default)
  - a contaminant contig in a host at 3% divergence        -> host contig removed (97% identity passes)
  - a fungal contig inside a bacterial genome (a "MAG")    -> bacterial contig removed
  - an 800 bp stretch shared by a host and a bacterial contig -> kept in both, listed
  - a 700 bp piece of that bacterial contig in a large host contig -> kept in both, listed (needs -p 0)
  - a 3 kb piece copied into a host at 12% divergence      -> ignored (below 95% identity)
A control run with minimap2's default -p must miss the two "needs -p 0" cases, which shows
that the test can tell the settings apart.
Organelle DNA (PROTOCOL.md section 10, entry 10; minimap2 -x asm20, 500 bp at 80% identity):
  - a bacterial contig that is a fungal mitochondrion at 10% divergence -> removed as organelle DNA
  - a contaminant contig that is a plant chloroplast at 15% divergence  -> removed as organelle DNA
  - a 1.5 kb mitochondrial piece at 4% divergence in a 20 kb contig     -> contig kept, stretch listed
  - a 1.5 kb mitochondrial piece at 30% divergence                       -> ignored (below 80% identity)
  - the fungal contig removed as fungal, half of it also mitochondrial   -> listed once, as fungal
  - a fungal mitochondrion whose length is not the one listed            -> exit 1, nothing written
"""
import contextlib
import csv
import hashlib
import io
import json
import os
import random
import shutil
import sys
import tempfile
import urllib.parse
import zipfile
from pathlib import Path

here = Path(__file__).parent
sys.path.insert(0, str(here.parent / "scripts"))
import clean_references as M  # noqa: E402

MM2 = os.environ.get("MINIMAP2") or shutil.which("minimap2")
if not MM2:
    print("SKIP: minimap2 not found (activate phase0-reads, or set MINIMAP2)")
    sys.exit(0)

fails = []

# ---- unit checks of the rule's arithmetic
if M.union_length([(0, 10), (5, 20), (30, 40), (40, 45)]) != 35:
    fails.append("union_length: overlapping and touching intervals")
for block, match, want in [(500, 475, True), (500, 474, False), (499, 499, False), (10000, 9500, True)]:
    if M.passes({"block": block, "match": match}) != want:
        fails.append(f"passes({block} bp, {match} matches) should be {want}")

# ---- synthetic genomes
rng = random.Random(20261120)


def rand(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq, rate):
    s = list(seq)
    for i in rng.sample(range(len(s)), int(len(s) * rate)):
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


H1_c1 = rand(60000)
B1_c1, B1_c2 = rand(12000), rand(8000)
B2_c1, B2_c2 = rand(20000), rand(15000)
C1_c1, C1_c2 = rand(15000), rand(10000)
GENOMES = {
    # hosts
    "GCA_900000001.1": [("H1c1.1", H1_c1),
                        ("H1c2.1", B1_c1),                                             # bacterial contig in host
                        ("H1c3.1", rand(10000) + B2_c1[5000:5800] + rand(29200)),      # 800 bp shared
                        ("H1c4.1", rand(10000) + mutate(B2_c2[2000:5000], 0.12) + rand(17000)),   # 12% divergent
                        ("H1c5.1", mutate(B1_c1, 0.03)),                               # 2nd bacterial copy, 3% divergent
                        ("H1c6.1", rand(20000) + B1_c1[3000:3700] + rand(20000))],     # 700 bp shared
    "GCA_900000002.1": [("H2c1.1", rand(30000)),
                        ("H2c2.1", mutate(C1_c1, 0.03))],                              # contaminant at 97%
    # endosymbionts
    "GCA_800000001.1": [("B1c1.1", B1_c1), ("B1c2.1", B1_c2)],
    "GCA_800000002.1": [("B2c1.1", B2_c1), ("B2c2.1", B2_c2),
                        ("B2c3.1", H1_c1[20000:29000])],                               # fungal contig in a MAG
    "GCA_800000003.1": [("B3c1.1", rand(6000))],
    # contaminant
    "GCF_700000001.1": [("C1c1.1", C1_c1), ("C1c2.1", C1_c2)],
}
ORGANELLES = {"NC_900001.1": rand(5000), "NC_900002.1": rand(4000)}
MITO = rand(4000) + H1_c1[21000:26000] + rand(3000)          # a fungal mitochondrion, with a host piece in it
ORGANELLES["MN900001.1"] = MITO
GENOMES["GCA_800000001.1"].append(("B1c3.1", mutate(MITO[1000:7000], 0.10)))               # 90%: organelle DNA
GENOMES["GCF_700000001.1"].append(("C1c3.1", mutate(ORGANELLES["NC_900001.1"][500:3500], 0.15)))  # 85%: chloroplast
GENOMES["GCA_800000002.1"].append(("B2c4.1", rand(9000) + mutate(MITO[2000:3500], 0.04) + rand(9500)))  # stretch
GENOMES["GCA_800000003.1"].append(("B3c2.1", rand(4000) + mutate(MITO[9500:11000], 0.30) + rand(4000)))  # 70%
H1, H2, B1, B2, B3, C1 = ("GCA_900000001.1", "GCA_900000002.1", "GCA_800000001.1", "GCA_800000002.1",
                          "GCA_800000003.1", "GCF_700000001.1")


def wrap(seq, width=80):
    return "".join(seq[i:i + width] + "\n" for i in range(0, len(seq), width))


BAD_MD5 = set()
calls = {"download": 0}


def fake_http_get(url, params=None, dest=None, tries=4):
    params = params or {}
    if url.endswith("/version"):
        return json.dumps({"version": "18.41.1"}).encode()
    if "/genome/accession/" in url and url.endswith("/download"):
        assert params.get("include_annotation_type") == "GENOME_FASTA", params
        acc = urllib.parse.unquote(url.split("/genome/accession/")[1].split("/")[0])
        fasta = "".join(f">{name} synthetic sequence of {acc}\n" + wrap(seq) for name, seq in GENOMES[acc])
        member = f"ncbi_dataset/data/{acc}/{acc}_TEST_genomic.fna"
        md5 = "0" * 32 if acc in BAD_MD5 else hashlib.md5(fasta.encode()).hexdigest()
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as z:
            z.writestr("README.md", "synthetic")
            z.writestr("ncbi_dataset/data/assembly_data_report.jsonl", "{}\n")
            z.writestr(member, fasta)
            z.writestr("md5sum.txt", f"{md5}  {member}\n")
        dest.write_bytes(buf.getvalue())
        calls["download"] += 1
        return dest
    if url.endswith("efetch.fcgi"):
        return "".join(f">{i} synthetic organelle\n" + wrap(ORGANELLES[i]) for i in params["id"].split(",")).encode()
    raise AssertionError(f"unexpected request: {url}")


M.http_get = fake_http_get


def write_lists(refdir):
    refdir.mkdir(parents=True)
    t = "\t"
    (refdir / "hosts.tsv").write_text(f"accession{t}organism{t}note\n{H1}{t}Testomyces alpha{t}host: Testomyces alpha\n"
                                      f"{H2}{t}Testomyces beta{t}host: Testomyces beta\n")
    (refdir / "endosymbionts.tsv").write_text(
        f"accession{t}organism{t}note\n{B1}{t}Candidatus Glomeribacter testii{t}\n"
        f"{B2}{t}Mycoavidus testii{t}\n{B3}{t}Candidatus Moeniiplasma testii{t}held_out_in_T2\n")
    (refdir / "contaminants.tsv").write_text(f"accession{t}organism{t}note\n{C1}{t}Ralstonia testii{t}type species\n")
    (refdir / "organelles.tsv").write_text(
        f"plant{t}organelle{t}accession{t}title{t}length_bp\n"
        f"Plantago testii{t}chloroplast{t}NC_900001.1{t}x{t}5000\nPlantago testii{t}mitochondrion{t}NC_900002.1{t}x{t}4000\n")
    (refdir / "fungal_mitochondria.tsv").write_text(
        f"fungus{t}organelle{t}accession{t}title{t}length_bp{t}source{t}query\n"
        f"Testomyces alpha{t}mitochondrion{t}MN900001.1{t}x{t}{len(MITO)}{t}GenBank{t}q\n")


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def names(path):
    return [line[1:].strip() for line in open(path) if line.startswith(">")]


def run(work, refdir):
    err = io.StringIO()
    for p in (refdir / "overlaps.tsv", refdir / "genome_files.tsv", refdir / "CLEANING.md"):
        p.unlink(missing_ok=True)
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
        code = M.main(["--workdir", str(work), "--refdir", str(refdir), "--threads", "2", "--minimap2", MM2])
    return code, err.getvalue()


with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    refdir, work = tmp / "references", tmp / "work"
    write_lists(refdir)
    code, err = run(work, refdir)
    if code != 0:
        fails.append(f"exit code {code}: {err}")
    else:
        rows = read(refdir / "overlaps.tsv")
        kinds = {}
        for r in rows:
            kinds.setdefault(r["kind"], []).append(r)
        got = sorted(r["host_contig"] for r in kinds.get("removed_host_contig", []))
        if got != [f"{H1}|H1c2.1", f"{H1}|H1c5.1", f"{H2}|H2c2.1"]:
            fails.append(f"removed host contigs: got {got}")
        matched = {r["host_contig"]: r["matched"] for r in kinds.get("removed_host_contig", [])}
        if (matched.get(f"{H1}|H1c2.1"), matched.get(f"{H1}|H1c5.1"), matched.get(f"{H2}|H2c2.1")) != (B1, B1, C1):
            fails.append(f"removed host contigs matched the wrong genomes: {matched}")
        got = [(r["bacterial_contig"], r["matched"]) for r in kinds.get("removed_bacterial_contig", [])]
        if got != [(f"{B2}|B2c3.1", H1)]:
            fails.append(f"removed bacterial contigs: got {got}")
        sh = {(r["host_contig"], r["bacterial_contig"]): r for r in kinds.get("shared_stretch", [])}
        want = {(f"{H1}|H1c3.1", f"{B2}|B2c1.1"): (800, 10000, 5000), (f"{H1}|H1c6.1", f"{B1}|B1c1.1"): (700, 20000, 3000)}
        if sorted(sh) != sorted(want) or len(kinds.get("shared_stretch", [])) != 2:
            fails.append(f"shared stretches: got {[(r['host_contig'], r['bacterial_contig']) for r in kinds.get('shared_stretch', [])]}")
        else:
            for k, (bp, hs, bs) in want.items():
                s = sh[k]
                if not (bp - 10 <= int(s["aligned_bp"]) <= bp + 30 and float(s["identity"]) >= 0.99
                        and abs(int(s["host_start"]) - hs) <= 15 and abs(int(s["bacterial_start"]) - bs) <= 15):
                    fails.append(f"shared stretch coordinates off: {s}")
        if any("B2c2.1" in r["bacterial_contig"] for r in rows):
            fails.append("the 12% divergent copy should not pass the 95% identity step")
        if any(B1 in r["bacterial_contig"] for r in kinds.get("removed_bacterial_contig", [])):
            fails.append("B1c1 was removed, but its only host match is a removed host contig")

        if names(work / "clean/hosts" / f"{H1}.fna") != [f"{H1}|H1c1.1", f"{H1}|H1c3.1", f"{H1}|H1c4.1", f"{H1}|H1c6.1"]:
            fails.append(f"cleaned host H1: {names(work / 'clean/hosts' / f'{H1}.fna')}")
        if names(work / "clean/bacteria" / f"{B2}.fna") != [f"{B2}|B2c1.1", f"{B2}|B2c2.1", f"{B2}|B2c4.1"]:
            fails.append(f"cleaned genome B2: {names(work / 'clean/bacteria' / f'{B2}.fna')}")
        if names(work / "clean/bacteria" / f"{B1}.fna") != [f"{B1}|B1c1.1", f"{B1}|B1c2.1"]:
            fails.append("cleaned genome B1 should be unchanged")
        if names(work / "clean/organelles.fna") != ["NC_900001.1|NC_900001.1", "NC_900002.1|NC_900002.1",
                                                    "MN900001.1|MN900001.1"]:
            fails.append(f"organelles: {names(work / 'clean/organelles.fna')}")

        # Organelle DNA (entry 10).
        got = sorted((r["bacterial_contig"], r["matched"]) for r in kinds.get("removed_as_organelle", []))
        if got != [(f"{B1}|B1c3.1", "MN900001.1"), (f"{C1}|C1c3.1", "NC_900001.1")]:
            fails.append(f"contigs removed as organelle DNA: got {got}")
        cov = {r["bacterial_contig"]: float(r["bacterial_covered"]) for r in kinds.get("removed_as_organelle", [])}
        if not all(v >= 0.9 for v in cov.values()):
            fails.append(f"contigs removed as organelle DNA should be nearly all covered: {cov}")
        st = kinds.get("organelle_stretch", [])
        if [(r["bacterial_contig"], r["host_assembly"], r["host_organism"]) for r in st] != \
                [(f"{B2}|B2c4.1", "MN900001.1", "Testomyces alpha mitochondrion")]:
            fails.append(f"organelle stretches: got {[(r['bacterial_contig'], r['host_contig']) for r in st]}")
        elif not (abs(int(st[0]["bacterial_start"]) - 9000) <= 30 and 1400 <= int(st[0]["aligned_bp"]) <= 1560
                  and abs(int(st[0]["host_start"]) - 2000) <= 30 and float(st[0]["identity"]) >= 0.94):
            fails.append(f"organelle stretch coordinates off: {st[0]}")
        if any("B3c2.1" in r["bacterial_contig"] for r in rows):
            fails.append("the 30% divergent mitochondrial piece should not pass the 80% identity step")
        if [r["kind"] for r in rows if r["bacterial_contig"] == f"{B2}|B2c3.1"] != ["removed_bacterial_contig"]:
            fails.append("B2c3.1 should be listed once, as removed for being fungal")
        if names(work / "clean/bacteria" / f"{C1}.fna") != [f"{C1}|C1c1.1", f"{C1}|C1c2.1"]:
            fails.append(f"cleaned genome C1: {names(work / 'clean/bacteria' / f'{C1}.fna')}")

        files = {r["accession"]: r for r in read(refdir / "genome_files.tsv")}
        if len(files) != 9:
            fails.append(f"genome_files.tsv: expected 9 rows, got {len(files)}")
        if (files[H1]["removed_sequences"], files[H1]["removed_bp"]) != ("2", "24000"):
            fails.append(f"genome_files.tsv H1 removed: {files[H1]['removed_sequences']}, {files[H1]['removed_bp']}")
        if (files[B2]["removed_sequences"], files[B2]["removed_bp"]) != ("1", "9000"):
            fails.append(f"genome_files.tsv B2 removed: {files[B2]['removed_sequences']}, {files[B2]['removed_bp']}")
        if (files[B1]["removed_sequences"], files[B1]["removed_bp"]) != ("1", "6000"):
            fails.append(f"genome_files.tsv B1 removed: {files[B1]['removed_sequences']}, {files[B1]['removed_bp']}")
        if (files["MN900001.1"]["role"], files["MN900001.1"]["organism"]) != ("organelle",
                                                                               "Testomyces alpha mitochondrion"):
            fails.append(f"genome_files.tsv fungal mitochondrion: {files['MN900001.1']}")
        if any(files[a]["ncbi_md5_checked"] != "yes" for a in GENOMES):
            fails.append("every synthetic download carries an md5 and should be checked")
        cleaning = (refdir / "CLEANING.md").read_text()
        for want in ("minimap2 2.31", "-x asm5 -c -p 0", "NCBI Datasets API 18.41.1", "section 10, entry 3",
                     "section 10, entry 10", "minimap2 -x asm20 -c -p 0", "paf/organelles.paf",
                     "| GCA_800000001.1 | Candidatus Glomeribacter testii | 3 | 1 | 6000 |",
                     "Organelle stretches left in the bacterial genomes (reads from them count for no lineage): 1,"):
            if want not in cleaning:
                fails.append(f"CLEANING.md does not mention {want!r}")

        # A rerun downloads nothing again and gives the same overlaps.
        first = (refdir / "overlaps.tsv").read_bytes()
        n_dl = calls["download"]
        code, err = run(work, refdir)
        if code != 0 or calls["download"] != n_dl:
            fails.append(f"rerun: exit {code}, {calls['download'] - n_dl} new downloads (expected 0)")
        if (refdir / "overlaps.tsv").read_bytes() != first:
            fails.append("rerun gave different overlaps.tsv")

        # Control: with minimap2's default -p, the two "needs -p 0" cases are missed.
        saved = M.MINIMAP2_ARGS
        M.MINIMAP2_ARGS = ["-x", "asm5", "-c"]
        code, err = run(work, refdir)
        M.MINIMAP2_ARGS = saved
        rows = read(refdir / "overlaps.tsv") if code == 0 else []
        hidden_copy = not any(r["host_contig"] == f"{H1}|H1c5.1" for r in rows)
        hidden_stretch = not any(r["host_contig"] == f"{H1}|H1c6.1" for r in rows)
        if code != 0 or not (hidden_copy and hidden_stretch):
            fails.append(f"control with default -p: exit {code}; divergent copy hidden: {hidden_copy}, "
                         f"700 bp stretch hidden: {hidden_stretch} (both should be hidden, or the test cannot "
                         "tell the settings apart)")

        # Control: with -x asm5 for the organelle step, the mitochondrion at 10% and the chloroplast at 15%
        # divergence are missed, which shows that the test can tell the presets apart.
        saved = M.ORGANELLE_ARGS
        M.ORGANELLE_ARGS = ["-x", "asm5", "-c", "-p", "0"]
        code, err = run(work, refdir)
        M.ORGANELLE_ARGS = saved
        rows = read(refdir / "overlaps.tsv") if code == 0 else []
        if code != 0 or any(r["kind"] == "removed_as_organelle" for r in rows):
            fails.append(f"control with -x asm5: exit {code}; removed as organelle DNA: "
                         f"{[r['bacterial_contig'] for r in rows if r['kind'] == 'removed_as_organelle']} (none should "
                         "be, or the test cannot tell the presets apart)")

# A download whose md5 does not match: exit 1, nothing written to refdir.
BAD_MD5.add(B2)
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    refdir, work = tmp / "references", tmp / "work"
    write_lists(refdir)
    code, err = run(work, refdir)
    if code != 1 or "md5" not in err:
        fails.append(f"bad md5: exit {code}, message {err!r}")
    for f in ("overlaps.tsv", "genome_files.tsv", "CLEANING.md"):
        if (refdir / f).exists():
            fails.append(f"bad md5: {f} was written anyway")
BAD_MD5.clear()

# A fungal mitochondrion whose length is not the listed one: exit 1, nothing written to refdir.
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    refdir, work = tmp / "references", tmp / "work"
    write_lists(refdir)
    lst = refdir / "fungal_mitochondria.tsv"
    lst.write_text(lst.read_text().replace(f"\t{len(MITO)}\t", f"\t{len(MITO) + 1}\t"))
    code, err = run(work, refdir)
    if code != 1 or "fungal_mitochondria.tsv says" not in err:
        fails.append(f"wrong mitochondrion length: exit {code}, message {err!r}")
    for f in ("overlaps.tsv", "genome_files.tsv", "CLEANING.md"):
        if (refdir / f).exists():
            fails.append(f"wrong mitochondrion length: {f} was written anyway")

if fails:
    print("FAIL\n" + "\n".join(fails))
    sys.exit(1)
print("PASS: bacterial contigs in hosts removed (exact, a second copy at 3% divergence, a contaminant at 3%), "
      "fungal contig in a bacterial genome removed, 800 bp and 700 bp shared stretches kept and listed, 12% "
      "divergent copy ignored, organelle DNA (mitochondrion at 90% and chloroplast at 85% removed, 1.5 kb stretch "
      "listed, 70% piece ignored, a contig removed as fungal listed once), cleaned files, md5 checks, rerun "
      "without new downloads, default -p misses the second copy and the 700 bp stretch (control), -x asm5 misses "
      "both organelle contigs (control), bad md5 and a "
      "wrong mitochondrion length write nothing")
