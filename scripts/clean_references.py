#!/usr/bin/env python3
"""Clean the Phase 0 reference genomes, once, before stage 0a (PROTOCOL.md 3.4).

Downloads every genome in the reference lists, aligns all bacterial genomes
(the endosymbionts, held-out ones included, and the kit contaminants) to each
host assembly with minimap2, and applies the protocol's rule:

  1. keep alignments of at least 500 bp at 95% identity or more;
  2. a host contig at least half covered by them is bacterial: it is removed
     from the host assembly;
  3. a bacterial contig at least half covered by alignments to the host
     contigs that remain is fungal: it is removed from its genome;
  4. anything else that is shared stays in both, and is listed.

The settings the protocol leaves open are fixed in PROTOCOL.md section 10,
entry 3: minimap2 runs with -x asm5 -c -p 0 (-c for exact identity, -p 0 so
that no secondary alignment is hidden); identity is matching bases divided by
alignment block length (PAF columns 10 and 11) and the 500 bp is the block
length; primary and secondary alignments both count; a contig is one FASTA
record, Ns included in its length; its coverage is the union of its aligned
intervals divided by its length; step 3 uses the remaining contigs of all host
assemblies together.

Writes under --workdir (large files, not committed):
  genomes/<accession>.fna          downloaded genomes, sequences renamed <accession>|<sequence id>
  genomes/organelles/<acc>.fna     downloaded organelle genomes, renamed the same way
  bacteria.fna                     all bacterial genomes together (the minimap2 query)
  paf/<host accession>.paf         minimap2 output, one file per host assembly
  clean/hosts/<accession>.fna      cleaned host assemblies
  clean/bacteria/<accession>.fna   cleaned bacterial genomes
  clean/organelles.fna             organelle genomes (not part of the cleaning; copied as they are)
and in --refdir (default references/; small files, committed before stage 0a):
  overlaps.tsv       every removed contig and every remaining shared stretch
  genome_files.tsv   sequences, length and md5 of every downloaded and cleaned file
  CLEANING.md        date, versions, commands and counts

Nothing is written to --refdir unless every step succeeds. Downloads that
finished are kept, so a rerun only fetches what is missing.

Needs minimap2 on PATH (activate the phase0-reads environment) and network
access to api.ncbi.nlm.nih.gov and eutils.ncbi.nlm.nih.gov. Setting
NCBI_API_KEY is optional.

Usage:
  conda activate phase0-reads
  python scripts/clean_references.py --workdir ~/phase0-data --threads 16
"""

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path

__version__ = "0.2"

DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY = os.environ.get("NCBI_API_KEY", "")

MINIMAP2_ARGS = ["-x", "asm5", "-c", "-p", "0"]   # PROTOCOL.md 3.4 step 1; -c and -p 0: section 10, entry 3
MIN_BLOCK = 500                         # bp, PAF column 11
MIN_IDENTITY_PCT = 95                   # PAF column 10 / column 11
MIN_COVERED = 0.5                       # half of the contig

OVERLAP_COLS = ["kind", "host_assembly", "host_organism", "host_contig", "host_length", "host_start", "host_end",
                "host_covered", "bacterial_genome", "bacterial_organism", "bacterial_contig", "bacterial_length",
                "bacterial_start", "bacterial_end", "bacterial_covered", "identity", "aligned_bp", "matched"]
FILE_COLS = ["accession", "role", "organism", "downloaded_sequences", "downloaded_bp", "downloaded_md5",
             "ncbi_md5_checked", "removed_sequences", "removed_bp", "clean_sequences", "clean_bp", "clean_md5"]


# ---------------------------------------------------------------- network ---

def http_get(url, params=None, dest=None, tries=4):
    """GET a URL. Returns the body as bytes, or streams it to dest (a Path). Tests replace this."""
    params = dict(params or {})
    headers = {"User-Agent": f"amf-endosym-phase0-clean/{__version__}"}
    if API_KEY:
        if url.startswith(EUTILS):
            params["api_key"] = API_KEY
        else:
            headers["api-key"] = API_KEY
    full = url + ("?" + urllib.parse.urlencode(params) if params else "")
    last = None
    for attempt in range(tries):
        if attempt:
            time.sleep(2 ** attempt)
        try:
            req = urllib.request.Request(full, headers=headers)
            with urllib.request.urlopen(req, timeout=300) as resp:
                if dest is None:
                    return resp.read()
                with open(dest, "wb") as out:
                    shutil.copyfileobj(resp, out, 1 << 20)
                return dest
        except urllib.error.HTTPError as e:
            last = f"HTTP {e.code} for {full}"
            if e.code not in (429, 500, 502, 503, 504):
                break
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last = f"{e} for {full}"
    raise RuntimeError(last)


def md5_listed(z, member):
    """The md5 that the NCBI download's md5sum.txt gives for member, or '' if it gives none."""
    if "md5sum.txt" not in z.namelist():
        return ""
    for line in z.read("md5sum.txt").decode("utf-8", "replace").splitlines():
        parts = line.split()
        if len(parts) == 2:
            path = parts[1].lstrip("*")
            if path == member or member.endswith("/" + path) or path.endswith(member):
                return parts[0].lower()
    return ""


def download_genome(acc, dest):
    """Fetch one assembly's genomic FASTA from NCBI Datasets; rename sequences <acc>|<id>.

    Returns True when the md5 of the FASTA inside the download was checked
    against the download's own md5sum.txt.
    """
    if dest.exists() and dest.stat().st_size > 0:
        return dest.with_name(dest.name + ".md5ok").exists()
    part = dest.with_name(dest.name + ".zip.part")
    http_get(f"{DATASETS}/genome/accession/{urllib.parse.quote(acc, safe='')}/download",
             {"include_annotation_type": "GENOME_FASTA", "hydrated": "FULLY_HYDRATED"}, dest=part)
    tmp = dest.with_name(dest.name + ".part")
    with zipfile.ZipFile(part) as z:
        members = [m for m in z.namelist()
                   if m.startswith(f"ncbi_dataset/data/{acc}/") and m.endswith((".fna", ".fa", ".fasta"))]
        if len(members) != 1:
            raise RuntimeError(f"{acc}: expected one genomic FASTA in the NCBI download, found {len(members)}")
        expected = md5_listed(z, members[0])
        h = hashlib.md5()
        with z.open(members[0]) as raw, open(tmp, "wb") as out:
            for line in raw:
                h.update(line)
                if line.startswith(b">"):
                    seqid = line[1:].split(None, 1)[0].decode("ascii")
                    out.write(f">{acc}|{seqid}\n".encode("ascii"))
                else:
                    out.write(line.rstrip(b"\r\n") + b"\n")
    if expected and h.hexdigest() != expected:
        tmp.unlink()
        raise RuntimeError(f"{acc}: the FASTA in the NCBI download does not match its md5sum.txt")
    os.replace(tmp, dest)
    part.unlink()
    if expected:
        dest.with_name(dest.name + ".md5ok").write_text(expected + "\n")
    return bool(expected)


def parse_fasta_text(text):
    """{first word of header: sequence} from FASTA text."""
    seqs, name, chunks = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks)
            name, chunks = line[1:].split(None, 1)[0], []
        elif name is not None:
            chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def download_organelles(rows, outdir):
    """Fetch the organelle genomes with E-utilities efetch; check each length against organelles.tsv."""
    outdir.mkdir(parents=True, exist_ok=True)
    todo = [r for r in rows if not (outdir / f"{r['accession']}.fna").exists()]
    if not todo:
        return
    text = http_get(f"{EUTILS}/efetch.fcgi", {"db": "nuccore", "id": ",".join(r["accession"] for r in todo),
                                               "rettype": "fasta", "retmode": "text"}).decode("ascii", "replace")
    seqs = parse_fasta_text(text)
    for r in todo:
        acc = r["accession"]
        seq = seqs.get(acc)
        if seq is None:
            raise RuntimeError(f"{acc}: efetch did not return this organelle genome")
        if r.get("length_bp") and len(seq) != int(r["length_bp"]):
            raise RuntimeError(f"{acc}: efetch returned {len(seq)} bp, organelles.tsv says {r['length_bp']}")
        tmp = outdir / f"{acc}.fna.part"
        with open(tmp, "w", encoding="ascii", newline="\n") as out:
            out.write(f">{acc}|{acc}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")
        os.replace(tmp, outdir / f"{acc}.fna")


# ------------------------------------------------------------------ FASTA ---

def fasta_stats(path):
    """(number of sequences, bp, md5 of the file)."""
    n = bp = 0
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for line in fh:
            h.update(line)
            if line.startswith(b">"):
                n += 1
            else:
                bp += len(line.strip())
    return n, bp, h.hexdigest()


def fasta_lengths(path):
    """{sequence name: length}."""
    lengths, name = {}, None
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                name = line[1:].split(None, 1)[0].decode("ascii")
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def write_without(src, dest, drop):
    """Copy a FASTA, leaving out the sequences named in drop."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    keep = True
    with open(src, "rb") as fh, open(tmp, "wb") as out:
        for line in fh:
            if line.startswith(b">"):
                keep = line[1:].split(None, 1)[0].decode("ascii") not in drop
            if keep:
                out.write(line)
    os.replace(tmp, dest)


# -------------------------------------------------------------- the rule ---

def read_paf(path):
    """Alignments from a PAF file; an alignment listed twice is kept once."""
    hits, seen = [], set()
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            key = (f[0], f[2], f[3], f[4], f[5], f[7], f[8])
            if key in seen:
                continue
            seen.add(key)
            hits.append({"q": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4],
                         "t": f[5], "tlen": int(f[6]), "ts": int(f[7]), "te": int(f[8]),
                         "match": int(f[9]), "block": int(f[10])})
    return hits


def passes(hit):
    """Step 1: at least 500 bp (alignment block) at 95% identity or more (integers, no rounding)."""
    return hit["block"] >= MIN_BLOCK and hit["match"] * 100 >= MIN_IDENTITY_PCT * hit["block"]


def union_length(intervals):
    """Total length covered by half-open intervals [start, end)."""
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                total += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s
    return total


def apply_rule(hits):
    """Steps 1-4. Returns kept hits, host coverage, removed host contigs, bacterial coverage,
    removed bacterial contigs and the shared stretches that stay in both."""
    kept = [h for h in hits if passes(h)]
    iv, length = defaultdict(list), {}
    for h in kept:
        iv[h["t"]].append((h["ts"], h["te"]))
        length[h["t"]] = h["tlen"]
    host_cov = {t: union_length(v) / length[t] for t, v in iv.items()}
    removed_hosts = {t for t, c in host_cov.items() if c >= MIN_COVERED}

    iv, length = defaultdict(list), {}
    for h in kept:
        if h["t"] in removed_hosts:
            continue
        iv[h["q"]].append((h["qs"], h["qe"]))
        length[h["q"]] = h["qlen"]
    bact_cov = {q: union_length(v) / length[q] for q, v in iv.items()}
    removed_bact = {q for q, c in bact_cov.items() if c >= MIN_COVERED}

    shared = [h for h in kept if h["t"] not in removed_hosts and h["q"] not in removed_bact]
    return kept, host_cov, removed_hosts, bact_cov, removed_bact, shared


def genome_of(name):
    return name.split("|", 1)[0]


# ------------------------------------------------------------------- main ---

def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, cols, rows):
    tmp = path.with_name(path.name + ".part")
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return tmp


def tool_version(minimap2):
    out = subprocess.run([minimap2, "--version"], capture_output=True, text=True, check=True)
    return out.stdout.strip()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True, help="where genomes and cleaned files go (a few GB)")
    ap.add_argument("--refdir", default="references", help="the reference lists; the outputs to commit go here too")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--minimap2", default="minimap2", help="minimap2 to run [minimap2 on PATH]")
    args = ap.parse_args(argv)
    work, refdir = Path(args.workdir).expanduser(), Path(args.refdir)

    try:
        endo = read_tsv(refdir / "endosymbionts.tsv")
        cont = read_tsv(refdir / "contaminants.tsv")
        hosts = read_tsv(refdir / "hosts.tsv")
        orgs = read_tsv(refdir / "organelles.tsv")
    except FileNotFoundError as e:
        print(f"Missing reference list: {e.filename}. Run this from the repository folder.", file=sys.stderr)
        return 1
    bacteria = [(r, "endosymbiont") for r in endo] + [(r, "contaminant") for r in cont]
    organism = {r["accession"]: r["organism"] for r, _ in bacteria}
    organism.update({r["accession"]: r["organism"] for r in hosts})

    try:
        mm2_version = tool_version(args.minimap2)
    except (OSError, subprocess.CalledProcessError) as e:
        print(f"Could not run minimap2 ({e}). Activate the phase0-reads environment first.", file=sys.stderr)
        return 1

    gdir, odir = work / "genomes", work / "genomes" / "organelles"
    gdir.mkdir(parents=True, exist_ok=True)
    md5_checked = {}
    try:
        api_version = json.loads(http_get(f"{DATASETS}/version")).get("version", "")
        todo = [r["accession"] for r, _ in bacteria] + [r["accession"] for r in hosts]
        for i, acc in enumerate(todo, 1):
            print(f"[{i}/{len(todo)}] {acc}", flush=True)
            md5_checked[acc] = download_genome(acc, gdir / f"{acc}.fna")
        print(f"organelle genomes ({len(orgs)})", flush=True)
        download_organelles(orgs, odir)
    except (RuntimeError, OSError, zipfile.BadZipFile, json.JSONDecodeError) as e:
        print(f"Download failed, nothing was written to {refdir}/: {e}", file=sys.stderr)
        print("Downloads that finished are kept; run the same command again.", file=sys.stderr)
        return 1

    # All bacterial genomes in one query file, held-out ones included (3.4 step 1).
    bact_fna = work / "bacteria.fna"
    with open(bact_fna.with_name("bacteria.fna.part"), "wb") as out:
        for r, _ in bacteria:
            with open(gdir / f"{r['accession']}.fna", "rb") as fh:
                shutil.copyfileobj(fh, out)
    os.replace(bact_fna.with_name("bacteria.fna.part"), bact_fna)

    pdir = work / "paf"
    pdir.mkdir(exist_ok=True)
    hits, commands = [], []
    for h in hosts:
        acc = h["accession"]
        paf = pdir / f"{acc}.paf"
        cmd = [args.minimap2, *MINIMAP2_ARGS, "-t", str(args.threads), str(gdir / f"{acc}.fna"), str(bact_fna)]
        print(f"aligning bacterial genomes to {acc} ({h['organism']})", flush=True)
        try:
            with open(paf.with_name(paf.name + ".part"), "w") as out, open(work / "minimap2.log", "a") as err:
                subprocess.run(cmd, stdout=out, stderr=err, check=True)
        except (OSError, subprocess.CalledProcessError) as e:
            print(f"minimap2 failed for {acc} ({e}); see {work / 'minimap2.log'}. Nothing was written to {refdir}/.",
                  file=sys.stderr)
            return 1
        os.replace(paf.with_name(paf.name + ".part"), paf)
        commands.append(" ".join(["minimap2", *MINIMAP2_ARGS, "-t", str(args.threads),
                                  f"genomes/{acc}.fna", "bacteria.fna", ">", f"paf/{acc}.paf"]))
        hits.extend(read_paf(paf))

    kept, host_cov, removed_hosts, bact_cov, removed_bact, shared = apply_rule(hits)

    # overlaps.tsv
    rows = []
    for t in sorted(removed_hosts):
        rows.append({"kind": "removed_host_contig", "host_assembly": genome_of(t),
                     "host_organism": organism.get(genome_of(t), ""), "host_contig": t,
                     "host_length": next(h["tlen"] for h in kept if h["t"] == t),
                     "host_covered": f"{host_cov[t]:.3f}",
                     "matched": ",".join(sorted({genome_of(h["q"]) for h in kept if h["t"] == t}))})
    for q in sorted(removed_bact):
        rows.append({"kind": "removed_bacterial_contig", "bacterial_genome": genome_of(q),
                     "bacterial_organism": organism.get(genome_of(q), ""), "bacterial_contig": q,
                     "bacterial_length": next(h["qlen"] for h in kept if h["q"] == q),
                     "bacterial_covered": f"{bact_cov[q]:.3f}",
                     "matched": ",".join(sorted({genome_of(h["t"]) for h in kept
                                                 if h["q"] == q and h["t"] not in removed_hosts}))})
    for h in sorted(shared, key=lambda x: (x["t"], x["ts"], x["q"], x["qs"])):
        rows.append({"kind": "shared_stretch", "host_assembly": genome_of(h["t"]),
                     "host_organism": organism.get(genome_of(h["t"]), ""), "host_contig": h["t"],
                     "host_length": h["tlen"], "host_start": h["ts"], "host_end": h["te"],
                     "host_covered": f"{host_cov.get(h['t'], 0):.3f}", "bacterial_genome": genome_of(h["q"]),
                     "bacterial_organism": organism.get(genome_of(h["q"]), ""), "bacterial_contig": h["q"],
                     "bacterial_length": h["qlen"], "bacterial_start": h["qs"], "bacterial_end": h["qe"],
                     "bacterial_covered": f"{bact_cov.get(h['q'], 0):.3f}",
                     "identity": f"{h['match'] / h['block']:.4f}", "aligned_bp": h["block"]})

    # Cleaned files and their stats.
    cdir = work / "clean"
    files = []
    removed = removed_hosts | removed_bact
    for (r, role) in [(h, "host") for h in hosts] + bacteria:
        acc = r["accession"]
        src = gdir / f"{acc}.fna"
        dest = cdir / ("hosts" if role == "host" else "bacteria") / f"{acc}.fna"
        drop = {n for n in removed if genome_of(n) == acc}
        lengths = fasta_lengths(src) if drop else {}
        write_without(src, dest, drop)
        n0, bp0, md0 = fasta_stats(src)
        n1, bp1, md1 = fasta_stats(dest)
        files.append({"accession": acc, "role": role, "organism": r["organism"], "downloaded_sequences": n0,
                      "downloaded_bp": bp0, "downloaded_md5": md0, "ncbi_md5_checked": "yes" if md5_checked[acc] else "no",
                      "removed_sequences": len(drop), "removed_bp": sum(lengths[n] for n in drop),
                      "clean_sequences": n1, "clean_bp": bp1, "clean_md5": md1})
    (cdir).mkdir(exist_ok=True)
    with open(cdir / "organelles.fna.part", "wb") as out:
        for r in orgs:
            src = odir / f"{r['accession']}.fna"
            with open(src, "rb") as fh:
                shutil.copyfileobj(fh, out)
            n, bp, md = fasta_stats(src)
            files.append({"accession": r["accession"], "role": "organelle",
                          "organism": f"{r['plant']} {r['organelle']}", "downloaded_sequences": n,
                          "downloaded_bp": bp, "downloaded_md5": md, "ncbi_md5_checked": "length checked",
                          "removed_sequences": 0, "removed_bp": 0, "clean_sequences": n, "clean_bp": bp,
                          "clean_md5": md})
    os.replace(cdir / "organelles.fna.part", cdir / "organelles.fna")

    # CLEANING.md
    now = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    by_acc = {f["accession"]: f for f in files}
    md = ["# Reference cleaning", "",
          f"Made by clean_references.py v{__version__} on {now}, with minimap2 {mm2_version}, "
          f"from genomes downloaded through NCBI Datasets API {api_version or 'version not reported'}.",
          "Rule: PROTOCOL.md section 3.4. Settings the protocol leaves open: section 10, entry 3.", "",
          "## Commands (in the work folder)", "", "```"] + commands + ["```", "",
          f"Alignments found: {len(hits)}; of these, {len(kept)} pass step 1 (at least {MIN_BLOCK} bp at "
          f"{MIN_IDENTITY_PCT}% identity or more).", "",
          "## Host assemblies", "", "| Assembly | Organism | Contigs | Removed as bacterial | Removed bp |",
          "|---|---|---|---|---|"]
    for h in hosts:
        f = by_acc[h["accession"]]
        md.append(f"| {f['accession']} | {f['organism']} | {f['downloaded_sequences']} | {f['removed_sequences']} | "
                  f"{f['removed_bp']} |")
    md += ["", "## Bacterial genomes with contigs removed as fungal", ""]
    changed = [by_acc[r["accession"]] for r, _ in bacteria if by_acc[r["accession"]]["removed_sequences"]]
    if changed:
        md += ["| Genome | Organism | Contigs | Removed as fungal | Removed bp |", "|---|---|---|---|---|"]
        md += [f"| {f['accession']} | {f['organism']} | {f['downloaded_sequences']} | {f['removed_sequences']} | "
               f"{f['removed_bp']} |" for f in changed]
    else:
        md.append("None.")
    md += ["", f"Shared stretches left in both (reads from them count for neither): {len(shared)}, "
               f"{sum(h['block'] for h in shared)} bp in all.", "",
           "The held-out MRE genomes took part in all of this, as section 3.4 says.", "",
           f"Downloads: {len(todo)} genomes, {sum(1 for a in todo if md5_checked[a])} of them checked against "
           f"the md5 in NCBI's own download; {len(orgs)} organelle genomes, each checked against its length "
           "in organelles.tsv.", "", "Details: `overlaps.tsv` (every removed contig and every shared stretch) and "
           "`genome_files.tsv` (every file, with md5)."]

    # Everything worked: now write the outputs to commit.
    parts = [write_tsv(refdir / "overlaps.tsv", OVERLAP_COLS, rows),
             write_tsv(refdir / "genome_files.tsv", FILE_COLS, files)]
    cmd_md = refdir / "CLEANING.md.part"
    cmd_md.write_text("\n".join(md) + "\n", encoding="utf-8")
    parts.append(cmd_md)
    for p in parts:
        os.replace(p, p.with_name(p.name[:-len(".part")]))

    print(f"host contigs removed as bacterial: {len(removed_hosts)}")
    print(f"bacterial contigs removed as fungal: {len(removed_bact)}")
    print(f"shared stretches left in both: {len(shared)}")
    print(f"-> {refdir}/overlaps.tsv, genome_files.tsv, CLEANING.md; cleaned genomes in {cdir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
