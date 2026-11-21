#!/usr/bin/env python3
"""Build the Phase 0 reference indexes from the cleaned genomes (PROTOCOL.md 3.4 step 5, 4.1).

Run once, after clean_references.py and before stage 0a. Makes two sets in --workdir/index/:
  full/   every cleaned genome: host assemblies, bacterial genomes (endosymbionts and contaminants)
          and plant organelle genomes
  t2/     the same without the 9 genomes marked held_out_in_T2 (test T2)
Each set holds:
  refs.fna         the genomes in one FASTA, in that order
  sequences.tsv    every sequence in refs.fna with its genome and lineage (read by assign_reads.py)
  refs.sr.mmi      minimap2 index for short reads (-x sr)
  refs.hifi.mmi    minimap2 index for HiFi reads (-x map-hifi)
  bacteria.syldb   sylph database of the cleaned bacterial genomes (the cross-check in 4.2)
  index.json       versions, md5s, counts and commands
Then, only when everything worked, in --refdir (commit these before stage 0a):
  lineages.tsv     every genome with its lineage, and whether it is in each set
  INDEX.md         versions, commands, md5s and counts

Lineages follow PROTOCOL.md 4.1: CaGg (Candidatus Glomeribacter genomes), CaMg/MRE (every genome
marked held_out_in_T2), Mycoavidus, each contaminant genus, the plant organelles, and each host
fungus. The settings the protocol leaves open are in PROTOCOL.md section 10, entry 4.

Needs minimap2 and sylph on PATH (activate the phase0-reads environment) and about 15 GB of disk.

Usage:
  conda activate phase0-reads
  python scripts/make_index.py --workdir ~/phase0-data --threads 16
"""

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

__version__ = "0.1"

MAX_INDEX_BP = 8_000_000_000   # minimap2 -I default (8G): a bigger set would be split into parts
SEQ_COLS = ["sequence", "accession", "kind", "lineage", "length"]
LINEAGE_COLS = ["accession", "organism", "kind", "lineage", "held_out_in_T2", "in_full", "in_t2", "sequences", "bp"]
ORGANELLES = "plant organelles"


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def genus(organism):
    words = organism.split()
    return words[1] if words and words[0] == "Candidatus" and len(words) > 1 else words[0]


def endosymbiont_lineage(row):
    """PROTOCOL.md 4.1: CaMg/MRE is every genome marked held_out_in_T2; then CaGg; then Mycoavidus."""
    if "held_out_in_T2" in row.get("note", ""):
        return "CaMg/MRE"
    if row["organism"].startswith("Candidatus Glomeribacter"):
        return "CaGg"
    if genus(row["organism"]) == "Mycoavidus":
        return "Mycoavidus"
    return None


def genomes_from_lists(refdir, workdir):
    """[(accession, organism, kind, lineage, held_out, path)] in index order, or raise ValueError."""
    clean = workdir / "clean"
    out, problems = [], []
    for r in read_tsv(refdir / "hosts.tsv"):
        out.append((r["accession"], r["organism"], "host fungus", r["organism"], False,
                    clean / "hosts" / f"{r['accession']}.fna"))
    for r in read_tsv(refdir / "endosymbionts.tsv"):
        lin = endosymbiont_lineage(r)
        if lin is None:
            problems.append(f"{r['accession']} ({r['organism']}) fits none of the endosymbiont lineages in 4.1")
        out.append((r["accession"], r["organism"], "endosymbiont", lin, "held_out_in_T2" in r.get("note", ""),
                    clean / "bacteria" / f"{r['accession']}.fna"))
    for r in read_tsv(refdir / "contaminants.tsv"):
        out.append((r["accession"], r["organism"], "contaminant", genus(r["organism"]), False,
                    clean / "bacteria" / f"{r['accession']}.fna"))
    org_file = clean / "organelles.fna"
    for r in read_tsv(refdir / "organelles.tsv"):
        out.append((r["accession"], f"{r['plant']} {r['organelle']}", "organelle", ORGANELLES, False, org_file))
    for acc, _, _, _, _, path in out:
        if not path.exists():
            problems.append(f"{acc}: missing cleaned file {path}")
    host_names = [g[3] for g in out if g[2] == "host fungus"]
    other = {g[3] for g in out if g[2] != "host fungus"}
    if len(set(host_names)) != len(host_names) or set(host_names) & other:
        problems.append("lineage names clash: each host fungus needs its own name, different from every other lineage")
    if problems:
        raise ValueError("\n".join(problems))
    return out


def fasta_records(path):
    """Yield (name, length) for each record."""
    name, n = None, 0
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    yield name, n
                name, n = line[1:].split(None, 1)[0].decode("ascii"), 0
            else:
                n += len(line.strip())
    if name is not None:
        yield name, n


def build_set(name, genomes, outdir, args, log):
    """Write refs.fna, sequences.tsv, both minimap2 indexes and the sylph database for one set."""
    outdir.mkdir(parents=True, exist_ok=True)
    refs = outdir / "refs.fna"
    seq_rows, seen, md5, total = [], set(), hashlib.md5(), 0
    organelle_accs = {g[0] for g in genomes if g[2] == "organelle"}
    per_genome = Counter()
    per_genome_bp = Counter()
    written_files = set()
    with open(refs.with_name("refs.fna.part"), "wb") as out:
        for acc, _, kind, lin, _, path in genomes:
            if path in written_files:          # organelles.fna holds all organelle genomes
                continue
            written_files.add(path)
            last = b"\n"
            with open(path, "rb") as fh:
                for line in fh:
                    out.write(line)
                    md5.update(line)
                    last = line
            if not last.endswith(b"\n"):      # never let a file run into the next one
                out.write(b"\n")
                md5.update(b"\n")
            for seq, n in fasta_records(path):
                if seq in seen:
                    raise ValueError(f"sequence name {seq} appears twice")
                seen.add(seq)
                seq_acc = seq.split("|", 1)[0]
                if kind == "organelle":
                    if seq_acc not in organelle_accs:
                        raise ValueError(f"{seq} in organelles.fna is not in organelles.tsv")
                    row_kind, row_lin = "organelle", ORGANELLES
                else:
                    if seq_acc != acc:
                        raise ValueError(f"{seq} is in the file for {acc}")
                    row_kind, row_lin = kind, lin
                seq_rows.append({"sequence": seq, "accession": seq_acc, "kind": row_kind, "lineage": row_lin,
                                 "length": n})
                per_genome[seq_acc] += 1
                per_genome_bp[seq_acc] += n
                total += n
    os.replace(refs.with_name("refs.fna.part"), refs)
    if total >= MAX_INDEX_BP:
        raise ValueError(f"{name}: {total} bp is too much for one minimap2 index part ({MAX_INDEX_BP})")
    missing = [g[0] for g in genomes if per_genome[g[0]] == 0]
    if missing:
        raise ValueError(f"{name}: no sequences for {', '.join(missing)}")

    tmp = outdir / "sequences.tsv.part"
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=SEQ_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(seq_rows)
    os.replace(tmp, outdir / "sequences.tsv")

    commands = []
    for preset, fname in (("sr", "refs.sr.mmi"), ("map-hifi", "refs.hifi.mmi")):
        dest = outdir / fname
        cmd = [args.minimap2, "-x", preset, "-t", str(args.threads), "-d", str(dest) + ".part", str(refs)]
        log(f"{name}: minimap2 index for -x {preset}")
        with open(outdir / "minimap2.log", "a") as err:
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=err, check=True)
        os.replace(str(dest) + ".part", dest)
        commands.append(f"minimap2 -x {preset} -t {args.threads} -d index/{name}/{fname} index/{name}/refs.fna")

    bact = [g[5] for g in genomes if g[2] in ("endosymbiont", "contaminant")]
    log(f"{name}: sylph database of {len(bact)} bacterial genomes")
    prefix = outdir / "bacteria"
    for old in (outdir / "bacteria.syldb",):
        if old.exists():
            old.unlink()
    with open(outdir / "sylph.log", "a") as err:
        subprocess.run([args.sylph, "sketch", "-g", *map(str, bact), "-o", str(prefix), "-t", str(args.threads)],
                       stdout=err, stderr=err, check=True)
    if not (outdir / "bacteria.syldb").exists():
        raise ValueError(f"{name}: sylph did not write bacteria.syldb")
    commands.append(f"sylph sketch -g clean/bacteria/<accession>.fna ({len(bact)} files) "
                    f"-o index/{name}/bacteria -t {args.threads}")

    info = {"set": name, "made_by": f"make_index.py v{__version__}",
            "made": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
            "minimap2": args.mm2_version, "sylph": args.sylph_version,
            "refs_md5": md5.hexdigest(), "sequences": len(seq_rows), "bp": total,
            "genomes": len({r["accession"] for r in seq_rows}),
            "lineages": dict(Counter(lin for _, _, _, lin, _, _ in genomes)),
            "syldb_md5": hashlib.md5((outdir / "bacteria.syldb").read_bytes()).hexdigest(),
            "commands": commands}
    (outdir / "index.json").write_text(json.dumps(info, indent=1) + "\n")
    return info, per_genome, per_genome_bp


def version_of(cmd):
    out = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return (out.stdout or out.stderr).strip()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True, help="the work folder used by clean_references.py")
    ap.add_argument("--refdir", default="references")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--sylph", default="sylph")
    args = ap.parse_args(argv)
    work, refdir = Path(args.workdir).expanduser(), Path(args.refdir)

    try:
        args.mm2_version = version_of([args.minimap2, "--version"])
        args.sylph_version = version_of([args.sylph, "--version"])
    except (OSError, subprocess.CalledProcessError) as e:
        print(f"Could not run minimap2 or sylph ({e}). Activate the phase0-reads environment first.", file=sys.stderr)
        return 1
    try:
        genomes = genomes_from_lists(refdir, work)
    except (OSError, ValueError, KeyError) as e:
        print(f"Cannot build the index, nothing was written:\n{e}", file=sys.stderr)
        return 1

    def log(msg):
        print(msg, flush=True)

    sets = {"full": genomes, "t2": [g for g in genomes if not g[4]]}
    infos, counts = {}, {}
    try:
        for name, gs in sets.items():
            infos[name], n_seq, n_bp = build_set(name, gs, work / "index" / name, args, log)
            counts[name] = (n_seq, n_bp)
    except (OSError, ValueError, subprocess.CalledProcessError) as e:
        print(f"Index build failed, nothing was written to {refdir}/: {e}", file=sys.stderr)
        return 1

    held = {g[0] for g in genomes if g[4]}
    rows = []
    for acc, org, kind, lin, ho, _ in genomes:
        rows.append({"accession": acc, "organism": org, "kind": kind, "lineage": lin,
                     "held_out_in_T2": "yes" if ho else "no", "in_full": "yes",
                     "in_t2": "no" if acc in held else "yes",
                     "sequences": counts["full"][0][acc], "bp": counts["full"][1][acc]})
    tmp = refdir / "lineages.tsv.part"
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=LINEAGE_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    lin_order = list(dict.fromkeys(r["lineage"] for r in rows))
    md = ["# Reference indexes", "",
          f"Made by make_index.py v{__version__} on {infos['full']['made']}, with minimap2 {args.mm2_version} "
          f"and {args.sylph_version}.",
          "From the cleaned genomes (`CLEANING.md`). Lineages: PROTOCOL.md 4.1. Settings the protocol leaves open: "
          "section 10, entry 4.", "",
          "| Set | Genomes | Sequences | bp | md5 of refs.fna | md5 of bacteria.syldb |", "|---|---|---|---|---|---|"]
    for name in ("full", "t2"):
        i = infos[name]
        md.append(f"| {name} | {i['genomes']} | {i['sequences']} | {i['bp']} | {i['refs_md5']} | {i['syldb_md5']} |")
    md += ["", "The t2 set leaves out the genomes marked held_out_in_T2: " + ", ".join(sorted(held)) + ".", "",
           "## Lineages", "", "| Lineage | Kind | Genomes | Sequences | bp | In t2 |", "|---|---|---|---|---|---|"]
    for lin in lin_order:
        rs = [r for r in rows if r["lineage"] == lin]
        md.append(f"| {lin} | {rs[0]['kind']} | {len(rs)} | {sum(r['sequences'] for r in rs)} | "
                  f"{sum(r['bp'] for r in rs)} | {'no' if all(r['in_t2'] == 'no' for r in rs) else 'yes'} |")
    md += ["", "## Commands (in the work folder, for each set)", "", "```"] + infos["full"]["commands"] + ["```", "",
           "Details: `lineages.tsv` (every genome and its lineage). Per-sequence tables are in the work folder "
           "(`index/<set>/sequences.tsv`)."]
    tmp_md = refdir / "INDEX.md.part"
    tmp_md.write_text("\n".join(md) + "\n", encoding="utf-8")
    os.replace(tmp, refdir / "lineages.tsv")
    os.replace(tmp_md, refdir / "INDEX.md")

    for name in ("full", "t2"):
        i = infos[name]
        print(f"{name}: {i['genomes']} genomes, {i['sequences']} sequences, {i['bp']} bp")
    print(f"-> {refdir}/lineages.tsv, {refdir}/INDEX.md; indexes in {work / 'index'}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
