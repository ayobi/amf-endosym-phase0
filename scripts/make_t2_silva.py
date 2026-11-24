#!/usr/bin/env python3
"""T2's SILVA: the removal list, and the database without it (PROTOCOL.md section 5, T2; section 10, entry 11).

Section 5 takes out of SILVA, for T2, every sequence whose taxonomy contains Mollicutes, Mycoplasmatota,
Tenericutes, Mycoplasmatales, Entomoplasmatales, Acholeplasmatales, Anaeroplasmatales or Haloplasmatales,
and every sequence whose name mentions Moeniiplasma, Mycoplasma-related endobacteria or MRE. As applied
here (entry 11), each SILVA header after the sequence id (the seven ranks and the organism name) is
searched for:
  - each of the eight taxon names, ignoring case;
  - "Moeniiplasma", and "Mycoplasma-related" (with a hyphen or a space) followed by any word, so
    endobacteria, endobacterium and endosymbiont alike, ignoring case;
  - "MRE" in capitals, not inside a longer word in capitals ("GvMRE" counts, "MREB" does not).
A sequence is removed if any of these is found. Searching the ranks for the names, and the names for the
taxa, only removes more, which makes T2 harder, never easier. As a check of the rule's reach, the report
lists the kept sequences whose genus or name is a Mollicutes genus (CHECK_GENERA); it removes nothing.

Reads WORKDIR/db/silva/ (get_silva.py). Writes:
  t2/silva_removed.tsv   every removed sequence: its id, what matched and its SILVA path
  t2/SILVA_T2.md         the counts and md5s
  (commit both before T2 runs)
  WORKDIR/db/silva_t2/   the T2 database: the FASTA without the removed sequences, taxonomy.tsv, and
                         manifest.json, written last, which run_t2.py checks against t2/silva_removed.tsv

Usage (phase0-reads environment):
  conda activate phase0-reads
  python scripts/make_t2_silva.py --workdir ~/phase0-data
"""

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
import get_silva as GS  # noqa: E402  the VSEARCH version check

__version__ = "0.1"

TAXA = ["Mollicutes", "Mycoplasmatota", "Tenericutes", "Mycoplasmatales", "Entomoplasmatales", "Acholeplasmatales",
        "Anaeroplasmatales", "Haloplasmatales"]                                     # section 5
TERMS = [(t, re.compile(re.escape(t), re.IGNORECASE)) for t in TAXA] + [
    ("Moeniiplasma", re.compile("moeniiplasma", re.IGNORECASE)),
    ("Mycoplasma-related", re.compile("mycoplasma[- ]related", re.IGNORECASE)),
    ("MRE", re.compile(r"(?<![A-Z])MRE(?![A-Za-z])"))]
REMOVED_COLS = ["sequence", "matched", "taxonomy"]
CHECK_GENERA = ["Acholeplasma", "Anaeroplasma", "Asteroleplasma", "Entomoplasma", "Haloplasma", "Hepatoplasma",
                "Malacoplasma", "Mesomycoplasma", "Mesoplasma", "Metamycoplasma", "Mycoplasma", "Mycoplasmoides",
                "Mycoplasmopsis", "Phytoplasma", "Spiroplasma", "Ureaplasma"]      # entry 11: a check, removes nothing
CHECK = re.compile(r"\b(" + "|".join(CHECK_GENERA) + r")\b", re.IGNORECASE)


def matches(path):
    """The terms found in a SILVA path (ranks and organism name)."""
    return [name for name, rx in TERMS if rx.search(path)]


def check_genus(path):
    """A Mollicutes genus named in the genus rank or the organism name (the report's check), or None."""
    hit = CHECK.search(";".join(path.split(";")[6:]))
    return hit.group(1) if hit else None


def md5_text(text):
    return hashlib.md5(text.encode("utf-8")).hexdigest()


def split_silva(fasta, out_fasta, out_tax):
    """Write the T2 FASTA and taxonomy table; return (total, removed rows, kept sequences the check flags,
    md5 of the input, md5 of the output)."""
    removed, flagged, total, keep = [], [], 0, True
    md_in, md_out = hashlib.md5(), hashlib.md5()
    tmp_fa, tmp_tax = out_fasta.with_name(out_fasta.name + ".part"), out_tax.with_name(out_tax.name + ".part")
    with open(fasta, "rb") as src, open(tmp_fa, "wb") as fa, open(tmp_tax, "w", encoding="ascii",
                                                                     errors="replace") as tax:
        tax.write("accession\ttaxonomy\n")
        for raw in src:
            md_in.update(raw)
            if raw.startswith(b">"):
                total += 1
                sid, _, path = raw[1:].decode("ascii", "replace").rstrip("\n").partition(" ")
                hit = matches(path)
                keep = not hit
                if hit:
                    removed.append({"sequence": sid, "matched": ",".join(hit), "taxonomy": path})
                else:
                    tax.write(f"{sid}\t{path}\n")
                    g = check_genus(path)
                    if g:
                        flagged.append((sid, g, path))
            if keep:
                fa.write(raw)
                md_out.update(raw)
    os.replace(tmp_fa, out_fasta)
    os.replace(tmp_tax, out_tax)
    return total, removed, flagged, md_in.hexdigest(), md_out.hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--t2-dir", default="t2", help="where the removal list goes [t2]")
    ap.add_argument("--vsearch", default="vsearch")
    args = ap.parse_args(argv)
    work, t2_dir = Path(args.workdir).expanduser(), Path(args.t2_dir)
    out = work / "db" / "silva_t2"
    try:
        silva = json.loads((work / "db" / "silva" / "manifest.json").read_text())
        vsearch = GS.vsearch_version(args.vsearch)
        if vsearch != silva["vsearch"]:
            raise RuntimeError(f"VSEARCH {vsearch} is not the {silva['vsearch']} get_silva.py checked")
        out.mkdir(parents=True, exist_ok=True)
        t2_dir.mkdir(parents=True, exist_ok=True)
        (out / "manifest.json").unlink(missing_ok=True)
        fasta = Path(silva["fasta"])
        t2_fasta, t2_tax = out / fasta.name.replace(".fasta", ".T2.fasta"), out / "taxonomy.tsv"
        print(f"reading {fasta.name}", flush=True)
        total, removed, flagged, md_in, md_out = split_silva(fasta, t2_fasta, t2_tax)
        if total != silva["sequences"]:
            raise RuntimeError(f"{fasta} holds {total:,} sequences; get_silva.py counted {silva['sequences']:,}")
        if not removed:
            raise RuntimeError("no sequence matched; has the SILVA header format changed?")
    except (OSError, ValueError, KeyError, RuntimeError) as e:
        print(f"Failed, so nothing was written to {t2_dir}/ and no manifest: {e}", file=sys.stderr)
        return 1

    lines = ["\t".join(REMOVED_COLS)] + [f"{r['sequence']}\t{r['matched']}\t{r['taxonomy']}" for r in removed]
    list_text = "\n".join(lines) + "\n"
    per_term = {name: sum(name in r["matched"].split(",") for r in removed) for name, _ in TERMS}
    made = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    md = ["# SILVA for T2", "",
          f"Made by make_t2_silva.py v{__version__} on {made}, from `{silva['file']}` ({total:,} sequences; the "
          f"download's md5 {silva['md5']}, as SILVA publishes it; the unpacked FASTA's md5 {md_in}).", "",
          "Rule: PROTOCOL.md section 5 (T2), as applied in section 10, entry 11. Each header after the sequence id "
          "(the seven ranks and the organism name) is searched for the eight taxon names, Moeniiplasma and "
          "Mycoplasma-related (followed by any word), ignoring case, and for MRE in capitals, not inside a longer "
          "word in capitals. A sequence is removed if any of them is found.", "",
          f"**Removed: {len(removed):,} sequences**; kept: {total - len(removed):,} (md5 of the T2 FASTA {md_out}).",
          "", "| Term | Sequences it matches |", "|---|---|"]
    md += [f"| {name} | {per_term[name]:,} |" for name, _ in TERMS]
    md += ["", "A sequence can match several terms. Every removed sequence, with what matched and its SILVA path: "
               "`silva_removed.tsv`. run_t2.py checks that the T2 database was made from this list."]
    md += ["", "## Check: Mollicutes genera left in", "",
           "A check of the rule's reach (entry 11); it removes nothing. Kept sequences whose genus rank or organism "
           "name holds one of " + ", ".join(CHECK_GENERA) + f": **{len(flagged):,}**."]
    if flagged:
        by_genus = {}
        for _, g, _ in flagged:
            by_genus[g.capitalize()] = by_genus.get(g.capitalize(), 0) + 1
        md += ["", "| Genus named | Sequences |", "|---|---|"]
        md += [f"| {g} | {n:,} |" for g, n in sorted(by_genus.items())]
        md += ["", "| Sequence | SILVA path |", "|---|---|"]
        md += [f"| {sid} | {path} |" for sid, _, path in flagged[:50]]
        if len(flagged) > 50:
            md.append(f"| ... | {len(flagged) - 50:,} more |")
    (t2_dir / "silva_removed.tsv").write_text(list_text, encoding="utf-8")
    (t2_dir / "SILVA_T2.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    manifest = {"made_by": f"make_t2_silva.py v{__version__}", "made": made, "source": silva["file"],
                "source_md5": silva["md5"], "source_fasta_md5": md_in, "sequences": total - len(removed),
                "removed": len(removed), "removed_list_md5": md5_text(list_text), "fasta_md5": md_out,
                "fasta": str(t2_fasta.resolve()), "taxonomy": str(t2_tax.resolve()), "vsearch": vsearch,
                "file": t2_fasta.name}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    print(f"\nremoved {len(removed):,} of {total:,} sequences: "
          + ", ".join(f"{name} {n:,}" for name, n in per_term.items() if n))
    if flagged:
        print(f"check: {len(flagged):,} kept sequences name a Mollicutes genus; see {t2_dir}/SILVA_T2.md")
    print(f"-> {t2_dir}/silva_removed.tsv, {t2_dir}/SILVA_T2.md; the T2 database in {out}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
