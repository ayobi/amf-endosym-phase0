#!/usr/bin/env python3
"""Get SILVA 144 SSU Ref NR99 for SSU placement (references/DATABASES.md; PROTOCOL.md 4.1 and section 10,
entry 7).

Downloads SILVA_144_SSURef_NR99_tax_silva_trunc.fasta.gz (release 144 offers SSU Ref NR99 only as
truncated files; this is the unaligned one), with the md5 file SILVA publishes beside it. It tries the release folder first and then SILVA's
"current" folder, which held release 144 on 2026-11-22. Checks:
  - the download matches SILVA's md5;
  - it holds 905,628 sequences, as references/DATABASES.md says;
  - VSEARCH is the version in env/phase0-reads.lock.yml.
Writes, in WORKDIR/db/silva/:
  SILVA_144_SSURef_NR99_tax_silva_trunc.fasta   the database VSEARCH reads (bases as SILVA gives them;
                                                VSEARCH reads U as T)
  taxonomy.tsv                                   each sequence's accession and SILVA taxonomy path
  manifest.json                                  written last, and only when every check passes
Reruns reuse the download.

Usage (phase0-reads environment):
  conda activate phase0-reads
  python scripts/get_silva.py --workdir ~/phase0-data
"""

import argparse
import datetime as dt
import gzip
import json
import os
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
import get_marker_databases as DB  # noqa: E402  fetch, md5_of, locked_version, tool_output

__version__ = "0.1"

FILE = "SILVA_144_SSURef_NR99_tax_silva_trunc.fasta.gz"
BASES = ["https://www.arb-silva.de/fileadmin/silva_databases/release_144/Exports/",
         "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/"]
EXPECTED_SEQUENCES = 905628                                   # references/DATABASES.md


def vsearch_version(vsearch):
    hit = re.search(r"vsearch v(\d+\.\d+\.\d+)", DB.tool_output([vsearch, "--version"]))
    got, want = (hit.group(1) if hit else "missing"), DB.locked_version("phase0-reads", "vsearch")
    if got != want:
        raise RuntimeError(f"VSEARCH is {got}, env/phase0-reads.lock.yml has {want}; activate phase0-reads")
    return got


def published_md5(sdir):
    """SILVA's md5 for FILE, and the folder it came from (the release folder first)."""
    last = None
    for base in BASES:
        try:
            with urllib.request.urlopen(urllib.request.Request(base + FILE + ".md5", headers=DB.HEADERS),
                                        timeout=60) as r:
                text = r.read().decode("ascii", "replace")
        except (urllib.error.URLError, OSError) as e:
            last = e
            continue
        hit = re.search(r"\b([0-9a-f]{32})\b", text)
        if not hit:
            raise RuntimeError(f"{base}{FILE}.md5 holds no md5")
        (sdir / (FILE + ".md5")).write_text(text)
        return hit.group(1), base
    raise RuntimeError(f"could not get {FILE}.md5 from SILVA: {last}")


def unpack(gz, fasta, taxonomy):
    """Write the plain FASTA and the taxonomy table; return the number of sequences."""
    n = 0
    tmp_fa, tmp_tax = fasta.with_name(fasta.name + ".part"), taxonomy.with_name(taxonomy.name + ".part")
    with gzip.open(gz, "rt", encoding="ascii", errors="replace") as src, open(tmp_fa, "w") as fa, \
            open(tmp_tax, "w") as tax:
        tax.write("accession\ttaxonomy\n")
        for line in src:
            fa.write(line)
            if line.startswith(">"):
                acc, _, path = line[1:].rstrip("\n").partition(" ")
                tax.write(f"{acc}\t{path}\n")
                n += 1
    os.replace(tmp_fa, fasta)
    os.replace(tmp_tax, taxonomy)
    return n


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--vsearch", default="vsearch")
    args = ap.parse_args(argv)
    sdir = Path(args.workdir).expanduser() / "db" / "silva"
    sdir.mkdir(parents=True, exist_ok=True)
    gz, fasta, taxonomy = sdir / FILE, sdir / FILE[:-3], sdir / "taxonomy.tsv"
    try:
        version = vsearch_version(args.vsearch)
        print(f"VSEARCH {version}: as locked", flush=True)
        md5, base = published_md5(sdir)
        if not (gz.exists() and DB.md5_of(gz) == md5):
            print(f"downloading {base}{FILE} (about 310 MB)", flush=True)
            gz.unlink(missing_ok=True)
            DB.fetch(base + FILE, gz)
            if DB.md5_of(gz) != md5:
                gz.unlink()
                raise RuntimeError(f"the download does not match SILVA's md5 ({md5})")
        print("unpacking and counting sequences", flush=True)
        n = unpack(gz, fasta, taxonomy)
        if n != EXPECTED_SEQUENCES:
            raise RuntimeError(f"{FILE} holds {n:,} sequences; references/DATABASES.md says {EXPECTED_SEQUENCES:,}")
    except (OSError, RuntimeError, ValueError, EOFError) as e:
        print(f"Failed, so no manifest was written: {e}", file=sys.stderr)
        return 1
    manifest = {"made_by": f"get_silva.py v{__version__}",
                "made": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
                "file": FILE, "url": base + FILE, "md5": md5, "sequences": n, "vsearch": version,
                "fasta": str(fasta.resolve()), "taxonomy": str(taxonomy.resolve())}
    (sdir / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    print(f"\n{FILE}: md5 {md5} as SILVA publishes it, {n:,} sequences as in DATABASES.md")
    print(f"-> {sdir / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
