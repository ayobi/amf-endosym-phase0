#!/usr/bin/env python3
"""Get the databases for the marker-gene and SSU route (references/DATABASES.md; PROTOCOL.md section 10,
entry 6).

  1. The SingleM metapackage S6.5.0.GTDB_r232.metapackage_20260319, from its own Zenodo record
     (doi:10.5281/zenodo.20150069). "singlem data" would take the newest metapackage in the series
     instead, because SingleM 0.21.4 does not pin the version. The download is checked against the
     md5 in references/DATABASES.md, and "singlem data --verify-only" then checks version 6.5.0 and
     every file in it.
  2. Metaxa2 2.2.3 from microbiology.se, and from it SSUplex's five SSU profile files, made exactly as
     SSUplex's dev/assemble_metaxa2_hmms.sh makes them (mitochondria = Metaxa2's M and N sets).
  3. Checks that SingleM, SSUplex and HMMER are the versions in env/*.lock.yml, and that SSUplex's
     quick start on Metaxa2's test.fasta gives what the SSUplex README says it gives with HMMER 3.4
     (12 bacteria, 10 archaea, 10 eukaryota, 9 mitochondria, 9 chloroplast, 50 unclassified).

Files go to WORKDIR/db/. Only when every check passes does it write WORKDIR/db/manifest.json, which
run_stage0a_m.py needs. Reruns reuse the downloads.

Usage (in the phase0-reads environment; SingleM is called through conda):
  conda activate phase0-reads
  python scripts/get_marker_databases.py --workdir ~/phase0-data
"""

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

__version__ = "0.1"

REPO = Path(__file__).resolve().parent.parent
METAPACKAGE = "S6.5.0.GTDB_r232.metapackage_20260319"          # references/DATABASES.md
METAPACKAGE_DOI = "10.5281/zenodo.20150069"
METAPACKAGE_MD5 = "580f9ae9ae6a0972f15665c7397c6369"
ZENODO_API = "https://zenodo.org/api/records/"
METAXA2_URL = "https://microbiology.se/sw/Metaxa2_2.2.3.tar.gz"
HMM_SETS = {"bacteria": ["B"], "archaea": ["A"], "eukaryota": ["E"], "mitochondria": ["M", "N"],
            "chloroplast": ["C"]}                                  # dev/assemble_metaxa2_hmms.sh
QUICKSTART = {"bacteria": 12, "archaea": 10, "eukaryota": 10, "mitochondria": 9, "chloroplast": 9,
              "unclassified": 50}                                  # SSUplex README, quick start
HEADERS = {"User-Agent": f"amf-endosym-phase0 get_marker_databases.py/{__version__}"}


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url, dest, want_md5=None):
    """Download url to dest through dest.part, with retries; check the md5 when one is given."""
    tmp = dest.with_name(dest.name + ".part")
    last = None
    for attempt in range(4):
        if attempt:
            time.sleep(10 * attempt)
        try:
            h = hashlib.md5()
            with urllib.request.urlopen(urllib.request.Request(url, headers=HEADERS), timeout=300) as resp, \
                    open(tmp, "wb") as out:
                for chunk in iter(lambda: resp.read(1 << 20), b""):
                    h.update(chunk)
                    out.write(chunk)
            break
        except OSError as e:
            last = e
    else:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"could not download {url}: {last}")
    if want_md5 and h.hexdigest() != want_md5:
        tmp.unlink()
        raise RuntimeError(f"{url} has md5 {h.hexdigest()}, but references/DATABASES.md gives {want_md5}")
    os.replace(tmp, dest)


def run_logged(cmd, log, what, env=None):
    with open(log, "a") as fh:
        fh.write("$ " + " ".join(map(str, cmd)) + "\n")
        fh.flush()
        p = subprocess.run(list(map(str, cmd)), stdout=fh, stderr=subprocess.STDOUT, env=env)
    if p.returncode != 0:
        raise RuntimeError(f"{what} failed (exit {p.returncode}); see {log}")


def tool_output(cmd):
    try:
        p = subprocess.run(list(map(str, cmd)), capture_output=True, text=True)
    except OSError as e:
        raise RuntimeError(f"cannot run {cmd[0]}: {e}") from None
    return (p.stdout or "") + (p.stderr or "")


def locked_version(env_name, package):
    text = (REPO / "env" / f"{env_name}.lock.yml").read_text()
    hit = re.search(rf"^\s*-\s*{re.escape(package)}=([^=\s]+)", text, flags=re.M)
    if not hit:
        raise RuntimeError(f"{package} is not in env/{env_name}.lock.yml")
    return hit.group(1)


def check_versions(args):
    """SSUplex, HMMER (nhmmer, which SSUplex calls) and SingleM here, checked against env/*.lock.yml."""
    if not args.conda:
        raise RuntimeError("conda was not found; give its path with --conda")
    words = tool_output([args.ssuplex, "--version"]).split()
    hmmer = re.search(r"HMMER (\S+)", tool_output(["nhmmer", "-h"]))
    singlem = tool_output([args.conda, "run", "-n", args.singlem_env, "singlem", "--version"]).split()
    got = {"ssuplex": words[-1] if words else "missing", "hmmer": hmmer.group(1) if hmmer else "missing",
           "singlem": singlem[-1] if singlem else "missing"}
    want = {"ssuplex": locked_version("phase0-reads", "ssuplex"), "hmmer": locked_version("phase0-reads", "hmmer"),
            "singlem": locked_version("phase0-markers", "singlem")}
    wrong = [f"{k} is {got[k]}, the lock file has {want[k]}" for k in want if got[k] != want[k]]
    if wrong:
        raise RuntimeError("; ".join(wrong) + f". Activate phase0-reads (SingleM runs from {args.singlem_env}).")
    return got


def parse_ssuplex_summary(path):
    """Total reads and reads per origin from an SSUplex {prefix}.summary.txt."""
    text = Path(path).read_text()
    total = re.search(r"^total input reads:\s+(\d+)\s*$", text, flags=re.M)
    if not total or "Per-origin counts:" not in text:
        raise RuntimeError(f"{path} is not an SSUplex 0.2.0 summary")
    per_origin = {}
    for line in text.split("Per-origin counts:", 1)[1].splitlines():
        f = line.split()
        if len(f) == 2 and f[1].isdigit():
            per_origin[f[0]] = int(f[1])
    return {"total": int(total.group(1)), "per_origin": per_origin}


def singlem_metapackage(args, dbdir, log):
    """Download the metapackage from its Zenodo record unless it is here; verify it either way."""
    sdir = dbdir / "singlem"
    found = sorted(p for p in sdir.glob(METAPACKAGE + "*") if p.is_dir()) if sdir.exists() else []
    if not found:
        sdir.mkdir(parents=True, exist_ok=True)
        record = METAPACKAGE_DOI.rsplit(".", 1)[-1]
        with urllib.request.urlopen(urllib.request.Request(ZENODO_API + record, headers=HEADERS), timeout=60) as r:
            meta = json.load(r)
        files = [f for f in meta.get("files", []) if f["key"].endswith(".tar.gz")]
        if len(files) != 1:
            raise RuntimeError(f"Zenodo record {record} does not hold exactly one .tar.gz file")
        archive = sdir / files[0]["key"]
        print(f"  downloading {archive.name} ({files[0].get('size', 0) / 1e9:.1f} GB) from Zenodo record {record}",
              flush=True)
        fetch(files[0]["links"]["self"], archive, METAPACKAGE_MD5)
        tmp = sdir / "extract.part"
        shutil.rmtree(tmp, ignore_errors=True)
        tmp.mkdir()
        with tarfile.open(archive) as tf:
            tf.extractall(tmp, filter="data")
        tops = [p for p in tmp.iterdir() if p.is_dir()]
        if len(tops) != 1 or not tops[0].name.startswith(METAPACKAGE):
            raise RuntimeError(f"{archive.name} does not hold one folder named {METAPACKAGE}*")
        os.replace(tops[0], sdir / tops[0].name)
        shutil.rmtree(tmp)
        archive.unlink()
        found = [sdir / tops[0].name]
    if len(found) != 1:
        raise RuntimeError(f"{sdir} holds more than one {METAPACKAGE} folder")
    print("  checking it with singlem data --verify-only", flush=True)
    run_logged([args.conda, "run", "-n", args.singlem_env, "--no-capture-output", "singlem", "data", "--verify-only"],
               log, "singlem data --verify-only", env={**os.environ, "SINGLEM_METAPACKAGE_PATH": str(found[0])})
    return found[0]


def ssuplex_profiles(dbdir):
    """Download Metaxa2 2.2.3 unless it is here and make SSUplex's five profile files from it."""
    mdir = dbdir / "metaxa2"
    mdir.mkdir(parents=True, exist_ok=True)
    tarball = mdir / "Metaxa2_2.2.3.tar.gz"
    if not tarball.exists():
        print(f"  downloading {METAXA2_URL}", flush=True)
        fetch(METAXA2_URL, tarball)
    parts, test = {}, None
    with tarfile.open(tarball, "r:gz") as tf:
        for m in tf.getmembers():
            if not m.isfile():
                continue
            hit = re.fullmatch(r"(?:\./)?Metaxa2_2\.2\.3/metaxa2_db/SSU/HMMs/([A-Z])\.hmm", m.name)
            if hit:
                parts[hit.group(1)] = tf.extractfile(m).read()
            elif re.fullmatch(r"(?:\./)?Metaxa2_2\.2\.3/test\.fasta", m.name):
                test = tf.extractfile(m).read()
    missing = [f"{c}.hmm" for codes in HMM_SETS.values() for c in codes if not parts.get(c)]
    if missing or not test:
        raise RuntimeError(f"{tarball.name} lacks " + ", ".join(missing + ([] if test else ["test.fasta"])))
    test_fa = mdir / "test.fasta"
    test_fa.write_bytes(test)
    hmm_dir, tmp = dbdir / "ssuplex_hmms", dbdir / "ssuplex_hmms.part"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir()
    for origin, codes in HMM_SETS.items():
        (tmp / f"{origin}.hmm").write_bytes(b"".join(parts[c] for c in codes))
    shutil.rmtree(hmm_dir, ignore_errors=True)
    os.replace(tmp, hmm_dir)
    return hmm_dir, tarball, test_fa


def quickstart(args, hmm_dir, test_fa, dbdir, log):
    out = dbdir / "ssuplex_quickstart" / "test"
    out.parent.mkdir(parents=True, exist_ok=True)
    run_logged([args.ssuplex, "-i", test_fa, "-o", out, "--hmm-dir", hmm_dir, "--threads", "4", "--force"], log,
               "SSUplex quick start")
    s = parse_ssuplex_summary(out.with_name("test.summary.txt"))
    if s["total"] != 100 or s["per_origin"] != QUICKSTART:
        raise RuntimeError(f"SSUplex quick start gave {s['per_origin']} for {s['total']} sequences; the SSUplex README "
                           f"gives {QUICKSTART} for 100 with HMMER 3.4")
    return s


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--ssuplex", default="ssuplex")
    ap.add_argument("--conda", default=os.environ.get("CONDA_EXE") or shutil.which("conda"))
    ap.add_argument("--singlem-env", default="phase0-markers")
    args = ap.parse_args(argv)
    dbdir = Path(args.workdir).expanduser() / "db"
    dbdir.mkdir(parents=True, exist_ok=True)
    log = dbdir / "get_marker_databases.log"
    try:
        versions = check_versions(args)
        print(f"SingleM {versions['singlem']}, SSUplex {versions['ssuplex']}, HMMER {versions['hmmer']}: as locked")
        print(f"SingleM metapackage {METAPACKAGE}", flush=True)
        package = singlem_metapackage(args, dbdir, log)
        print("Metaxa2 2.2.3 and SSUplex's profile files", flush=True)
        hmm_dir, tarball, test_fa = ssuplex_profiles(dbdir)
        print("SSUplex quick start on Metaxa2's test.fasta", flush=True)
        qs = quickstart(args, hmm_dir, test_fa, dbdir, log)
    except (OSError, RuntimeError, ValueError, KeyError, tarfile.TarError) as e:
        print(f"Failed, so no manifest was written: {e}", file=sys.stderr)
        return 1
    hmms = {}
    for origin, codes in HMM_SETS.items():
        path = hmm_dir / f"{origin}.hmm"
        hmms[origin] = {"from": " + ".join(f"{c}.hmm" for c in codes), "md5": md5_of(path),
                        "profiles": sum(1 for line in path.read_text().splitlines() if line.startswith("NAME"))}
    manifest = {
        "made_by": f"get_marker_databases.py v{__version__}",
        "made": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "versions": versions,
        "singlem_metapackage": str(package.resolve()), "singlem_metapackage_name": package.name,
        "singlem_metapackage_doi": METAPACKAGE_DOI, "singlem_metapackage_md5": METAPACKAGE_MD5,
        "metaxa2_url": METAXA2_URL, "metaxa2_md5": md5_of(tarball),
        "ssuplex_hmm_dir": str(hmm_dir.resolve()), "ssuplex_hmms": hmms, "ssuplex_quickstart": qs,
    }
    (dbdir / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    print(f"\nmetapackage {package.name}: md5 as in DATABASES.md, verified by SingleM")
    print("SSUplex profiles: " + ", ".join(f"{o} {v['profiles']}" for o, v in hmms.items()))
    print(f"quick start: {qs['total'] - qs['per_origin']['unclassified']} of {qs['total']} classified, as in the README")
    print(f"-> {dbdir / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
