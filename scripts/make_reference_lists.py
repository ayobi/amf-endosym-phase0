#!/usr/bin/env python3
"""Make the Phase 0 reference lists from NCBI. Run it on the lock date.

Writes to --outdir (default: references/):
  endosymbionts.tsv  every GenBank assembly under Candidatus Moeniiplasma,
                     Candidatus Glomeribacter and Mycoavidus, plus the bacterial
                     assemblies in the control BioProjects (MRE genomes filed under
                     other names, such as GvMRE). MRE genomes are marked
                     "held_out_in_T2".
  contaminants.tsv   the RefSeq reference genome of each kit-contaminant type species
  organelles.tsv     RefSeq plastid and mitochondrial genomes of the host plants;
                     GenBank complete genomes where RefSeq has none
  hosts.tsv          one assembly per control host (rule in pick_host below)
  QUERIES.md         the date, the NCBI Datasets version, every query and the counts

Nothing is written unless every query succeeds. It first checks that both NCBI
services answer (exit code 2 if not), and it writes no files if any query fails
(exit code 1), so a bad run never leaves short lists behind.

Standard library only. Needs network access to api.ncbi.nlm.nih.gov and
eutils.ncbi.nlm.nih.gov. Set NCBI_API_KEY to go faster (optional).

Usage:
  python make_reference_lists.py -o references
"""

import argparse
import csv
import datetime as dt
import json
import os
import re
import socket
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

__version__ = "0.2"

DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

ENDO_TAXA = ["Candidatus Moeniiplasma", "Candidatus Glomeribacter", "Mycoavidus"]
# Control BioProjects that hold endosymbiont genomes under other names.
ENDO_BIOPROJECTS = ["PRJNA407498", "PRJNA670049", "PRJNA1176334"]
ENDO_NAME = re.compile(r"moeniiplasma|glomeribacter|mycoavidus|mycoplasma|mollicutes|"
                       r"mre\b|endosymbiont|endobacteri", re.I)
MRE_NAME = re.compile(r"moeniiplasma|mycoplasma|mollicutes|mre\b", re.I)
# MRE genomes are often filed under loose names ("Mycoplasmataceae bacterium",
# "endosymbiont GvMRE of ..."), so the taxon queries miss them. This Assembly search
# catches them; everything it finds is held out in T2.
MRE_CROSSCHECK = ('("Mollicutes"[Organism] OR "Mycoplasmatales"[Organism]) AND '
                  '(endosymbiont[All Fields] OR MRE[All Fields] OR Moeniiplasma[All Fields] '
                  'OR Glomeromycota[All Fields])')

# Kit-contaminant genera (Salter et al. 2014) and their type species.
CONTAMINANT_SPECIES = {
    "Ralstonia": "Ralstonia pickettii",
    "Bradyrhizobium": "Bradyrhizobium japonicum",
    "Burkholderia": "Burkholderia cepacia",
    "Herbaspirillum": "Herbaspirillum seropedicae",
    "Pseudomonas": "Pseudomonas aeruginosa",
    "Sphingomonas": "Sphingomonas paucimobilis",
    "Stenotrophomonas": "Stenotrophomonas maltophilia",
}

PLANTS = ["Daucus carota", "Lotus japonicus", "Medicago truncatula", "Trifolium repens"]
ORGANELLES = ["chloroplast", "mitochondrion"]

# Control hosts: (label, query kind, query, organism-name filter or None)
HOSTS = [
    ("Gigaspora margarita", "taxon", "Gigaspora margarita", None),
    ("Rhizophagus irregularis", "taxon", "Rhizophagus irregularis", None),
    ("Rhizophagus sp. LPA8", "bioproject", "PRJNA670049", re.compile(r"LPA8", re.I)),
    ("Diversispora epigaea", "taxon", "Diversispora epigaea", None),
    ("Linnemannia elongata (T5 only)", "taxon", "Linnemannia elongata", None),
]
LEVEL_RANK = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}

ASSEMBLY_COLS = ["accession", "organism", "strain_or_isolate", "tax_id", "assembly_level",
                 "size_bp", "gc_percent", "release_date", "bioproject", "refseq_category",
                 "checkm_completeness", "checkm_contamination", "query", "note"]
ORGANELLE_COLS = ["plant", "organelle", "accession", "title", "length_bp", "source", "query"]

API_KEY = os.environ.get("NCBI_API_KEY", "")
PAUSE = 0.12 if API_KEY else 0.4   # NCBI allows 10 requests/s with a key, 3 without


# ------------------------------------------------------------------ network ---

def get_json(url, params=None):
    """GET a JSON document, with a pause and a few retries. Tests replace this."""
    params = dict(params or {})
    headers = {"User-Agent": f"amf-endosym-phase0/{__version__}", "Accept": "application/json"}
    if API_KEY:
        if url.startswith(EUTILS):
            params["api_key"] = API_KEY
        else:
            headers["api-key"] = API_KEY
    full = url + ("?" + urllib.parse.urlencode(params) if params else "")
    last, cause = None, None
    for attempt in range(4):
        time.sleep(PAUSE * (2 ** attempt))
        try:
            req = urllib.request.Request(full, headers=headers)
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            last, cause = f"HTTP {e.code} for {full}", e
            if e.code not in (429, 500, 502, 503, 504):
                break
        except (urllib.error.URLError, TimeoutError) as e:
            last, cause = f"{e} for {full}", e
    raise RuntimeError(last) from cause


def check_ncbi():
    """Ask each NCBI service for a tiny answer before any query is run.

    Returns (problems, datasets_version). A failed host-name lookup gets a plain
    message, because that is the local network's DNS, not NCBI.
    """
    problems, version = [], ""
    for name, url, params in (("NCBI Datasets", f"{DATASETS}/version", None),
                              ("NCBI E-utilities", f"{EUTILS}/einfo.fcgi", {"retmode": "json"})):
        host = urllib.parse.urlsplit(url).hostname
        try:
            data = get_json(url, params)
        except RuntimeError as e:
            if isinstance(getattr(e.__cause__, "reason", None), socket.gaierror):
                problems.append(f"{name}: this computer could not look up {host}. That is a DNS problem on "
                                f"this network, not at NCBI. Try another network, or set your DNS servers "
                                f"to 1.1.1.1 and 8.8.8.8, then run again.")
            else:
                problems.append(f"{name}: {host} did not answer ({e}).")
            continue
        if name == "NCBI Datasets":
            version = str(data.get("version", ""))
    return problems, version


def dataset_reports(kind, value, **filters):
    """All assembly reports for a taxon, BioProject or accession list (follows pages)."""
    url = f"{DATASETS}/genome/{kind}/{urllib.parse.quote(value, safe='')}/dataset_report"
    params = {"page_size": 1000}
    params.update({f"filters.{k}": v for k, v in filters.items()})
    reports, token = [], None
    while True:
        if token:
            params["page_token"] = token
        data = get_json(url, params)
        reports.extend(data.get("reports", []))
        token = data.get("next_page_token")
        if not token:
            return reports


# ------------------------------------------------------------------ parsing ---

def assembly_row(rep, query, note=""):
    org = rep.get("organism", {})
    names = org.get("infraspecific_names", {}) or {}
    info = rep.get("assembly_info", {})
    stats = rep.get("assembly_stats", {})
    checkm = rep.get("checkm_info", {}) or {}
    return {
        "accession": rep.get("accession", ""),
        "organism": org.get("organism_name", ""),
        "strain_or_isolate": names.get("strain") or names.get("isolate") or "",
        "tax_id": str(org.get("tax_id", "")),
        "assembly_level": info.get("assembly_level", ""),
        "size_bp": str(stats.get("total_sequence_length", "")),
        "gc_percent": str(stats.get("gc_percent", "")),
        "release_date": info.get("release_date", ""),
        "bioproject": info.get("bioproject_accession", ""),
        "refseq_category": info.get("refseq_category", ""),
        "checkm_completeness": str(checkm.get("completeness", "")),
        "checkm_contamination": str(checkm.get("contamination", "")),
        "query": query,
        "note": note,
    }


def pick_host(reports):
    """One assembly per host: the NCBI reference genome if there is one; otherwise
    the highest assembly level, then the largest contig N50, then the latest release."""
    if not reports:
        return None
    refs = [r for r in reports if r.get("assembly_info", {}).get("refseq_category") == "reference genome"]
    pool = refs or reports

    def key(r):
        info, stats = r.get("assembly_info", {}), r.get("assembly_stats", {})
        return (LEVEL_RANK.get(info.get("assembly_level", ""), 0),
                int(stats.get("contig_n50", 0) or 0),
                info.get("release_date", ""))
    return max(pool, key=key)


# ---------------------------------------------------------------- the lists ---

def endosymbionts(log):
    rows, seen = [], set()
    queries = [("taxon", t) for t in ENDO_TAXA] + [("bioproject", b) for b in ENDO_BIOPROJECTS]
    for kind, value in queries:
        q = f"{kind}: {value} (GenBank)"
        try:
            reps = dataset_reports(kind, value, assembly_source="genbank")
        except RuntimeError as e:
            log["errors"].append(f"{q}: {e}")
            continue
        kept = 0
        for rep in reps:
            name = rep.get("organism", {}).get("organism_name", "")
            if kind == "bioproject" and not ENDO_NAME.search(name):
                continue  # host fungal assemblies in the same BioProject
            acc = rep.get("accession", "")
            if acc in seen:
                continue
            seen.add(acc)
            note = "held_out_in_T2" if MRE_NAME.search(name) else ""
            rows.append(assembly_row(rep, q, note))
            kept += 1
        log["queries"].append((q, len(reps), kept))

    q = f"assembly: {MRE_CROSSCHECK}"
    try:
        accs = assembly_accessions(MRE_CROSSCHECK)
        new = [a for a in accs if a not in seen]
        reps = dataset_reports("accession", ",".join(new)) if new else []
        for rep in reps:
            seen.add(rep.get("accession", ""))
            rows.append(assembly_row(rep, "MRE cross-check (NCBI Assembly search)", "held_out_in_T2"))
        log["queries"].append((q, len(accs), len(reps)))
    except RuntimeError as e:
        log["errors"].append(f"{q}: {e}")
    return rows


def assembly_accessions(term):
    """GenBank accessions from an NCBI Assembly (E-utilities) search."""
    ids = get_json(f"{EUTILS}/esearch.fcgi",
                   {"db": "assembly", "retmode": "json", "retmax": 500, "term": term}
                   ).get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    res = get_json(f"{EUTILS}/esummary.fcgi",
                   {"db": "assembly", "retmode": "json", "id": ",".join(ids)}).get("result", {})
    accs = []
    for uid in res.get("uids", []):
        rec = res.get(uid, {})
        acc = rec.get("synonym", {}).get("genbank") or rec.get("assemblyaccession", "")
        if acc:
            accs.append(acc)
    return accs


def contaminants(log):
    rows = []
    for genus, species in CONTAMINANT_SPECIES.items():
        q = f"taxon: {species} (RefSeq reference genome)"
        try:
            reps = dataset_reports("taxon", species, reference_only="true", assembly_source="refseq")
        except RuntimeError as e:
            log["errors"].append(f"{q}: {e}")
            continue
        log["queries"].append((q, len(reps), min(len(reps), 1)))
        if not reps:
            log["errors"].append(f"{q}: no reference genome found")
            continue
        rep = sorted(reps, key=lambda r: r.get("accession", ""))[0]
        rows.append(assembly_row(rep, q, f"type species of {genus}"))
    return rows


def esummaries(term):
    ids = get_json(f"{EUTILS}/esearch.fcgi",
                   {"db": "nuccore", "retmode": "json", "retmax": 50, "term": term}
                   ).get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    res = get_json(f"{EUTILS}/esummary.fcgi",
                   {"db": "nuccore", "retmode": "json", "id": ",".join(ids)}).get("result", {})
    return [res[i] for i in res.get("uids", []) if i in res]


def organelles(log):
    rows = []
    for plant in PLANTS:
        for org in ORGANELLES:
            hits, source = [], ""
            for src, term in [
                ("RefSeq", f'"{plant}"[Organism] AND refseq[filter] AND {org}[filter]'),
                ("GenBank", f'"{plant}"[Organism] AND {org}[filter] AND "complete genome"[Title]'),
            ]:
                try:
                    found = [h for h in esummaries(term)
                             if h.get("genome") == org and h.get("completeness", "complete") == "complete"]
                except RuntimeError as e:
                    log["errors"].append(f"{term}: {e}")
                    found = []
                log["queries"].append((f"nuccore: {term}", len(found), len(found)))
                if found:
                    hits, source = found, src
                    break
            if not hits:
                log["missing"].append(f"{plant} {org}: none in RefSeq or GenBank")
            for h in sorted(hits, key=lambda h: h.get("accessionversion", "")):
                rows.append({"plant": plant, "organelle": org,
                             "accession": h.get("accessionversion", ""),
                             "title": h.get("title", ""), "length_bp": str(h.get("slen", "")),
                             "source": source, "query": term})
    return rows


def hosts(log):
    rows = []
    for label, kind, value, name_filter in HOSTS:
        q = f"{kind}: {value} (GenBank)"
        try:
            reps = dataset_reports(kind, value, assembly_source="genbank")
        except RuntimeError as e:
            log["errors"].append(f"{q}: {e}")
            continue
        if name_filter:
            reps = [r for r in reps if name_filter.search(r.get("organism", {}).get("organism_name", ""))]
        chosen = pick_host(reps)
        log["queries"].append((q, len(reps), 1 if chosen else 0))
        if chosen is None:
            log["missing"].append(f"host {label}: no assembly found")
            continue
        rows.append(assembly_row(chosen, q, f"host: {label}; picked from {len(reps)} assemblies"))
    return rows


# ------------------------------------------------------------------- output ---

def write_tsv(path, cols, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", default="references")
    args = ap.parse_args(argv)
    out = Path(args.outdir)

    problems, api_version = check_ncbi()
    if problems:
        print("Could not reach NCBI, so no queries were run and nothing was written:",
              *problems, sep="\n  ", file=sys.stderr)
        return 2

    log = {"queries": [], "errors": [], "missing": []}
    lists = {
        "endosymbionts.tsv": (ASSEMBLY_COLS, endosymbionts(log)),
        "contaminants.tsv": (ASSEMBLY_COLS, contaminants(log)),
        "organelles.tsv": (ORGANELLE_COLS, organelles(log)),
        "hosts.tsv": (ASSEMBLY_COLS, hosts(log)),
    }
    if log["errors"]:
        print(f"Some NCBI queries failed, so nothing was written to {out}/. Run it again later:",
              *log["errors"], sep="\n  ", file=sys.stderr)
        return 1

    out.mkdir(parents=True, exist_ok=True)
    for name, (cols, rows) in lists.items():
        write_tsv(out / name, cols, rows)

    now = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    held = sum(1 for r in lists["endosymbionts.tsv"][1] if r["note"] == "held_out_in_T2")
    md = [f"# Reference lists", "",
          f"Made by make_reference_lists.py v{__version__} on {now}, "
          f"against NCBI Datasets API {api_version or 'version not reported'}.", "",
          "| File | Entries |", "|---|---|"]
    md += [f"| {name} | {len(rows)} |" for name, (_, rows) in lists.items()]
    md += ["", f"MRE genomes marked held_out_in_T2: {held}", "", "## Queries", "",
           "| Query | Hits | Kept |", "|---|---|---|"]
    md += [f"| {q} | {n} | {k} |" for q, n, k in log["queries"]]
    if log["missing"]:
        md += ["", "## Not available"] + [f"- {m}" for m in log["missing"]]
    (out / "QUERIES.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    for name, (_, rows) in lists.items():
        print(f"{name}: {len(rows)}")
    print(f"MRE genomes held out in T2: {held}")
    for m in log["missing"]:
        print(f"not available: {m}")
    print(f"-> {out}/ (see QUERIES.md)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
