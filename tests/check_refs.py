"""Check make_reference_lists.py offline, with canned NCBI responses.

The response shapes copy real NCBI Datasets v2 and E-utilities answers seen on
2026-11-20 (trimmed to the fields the script reads).
"""
import contextlib
import csv
import io
import os
import socket
import sys
import tempfile
import urllib.error
import urllib.parse
from pathlib import Path

here = Path(__file__).parent
sys.path.insert(0, str(here.parent / "scripts"))
import make_reference_lists as M  # noqa: E402

fails = []


def rep(acc, name, level="Contig", n50=10000, refcat="", strain=None, isolate=None,
        size="1000", bioproject="PRJNA0", date="2020-01-01"):
    names = {}
    if strain:
        names["strain"] = strain
    if isolate:
        names["isolate"] = isolate
    info = {"assembly_level": level, "release_date": date, "bioproject_accession": bioproject}
    if refcat:
        info["refseq_category"] = refcat
    return {"accession": acc, "organism": {"tax_id": 1, "organism_name": name, "infraspecific_names": names},
            "assembly_info": info, "assembly_stats": {"total_sequence_length": size, "contig_n50": n50,
                                                      "gc_percent": 33}}


DATASETS = {
    ("taxon", "Candidatus Moeniiplasma"): [[rep("GCA_021409705.1", "Candidatus Moeniiplasma glomeromycotorum", isolate="A2_WG"),
                                            rep("GCA_021409605.1", "Candidatus Moeniiplasma glomeromycotorum", isolate="LPA8_4")]],
    ("taxon", "Candidatus Glomeribacter"): [[rep("GCA_000000001.1", "Candidatus Glomeribacter gigasporarum BEG34")]],
    # two pages, to test next_page_token
    ("taxon", "Mycoavidus"): [[rep("GCA_003966915.1", "Mycoavidus cysteinexigens", strain="B1-EB")],
                              [rep("GCA_014218255.1", "Mycoavidus sp. B2-EB", strain="B2-EB")]],
    ("bioproject", "PRJNA407498"): [[rep("GCA_003546845.1", "endosymbiont GvMRE of Glomus versiforme", strain="GvMRE")]],
    ("bioproject", "PRJNA670049"): [[rep("GCA_021464345.1", "Rhizophagus irregularis", isolate="A2_WG_AMF"),
                                     rep("GCA_021409705.1", "Candidatus Moeniiplasma glomeromycotorum", isolate="A2_WG"),
                                     rep("GCA_021464375.1", "Rhizophagus sp. LPA8_1_AMF", n50=14835),
                                     rep("GCA_021464405.1", "Rhizophagus sp. LPA8_3_AMF", n50=14792)]],
    ("bioproject", "PRJNA1176334"): [[rep("GCA_044000001.1", "Candidatus Moeniiplasma glomeromycotorum", isolate="CaMg-9")]],
    ("taxon", "Gigaspora margarita"): [[rep("GCA_100000001.1", "Gigaspora margarita", level="Scaffold", n50=100000),
                                        rep("GCA_100000002.1", "Gigaspora margarita", level="Chromosome", n50=50000)]],
    ("taxon", "Rhizophagus irregularis"): [[rep("GCA_200000001.1", "Rhizophagus irregularis", level="Scaffold",
                                                refcat="reference genome", strain="DAOM 197198"),
                                            rep("GCA_200000002.1", "Rhizophagus irregularis", level="Chromosome")]],
    ("taxon", "Diversispora epigaea"): [[rep("GCA_003550325.1", "Diversispora epigaea", strain="IT104")]],
    ("taxon", "Linnemannia elongata"): [[]],
}
for species, acc in [("Ralstonia pickettii", "GCF_014142625.1"), ("Bradyrhizobium japonicum", "GCF_030584785.1"),
                     ("Burkholderia cepacia", "GCF_009586235.1"), ("Herbaspirillum seropedicae", "GCF_001040945.1"),
                     ("Pseudomonas aeruginosa", "GCF_000006765.1"), ("Sphingomonas paucimobilis", "GCF_016027095.1"),
                     ("Stenotrophomonas maltophilia", "GCF_020641395.2")]:
    DATASETS[("taxon", species)] = [[rep(acc, species, level="Complete Genome", refcat="reference genome")]]

ORG = {  # (plant, organelle, refseq?) -> esummary records
    ("Daucus carota", "chloroplast", True): [("114107112", "NC_008325.1", 155911)],
    ("Daucus carota", "mitochondrion", True): [("2692969881", "NC_087060.1", 250368)],
    ("Lotus japonicus", "chloroplast", True): [("13518417", "NC_002694.1", 150519)],
    ("Lotus japonicus", "mitochondrion", True): [("387866040", "NC_016743.2", 380861)],
    ("Medicago truncatula", "chloroplast", True): [("1253235908", "NC_003119.8", 124033)],
    ("Medicago truncatula", "mitochondrion", True): [("1003725997", "NC_029641.1", 271618)],
    ("Trifolium repens", "chloroplast", False): [("2531238564", "CP125853.1", 126579), ("2497507300", "OQ694001.1", 132344),
                                                 ("1899342395", "MT506238.1", 132429), ("1839353213", "MT120812.1", 132328)],
}
SUMMARY = {uid: (acc, slen, key) for key, recs in ORG.items() for uid, acc, slen in recs}
calls = {"page": {}}
FAIL_ON = set()


ASSEMBLY_SEARCH = {"11732111": "GCA_021409705.1", "18653271": "GCA_031675535.1"}
DATASETS[("accession", "GCA_031675535.1")] = [[rep("GCA_031675535.1", "Mycoplasmataceae bacterium",
                                                   level="Complete Genome", isolate="MRE")]]


DNS_FAIL = set()  # hosts whose names "cannot be looked up"


def fake_get_json(url, params=None):
    params = params or {}
    host = urllib.parse.urlsplit(url).hostname
    if host in DNS_FAIL:  # what the real get_json raises when the lookup fails (Errno 8 on macOS)
        err = urllib.error.URLError(socket.gaierror(8, "nodename nor servname provided, or not known"))
        raise RuntimeError(f"{err} for {url}") from err
    if url.endswith("/version"):
        return {"version": "18.41.1", "major_ver": 18, "minor_ver": 41, "patch_ver": 1}
    if url.endswith("einfo.fcgi"):
        return {"einforesult": {"dblist": ["assembly", "nuccore"]}}
    if params.get("db") == "assembly" and url.endswith("esearch.fcgi"):
        return {"esearchresult": {"idlist": list(ASSEMBLY_SEARCH)}}
    if params.get("db") == "assembly" and url.endswith("esummary.fcgi"):
        uids = params["id"].split(",")
        return {"result": {"uids": uids, **{u: {"assemblyaccession": ASSEMBLY_SEARCH[u],
                                                 "synonym": {"genbank": ASSEMBLY_SEARCH[u]}} for u in uids}}}
    if url.endswith("esearch.fcgi"):
        term = params["term"]
        plant = term.split('"')[1]
        org = "chloroplast" if "chloroplast[filter]" in term else "mitochondrion"
        refseq = "refseq[filter]" in term
        ids = [u for u, _, _ in ORG.get((plant, org, refseq), [])]
        return {"esearchresult": {"idlist": ids}}
    if url.endswith("esummary.fcgi"):
        uids = params["id"].split(",")
        res = {"uids": uids}
        for u in uids:
            acc, slen, (plant, org, _) = SUMMARY[u]
            res[u] = {"accessionversion": acc, "title": f"{plant} {org}, complete genome",
                      "slen": slen, "genome": org, "completeness": "complete"}
        return {"result": res}
    parts = url.split("/genome/")[1].split("/")
    kind, value = parts[0], urllib.parse.unquote(parts[1])
    if (kind, value) in FAIL_ON:
        raise RuntimeError(f"HTTP 500 for {url}")
    pages = DATASETS[(kind, value)]
    i = 1 if params.get("page_token") == "p2" else 0
    out = {"reports": pages[i]}
    if i + 1 < len(pages):
        out["next_page_token"] = "p2"
    return out


real_get_json = M.get_json
M.get_json = fake_get_json


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


with tempfile.TemporaryDirectory() as tmp:
    code = M.main(["-o", tmp])
    if code != 0:
        fails.append(f"exit code {code}, expected 0")
    endo = read(Path(tmp) / "endosymbionts.tsv")
    accs = [r["accession"] for r in endo]
    if sorted(accs) != sorted(["GCA_021409705.1", "GCA_021409605.1", "GCA_000000001.1", "GCA_003966915.1",
                               "GCA_014218255.1", "GCA_003546845.1", "GCA_044000001.1",
                               "GCA_031675535.1"]):
        fails.append(f"endosymbionts: got {accs}")
    held = sorted(r["accession"] for r in endo if r["note"] == "held_out_in_T2")
    if held != sorted(["GCA_021409705.1", "GCA_021409605.1", "GCA_003546845.1", "GCA_044000001.1",
                       "GCA_031675535.1"]):
        fails.append(f"held_out_in_T2: got {held}")
    if "GCA_021464345.1" in accs:
        fails.append("a host fungal assembly from PRJNA670049 leaked into endosymbionts")

    cont = read(Path(tmp) / "contaminants.tsv")
    if len(cont) != 7 or "GCF_000006765.1" not in [r["accession"] for r in cont]:
        fails.append(f"contaminants: got {[r['accession'] for r in cont]}")

    orgs = read(Path(tmp) / "organelles.tsv")
    if len(orgs) != 10:
        fails.append(f"organelles: expected 10 rows, got {len(orgs)}")
    tr = [r for r in orgs if r["plant"] == "Trifolium repens"]
    if len(tr) != 4 or any(r["source"] != "GenBank" for r in tr):
        fails.append(f"Trifolium repens plastids should come from GenBank (4): got {tr}")
    queries = (Path(tmp) / "QUERIES.md").read_text()
    if "Trifolium repens mitochondrion: none in RefSeq or GenBank" not in queries:
        fails.append("QUERIES.md does not list the missing T. repens mitochondrion")
    if "NCBI Datasets API 18.41.1" not in queries:
        fails.append("QUERIES.md does not record the NCBI Datasets version")

    hosts = {r["note"].split(";")[0]: r["accession"] for r in read(Path(tmp) / "hosts.tsv")}
    want = {"host: Gigaspora margarita": "GCA_100000002.1",       # chromosome beats scaffold
            "host: Rhizophagus irregularis": "GCA_200000001.1",    # reference genome wins
            "host: Rhizophagus sp. LPA8": "GCA_021464375.1",       # larger contig N50
            "host: Diversispora epigaea": "GCA_003550325.1"}
    if hosts != want:
        fails.append(f"hosts: got {hosts}")


def run_quietly(argv):
    """main() with stderr captured, so the expected failures below don't clutter the output."""
    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        code = M.main(argv)
    return code, err.getvalue()


# A failing query must give exit code 1, name the query, and write nothing (no short lists).
FAIL_ON.add(("taxon", "Mycoavidus"))
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "refs"
    code, err = run_quietly(["-o", str(out)])
    if code != 1:
        fails.append(f"with a failing query: exit code {code}, expected 1")
    if "taxon: Mycoavidus" not in err:
        fails.append("with a failing query: the error output does not name the query")
    if out.exists():
        fails.append(f"with a failing query: files were written: {sorted(p.name for p in out.iterdir())}")
FAIL_ON.clear()

# If the Datasets host name can't be looked up (Aaron's 2026-11-20 run): exit code 2, a plain
# DNS message, no queries run, nothing written.
DNS_FAIL.add("api.ncbi.nlm.nih.gov")
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "refs"
    code, err = run_quietly(["-o", str(out)])
    if code != 2:
        fails.append(f"with a DNS failure: exit code {code}, expected 2")
    if "could not look up api.ncbi.nlm.nih.gov" not in err or "DNS problem" not in err:
        fails.append(f"with a DNS failure: no plain DNS message, got: {err!r}")
    if "eutils" in err:
        fails.append("with a DNS failure: E-utilities was reported as failing too")
    if out.exists():
        fails.append("with a DNS failure: files were written")
DNS_FAIL.clear()

# The same message through real urllib: a name that can never resolve (.invalid, RFC 2606).
# Proxies are switched off for it, so the lookup happens here, as it would on a laptop.
saved = {k: os.environ.get(k) for k in ("no_proxy", "NO_PROXY")}
os.environ["no_proxy"] = os.environ["NO_PROXY"] = "*"
M.get_json, M.PAUSE, M.DATASETS = real_get_json, 0, "https://ncbi-check.invalid/datasets/v2"
try:
    problems, _ = M.check_ncbi()
finally:
    M.get_json, M.DATASETS = fake_get_json, "https://api.ncbi.nlm.nih.gov/datasets/v2"
    for k, v in saved.items():
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v
if not any("could not look up ncbi-check.invalid" in p for p in problems):
    fails.append(f"real urllib lookup failure did not give the DNS message: {problems}")

if fails:
    print("FAIL\n" + "\n".join(fails))
    sys.exit(1)
print("PASS: endosymbionts (8, 5 held out, incl. MRE cross-check), contaminants (7), organelles (10, GenBank fallback), "
      "hosts (4), pagination, Datasets version recorded, failed query writes nothing (exit 1), "
      "DNS failure gives a plain message and writes nothing (exit 2)")
