# AMF endosymbiont caller: Phase 0

Can we tell real endosymbionts from contaminants in public AMF sequencing data? Phase 0 tests that on known cases before any survey.

| File | What it is |
|---|---|
| `PROTOCOL.md` | The plan: calls, evidence, tests, kill criteria. Locked before any control is run. |
| `controls.tsv` | Known positives, negatives and backgrounds, with run accessions and how each was checked. |
| `overrides.tsv` | Hand fixes to the run labels, one row per fix, each with a reason. |
| `scripts/tag_runs.py` | Labels SRA runs from metadata and writes the inventory. |
| `scripts/make_reference_lists.py` | Makes the reference genome lists from NCBI (run on the lock date). |
| `references/` | The reference lists, `QUERIES.md` (how they were made) and `DATABASES.md` (database releases). |
| `env/` | Conda environments for the toolchain. |
| `tests/` | Offline tests for both scripts. |

## Make the inventory

1. In NCBI SRA, run each query from `PROTOCOL.md` section 3.1 (A, B and C).
2. For each: "Send results to Run selector", then "Metadata" in the Download column of the Total row. Rename the files `SraRunTable_A.csv`, `SraRunTable_B.csv`, `SraRunTable_C.csv`.
3. Run:

```bash
python scripts/tag_runs.py SraRunTable_A.csv SraRunTable_B.csv SraRunTable_C.csv \
    --overrides overrides.tsv -o inventory
```

## Make the reference lists (lock date)

Standard-library Python; needs internet access to NCBI. Setting `NCBI_API_KEY` makes it faster but is optional.

```bash
python scripts/make_reference_lists.py -o references
```

It prints how many entries each list has and writes `references/QUERIES.md`, which records the date and the NCBI Datasets version. If NCBI can't be reached or any query fails, it writes nothing and exits non-zero, so a bad run never leaves short lists behind. Just run it again. If it says your computer could not look up a host name, the problem is your network's DNS, not NCBI: switch networks, or set your DNS servers to 1.1.1.1 and 8.8.8.8.

## Software environments

Create them on the Linux machine that will run Phase 0, not on a laptop: lock files are specific to one platform, and some of the tools (GTDB-Tk among them) don't install on macOS.

```bash
for f in env/phase0-reads.yml env/phase0-markers.yml env/phase0-assembly.yml; do
    conda env create -f "$f"
done
for e in phase0-reads phase0-markers phase0-assembly; do
    conda env export -n "$e" --no-builds | grep -v '^prefix:' > "env/$e.lock.yml"
done
```

The `.lock.yml` files hold the exact versions; commit them. (The `grep` drops the line with your local install path.) Then fill in `references/DATABASES.md`. `env/phase0-ssu-fallback.yml` is only needed if stage 0a shows that SSUplex misses short-read SSU.

## Lock

Once every box in `PROTOCOL.md` section 9 is ticked except the last. If the folder is not a git repository yet, first run `git init -b main`, make an empty repository on GitHub, and `git remote add origin <its URL>`. Then:

```bash
git add -A && git commit -m "Phase 0 protocol v1"
git tag -a protocol-v1-lock.1 -m "Phase 0 protocol locked"
git push origin main protocol-v1-lock.1
```

With the repository enabled in Zenodo's GitHub settings, a GitHub release made from the tag gets archived with a DOI.

## Tests

```bash
python tests/check_tags.py
python tests/check_refs.py
```
