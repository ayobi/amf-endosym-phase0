# AMF endosymbiont caller: Phase 0

Can we tell real endosymbionts from contaminants in public AMF sequencing data? Phase 0 tests that on known cases before any survey.

| File | What it is |
|---|---|
| `PROTOCOL.md` | The plan: calls, evidence, tests, kill criteria. Locked before any control is run. |
| `controls.tsv` | Known positives, negatives and backgrounds, with run accessions and how each was checked. |
| `overrides.tsv` | Hand fixes to the run labels, one row per fix, each with a reason. |
| `inventory/` | The three SRA Run Selector exports of 2026-11-20 and the labelled inventory made from them (`runs_tagged.tsv`, `inventory_summary.md`). |
| `scripts/tag_runs.py` | Labels SRA runs from metadata and writes the inventory. |
| `scripts/make_reference_lists.py` | Makes the reference genome lists from NCBI (run on the lock date). |
| `references/` | The reference lists, `QUERIES.md` (how they were made) and `DATABASES.md` (database releases). |
| `env/` | Conda environments for the toolchain. |
| `tests/` | Offline tests for both scripts. |

## Make the inventory

1. In NCBI SRA, run each query from `PROTOCOL.md` section 3.1 (A, B and C).
2. For each: "Send results to Run selector", then "Metadata" in the Download column of the Total row. Save the files as `inventory/SraRunTable_A.csv`, `inventory/SraRunTable_B.csv` and `inventory/SraRunTable_C.csv`.
3. Run:

```bash
python scripts/tag_runs.py inventory/SraRunTable_A.csv inventory/SraRunTable_B.csv \
    inventory/SraRunTable_C.csv --overrides overrides.tsv -o inventory
```

The committed `inventory/` holds the exports of 2026-11-20 and the outputs of the final rerun. SRA grows every day, so new exports will differ: the committed files are the Phase 0 run set.

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
    CONDA_CHANNEL_PRIORITY=strict conda env create -f "$f"
done
for e in phase0-reads phase0-markers phase0-assembly; do
    conda env export -n "$e" --no-builds | grep -v '^prefix:' > "env/$e.lock.yml"
done
```

The env files list `nodefaults`, and `CONDA_CHANNEL_PRIORITY=strict` applies bioconda's recommended channel order for this command only, without changing your conda settings. The `.lock.yml` files hold the exact versions; commit them. (The `grep` drops the line with your local install path.) Then fill in `references/DATABASES.md`. `env/phase0-ssu-fallback.yml` is only needed if stage 0a shows that SSUplex misses short-read SSU.

The committed lock files were exported on the CSB server (linux-64, conda 26.3.2) on 2026-11-20. To rebuild the exact environments elsewhere:

```bash
for e in phase0-reads phase0-markers phase0-assembly; do
    CONDA_CHANNEL_PRIORITY=strict conda env create -f "env/$e.lock.yml"
done
```

Two lines describe that machine rather than a tool: `_x86_64-microarch-level=4` (its CPU supports AVX-512) in `phase0-assembly.lock.yml`, and `sysroot_linux-64=2.39` in `phase0-reads.lock.yml`. If conda says one of them cannot be installed on another machine, delete that line and run the command again; the tool versions are what the lock fixes.

## Lock

`PROTOCOL.md` v1 is the locked version, at the tag `protocol-v1-lock.2`, archived on Zenodo under concept DOI 10.5281/zenodo.24406915. (The earlier tag `protocol-v1-lock.1` was made from the v0.5 working folder by mistake; `PROTOCOL.md` section 10, entry 1, says what differs.) The lock was made like this:

```bash
git add -A && git commit -m "Phase 0 protocol v1, lock 2"
git tag -a protocol-v1-lock.2 -m "Phase 0 protocol locked (replaces lock.1)"
git push origin main protocol-v1-lock.2
```

For the Zenodo archive, the repository is switched on in Zenodo's GitHub settings, and each GitHub release made from a tag becomes a new version of the same Zenodo record, with its own DOI.

After the lock, any change to the protocol, the controls, the reference lists or the environments needs an entry in the deviations log (`PROTOCOL.md` section 10). Planned outputs, such as `references/overlaps.tsv` and the T2 removal lists, are simply added.

## Tests

```bash
python tests/check_tags.py
python tests/check_refs.py
```
