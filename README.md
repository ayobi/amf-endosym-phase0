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
| `scripts/clean_references.py` | Stage 0: downloads the reference genomes and cleans them (`PROTOCOL.md` 3.4). |
| `scripts/make_index.py` | Stage 0: builds the lineage table, the minimap2 indexes and the sylph database from the cleaned genomes (full and T2). |
| `scripts/run_stage0a.py` | Stage 0a, step 1: fetches the first 20 million spots of each positive run, assigns them and works out the reference titer t. |
| `scripts/get_marker_databases.py` | Stage 0a, step 2: downloads and checks the SingleM metapackage and SSUplex's profile files. |
| `scripts/run_stage0a_m.py` | Stage 0a, step 2: runs SSUplex and SingleM on the step 1 reads and works out the marker share m. |
| `scripts/get_silva.py` | Downloads and checks SILVA 144 SSU Ref NR99 for SSU placement. |
| `scripts/rescore_stage0a_m.py` | Stage 0a, step 2 rescored under entry 7's route rules: places the symbiont's SSU reads with VSEARCH and works out m again. |
| `scripts/assign_reads.py` | Counts the reads of one run per lineage (`PROTOCOL.md` 4.1) and writes the reads that go on to SingleM and SSUplex. |
| `references/` | The reference lists, `QUERIES.md` (how they were made) and `DATABASES.md` (database releases). |
| `env/` | Conda environments for the toolchain. |
| `stage0a/` | Stage 0a: the positive runs (`runs.tsv`), the fractions and t (step 1), and the marker shares and m (step 2). |
| `tests/` | Tests for the scripts (`check_clean.py`, `check_assign.py`, `check_stage0a.py` and `check_stage0a_m.py` need minimap2; `check_rescore.py` needs VSEARCH). |

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

## Clean the references (stage 0, before stage 0a)

The settings are in `PROTOCOL.md` section 10, entry 3, which must be committed before this runs. Needs the `phase0-reads` environment (for minimap2), internet access to NCBI, and about 4 GB of disk for the work folder. Keep the work folder outside the repository.

```bash
conda activate phase0-reads
python tests/check_clean.py                       # synthetic genomes, real minimap2
python scripts/clean_references.py --workdir ~/phase0-data --threads 16
```

It writes `references/overlaps.tsv` (every removed contig and every shared stretch), `references/genome_files.tsv` (every downloaded and cleaned file, with md5) and `references/CLEANING.md` (versions, commands, counts). Commit those three before stage 0a. The cleaned genomes stay in `~/phase0-data/clean/`; every later step uses them. If a download fails, nothing is written to `references/`; run the same command again and it fetches only what is missing.

## Build the index (stage 0, after cleaning, before stage 0a)

The settings are in `PROTOCOL.md` section 10, entry 4, which must be committed before this runs. Needs the `phase0-reads` environment (minimap2 and sylph) and about 15 GB more in the work folder.

```bash
conda activate phase0-reads
python tests/check_assign.py                      # synthetic genomes and reads, real minimap2
python scripts/make_index.py --workdir ~/phase0-data --threads 16
```

It writes `references/lineages.tsv` (every genome and its lineage) and `references/INDEX.md` (versions, commands, md5s). Commit both before stage 0a. The indexes go to `~/phase0-data/index/full/` and `~/phase0-data/index/t2/` (the T2 set leaves out the 9 held-out MRE genomes).

Reads are then assigned one run at a time; give both files of a paired run, and each mate is counted on its own:

```bash
python scripts/assign_reads.py --index ~/phase0-data/index/full --platform short \
    --reads SRR123_1.fastq.gz SRR123_2.fastq.gz --out ~/phase0-data/runs/SRR123/full --threads 16
```

This writes `counts.tsv` (reads per lineage, the fungal read count for E2, and the reads counted for none or routed, by reason), `routed.fastq.gz` (the reads that go on to SingleM and SSUplex) and `assign.json`. Use `--platform hifi` for HiFi runs and `--index ~/phase0-data/index/t2` for T2.

## Stage 0a, step 1: reference titer t

The settings are in `PROTOCOL.md` section 10, entry 5, which must be committed before this runs. It fetches the first 20 million spots of each of the 31 runs in `stage0a/runs.tsv` (3 at a time), assigns them against the full index, and works out t. It keeps the reads for step 2, about 80 GB in `~/phase0-data/stage0a/reads/`, and needs internet access to SRA.

```bash
conda activate phase0-reads
python tests/check_stage0a.py
python scripts/run_stage0a.py --workdir ~/phase0-data --threads 16
```

It takes many hours, so run it inside `tmux` or `screen`, or with `nohup` (see below). If it stops, or a run fails to download, run the same command again: finished runs are skipped. When every run is done it writes `stage0a/fractions.tsv` and `stage0a/STAGE0A.md`; commit both before step 2.

```bash
nohup python scripts/run_stage0a.py --workdir ~/phase0-data --threads 16 > ~/stage0a.log 2>&1 &
tail -f ~/stage0a.log      # Ctrl-C stops watching, not the run
```

## Stage 0a, step 2: marker share m

The settings are in `PROTOCOL.md` section 10, entry 6, which must be committed before this runs. First get the databases: the SingleM metapackage (about 2.2 GB, from Zenodo) and Metaxa2 2.2.3 for SSUplex's profile files (from microbiology.se).

```bash
conda activate phase0-reads
python tests/check_stage0a_m.py
python scripts/get_marker_databases.py --workdir ~/phase0-data
```

It checks the tool versions against the lock files, the metapackage against the md5 in `references/DATABASES.md`, and SSUplex's quick start against its README. Only if all of that passes does it write `~/phase0-data/db/manifest.json`, which step 2 needs. Then run step 2 on the reads kept from step 1:

```bash
nohup python scripts/run_stage0a_m.py --workdir ~/phase0-data --threads 16 > ~/stage0a_m.log 2>&1 &
tail -f ~/stage0a_m.log
```

It runs SSUplex and SingleM on every read of the 31 runs, DNA positives first, and aligns the reads they find again, exactly as in step 1. Expect about an hour per run with 16 threads, so a day or two in all. SSUplex holds one mate file in memory (around 10 GB), and a run needs about 10 GB of temporary disk in `~/phase0-data/stage0a/route/`. If it stops, run the same command again: finished work is skipped. When every run is done it writes `stage0a/marker_share.tsv` and `stage0a/STAGE0A_M.md`; commit both.

## Stage 0a, step 2 rescored under entry 7

Entry 7 (`PROTOCOL.md` section 10) narrows what the route counts as found, and m is worked out again from step 2's work files, without running SSUplex or SingleM again. Commit entry 7 first. Then get SILVA (about 310 MB) and run the rescore:

```bash
conda activate phase0-reads
python tests/check_rescore.py
python scripts/get_silva.py --workdir ~/phase0-data
nohup python scripts/rescore_stage0a_m.py --workdir ~/phase0-data --threads 30 > ~/rescore.log 2>&1 &
tail -f ~/rescore.log
```

The search against SILVA is exhaustive (section 4.1), so it is slow: expect a day or more for the roughly 27,000 SSU reads. Use as many threads as the server can spare. The log prints a time estimate after each batch of 500 reads, and gives m as soon as the DNA runs are done. If it stops, run the same command again: finished batches are kept. When it is done it writes `stage0a/marker_share_entry7.tsv` and `stage0a/STAGE0A_M_ENTRY7.md`; commit both.

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

`PROTOCOL.md` v1 is the locked version, at the tag `protocol-v1-lock.3`, archived on Zenodo under concept DOI 10.5281/zenodo.24406915. (The earlier tags `protocol-v1-lock.1` and `protocol-v1-lock.2` were incomplete; `PROTOCOL.md` section 10, entries 1 and 2, says how.) The lock was made like this:

```bash
git add -A && git commit -m "Phase 0 protocol v1, lock 3"
git tag -a protocol-v1-lock.3 -m "Phase 0 protocol locked, with inventory (replaces lock.2)"
git push origin main protocol-v1-lock.3
```

For the Zenodo archive, the repository is switched on in Zenodo's GitHub settings, and each GitHub release made from a tag becomes a new version of the same Zenodo record, with its own DOI.

After the lock, any change to the protocol, the controls, the reference lists or the environments needs an entry in the deviations log (`PROTOCOL.md` section 10). Planned outputs, such as `references/overlaps.tsv` and the T2 removal lists, are simply added.

## Tests

```bash
python tests/check_tags.py
python tests/check_refs.py
python tests/check_clean.py    # needs minimap2: activate phase0-reads first
python tests/check_assign.py   # needs minimap2 too
python tests/check_stage0a.py  # needs minimap2 too
python tests/check_stage0a_m.py  # needs minimap2 too
python tests/check_rescore.py   # needs VSEARCH
```
