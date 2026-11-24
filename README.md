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
| `scripts/clean_references.py` | Stage 0: downloads the reference genomes and cleans them (`PROTOCOL.md` 3.4; from v0.3 also against the organelle genomes, entry 10). |
| `scripts/make_index.py` | Stage 0: builds the lineage table, the minimap2 indexes and the sylph database from the cleaned genomes (full and T2). |
| `scripts/run_stage0a.py` | Stage 0a, step 1: fetches the first 20 million spots of each positive run, assigns them and works out the reference titer t (from v0.2 also again, with a new index). |
| `scripts/get_marker_databases.py` | Stage 0a, step 2: downloads and checks the SingleM metapackage and SSUplex's profile files. |
| `scripts/run_stage0a_m.py` | Stage 0a, step 2: runs SSUplex and SingleM on the step 1 reads and works out the marker share m. |
| `scripts/get_silva.py` | Downloads and checks SILVA 144 SSU Ref NR99 for SSU placement. |
| `scripts/rescore_stage0a_m.py` | Stage 0a, step 2 rescored under the route rules of entries 7 and 8: places the symbiont's SSU reads with VSEARCH and works out m again. |
| `scripts/rescore_m_entry10.py` | Stage 0a, step 2 under entry 10: aligns step 2's found reads again with the entry 10 index, applies the route rules of entries 7, 8 and 10, and works out m again. |
| `scripts/assign_reads.py` | Counts the reads of one run per lineage (`PROTOCOL.md` 4.1), writes the reads that go on to SingleM and SSUplex, and (from v0.2) how the mapped reads spread over each bacterial genome. |
| `scripts/call_lineages.py` | The caller: route-read labels and lineages, levels, cross-talk and the spread check, evidence and calls (`PROTOCOL.md` 4.1, entries 7 to 10). |
| `scripts/run_t1.py` | Test T1: fetches and assigns the runs, runs the route, places the SSU reads, runs sylph, makes the calls and checks K1 (v0.2: again under entry 10). |
| `references/` | The reference lists (with the fungal mitochondria of entry 10 in `fungal_mitochondria.tsv`), `QUERIES.md` (how they were made) and `DATABASES.md` (database releases). |
| `env/` | Conda environments for the toolchain. |
| `stage0a/` | Stage 0a: the positive runs (`runs.tsv`), the fractions and t (step 1), and the marker shares and m (step 2); each measured again under later entries (`_entry7`, `_entry8`, `_entry10`). |
| `t1/` | Test T1: the runs (`runs.tsv`) and, once T1 has run, the levels, the calls and `T1.md`; under entry 10, `levels_entry10.tsv`, `calls_entry10.tsv` and `T1_ENTRY10.md`. |
| `tests/` | Tests for the scripts (`check_clean.py`, `check_assign.py`, `check_stage0a.py` and `check_stage0a_m.py` need minimap2; `check_rescore.py` needs VSEARCH; `check_caller.py` and `check_rescore_entry10.py` need both). |

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

## Stage 0a, step 2 rescored (entries 7 and 8)

Entry 7 (`PROTOCOL.md` section 10) narrows what the route counts as found, and entry 8 fixes how VSEARCH scores a hit (`--iddef 0`). m is worked out again from step 2's work files, without running SSUplex or SingleM again. Version 0.1 of the script ran under entry 7 and wrote `stage0a/marker_share_entry7.tsv` and `stage0a/STAGE0A_M_ENTRY7.md`, which stay in the repository. Version 0.2 runs under entry 8. Commit entry 8 first; SILVA is already in `~/phase0-data/db/silva/` (`scripts/get_silva.py` fetched it for entry 7). Then:

```bash
conda activate phase0-reads
python tests/check_rescore.py
nohup python scripts/rescore_stage0a_m.py --workdir ~/phase0-data --threads 30 > ~/rescore8.log 2>&1 &
tail -f ~/rescore8.log
```

The search against SILVA is exhaustive (section 4.1). Under entry 7 it took about 1.7 hours with 30 threads for the 5,269 distinct SSU sequences. The log prints a time estimate after each batch of 500, and gives m as soon as the DNA runs are done. Work files go to `~/phase0-data/stage0a/rescore_entry8/`, so nothing from version 0.1 is reused. If it stops, run the same command again: finished batches are kept. When it is done it writes `stage0a/marker_share_entry8.tsv` and `stage0a/STAGE0A_M_ENTRY8.md`; commit both.

## Test T1 (entry 9)

T1 needs entry 8's m, so commit `stage0a/marker_share_entry8.tsv` and `stage0a/STAGE0A_M_ENTRY8.md` first, then entry 9 with `t1/runs.tsv`. The 18 B+ runs reuse stage 0a's reads, counts and route; the other 25 runs are fetched (the first 20 million spots; the HiFi run, its first million reads), which takes about 80 to 100 GB of disk in `~/phase0-data/t1/`. Check the free space first:

```bash
conda activate phase0-reads
python tests/check_caller.py
df -h ~/phase0-data
nohup python scripts/run_t1.py --workdir ~/phase0-data --threads 30 > ~/t1.log 2>&1 &
tail -f ~/t1.log
```

It fetches and assigns the new runs, runs SSUplex (and SingleM, for the B− Hi-C run) on the 14 new scored runs, about an hour each, then places every routed SSU read with VSEARCH (the log prints a time estimate after each batch of 500) and runs sylph on all 43 runs. Expect a day or two. The reported runs (HiFi and JGI) are mapped only. If it stops, run the same command again: finished work is kept. When every run is done it writes `t1/levels.tsv`, `t1/calls.tsv` and `t1/T1.md`, with K1 at the top of `T1.md`; commit all three.

## Entry 10: references, t, m and T1 again

Entry 10 (`PROTOCOL.md` section 10) adds the fungal mitochondria to the references, cleans the bacterial genomes against all organelle genomes, and changes the route and the levels. The references and the index are made again, then t, m and T1. Commit entry 10 and the code first. Then, in the `phase0-reads` environment:

```bash
conda activate phase0-reads
python tests/check_clean.py && python tests/check_assign.py && python tests/check_stage0a.py
python tests/check_rescore_entry10.py && python tests/check_caller.py
python scripts/clean_references.py --workdir ~/phase0-data --threads 30
python scripts/make_index.py --workdir ~/phase0-data --threads 30
```

The cleaning downloads only the 25 mitochondrial genomes (the other genomes are already there) and aligns everything again. It replaces `~/phase0-data/clean/` and the index replaces `~/phase0-data/index/`, so the entry 9 versions are gone afterwards; the committed md5s still describe them. Commit `references/overlaps.tsv`, `genome_files.tsv`, `CLEANING.md`, `lineages.tsv` and `INDEX.md`. Then measure t again, on the reads stage 0a already has (about 8 hours):

```bash
nohup python scripts/run_stage0a.py --workdir ~/phase0-data --threads 30 \
    --work-name stage0a_entry10 --suffix _entry10 > ~/stage0a_e10.log 2>&1 &
```

It assigns the 31 runs into `~/phase0-data/stage0a_entry10/assign/` and writes `stage0a/fractions_entry10.tsv` and `stage0a/STAGE0A_ENTRY10.md`; the earlier outputs stay. Commit both. Then m:

```bash
nohup python scripts/rescore_m_entry10.py --workdir ~/phase0-data --threads 30 > ~/rescore10.log 2>&1 &
```

It aligns step 2's found reads again and places only the SSU pieces entry 8 did not place, so it should take well under an hour. It writes `stage0a/marker_share_entry10.tsv` and `stage0a/STAGE0A_M_ENTRY10.md`; commit both. Then T1:

```bash
nohup python scripts/run_t1.py --workdir ~/phase0-data --threads 30 > ~/t1_e10.log 2>&1 &
```

It reuses the reads in `~/phase0-data/t1/reads/` and the outputs of SSUplex and SingleM (linked, not copied). The 18 B+ runs take stage 0a's new counts; the others are assigned again, and sylph runs again on all of them. Expect most of a day. Work files go to `~/phase0-data/t1_entry10/`. It writes `t1/levels_entry10.tsv`, `t1/calls_entry10.tsv` and `t1/T1_ENTRY10.md`; the entry 9 results stay. Commit all three. Each script can be run again if it stops: finished work is kept.

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
python tests/check_caller.py    # needs minimap2 and VSEARCH
python tests/check_rescore_entry10.py  # needs minimap2 and VSEARCH
```
