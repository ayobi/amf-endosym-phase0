# Databases used in Phase 0

Fixed on 2026-11-20 from the tool versions in `env/*.lock.yml`: each release is the one that version of the tool downloads or requires. Download them in stage 0a, not before.

| Database | Used by | Release / version | Source |
|---|---|---|---|
| GTDB | GTDB-Tk 2.8.0 (E6 placement) | R232 (232.0), the release GTDB-Tk 2.8.0 requires; about 56 GB | https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz |
| SingleM metapackage | SingleM 0.21.4 (marker genes) | `S6.5.0.GTDB_r232.metapackage_20260319` (version 6.5.0, built from GTDB R232; 2.2 GB; md5 `580f9ae9ae6a0972f15665c7397c6369`), the version SingleM 0.21.4 asks for | Zenodo, doi:10.5281/zenodo.20150069 (`singlem data` downloads it) |
| CheckM2 database | CheckM2 1.1.0 (E4) | version 3 (`CheckM2_database/uniref100.KO.1.dmnd`), the highest version CheckM2 1.1.0 accepts | Zenodo, doi:10.5281/zenodo.14897628 (`checkm2 database --download` fetches it) |
| SILVA SSU Ref NR99 | VSEARCH 2.32.0 (SSU placement, section 4.1) | release 144 (14 August 2026; 905,628 sequences) | https://www.arb-silva.de/ |
| Metaxa2 SSU profile database | SSUplex 0.2.0 (SSU origin sorting) | the SSU profiles in Metaxa2 2.2.3 (`metaxa2_db/SSU/HMMs`), assembled into SSUplex's five files with its `dev/assemble_metaxa2_hmms.sh` | https://microbiology.se/sw/Metaxa2_2.2.3.tar.gz |

GTDB-Tk and the SingleM metapackage both use GTDB R232, so genome placement and marker genes share one taxonomy.

None of the 9 MRE genomes held out in T2 passes GTDB's quality checks (the 7 CaMg genomes and GvMRE fail them outright, and GCA_031675535.1 has only 43 of the 120 bacterial marker genes; checked in the GTDB API on 2026-11-20), so GTDB R232 and the SingleM metapackage hold no MRE genome and are used unchanged in T2 (section 5).

The reference genome lists (endosymbionts, contaminants, organelles, hosts) come from `scripts/make_reference_lists.py`; see `QUERIES.md` in this folder for the date and the queries.
