# Phase 0 protocol: calling endosymbionts in public AMF sequencing data

**Status:** LOCKED v1 (2026-11-20), git tag `protocol-v1-lock.3`, archived on Zenodo from the GitHub release (concept DOI 10.5281/zenodo.24406915, which resolves to the latest version). The earlier tags `protocol-v1-lock.1` and `protocol-v1-lock.2` were incomplete; see the deviations log (section 10, entries 1 and 2).
**After locking:** any change goes in the deviations log (section 10), with the date and the reason, before the step it affects runs.

**Changes in v1 (lock):** settings the tools would otherwise leave open are fixed: tool defaults (4.2), SSU placement and marker-gene and SSU lineages (4.1), the spike error model (4.2) and the T2 removal lists, including SILVA (section 5); database releases recorded (`references/DATABASES.md`); environments created and exported (`env/*.lock.yml`).

**Changes in v0.5:** reference lists made (section 3.4: result); reference cleaning added (3.4), because three of the five host assemblies were made from material that carries a symbiont; read assignment added (4.1), because several genomes per lineage sit in one mapping index; T5 noted as a known-reference pair; CaMg-9 (PRJNA1176334) noted as not yet public.

**Changes in v0.4:** call rule (4.1) and toolchain (4.2) agreed; SSU origin sorting with SSUplex added to the toolchain; reference rules made exact and scripted (`scripts/make_reference_lists.py`), including a search that catches MRE genomes filed under loose names; host references limited to the control hosts; final inventory rerun done (12 runs excluded).

**Changes in v0.3:** numbers fixed (section 11); K2 and K3 rewritten; a "lab kit" spike added to T3; E1 coverage measured from the inventory; PRJEB55531 confirmed single-nucleus; *D. epigaea* reads confirmed; overrides for PRJNA1208279 and PRJNA1036787.

## 1. Question

Can a calibrated caller tell resident endosymbionts from contaminants in public AMF sequencing data, and say "unresolved" when the data cannot decide, including for bacterial lineages that are missing from the reference set?

## 2. What gets called

Step 1, per run: **informative** or **uninformative**. A run is uninformative when it could not have detected a symbiont even if one were there (too few reads, single-nucleus/WGA data, small-RNA libraries, fungal-marker amplicons). Uninformative runs never produce "absent".

Step 2, per isolate or line and bacterial lineage (B+ and B- count as different lines):

| Call | Meaning |
|---|---|
| resident | evidence says the lineage lives in the fungus |
| contaminant | evidence says it came from outside (kit, lab, plant, soil, hyphosphere) |
| unresolved | signal is there, but the evidence cannot decide |
| absent | too little signal to matter (under a tenth of what the reference titer would give), in runs with enough power to see it |

## 3. Data

**3.1 Inventory.** Three NCBI SRA queries, each exported with Run Selector "Metadata" on **2026-11-20**, then labelled together with `scripts/tag_runs.py` v0.3 and `overrides.tsv`:

- **A** `txid214504[Organism:exp]` (Glomeromycota and everything under it).
- **B** `"fungus metagenome"[Organism] AND (Gigaspora OR Rhizophagus OR Glomus OR Diversispora OR Funneliformis OR Claroideoglomus OR Entrophospora OR Paraglomus OR Acaulospora OR Scutellospora OR Racocetra OR Dentiscutata OR Cetraspora OR Septoglomus OR Ambispora OR Geosiphon)`. Catches AMF data filed as "fungus metagenome", such as the JGI metatranscriptomes and spore 16S amplicons.
- **C** `"Candidatus Moeniiplasma glomeromycotorum"[Organism] OR "Candidatus Glomeribacter"[Organism]`. Catches AMF data filed under the endobacterium's name: 7 runs, all T2 positives (*R. irregularis* A2 and *Rhizophagus* sp. LPA8, PRJNA670049).

Result: 4,919 unique runs. The three exports and the labelled inventory from the final rerun (`runs_tagged.tsv`, `inventory_summary.md`) are committed in `inventory/` with the lock, so the run set stays fixed as SRA grows. Not covered on purpose: about 54 Mortierellaceae runs that a looser text search adds (not AMF), and host-plant RNA-seq filed under the plant (the T1 *Lotus* controls are fetched by accession).

**3.2 Labels and overrides.** `tag_runs.py` labels each run from metadata and flags what metadata cannot settle. Hand corrections live in `overrides.tsv` (10 rows on 2026-11-20), each with a reason. Every flagged BioProject has a decision.

**3.3 Controls.** `controls.tsv` (28 rows), frozen at lock. "Deferred" and "dropped" rows stay in the file with the reason.

**3.4 References.** Made on the lock date by `scripts/make_reference_lists.py`, which writes `references/*.tsv` and logs every query in `references/QUERIES.md`. Database releases go in `references/DATABASES.md`.
- *Endosymbionts:* all GenBank assemblies under *Candidatus* Moeniiplasma, *Candidatus* Glomeribacter and *Mycoavidus*; the bacterial assemblies in the control BioProjects PRJNA407498, PRJNA670049 and PRJNA1176334; and anything an NCBI Assembly search for Mollicutes or Mycoplasmatales genomes that mention "endosymbiont", "MRE", *Moeniiplasma* or Glomeromycota turns up. That last search exists because MRE genomes are filed under loose names (for example "endosymbiont GvMRE of Glomus versiforme" and "Mycoplasmataceae bacterium"). Every MRE genome is marked `held_out_in_T2`.
- *Kit contaminants:* the RefSeq reference genome of the type species of *Ralstonia*, *Bradyrhizobium*, *Burkholderia*, *Herbaspirillum*, *Pseudomonas*, *Sphingomonas* and *Stenotrophomonas* (genera from Salter et al. 2014).
- *Host plant organelles:* RefSeq plastid and mitochondrial genomes of *Daucus carota*, *Lotus japonicus*, *Medicago truncatula* and *Trifolium repens*; where RefSeq has none, all GenBank complete genomes for that organelle. (On 2026-11-20 RefSeq had no *T. repens* organelle genomes, and no *T. repens* mitochondrial genome existed in GenBank either.)
- *AMF and other fungal hosts:* one assembly per control host (*G. margarita*, *R. irregularis*, *Rhizophagus* sp. LPA8, *D. epigaea*, and *L. elongata* for T5): the NCBI reference genome if there is one, otherwise the highest assembly level, then the largest contig N50, then the latest release. Used to count fungal reads (E2) and to keep fungal reads out of the bacterial calls.

*Result on 2026-11-20* (`make_reference_lists.py` v0.2, NCBI Datasets API 18.41.1, no query errors): 33 endosymbiont genomes (7 CaMg, all from PRJNA670049; 5 CaGg; 19 *Mycoavidus*; GvMRE; and an MRE filed as "Mycoplasmataceae bacterium", GCA_031675535.1), of which the 9 MRE genomes are held out in T2; 7 contaminant genomes; 10 organelle genomes; 5 host assemblies. The CaMg-9 project (PRJNA1176334, Nanjing Agricultural University, registered October 2024) had no public assembly yet, so it adds nothing; its query stays in the script.

*Reference cleaning* (once, before stage 0a). Three of the five host assemblies were made from material that carries a symbiont: *G. margarita* BEG34 from B+ spores (its BioSample, SAMN12878420, says so), *D. epigaea* IT104 from sporocarps of the isolate that carries GvMRE (same lab, both collected in 2013), and *Rhizophagus* sp. LPA8 from the same plate (LPA8_1) as one of its CaMg genomes. Bacterial contigs left in a host assembly would take the symbiont's reads as "fungal". The other direction can happen too: most of the bacterial genomes were binned from fungal material (the CaMg genomes, GvMRE and most of the *Mycoavidus* genomes), and binned genomes can carry host contigs. (GvMRE is also 2.5 Mb with a CheckM contamination of 114% in NCBI's report, so it holds more than one MRE genome; that does not matter here, since all MRE genomes form one lineage.) So:
1. Align every bacterial reference (all of `endosymbionts.tsv`, held-out genomes included, and `contaminants.tsv`) to every host assembly with minimap2 (`-x asm5`). Keep alignments of at least 500 bp at 95% identity or more.
2. A host contig that is at least half covered by these alignments is bacterial: remove it from the host assembly.
3. A bacterial contig that is at least half covered by alignments to the host contigs that remain is fungal: remove it from that bacterial reference.
4. Anything else that is shared stays in both; reads from it count for neither (read assignment, 4.1).
5. Every removed contig and every remaining shared stretch goes in `references/overlaps.tsv`, committed before stage 0a. The cleaned files are what every later step uses (the minimap2 index and the sylph database).

The held-out genomes take part in steps 1 and 2 even for T2. Those steps only ever remove host sequence, so they cannot create an MRE call; leaving MRE contigs inside the *D. epigaea* and LPA8 assemblies would instead hide the MRE.

**3.5 Inputs per run.** Raw reads for all controls. Logan unitigs for the same runs (test T4).

## 4. Evidence

- **E1 Replication.** The lineage is present in at least 2 independent runs of the same isolate or line: different BioProjects or different years. Replicate plates or libraries in one BioProject do not count as independent.
- **E2 Tracking.** Across the runs used for E1, the ratio of lineage reads to fungal reads stays within 10-fold.
- **E3 Activity.** The lineage is present in a total-RNA or metatranscriptome run of the same isolate or line. Poly-A runs can support presence but never count as "absent".
- **E4 Genome traits.** Where reads assemble: reduced genome, lost DNA-repair or biosynthesis genes, genes of fungal origin. Reported; not required for any call.
- **E5 Background.** The lineage is not present at a similar level (within 10-fold) in 2 or more unrelated isolates from the same lab; it is not organelle DNA; and the signal is not cross-talk (section 4.1).
- **E6 Placement.** One of: known endosymbiont clade (CaGg, CaMg/MRE, *Mycoavidus*); kit-contaminant genus; other known genus; novel (no genus-level match).
- **Genetic code.** Before gene calling, check which translation table MRE use. *Mycoplasma* and *Spiroplasma* read UGA as tryptophan; *Acholeplasma* does not. A wrong code silently cuts genes short. (NCBI annotates GvMRE and the complete MRE genome GCA_031675535.1 with table 4, UGA = Trp.)

### 4.1 Call rule (agreed 2026-11-20)

The rule is categorical; there is no score to tune. Its three-way logic (a call, "absent", or "unresolved" when the data cannot decide) follows O'Brien (in review, *IMA Fungus*), "Identified, absent or unresolved: error-controlled identification of fungal ITS2 sequences against a reference database", which applies the same logic to fungal ITS2 identification.

**Per run.**
- *Reference titer t.* Set once, in stage 0a (section 5), as the lowest lineage read fraction among the T1 and T2 positive runs, measured with the full references, divided by 10.
- *Marker share m.* Also set in stage 0a: the share of a known symbiont's reads that the marker-gene and SSU route recovers (reads found that way ÷ reads found by mapping to its reference), with that route run on all reads of the positive runs.
- *Read assignment* (added in v0.5, agreed 2026-11-20). One minimap2 index holds the cleaned host assemblies, the bacterial references and the organelle genomes (`-x sr` for short reads, `-x map-hifi` for HiFi), with secondary alignments kept (`--secondary=yes -N 50`). Each read is counted on its own, mates separately. A read counts for a lineage only if its best alignment has at least 95% identity over at least 80% of the read, and every alignment with the same best score falls in that lineage; otherwise it counts for none. The lineages are: CaGg (*Candidatus* Glomeribacter genomes), CaMg/MRE (every genome marked `held_out_in_T2`), *Mycoavidus*, each contaminant genus, the plant organelles, and each host fungus (whose count is the fungal read count for E2). Judging ties by lineage, not by genome, matters because one lineage has up to 19 genomes in the index: a read that matches several *Mycoavidus* genomes equally is still a *Mycoavidus* read.
- *Marker-gene and SSU lineages* (v1). Reads with no alignment that passes the read-assignment filter go on to SingleM and SSUplex; reads that match a reference, tied ones included, do not, so a symbiont with a reference genome is never counted a second time as a novel lineage. Marker-gene and SSU reads form lineages by the genus they are placed in (GTDB R232 names for SingleM, SILVA 144 names for SSU, matched across the two ignoring "Candidatus" and GTDB suffixes such as "_A"). Reads not placed in a genus form one novel lineage per route, under the lowest rank they are placed in.
- *SSU placement* (v1). Bacterial SSU reads from SSUplex are searched against SILVA SSU Ref NR99 with VSEARCH (`--usearch_global --id 0.75 --maxaccepts 0 --maxrejects 0 --top_hits_only`: an exhaustive search that keeps every equally best hit). A read is placed at the lowest rank all its top hits share. It is a genus-level match if they share a genus and the top identity is at least 94.5%, the genus threshold of Yarza et al. (2014); otherwise it is novel at genus level. A read with no hit at 75% or more is novel with no placement.
- *Power.* For lineages found by mapping to a reference, N = t × total reads. For lineages found only through marker genes or SSU (every novel lineage), N = t × m × total reads. A run is informative for "absent" if N ≥ 10. Poly-A RNA, single-nucleus/WGA, small-RNA and fungal-marker amplicon runs are never informative for "absent".
- *Level.* With n = reads assigned to the lineage: **present** if n ≥ N; **low** if N/10 ≤ n < N; **absent** if n < N/10.
- *Cross-talk.* Within one BioProject, if the lineage's read fraction in a run is below 2% of its highest fraction in another run of that BioProject, the run counts as **low**, never **present**. (Index hopping and barcode bleed sit around this level.)

**Per isolate or line**, pooling its runs, in this order:
1. **absent**: at least one informative run, and every informative run is "absent".
2. **contaminant**: the lineage is present somewhere and E5 fails; or it sits in a kit-contaminant genus and E1 and E2 do not both hold.
3. **resident**: step 2 did not apply, and
   - known endosymbiont clade (E6): present in 2 or more runs of the isolate or line (any BioProjects; RNA runs count); or
   - any other lineage, including novel: E1 and E2 both hold. Novel lineages are reported as "resident (novel)".
4. **unresolved**: everything else with signal.

No informative run and no signal: no call (uninformative). Organelle reads are removed before any of this and never become a lineage.

What this means for the tests: held-out MRE in T2 can reach "unresolved" but not "resident (novel)", because each positive isolate's runs sit in one BioProject. T2 therefore tests safety (never "absent"), not discovery power. Discovery power can only be tested later, on an isolate with runs in two or more BioProjects.

### 4.2 Toolchain (agreed 2026-11-20; environments in `env/`, exact versions in the exported `env/*.lock.yml` files, database releases in `references/DATABASES.md`)

| Step | Tool | Note |
|---|---|---|
| Get reads | SRA Toolkit (`fasterq-dump`) or the SRA cloud copies | Logan unitigs for T4 |
| Known lineages and contaminants | sylph against the locked references | fast; works at low coverage; reported as a cross-check, does not change calls |
| Organelle removal | competitive mapping to the organelle list (before SingleM), and SSUplex origin sorting (below) | SSUplex sorts by profile, so it also catches organelles of plants not on the list |
| Any bacterium, including novel: marker genes | SingleM, on reads that match no reference (4.1) | check that MRE reads are seen (UGA = Trp); if not, MRE rest on SSU |
| Any bacterium, including novel: SSU rRNA | SSUplex (O'Brien 2026, *Bioinformatics Advances*, doi:10.1093/bioadv/vbag318), on reads that match no reference, extracts SSU reads and sorts them into bacterial, archaeal, eukaryotic, mitochondrial and chloroplast; only bacterial SSU goes on, to placement against SILVA with VSEARCH (rule in 4.1) | stage 0a checks that SSUplex finds CaGg and MRE SSU in short-read positives; if not, phyloFlash (`env/phase0-ssu-fallback.yml`) replaces it |
| Reference cleaning | minimap2 (`-x asm5`), once before stage 0a | rule in section 3.4; output `references/overlaps.tsv` |
| Read counts n per lineage | competitive mapping (minimap2) to the cleaned endosymbionts, contaminants, organelles and AMF hosts together | ties judged by lineage (read assignment, 4.1); gives n and the fungal read count for E2 |
| Assembly, only where there is signal | metaSPAdes (short reads), metaFlye (long reads) | |
| Genome checks (E4, E6) | CheckM2, GTDB-Tk | GTDB R232 holds none of the held-out MRE genomes, so T2 uses it unchanged (section 5) |
| Spikes (T3) | InSilicoSeq (`--mode kde` with the built-in error model for the background run's instrument: HiSeq, NextSeq, NovaSeq or MiSeq), seeds 1, 2 and 3 | spiked reads = spike level × the background run's reads; Illumina backgrounds only in Phase 0 |

Every tool runs with the defaults of the version in `env/*.lock.yml` unless this protocol sets a value. A setting this protocol does not fix is chosen before the step that needs it runs, and written in the deviations log with the reason.

## 5. Tests

- **Stage 0a, calibration (not scored).** Run the T1 and T2 positives with the full references, only to measure the reference titer t and the marker share m. Nothing else from stage 0a is used to set rules.
- **T1 Easy pairs.** *Gigaspora margarita* BEG34 with CaGg (B+) against the cured line (B-).
  - DNA: Hi-C B+ (SRR36098601) against Hi-C B- (SRR36103000), same lab and depth (PRJNA1364746). Also the B+ genome reads (PRJNA575165).
  - RNA: spores, B+ against B- (PRJNA267628); *Lotus* roots, B+ against B- (PRJNA751155).
  - Plant-only negative: non-mycorrhizal *Lotus* roots from the same experiment (PRJNA751155).
  - Expected: B+ line resident for CaGg; B- line not resident; plant-only roots resident for nothing.
  - Reported, not scored: the JGI BEG34 metatranscriptomes (activity), and a blind call on the BEG34 HiFi run (SRR36031820), whose line is not stated.
- **T2 Unknown lineages.** Remove every CaMg/MRE genome from every database the caller uses, then rerun the MRE positives: *D. epigaea* (PRJNA407335), *R. irregularis* A2 (4 plates) and *Rhizophagus* sp. LPA8 (3 plates). What comes out:
  - the 9 genomes marked `held_out_in_T2`, including those assembled from these same reads (PRJNA670049) and GvMRE (PRJNA407498), from the minimap2 index and the sylph database;
  - nothing from GTDB-Tk or SingleM data: none of the 9 passes GTDB's quality checks, so GTDB R232 and the SingleM metapackage built from it hold none of them (checked 2026-11-20; see `references/DATABASES.md`);
  - from SILVA, every sequence whose taxonomy contains Mollicutes, Mycoplasmatota, Tenericutes, Mycoplasmatales, Entomoplasmatales, Acholeplasmatales, Anaeroplasmatales or Haloplasmatales, and every sequence whose name mentions *Moeniiplasma*, Mycoplasma-related endobacteria or MRE. MRE 16S sequences are filed under loose names, as the genomes were, so the whole class comes out; that makes T2 harder, never easier.

  The removal lists are made by these rules and committed before any T2 run. Expected: "unresolved" (see section 4.1); never "absent". The same runs with the references kept in are reported too, but they are easy by construction and do not count.
- **T3 Contaminant challenge.**
  - *Single-run spikes.* Simulated reads from each kit-contaminant genome, spiked into informative Illumina negatives (the B- Hi-C run, and the Illumina and Hi-C runs in PRJNA754762 and PRJNA922099) at 0.001%, 0.01%, 0.1% and 1% of reads, 3 replicates each (seeds 1, 2, 3). Expected: never "resident".
  - *Lab kit spike.* One contaminant genome spiked at 0.1% into every Illumina run from the Ottawa projects (PRJNA754762, PRJNA922099, PRJNA277029, PRJNA287285), across all isolates, to mimic a contaminant that follows a lab. This is the case E1 alone cannot catch. Expected: "contaminant" through E5.
  - *Outside the fungus.* The hyphosphere bacteria of MUCL 43194 (PRJNA943439) and the co-cultured bacterium (PRJNA1203428), as they are. Expected: not "resident".
- **T4 Substrate.** Repeat T1 and T2 on Logan unitigs and compare with raw reads.
- **T5 External pair** (reported, not a kill test). *Linnemannia elongata* AD073 wild type against the cured line. The AD073 *Mycoavidus* genome (GCA_031963465.1) is in the references, so T5, like T1, is a pair with a known reference.

## 6. Kill criteria

- **K1.** In T1: the B- line called "resident" for CaGg; the plant-only roots called "resident" for anything; or the B+ line not called "resident", or its Hi-C or genome runs not at least "present". Any of these: stop.
- **K2.** In T2: any of the 3 held-out isolates (*D. epigaea*, A2, LPA8) called "absent" for the held-out MRE, or no lineage reported where the MRE should be (a miss counts as "absent"): drop the unknown-lineage claim and rescope to known lineages. "Unresolved" passes.
- **K3.** In T3: any single-run spike called "resident", any isolate called "resident" in the lab kit spike, or either "outside the fungus" set called "resident": stop and rethink the evidence.

Not kill criteria, but they change the plan:
- If Logan misses more than 10% of control positives found in raw reads, the survey runs on raw reads.

## 7. Out of scope for Phase 0

Viruses (a later module: reuse LucaProt or RdRp profiles, and translate mitoviruses with UGA as tryptophan). A host identity check from rRNA reads, for the survey, where organism labels cannot be trusted (PRJNA1036787 is filed as *Glomus* but is fecal samples): SSUplex's eukaryotic SSU, ITSxRust for ITS, and AMFsplit on long reads that span the rRNA operon. The ddRAD panel test of "absent" against "unresolved" on sparse data (deferred; see `controls.tsv`). Spikes into long-read backgrounds. The full survey of all AMF data. Host-specificity or ecology questions.

## 8. Reporting

For every control: the call, the evidence behind it, the power (N) and the level (n). Confusion tables for T1-T3. Contaminant calls against spike level. Logan against raw reads. All deviations.

## 9. Lock checklist

- [x] Queries A, B and C exported (2026-11-20)
- [x] Final rerun of `tag_runs.py` with all 10 overrides: 4,919 runs, 12 excluded (2026-11-20); exports and outputs committed in `inventory/`
- [x] Every flagged BioProject decided (overrides for PRJEB45340, PRJEB55531, PRJNA754769, PRJNA1284771, PRJNA1036787)
- [x] PRJEB55531 confirmed single-nucleus: Sahraei et al. 2022 assembled two spores from 19 and 22 nuclei, matching the 41 runs
- [x] *D. epigaea* reads confirmed: BioSample SAMN07646446 links PRJNA407335 to the GvMRE BioProject PRJNA407498 (Fei lab, Cornell)
- [x] Every number fixed (section 11)
- [x] Call rule (4.1) and toolchain (4.2) agreed (2026-11-20)
- [x] `scripts/make_reference_lists.py` v0.2 run (2026-11-20, no query errors; result in section 3.4); `references/` goes in the lock commit
- [x] Reference cleaning (3.4) and read assignment (4.1) agreed (2026-11-20)
- [x] Environments created from `env/*.yml` on the CSB server (linux-64, conda 26.3.2, 2026-11-20); `env/*.lock.yml` exported; `references/DATABASES.md` filled; seeds fixed (1, 2, 3)
- [x] Settings the tools would leave open fixed (v1: tool defaults, SSU placement, marker-gene and SSU lineages, spike error model, T2 removal lists)
- [x] Tag `protocol-v1-lock.3` on this commit; Zenodo archive from the GitHub release (replaces `protocol-v1-lock.1` and `protocol-v1-lock.2`; section 10, entries 1 and 2)

Not a lock item: confirm with JGI or the authors that the BEG34 metatranscriptomes are B+ before any write-up (they are reported, not scored).

## 10. Deviations log

1. **2026-11-20, before any data were analysed: lock tag replaced.** `protocol-v1-lock.1` (commit 6e6ad92; Zenodo doi:10.5281/zenodo.24406916, published 20:54 UTC) was made from the v0.5 working folder by mistake, so its `PROTOCOL.md` is the v0.5 draft and its `references/DATABASES.md` is blank. `protocol-v1-lock.2` holds the intended v1 files. It differs from lock.1 in three files, `PROTOCOL.md` (the v1 changes listed at the top, and this entry), `README.md` and `references/DATABASES.md` (filled in), and adds `inventory/` (the SRA exports and the labelled inventory). The other 23 files, including the reference lists, the controls, the overrides, the environment and lock files and the scripts, are byte-identical in the two tags. Both Zenodo archives are versions of the same record.
2. **2026-11-20, before any data were analysed: inventory added.** `protocol-v1-lock.2` (commit 01c6c00; Zenodo doi:10.5281/zenodo.24407042, published 21:03 UTC) was committed before the inventory files had been copied into `inventory/`, so, contrary to entry 1, it does not contain them. `protocol-v1-lock.3` adds them: the three SRA Run Selector exports of 2026-11-20, and the labelled inventory made from them with `tag_runs.py` v0.3 and the 10 overrides (4,919 runs, 12 excluded). It also changes, in `PROTOCOL.md` and `README.md`, only the tag named as the lock and this entry. All other files are as in lock.2.

## 11. Fixed numbers (agreed 2026-11-20)

| Where | Value |
|---|---|
| Detection power | A run is informative for "absent" if N ≥ 10 expected reads at the reference titer t; t = lowest lineage read fraction among the T1/T2 positive runs (stage 0a) ÷ 10; for marker/SSU-only lineages N also includes the marker share m (section 4.1) |
| K2 | Any of the 3 held-out isolates called "absent" |
| K3 | Zero "resident" calls among spikes, in the lab kit spike, and in the "outside the fungus" sets |
| Logan switch | More than 10% of control positives missed |
| Read assignment | Best alignment at least 95% identity over at least 80% of the read; equal-best alignments in more than one lineage count for none (section 4.1) |
| Reference cleaning | Alignments of at least 500 bp at 95% identity or more; a contig at least half covered is removed (section 3.4) |
| SSU placement | Exhaustive VSEARCH search of SILVA 144, hits at 75% identity or more; genus-level match when all top hits share a genus at 94.5% identity or more (Yarza et al. 2014) (section 4.1) |
| Spikes | InSilicoSeq `--mode kde`, built-in error model for the background run's instrument; seeds 1, 2 and 3 (section 4.2) |
| T2 removal from SILVA | All Mollicutes, plus any sequence named as an MRE (section 5) |
| E1 coverage | Measured, not a threshold: about 8 isolates have informative runs in 2+ BioProjects (DAOM197198, A1, A4, A5, B3, C2, BEG34, *R. clarus* HR1); only DAOM197198 and BEG34 have them from 2+ labs. Counts come from metadata isolate fields and may be off by one or two. |
