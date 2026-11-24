# Reference cleaning

Made by clean_references.py v0.3 on 2026-11-24 14:40 UTC, with minimap2 2.31-r1302, from genomes downloaded through NCBI Datasets API 18.41.1.
Rule: PROTOCOL.md section 3.4. Settings the protocol leaves open: section 10, entry 3.

## Commands (in the work folder)

```
minimap2 -x asm5 -c -p 0 -t 30 genomes/GCA_009809945.1.fna bacteria.fna > paf/GCA_009809945.1.paf
minimap2 -x asm5 -c -p 0 -t 30 genomes/GCA_026210795.1.fna bacteria.fna > paf/GCA_026210795.1.paf
minimap2 -x asm5 -c -p 0 -t 30 genomes/GCA_021464375.1.fna bacteria.fna > paf/GCA_021464375.1.paf
minimap2 -x asm5 -c -p 0 -t 30 genomes/GCA_003547095.1.fna bacteria.fna > paf/GCA_003547095.1.paf
minimap2 -x asm5 -c -p 0 -t 30 genomes/GCA_036320915.1.fna bacteria.fna > paf/GCA_036320915.1.paf
minimap2 -x asm20 -c -p 0 -t 30 clean/organelles.fna bacteria.fna > paf/organelles.paf
```

Alignments found: 130; of these, 110 pass step 1 (at least 500 bp at 95% identity or more).

## Host assemblies

| Assembly | Organism | Contigs | Removed as bacterial | Removed bp |
|---|---|---|---|---|
| GCA_009809945.1 | Gigaspora margarita | 6489 | 5 | 17316 |
| GCA_026210795.1 | Rhizophagus irregularis | 32 | 0 | 0 |
| GCA_021464375.1 | Rhizophagus sp. LPA8_1_AMF | 15685 | 32 | 228341 |
| GCA_003547095.1 | Diversispora epigaea | 731 | 0 | 0 |
| GCA_036320915.1 | Linnemannia elongata | 33 | 0 | 0 |

## Bacterial genomes with contigs removed as fungal

| Genome | Organism | Contigs | Removed as fungal | Removed bp |
|---|---|---|---|---|
| GCA_021409605.1 | Candidatus Moeniiplasma glomeromycotorum | 53 | 2 | 8673 |
| GCA_021409645.1 | Candidatus Moeniiplasma glomeromycotorum | 54 | 1 | 4581 |

Shared stretches left in both (reads from them count for neither): 24, 34442 bp in all.

## Organelle DNA in the bacterial genomes (section 10, entry 10)

The bacterial genomes against the 35 organelle genomes (10 plant, from `organelles.tsv`; 25 fungal mitochondria, from `fungal_mitochondria.tsv`), with minimap2 -x asm20 -c -p 0. Alignments found: 158; of these, 96 have at least 500 bp at 80% identity or more. A contig at least half covered by them is organelle DNA and is removed (contigs already removed as fungal are left out of this).

| Genome | Organism | Contigs | Removed as organelle DNA | Removed bp |
|---|---|---|---|---|
| GCA_001684025.1 | Candidatus Glomeribacter gigasporarum | 851 | 29 | 50946 |
| GCA_001684155.1 | Candidatus Glomeribacter gigasporarum | 1172 | 1 | 968 |

Organelle stretches left in the bacterial genomes (reads from them count for no lineage): 3, 1720 bp in all.

The held-out MRE genomes took part in all of this, as section 3.4 says.

Downloads: 45 genomes, 45 of them checked against the md5 in NCBI's own download; 35 organelle genomes, each checked against its length in `organelles.tsv` or `fungal_mitochondria.tsv`.

Details: `overlaps.tsv` (every removed contig, every shared stretch and every organelle stretch; in the organelle rows the host columns hold the organelle genome) and `genome_files.tsv` (every file, with md5).
