# Reference lists

Made by make_reference_lists.py v0.2 on 2026-11-20 19:21 UTC, against NCBI Datasets API 18.41.1.

| File | Entries |
|---|---|
| endosymbionts.tsv | 33 |
| contaminants.tsv | 7 |
| organelles.tsv | 10 |
| hosts.tsv | 5 |

MRE genomes marked held_out_in_T2: 9

## Queries

| Query | Hits | Kept |
|---|---|---|
| taxon: Candidatus Moeniiplasma (GenBank) | 7 | 7 |
| taxon: Candidatus Glomeribacter (GenBank) | 5 | 5 |
| taxon: Mycoavidus (GenBank) | 19 | 19 |
| bioproject: PRJNA407498 (GenBank) | 1 | 1 |
| bioproject: PRJNA670049 (GenBank) | 14 | 0 |
| bioproject: PRJNA1176334 (GenBank) | 0 | 0 |
| assembly: ("Mollicutes"[Organism] OR "Mycoplasmatales"[Organism]) AND (endosymbiont[All Fields] OR MRE[All Fields] OR Moeniiplasma[All Fields] OR Glomeromycota[All Fields]) | 8 | 1 |
| taxon: Ralstonia pickettii (RefSeq reference genome) | 1 | 1 |
| taxon: Bradyrhizobium japonicum (RefSeq reference genome) | 1 | 1 |
| taxon: Burkholderia cepacia (RefSeq reference genome) | 1 | 1 |
| taxon: Herbaspirillum seropedicae (RefSeq reference genome) | 1 | 1 |
| taxon: Pseudomonas aeruginosa (RefSeq reference genome) | 1 | 1 |
| taxon: Sphingomonas paucimobilis (RefSeq reference genome) | 1 | 1 |
| taxon: Stenotrophomonas maltophilia (RefSeq reference genome) | 1 | 1 |
| nuccore: "Daucus carota"[Organism] AND refseq[filter] AND chloroplast[filter] | 1 | 1 |
| nuccore: "Daucus carota"[Organism] AND refseq[filter] AND mitochondrion[filter] | 1 | 1 |
| nuccore: "Lotus japonicus"[Organism] AND refseq[filter] AND chloroplast[filter] | 1 | 1 |
| nuccore: "Lotus japonicus"[Organism] AND refseq[filter] AND mitochondrion[filter] | 1 | 1 |
| nuccore: "Medicago truncatula"[Organism] AND refseq[filter] AND chloroplast[filter] | 1 | 1 |
| nuccore: "Medicago truncatula"[Organism] AND refseq[filter] AND mitochondrion[filter] | 1 | 1 |
| nuccore: "Trifolium repens"[Organism] AND refseq[filter] AND chloroplast[filter] | 0 | 0 |
| nuccore: "Trifolium repens"[Organism] AND chloroplast[filter] AND "complete genome"[Title] | 4 | 4 |
| nuccore: "Trifolium repens"[Organism] AND refseq[filter] AND mitochondrion[filter] | 0 | 0 |
| nuccore: "Trifolium repens"[Organism] AND mitochondrion[filter] AND "complete genome"[Title] | 0 | 0 |
| taxon: Gigaspora margarita (GenBank) | 2 | 1 |
| taxon: Rhizophagus irregularis (GenBank) | 34 | 1 |
| bioproject: PRJNA670049 (GenBank) | 3 | 1 |
| taxon: Diversispora epigaea (GenBank) | 1 | 1 |
| taxon: Linnemannia elongata (GenBank) | 10 | 1 |

## Not available
- Trifolium repens mitochondrion: none in RefSeq or GenBank
