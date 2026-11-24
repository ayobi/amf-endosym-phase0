# Reference indexes

Made by make_index.py v0.2 on 2026-11-24 14:40 UTC, with minimap2 2.31-r1302 and sylph 1.0.0.
From the cleaned genomes (`CLEANING.md`). Lineages: PROTOCOL.md 4.1. Settings the protocol leaves open: section 10, entry 4. The organelle genomes, plant and fungal, form one lineage, "organelles" (section 10, entry 10); each set also holds them alone in `organelles.fna`, for the organelle check on routed reads.

| Set | Genomes | Sequences | bp | md5 of refs.fna | md5 of bacteria.syldb |
|---|---|---|---|---|---|
| full | 80 | 27099 | 1359069380 | a2788e50287328b5ff926d49191ebe3e | 42bfd3f5c27a2c3ac81ec941e7456996 |
| t2 | 71 | 26910 | 1349415071 | 6452333ffcd096b16b746c5fa5c4791a | d7aa1182b75b91bd556d4f26142270b9 |

The t2 set leaves out the genomes marked held_out_in_T2: GCA_003546845.1, GCA_021409605.1, GCA_021409645.1, GCA_021409665.1, GCA_021409685.1, GCA_021409705.1, GCA_021409725.1, GCA_021409745.1, GCA_031675535.1.

## Lineages

| Lineage | Kind | Genomes | Sequences | bp | In t2 |
|---|---|---|---|---|---|
| Gigaspora margarita | host fungus | 1 | 6484 | 773081230 | yes |
| Rhizophagus irregularis | host fungus | 1 | 32 | 146773001 | yes |
| Rhizophagus sp. LPA8_1_AMF | host fungus | 1 | 15653 | 116967446 | yes |
| Diversispora epigaea | host fungus | 1 | 731 | 156576535 | yes |
| Linnemannia elongata | host fungus | 1 | 33 | 51153770 | yes |
| CaMg/MRE | endosymbiont | 9 | 189 | 9654309 | no |
| CaGg | endosymbiont | 5 | 2725 | 9630070 | yes |
| Mycoavidus | endosymbiont | 19 | 1153 | 47288638 | yes |
| Ralstonia | contaminant | 1 | 52 | 5430218 | yes |
| Bradyrhizobium | contaminant | 1 | 1 | 9677243 | yes |
| Burkholderia | contaminant | 1 | 3 | 8366868 | yes |
| Herbaspirillum | contaminant | 1 | 1 | 5509723 | yes |
| Pseudomonas | contaminant | 1 | 1 | 6264404 | yes |
| Sphingomonas | contaminant | 1 | 4 | 4367078 | yes |
| Stenotrophomonas | contaminant | 1 | 2 | 4646175 | yes |
| organelles | organelle | 35 | 35 | 3682672 | yes |

## Commands (in the work folder, for each set)

```
minimap2 -x sr -t 30 -d index/full/refs.sr.mmi index/full/refs.fna
minimap2 -x map-hifi -t 30 -d index/full/refs.hifi.mmi index/full/refs.fna
sylph sketch -g clean/bacteria/<accession>.fna (40 files) -o index/full/bacteria -t 30
```

Details: `lineages.tsv` (every genome and its lineage). Per-sequence tables are in the work folder (`index/<set>/sequences.tsv`).
