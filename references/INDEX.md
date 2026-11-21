# Reference indexes

Made by make_index.py v0.1 on 2026-11-21 00:37 UTC, with minimap2 2.31-r1302 and sylph 1.0.0.
From the cleaned genomes (`CLEANING.md`). Lineages: PROTOCOL.md 4.1. Settings the protocol leaves open: section 10, entry 4.

| Set | Genomes | Sequences | bp | md5 of refs.fna | md5 of bacteria.syldb |
|---|---|---|---|---|---|
| full | 55 | 27104 | 1357295612 | ecbdc1bd4f9f45a45af3d49dd6585468 | 4670b2605d8d0242ea79f5673bd747df |
| t2 | 46 | 26915 | 1347641303 | 27c4b2f6b53f5b11c97763ccc6c9a1ed | dbc6b819caa31c889a978047338e89e2 |

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
| CaGg | endosymbiont | 5 | 2755 | 9681984 | yes |
| Mycoavidus | endosymbiont | 19 | 1153 | 47288638 | yes |
| Ralstonia | contaminant | 1 | 52 | 5430218 | yes |
| Bradyrhizobium | contaminant | 1 | 1 | 9677243 | yes |
| Burkholderia | contaminant | 1 | 3 | 8366868 | yes |
| Herbaspirillum | contaminant | 1 | 1 | 5509723 | yes |
| Pseudomonas | contaminant | 1 | 1 | 6264404 | yes |
| Sphingomonas | contaminant | 1 | 4 | 4367078 | yes |
| Stenotrophomonas | contaminant | 1 | 2 | 4646175 | yes |
| plant organelles | organelle | 10 | 10 | 1856990 | yes |

## Commands (in the work folder, for each set)

```
minimap2 -x sr -t 16 -d index/full/refs.sr.mmi index/full/refs.fna
minimap2 -x map-hifi -t 16 -d index/full/refs.hifi.mmi index/full/refs.fna
sylph sketch -g clean/bacteria/<accession>.fna (40 files) -o index/full/bacteria -t 16
```

Details: `lineages.tsv` (every genome and its lineage). Per-sequence tables are in the work folder (`index/<set>/sequences.tsv`).
