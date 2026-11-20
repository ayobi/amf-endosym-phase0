# AMF run inventory

Made by tag_runs.py v0.3 on 2026-11-20.

## Inputs
- SraRunTable_A.csv: 4637 rows, 4637 unique runs; 282 runs found in other inputs but not in this one
- SraRunTable_B.csv: 275 rows, 275 unique runs; 4644 runs found in other inputs but not in this one
- SraRunTable_C.csv: 7 rows, 7 unique runs; 4912 runs found in other inputs but not in this one
- Unique runs after merging: 4919; BioProjects: 148; organisms: 128

## Runs missing from some inputs (first 25 per file)
Expected when the inputs come from different queries; worth a look when they are two exports of the same query.
- not in SraRunTable_A.csv: SRR12881897, SRR12881898, SRR12881899, SRR12881900, SRR12881901, SRR12881902, SRR12881903, SRR15365641, SRR15365642, SRR15365643, SRR15365644, SRR15365645, SRR15365646, SRR15365647, SRR15365648, SRR15365649, SRR15365650, SRR15365651, SRR15365652, SRR15365653, SRR15365654, SRR15365655, SRR15365656, SRR15365657, SRR15365658 ...
- not in SraRunTable_B.csv: DRR017626, DRR017627, DRR017628, DRR017629, DRR017630, DRR017631, DRR023818, DRR023819, DRR023820, DRR023821, DRR023822, DRR023823, DRR023824, DRR023825, DRR023826, DRR023827, DRR023828, DRR023829, DRR023830, DRR023831, DRR023832, DRR023833, DRR023834, DRR023835, DRR023836 ...
- not in SraRunTable_C.csv: DRR017626, DRR017627, DRR017628, DRR017629, DRR017630, DRR017631, DRR023818, DRR023819, DRR023820, DRR023821, DRR023822, DRR023823, DRR023824, DRR023825, DRR023826, DRR023827, DRR023828, DRR023829, DRR023830, DRR023831, DRR023832, DRR023833, DRR023834, DRR023835, DRR023836 ...

## Roles (set in the overrides file)
| role | runs |
|---|---|
| screen | 4907 |
| excluded | 12 |

Everything below counts only runs with role 'screen', unless it says otherwise.

## Runs by data class and platform
| data class | long | other | short | total |
|---|---|---|---|---|
| amplicon | 223 | 0 | 2451 | 2674 |
| dna_metagenome | 0 | 0 | 32 | 32 |
| hic | 0 | 0 | 11 | 11 |
| rad_reduced | 0 | 0 | 605 | 605 |
| rna | 1 | 1 | 725 | 727 |
| single_nucleus_wga | 47 | 0 | 621 | 668 |
| wgs_long | 107 | 0 | 0 | 107 |
| wgs_short | 0 | 0 | 83 | 83 |
| total | 378 | 1 | 4528 | 4907 |

## Screening value
| screen_value | runs | gbp |
|---|---|---|
| not_applicable | 2363 | 146.39 |
| low | 677 | 528.24 |
| uninformative | 668 | 4404.33 |
| check_prep | 617 | 3224.08 |
| relevant_16S | 311 | 16.79 |
| high | 155 | 1608.37 |
| medium | 90 | 670.72 |
| very_low | 26 | 20.03 |

## RNA library prep (from metadata)
| rna_prep | runs |
|---|---|
| unknown | 627 |
| polyA | 72 |
| small_rna | 26 |
| rRNA_depleted | 2 |

## Amplicon targets (guessed from text)
| target | runs |
|---|---|
| unknown | 1861 |
| fungal_marker | 502 |
| bacterial_16S | 311 |

## Top organisms
| organism | runs |
|---|---|
| Rhizophagus irregularis | 1458 |
| Glomeromycotina | 831 |
| uncultured Glomeromycotina | 622 |
| fungus metagenome | 275 |
| Funneliformis mosseae | 212 |
| Rhizophagus irregularis DAOM 181602=DAOM 197198 | 163 |
| Paraglomus sp. | 122 |
| Rhizophagus clarus | 84 |
| Gigaspora margarita | 76 |
| Glomeromycetes | 67 |
| Cetraspora pellucida | 52 |
| Gigaspora rosea | 50 |
| Entrophospora candida | 48 |
| Funneliformis geosporum | 43 |
| Ambispora callosa | 38 |

## BioProjects most worth screening (high / medium / 16S)
| bioproject | runs | gbp | organisms | classes | culture | isolates | pair_hints | flags |
|---|---|---|---|---|---|---|---|---|
| PRJNA1276452 | 189 | 9.93 | fungus metagenome (189) | amplicon (189) | spores (189); field (189) | INVAM Facility (189) |  |  |
| PRJNA943439 | 94 | 28.42 | Rhizophagus irregularis (94) | amplicon (90); wgs_short (2); wgs_long (2) | in_vitro_root_culture (94); field (94) | Rhizophagus irregularis MUCL 43194 (94) |  |  |
| PRJEB33553 | 66 | 214.95 | Rhizophagus irregularis (66) | wgs_long (66) |  | CHC3 (24); CHC2 (12); CHC5 (11) |  |  |
| PRJNA277029 | 23 | 0.0 | Rhizophagus irregularis (23) | wgs_short (23) | spores (23) | A1 (4); A4 (4); A5 (4) |  |  |
|  | 22 | 0.16 | uncultured Glomeromycotina (22) | dna_metagenome (22) |  |  |  |  |
| PRJNA754762 | 19 | 249.87 | Rhizophagus irregularis (19) | wgs_long (13); hic (5); wgs_short (1) | in_vitro_root_culture (19) | C2 (11); 4401 (3); B3 (2) |  |  |
| PRJNA287285 | 17 | 59.25 | Rhizophagus irregularis (17) | wgs_short (17) | in_planta (4); field (3) | A4 (4); C2 (4); A1 (3) |  |  |
| PRJNA1147075 | 17 | 0.58 | fungus metagenome (17) | amplicon (17) |  |  |  |  |
| PRJDB4945 | 13 | 73.98 | Rhizophagus irregularis DAOM 181602=DAOM 197198 (13) | wgs_long (11); wgs_short (2) |  | DAOM 181602 (13) |  |  |
| PRJNA922099 | 12 | 298.46 | Rhizophagus irregularis (12) | wgs_long (8); hic (4) | in_vitro_root_culture (12) | G1 (3); A4 (3); A5 (3) |  |  |
| PRJNA693419 | 8 | 0.65 | Rhizophagus irregularis (3); Rhizophagus prolifer (2); Septoglomus deserticola (1) | amplicon (8) | in_vitro_root_culture (8); spores (8); field (8); in_planta (2) | root organ culture of CMCC/AM-2901 (1); root organ culture of CMCC/AM-1503 (1); root organ culture of CMCC ROC 38 (1) |  |  |
| PRJNA670049 | 7 | 54.68 | Candidatus Moeniiplasma glomeromycotorum (7) | dna_metagenome (7) | in_vitro_root_culture (7) | LPA8_4 (1); LPA8_3 (1); LPA8_1 (1) |  |  |
| PRJNA407335 | 6 | 254.25 | Diversispora versiformis (6) | wgs_short (6) | in_planta (6); spores (6) | S1 (6) |  |  |
| PRJNA477348 | 6 | 92.91 | Glomus cerebriforme (2); Oehlia diaphana (2); Rhizophagus irregularis (2) | wgs_short (6) | spores (6); in_planta (2) | JS1 (2); SD1 (2); DAOM 240409 (2) |  |  |
| PRJNA230015 | 6 | 42.55 | Rhizophagus irregularis DAOM 197198w (6) | wgs_short (6) |  | DAOM197198 (6) |  |  |
| PRJNA1208279 | 6 | 0.24 | fungus metagenome (6) | amplicon (6) | spores (6) |  |  |  |
| PRJNA174749 | 5 | 1.31 | Rhizophagus irregularis DAOM 181602=DAOM 197198 (2); Rhizophagus cf. irregularis MUCL 46240 (2); Oehlia diaphana (1) | wgs_short (5) |  | DAOM 181602 (2); MUCL 46240 (2) |  |  |
| PRJNA357003 | 4 | 161.31 | Rhizophagus irregularis (4) | wgs_short (4) | in_vitro_root_culture (4); field (4) | multiple set 4 (1); multiple set 3 (1); multiple set 1 (1) |  |  |
| PRJNA1364746 | 3 | 184.39 | Gigaspora margarita (3) | hic (2); wgs_long (1) |  | BEG 34 (3) | Bminus (1); Bplus (1) |  |
| PRJNA575165 | 3 | 101.09 | Gigaspora margarita (3) | wgs_short (3) | in_planta (3); spores (3) | BEG34 (3) | B+ (3) |  |
| PRJNA885267 | 3 | 12.63 | Rhizophagus irregularis (3) | wgs_long (3) | spores (3) | DAOM-197198 (3) |  |  |
| PRJNA1285418 | 3 | 0.14 | Funneliformis mosseae (3) | dna_metagenome (3) |  | fungi (3) |  |  |
| PRJNA268659 | 2 | 57.38 | Rhizophagus irregularis (2) | wgs_short (2) | in_vitro_root_culture (2); field (2) | multiple (2) |  |  |
| PRJDB6444 | 2 | 51.98 | Rhizophagus clarus (2) | wgs_short (2) | field (2) | HR1 (2) |  |  |
| PRJNA368420 | 2 | 46.95 | fungus metagenome (2) | rna (2) | in_planta (2); field (2) |  |  |  |
| PRJNA368421 | 2 | 40.75 | fungus metagenome (2) | rna (2) | in_planta (2); field (2) |  |  |  |
| PRJNA368422 | 2 | 39.45 | fungus metagenome (2) | rna (2) | in_planta (2); field (2) |  |  |  |
| PRJNA747641 | 2 | 39.33 | Rhizophagus irregularis (2) | wgs_short (1); wgs_long (1) | in_vitro_root_culture (2); field (2) | C3 (2) |  |  |
| PRJNA368395 | 1 | 31.93 | fungus metagenome (1) | rna (1) | in_planta (1); field (1) |  |  |  |
| PRJNA535324 | 1 | 23.77 | fungus metagenome (1) | rna (1) |  |  |  |  |
| PRJNA535322 | 1 | 22.58 | fungus metagenome (1) | rna (1) |  |  |  |  |
| PRJNA535319 | 1 | 21.48 | fungus metagenome (1) | rna (1) |  |  |  |  |
| PRJDB8680 | 1 | 19.72 | Rhizophagus clarus (1) | wgs_long (1) | field (1) | HR1 (1) |  |  |
| PRJNA1097788 | 1 | 18.04 | Rhizophagus irregularis (1) | rna (1) |  | fungi (1) |  |  |
| PRJNA610605 | 1 | 15.73 | Geosiphon pyriformis (1) | wgs_short (1) | field (1) | uo1 (1) |  |  |
| PRJNA1301529 | 1 | 11.11 | Rhizophagus irregularis (1) | rna (1) |  | fungi (1) |  |  |
| PRJNA767608 | 1 | 6.99 | Paraglomus occultum (1) | wgs_short (1) |  | Uo1 (1) |  |  |
| PRJNA722321 | 1 | 6.8 | Rhizophagus irregularis (1) | wgs_long (1) | spores (1) | DAOM197198 (1) |  |  |
| PRJEB19799 | 1 | 0.07 | Glomus sp. acp3 (1) | wgs_short (1) |  |  |  |  |
| PRJNA1479122 | 1 | 0.07 | Glomeromycotina (1) | amplicon (1) | in_planta (1); field (1) | 1-127-341F_S58 (1) |  |  |

## Possible with/without-endosymbiont designs (text mentions B+/B-, cured, WT)
| bioproject | runs | organisms | classes | hints |
|---|---|---|---|---|
| PRJNA267628 | 18 | Gigaspora margarita (18) | rna (18) | B+ (11); wild type (11); cured (7); B- (7) |
| PRJNA575165 | 3 | Gigaspora margarita (3) | wgs_short (3) | B+ (3) |
| PRJNA1364746 | 2 | Gigaspora margarita (2) | hic (2) | Bminus (1); Bplus (1) |

## Review flags (all roles; labels are not changed by flags)
| bioproject | runs | gbp | organisms | classes | flags | role | overridden |
|---|---|---|---|---|---|---|---|
| PRJEB45340 | 431 | 2961.6 | Entrophospora candida (48); Cetraspora pellucida (41); Ambispora gerdemannii (24) | single_nucleus_wga (431) | n_libraries_one_biosample | screen (431) | yes |
| PRJEB55531 | 41 | 269.21 | Funneliformis geosporum (41) | single_nucleus_wga (41) | n_libraries_one_biosample | screen (41) | yes |
| PRJNA1036787 | 10 | 0.4 | Glomus (10) | amplicon (10) | metadata mentions human, patient, gut or mouse | excluded (10) | yes |
| PRJNA1284771 | 1 | 7.23 | Rhizophagus irregularis (1) | dna_metagenome (1) | metadata mentions human, patient, gut or mouse | excluded (1) | yes |
| PRJNA754769 | 1 | 0.64 | Oehlia diaphana (1) | wgs_short (1) | metadata mentions human, patient, gut or mouse | excluded (1) | yes |

A flagged BioProject with no override still needs a decision: add a row to the overrides file, or leave it and say why.

## Overrides applied
| accession | only_class | field | value | runs | note |
|---|---|---|---|---|---|
| PRJEB45340 | - | data_class | single_nucleus_wga | 438 | Uppsala single-nucleus genome project, 19 AMF species: about 24 libraries per BioSample, one per sorted nucleus. 438 runs, 3.1 Tbp. Was labelled wgs_short/high. |
| PRJEB55531 | - | data_class | single_nucleus_wga | 41 | Uppsala, Funneliformis geosporum: confirmed single-nucleus. Sahraei et al. 2022 Mycorrhiza 32:361-371 assembled two spores from 19 and 22 individually amplified nuclei, which matches the 41 runs on 2 BioSamples. |
| PRJNA1276452 | amplicon | amplicon_target | bacterial_16S | 189 | AMF Sporosphere Bacteria (William & Mary): 515f-806r 16S amplicons from spore DNA, per the SRA experiment design text (checked 2026-11-20). |
| PRJNA754769 | - | role | excluded | 1 | Stanford School of Medicine sample CHB0028 filed as Oehlia diaphana, with isolate 'pediatric patient'. Not an AMF culture. |
| PRJNA1284771 | - | role | excluded | 1 | Hunan Agricultural University sample UH2-Pardosa pseudoannulata (a wolf spider) filed as Rhizophagus irregularis, with isolate 'gut'. Not AMF data. |
| PRJNA1203428 | - | data_class | rna | 16 | SRA title 'RNA-Seq of bacterium', design 'Transcriptome sequencing of bacteria' (UCLouvain, AMF-bacteria interaction study, R. irregularis MUCL 43194; checked 2026-11-20). Was labelled single_nucleus_wga from free text. |
| PRJNA1147075 | amplicon | amplicon_target | bacterial_16S | 17 | BioProject title: 16S metabarcoding of bacterial communities in a rock phosphate experiment on leek (UM6P; checked 2026-11-20). Root and soil bacteria, not spore-internal. |
| PRJNA752715 | amplicon | amplicon_target | fungal_marker | 50 | BioProject title: LSU Acaulosporaceae spore sequences (University of Kansas; checked 2026-11-20). |
| PRJNA1208279 | amplicon | amplicon_target | bacterial_16S | 6 | BioProject 'Spore-associated bacteria of arbuscular mycorrhizal fungi' (University of Pisa): Racocetra persica and R. fulgida spores, DNeasy PowerSoil Pro, MiSeq (checked 2026-11-20). v0.2 had called these fungal because of the word AMF. |
| PRJNA1036787 | - | role | excluded | 10 | South China Agricultural University: every run is titled 'Fecal microorganism' (about 40 Mbp each); 10 of its 32 runs are filed under the organism Glomus, which is why query A caught them (checked 2026-11-20). Not AMF data. |

## Caveats
- Labels come from metadata and free text only. 'check_prep' and 'check' mean read the paper, not that the data are unusable.
- 'uninformative' runs cannot show absence: no endosymbiont call should be made from them.
- AMF data can be filed under 'fungus metagenome', plant taxa, or the endobacterium's own name; an organism=Glomeromycota query misses all three. See PROTOCOL.md section 3.1 for the extra queries.
