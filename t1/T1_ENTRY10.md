# Test T1 under entry 10: the paired positive and negative

Made by run_t1.py v0.2 (call_lineages.py v0.2) on 2026-11-25 20:15 UTC, with assign_reads.py v0.2, minimap2 2.31-r1302 and the full index (refs.fna md5 a2788e50287328b5ff926d49191ebe3e); SSUplex, SingleM and HMMER hmmer 3.4, singlem 0.21.4, ssuplex 0.2.0; VSEARCH 2.32.0 against SILVA `SILVA_144_SSURef_NR99_tax_silva_trunc.fasta.gz`; sylph 1.0.0.
Rule: PROTOCOL.md 4.1, section 5 (T1) and section 6 (K1), with section 10, entries 4 to 10. t = 0.0009396 and m = 0.001403, both measured again under entry 10 (`stage0a/STAGE0A_ENTRY10.md`, `stage0a/STAGE0A_M_ENTRY10.md`). The entry 9 results are in `T1.md`.

## K1: passed

None of the K1 conditions holds: the B- line is not called resident for CaGg, the plant-only roots are called resident for nothing, and the B+ line is resident for CaGg, present in its Hi-C and genome runs.

## Against controls.tsv

| Line | Scored | Expected | CaGg | Basis | Lineages called resident |
|---|---|---|---|---|---|
| BEG34 B+ | yes | resident | resident | known endosymbiont clade, present in 7 runs | CaGg, Glomeribacter (route) |
| BEG34 B- | yes | absent | absent | all 1 informative runs absent | none |
| Lotus plant-only roots | yes | absent (any bacterial lineage) | no call | no informative run and no signal | none |
| BEG34 HiFi | no | unknown until checked | unresolved | signal, but no rule above applies | none |
| BEG34 JGI MT | no | resident (expected) | resident | known endosymbiont clade, present in 3 runs | CaGg |
| BEG34 JGI EM | no | resident (expected) | unresolved | signal, but no rule above applies | none |

## Calls other than absent and no call

| Line | Lineage | Category | Call | Basis | Present / low runs | E1 | E2 | E5 |
|---|---|---|---|---|---|---|---|---|
| BEG34 B+ | CaGg | known endosymbiont clade | resident | known endosymbiont clade, present in 7 runs | 7 / 5 | yes | yes | yes |
| BEG34 HiFi | CaGg | known endosymbiont clade | unresolved | signal, but no rule above applies | 0 / 1 | no | no | yes |
| BEG34 JGI MT | CaGg | known endosymbiont clade | resident | known endosymbiont clade, present in 3 runs | 3 / 1 | no | no | yes |
| BEG34 JGI EM | CaGg | known endosymbiont clade | unresolved | signal, but no rule above applies | 0 / 3 | no | no | yes |
| BEG34 JGI EM | Bradyrhizobium | kit-contaminant genus | contaminant | kit-contaminant genus without both E1 and E2 | 0 / 1 | no | no | yes |
| BEG34 JGI EM | Herbaspirillum | kit-contaminant genus | contaminant | kit-contaminant genus without both E1 and E2 | 0 / 1 | no | no | yes |
| BEG34 JGI MT | Pseudomonas | kit-contaminant genus | contaminant | kit-contaminant genus without both E1 and E2 | 2 / 1 | no | no | yes |
| BEG34 JGI EM | Sphingomonas | kit-contaminant genus | contaminant | kit-contaminant genus without both E1 and E2 | 0 / 1 | no | no | yes |
| BEG34 JGI MT | Stenotrophomonas | kit-contaminant genus | contaminant | kit-contaminant genus without both E1 and E2 | 0 / 3 | no | no | yes |
| BEG34 B+ | Aquipseudomonas (route) | other known genus | unresolved | signal, but no rule above applies | 0 / 2 | no | no | yes |
| BEG34 B- | Aquipseudomonas (route) | other known genus | unresolved | signal, but no rule above applies | 0 / 1 | no | no | yes |
| Lotus plant-only roots | Brevibacillus (route) | other known genus | unresolved | signal, but no rule above applies | 0 / 1 | no | no | yes |
| BEG34 B+ | CGIF01 (route) | other known genus | unresolved | signal, but no rule above applies | 0 / 2 | no | no | yes |
| Lotus plant-only roots | Cohnella (route) | other known genus | unresolved | signal, but no rule above applies | 0 / 3 | no | no | yes |
| BEG34 B+ | Fusicatenibacter (route) | other known genus | unresolved | signal, but no rule above applies | 0 / 1 | no | no | yes |
| BEG34 B+ | Glomeribacter (route) | known endosymbiont clade | resident | known endosymbiont clade, present in 4 runs | 4 / 5 | yes | yes | yes |
| BEG34 B+ | Lentzea (route) | other known genus | unresolved | signal, but no rule above applies | 0 / 1 | no | no | yes |
| Lotus plant-only roots | Streptococcus (route) | other known genus | unresolved | signal, but no rule above applies | 0 / 1 | no | no | yes |
| BEG34 B+ | Thauera (route) | other known genus | unresolved | signal, but no rule above applies | 1 / 0 | no | yes | yes |
| BEG34 B+ | novel (SingleM) | novel | unresolved | signal, but no rule above applies | 1 / 3 | no | yes | yes |
| BEG34 B- | novel (SingleM) | novel | unresolved | signal, but no rule above applies | 0 / 1 | no | no | yes |

## CaGg per run

| Run | Line | Kind | Reads | n | N | Level | Informative | Spread | sylph ANI |
|---|---|---|---|---|---|---|---|---|---|
| SRR36098601 | BEG34 B+ | DNA  | 40000000 | 2079811 | 3.758e+04 | present | yes | 0.99 | 99.99 |
| SRR10215375 | BEG34 B+ | DNA  | 40000000 | 379799 | 3.758e+04 | present | yes | 0.99 | 99.99 |
| SRR10215376 | BEG34 B+ | DNA  | 40000000 | 375821 | 3.758e+04 | present | yes | 0.99 | 99.99 |
| SRR10215377 | BEG34 B+ | DNA  | 40000000 | 485853 | 3.758e+04 | present | yes | 0.99 | 99.99 |
| SRR1659836 | BEG34 B+ | RNA poly-A | 32562116 | 65154 | 3.059e+04 | present | no | - | 95.79 |
| SRR1659837 | BEG34 B+ | RNA poly-A | 31334876 | 22666 | 2.944e+04 | low | no | - | - |
| SRR1659838 | BEG34 B+ | RNA poly-A | 35025904 | 30286 | 3.291e+04 | low | no | - | - |
| SRR1659839 | BEG34 B+ | RNA poly-A | 31780558 | 2980 | 2.986e+04 | absent | no | - | - |
| SRR1659840 | BEG34 B+ | RNA poly-A | 20000000 | 9618 | 1.879e+04 | low | no | - | - |
| SRR1659841 | BEG34 B+ | RNA poly-A | 20000000 | 11505 | 1.879e+04 | low | no | - | - |
| SRR1659842 | BEG34 B+ | RNA poly-A | 20000000 | 23552 | 1.879e+04 | present | no | - | - |
| SRR1659846 | BEG34 B+ | RNA poly-A | 16947449 | 12269 | 1.592e+04 | low | no | - | - |
| SRR1659847 | BEG34 B+ | RNA poly-A | 19541952 | 21217 | 1.836e+04 | present | no | - | - |
| SRR1659850 | BEG34 B+ | RNA poly-A | 20000000 | 27 | 1.879e+04 | absent | no | - | - |
| SRR1659851 | BEG34 B+ | RNA poly-A | 20000000 | 54 | 1.879e+04 | absent | no | - | - |
| SRR15311774 | BEG34 B+ | RNA poly-A | 40000000 | 322 | 3.758e+04 | absent | no | - | - |
| SRR15311775 | BEG34 B+ | RNA poly-A | 40000000 | 168 | 3.758e+04 | absent | no | - | - |
| SRR15311776 | BEG34 B+ | RNA poly-A | 40000000 | 1935 | 3.758e+04 | absent | no | - | - |
| SRR36103000 | BEG34 B- | DNA  | 40000000 | 17 | 3.758e+04 | absent | yes | 1.00 | - |
| SRR1659843 | BEG34 B- | RNA poly-A | 20000000 | 7 | 1.879e+04 | absent | no | - | - |
| SRR1659844 | BEG34 B- | RNA poly-A | 20000000 | 2 | 1.879e+04 | absent | no | - | - |
| SRR1659845 | BEG34 B- | RNA poly-A | 20000000 | 22 | 1.879e+04 | absent | no | - | - |
| SRR1659848 | BEG34 B- | RNA poly-A | 20000000 | 31 | 1.879e+04 | absent | no | - | - |
| SRR1659849 | BEG34 B- | RNA poly-A | 20000000 | 4 | 1.879e+04 | absent | no | - | - |
| SRR1659852 | BEG34 B- | RNA poly-A | 20000000 | 5 | 1.879e+04 | absent | no | - | - |
| SRR1659853 | BEG34 B- | RNA poly-A | 20000000 | 9 | 1.879e+04 | absent | no | - | - |
| SRR15311773 | BEG34 B- | RNA poly-A | 40000000 | 1 | 3.758e+04 | absent | no | - | - |
| SRR15311772 | BEG34 B- | RNA poly-A | 40000000 | 3 | 3.758e+04 | absent | no | - | - |
| SRR15311771 | BEG34 B- | RNA poly-A | 40000000 | 67 | 3.758e+04 | absent | no | - | - |
| SRR15311779 | Lotus plant-only roots | RNA poly-A | 40000000 | 0 | 3.758e+04 | absent | no | - | - |
| SRR15311778 | Lotus plant-only roots | RNA poly-A | 40000000 | 5 | 3.758e+04 | absent | no | - | - |
| SRR15311777 | Lotus plant-only roots | RNA poly-A | 40000000 | 2 | 3.758e+04 | absent | no | - | - |
| SRR36031820 | BEG34 HiFi | DNA  | 1000000 | 764 | 939.6 | low | yes | 0.94 | 99.91 |
| SRR6249698 | BEG34 JGI MT | RNA unstated | 40000000 | 7733 | 3.758e+04 | low | no | - | 98.00 |
| SRR6249699 | BEG34 JGI MT | RNA unstated | 40000000 | 45578 | 3.758e+04 | present | no | - | 98.98 |
| SRR6249700 | BEG34 JGI MT | RNA unstated | 40000000 | 49030 | 3.758e+04 | present | no | - | 99.05 |
| SRR6249701 | BEG34 JGI MT | RNA unstated | 40000000 | 46852 | 3.758e+04 | present | no | - | 99.21 |
| SRR9002094 | BEG34 JGI MT | RNA unstated | 40000000 | 2135 | 3.758e+04 | absent | no | - | - |
| SRR9002095 | BEG34 JGI MT | RNA unstated | 40000000 | 2099 | 3.758e+04 | absent | no | - | - |
| SRR9002098 | BEG34 JGI MT | RNA unstated | 40000000 | 2749 | 3.758e+04 | absent | no | - | - |
| SRR9002096 | BEG34 JGI EM | RNA unstated | 40000000 | 36345 | 3.758e+04 | low | no | - | 97.63 |
| SRR9002097 | BEG34 JGI EM | RNA unstated | 40000000 | 27546 | 3.758e+04 | low | no | - | 97.67 |
| SRR9002099 | BEG34 JGI EM | RNA unstated | 40000000 | 20815 | 3.758e+04 | low | no | - | 98.23 |

## Spread check (entry 10)

In DNA runs, a mapped lineage none of whose genomes has a spread ratio of 0.5 or more is at most low: its reads sit in a few places, not across a genome (`levels_entry10.tsv`, column spread).

It lowered no level.

## The route per run

| Run | Routed reads found | SSU kept | SSU dropped: outside Bacteria / organelle | SSU unplaced: several genera / under 60 bp | SingleM kept | Organelle DNA | Split reads: bacterial / host |
|---|---|---|---|---|---|---|---|
| SRR36098601 | 2017 | 176 | 41 / 82 | 243 / 197 | 1237 | 5 | 936 / 7 |
| SRR10215375 | 723 | 47 | 30 / 11 | 87 / 167 | 343 | 3 | 155 / 21 |
| SRR10215376 | 759 | 45 | 26 / 11 | 138 / 194 | 302 | 2 | 131 / 40 |
| SRR10215377 | 199 | 9 | 8 / 6 | 25 / 16 | 123 | 1 | 0 / 0 |
| SRR1659836 | 376 | 46 | 81 / 18 | 135 / 96 | n/a | 4 | 10 / 0 |
| SRR1659837 | 154 | 20 | 27 / 12 | 60 / 35 | n/a | 2 | 2 / 0 |
| SRR1659838 | 1901 | 143 | 95 / 35 | 1475 / 153 | n/a | 17 | 3 / 0 |
| SRR1659839 | 456 | 25 | 19 / 36 | 344 / 32 | n/a | 0 | 1 / 0 |
| SRR1659840 | 304 | 15 | 66 / 58 | 161 / 4 | n/a | 2 | 3 / 0 |
| SRR1659841 | 625 | 33 | 112 / 93 | 383 / 4 | n/a | 4 | 2 / 0 |
| SRR1659842 | 337 | 23 | 61 / 47 | 202 / 4 | n/a | 1 | 0 / 0 |
| SRR1659846 | 701 | 46 | 117 / 92 | 438 / 8 | n/a | 4 | 3 / 0 |
| SRR1659847 | 8752 | 305 | 1543 / 2207 | 4665 / 32 | n/a | 23 | 3 / 0 |
| SRR1659850 | 253 | 26 | 36 / 46 | 143 / 2 | n/a | 3 | 0 / 0 |
| SRR1659851 | 229 | 12 | 41 / 54 | 120 / 2 | n/a | 0 | 0 / 0 |
| SRR15311774 | 2424 | 59 | 315 / 515 | 1527 / 8 | n/a | 2 | 3 / 0 |
| SRR15311775 | 1222 | 51 | 232 / 301 | 632 / 6 | n/a | 1 | 0 / 0 |
| SRR15311776 | 10766 | 714 | 1835 / 1617 | 6475 / 125 | n/a | 17 | 3 / 0 |
| SRR36103000 | 228 | 15 | 11 / 7 | 27 / 47 | 58 | 0 | 0 / 16 |
| SRR1659843 | 432 | 28 | 57 / 59 | 283 / 5 | n/a | 1 | 0 / 0 |
| SRR1659844 | 141 | 4 | 28 / 21 | 86 / 2 | n/a | 0 | 0 / 0 |
| SRR1659845 | 248 | 7 | 53 / 52 | 133 / 3 | n/a | 0 | 0 / 0 |
| SRR1659848 | 881 | 22 | 171 / 249 | 438 / 1 | n/a | 0 | 0 / 0 |
| SRR1659849 | 126 | 7 | 23 / 27 | 66 / 3 | n/a | 0 | 0 / 0 |
| SRR1659852 | 187 | 14 | 26 / 51 | 96 / 0 | n/a | 1 | 0 / 0 |
| SRR1659853 | 200 | 16 | 29 / 39 | 115 / 1 | n/a | 0 | 0 / 0 |
| SRR15311773 | 350 | 16 | 55 / 106 | 169 / 4 | n/a | 2 | 0 / 0 |
| SRR15311772 | 341 | 16 | 46 / 85 | 191 / 3 | n/a | 2 | 0 / 0 |
| SRR15311771 | 1515 | 117 | 237 / 243 | 904 / 14 | n/a | 2 | 0 / 0 |
| SRR15311779 | 317 | 37 | 28 / 134 | 116 / 2 | n/a | 5 | 0 / 0 |
| SRR15311778 | 482 | 27 | 35 / 157 | 258 / 5 | n/a | 1 | 0 / 0 |
| SRR15311777 | 471 | 50 | 25 / 137 | 257 / 2 | n/a | 1 | 0 / 0 |

Reads are counted per read (mates separately). `levels_entry10.tsv` gives every lineage in every run; `calls_entry10.tsv` gives every call with its evidence; novel lineages list their most common placements. SSU columns count the found reads by their SSU label; organelle DNA and split reads are found reads with a label that leave the route. sylph is a cross-check and does not change calls (4.2).
