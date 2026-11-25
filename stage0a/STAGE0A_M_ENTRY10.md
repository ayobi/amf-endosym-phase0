# Stage 0a, step 2 under entry 10: marker share m

Made by rescore_m_entry10.py v0.1 on 2026-11-25 16:29 UTC (call_lineages.py v0.2), from step 2's work files, with minimap2 2.31-r1302 and the entry 10 index (refs.fna md5 a2788e50287328b5ff926d49191ebe3e; `references/INDEX.md`), step 1 again with that index (`fractions_entry10.tsv`), and VSEARCH 2.32.0 against SILVA `SILVA_144_SSURef_NR99_tax_silva_trunc.fasta.gz` (md5 ccead7bba52075f5d8ba2bfab5adc4c1).

Rule: PROTOCOL.md 4.1 and 4.2, with section 10, entries 6, 7, 8 and 10. Step 2's found reads were aligned again with the entry 10 index; a read counts for the symbiont if that alignment counts it so. It is found by the route if SingleM places it in Bacteria (DNA runs only), or if SSUplex calls it bacterial and T1 would keep its SSU piece: not outside Bacteria or in an organelle, not in several genera at 94.5% or more, and at least 60 bp unless it matches a genus. A read whose alignments to the organelle genomes alone, at 90% identity or more, cover half of it is organelle DNA and is not found. m is the lowest share among the DNA positives, and it replaces entry 8's m.

**m = 0.001403** (entry 8: 0.001454): the lowest DNA-positive share, in SRR6039993 (D. epigaea S1, genome reads, mate-pair 3 kb (MP2_3K)).

Checks (section 4.2), under the same rules. A tool finds a symbiont when at least 10 of the reads it found count for that symbiont, in every DNA positive (section 10, entry 6):
- SSUplex finds CaGg and MRE SSU: **passed**. Fewest kept symbiont SSU reads in a DNA positive: 90, in SRR10215377.
- SingleM sees MRE: **yes**. Fewest MRE reads SingleM places in Bacteria, in a T2 positive: 534, in SRR12881899.

| Run | Kind | Symbiont | Sets m | Symbiont reads | Found | Organelle DNA | SSU kept | SSU unplaced: several genera / under 60 bp | SingleM kept | Share (entry 8) | Share |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SRR36098601 | DNA | CaGg | yes | 2079811 | 3932 | 54 | 419 | 269 / 28 | 2999 | 0.001753 | 0.001643 |
| SRR10215375 | DNA | CaGg | yes | 379799 | 849 | 22 | 96 | 60 / 3 | 640 | 0.00204 | 0.001938 |
| SRR10215376 | DNA | CaGg | yes | 375821 | 774 | 25 | 92 | 61 / 2 | 565 | 0.001823 | 0.001748 |
| SRR10215377 | DNA | CaGg | yes | 485853 | 1116 | 26 | 90 | 68 / 1 | 894 | 0.002057 | 0.002025 |
| SRR1659836 | RNA | CaGg | no | 65154 | 1853 | 248 | 104 | 241 / 4 | n/a | 0.01053 | 0.001596 |
| SRR1659837 | RNA | CaGg | no | 22666 | 722 | 93 | 34 | 232 / 8 | n/a | 0.009323 | 0.0015 |
| SRR1659838 | RNA | CaGg | no | 30286 | 1289 | 127 | 36 | 257 / 8 | n/a | 0.01176 | 0.001189 |
| SRR1659839 | RNA | CaGg | no | 2980 | 89 | 12 | 6 | 51 / 0 | n/a | 0.002478 | 0.002013 |
| SRR1659840 | RNA | CaGg | no | 9618 | 2061 | 321 | 491 | 232 / 1 | n/a | 0.168 | 0.05105 |
| SRR1659841 | RNA | CaGg | no | 11505 | 2314 | 360 | 521 | 219 / 0 | n/a | 0.1546 | 0.04528 |
| SRR1659842 | RNA | CaGg | no | 23552 | 4761 | 792 | 1025 | 237 / 0 | n/a | 0.1612 | 0.04352 |
| SRR1659846 | RNA | CaGg | no | 12269 | 2583 | 382 | 537 | 254 / 0 | n/a | 0.1678 | 0.04377 |
| SRR1659847 | RNA | CaGg | no | 21217 | 4589 | 800 | 934 | 223 / 1 | n/a | 0.1592 | 0.04402 |
| SRR1659850 | RNA | CaGg | no | 27 | 4 | 1 | 0 | 2 / 0 | n/a | 0.03061 | 0 |
| SRR1659851 | RNA | CaGg | no | 54 | 6 | 0 | 2 | 4 / 0 | n/a | 0.048 | 0.03704 |
| SRR15311774 | RNA | CaGg | no | 322 | 77 | 8 | 17 | 40 / 0 | n/a | 0.08344 | 0.0528 |
| SRR15311775 | RNA | CaGg | no | 168 | 36 | 2 | 9 | 19 / 0 | n/a | 0.1099 | 0.05357 |
| SRR15311776 | RNA | CaGg | no | 1935 | 352 | 44 | 48 | 173 / 2 | n/a | 0.1002 | 0.02481 |
| SRR6039993 | DNA | CaMg/MRE | yes | 3413545 | 5412 | 14 | 391 | 159 / 1 | 4399 | 0.001454 | 0.001403 |
| SRR6039994 | DNA | CaMg/MRE | yes | 1853553 | 10312 | 35 | 2278 | 278 / 3 | 7140 | 0.005252 | 0.005081 |
| SRR6039995 | DNA | CaMg/MRE | yes | 3132450 | 28469 | 0 | 3396 | 58 / 25 | 23814 | 0.008713 | 0.008686 |
| SRR6039996 | DNA | CaMg/MRE | yes | 2995463 | 7167 | 53 | 1666 | 503 / 3 | 4205 | 0.002146 | 0.00196 |
| SRR6039997 | DNA | CaMg/MRE | yes | 2182017 | 11943 | 44 | 2440 | 329 / 7 | 8555 | 0.005213 | 0.005039 |
| SRR6039998 | DNA | CaMg/MRE | yes | 2591996 | 14637 | 38 | 3208 | 420 / 6 | 10144 | 0.00533 | 0.005151 |
| SRR12881897 | DNA | CaMg/MRE | yes | 489398 | 1095 | 0 | 269 | 42 / 0 | 704 | 0.002074 | 0.001988 |
| SRR12881898 | DNA | CaMg/MRE | yes | 503527 | 1035 | 0 | 262 | 43 / 0 | 655 | 0.001907 | 0.001821 |
| SRR12881899 | DNA | CaMg/MRE | yes | 393098 | 837 | 0 | 210 | 42 / 0 | 534 | 0.002 | 0.001893 |
| SRR12881900 | DNA | CaMg/MRE | yes | 972422 | 2035 | 0 | 485 | 38 / 1 | 1413 | 0.001992 | 0.001952 |
| SRR12881901 | DNA | CaMg/MRE | yes | 972034 | 1967 | 0 | 482 | 45 / 1 | 1339 | 0.001921 | 0.001873 |
| SRR12881902 | DNA | CaMg/MRE | yes | 681857 | 1327 | 0 | 314 | 23 / 1 | 929 | 0.001858 | 0.001823 |
| SRR12881903 | DNA | CaMg/MRE | yes | 1708551 | 3391 | 0 | 768 | 55 / 1 | 2421 | 0.001899 | 0.001866 |

Symbiont reads: step 1 with the entry 10 index. Found: step 2's found reads (by SSUplex or SingleM) that the entry 10 index counts for the symbiont. Organelle DNA: those left out by the organelle check. Per run, with the SSU pieces dropped by reason: `marker_share_entry10.tsv`. SingleM is not used in RNA runs (n/a).
RNA runs where the reads entry 8 placed stand for the rest (--rna-sample 500): SRR1659836, SRR1659837, SRR1659838, SRR1659840, SRR1659841, SRR1659842, SRR1659846, SRR1659847. There SSU kept applies the kept share of the placed reads to all of the run's symbiont SSU reads; the other SSU columns count the placed reads. RNA runs are reported only.
