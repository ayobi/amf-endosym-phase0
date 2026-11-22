# Stage 0a, step 2 rescored under entry 7: marker share m

Made by rescore_stage0a_m.py v0.1 on 2026-11-22 21:01 UTC, from step 2's work files (run_stage0a_m.py v0.1), with VSEARCH 2.32.0 against SILVA `SILVA_144_SSURef_NR99_tax_silva_trunc.fasta.gz` (md5 ccead7bba52075f5d8ba2bfab5adc4c1, 905,628 sequences).

Rule: PROTOCOL.md 4.1 and 4.2, with section 10, entries 6 and 7. A symbiont read counts as found by the route if SingleM places it in Bacteria (DNA runs only), or if SSUplex calls it bacterial and none of its top SILVA hits lies outside Bacteria or in chloroplasts or mitochondria (reads with no hit at 75% or more stay). The reads and their mapping outcomes are step 2's. m is the lowest share among the DNA positives, and it replaces step 2's m.

**m = 0.001454** (step 2: 0.001585): the lowest DNA-positive share, in SRR6039993 (D. epigaea S1, genome reads, mate-pair 3 kb (MP2_3K)).

Checks (section 4.2), under the same rules. A tool finds a symbiont when at least 10 of the reads it found count for that symbiont, in every DNA positive (section 10, entry 6):
- SSUplex finds CaGg and MRE SSU: **passed**. Fewest kept symbiont SSU reads in a DNA positive: 142, in SRR10215375.
- SingleM sees MRE: **yes**. Fewest MRE reads SingleM places in Bacteria, in a T2 positive: 534, in SRR12881899.

| Run | Kind | Symbiont | Sets m | Symbiont reads | SSU (step 2) | SSU kept | SingleM (step 2) | SingleM kept | Share (step 2) | Share |
|---|---|---|---|---|---|---|---|---|---|---|
| SRR36098601 | DNA | CaGg | yes | 2150922 | 811 | 658 | 3122 | 2999 | 0.001829 | 0.0017 |
| SRR10215375 | DNA | CaGg | yes | 402536 | 196 | 142 | 653 | 640 | 0.002109 | 0.001943 |
| SRR10215376 | DNA | CaGg | yes | 409195 | 205 | 154 | 570 | 565 | 0.001894 | 0.001757 |
| SRR10215377 | DNA | CaGg | yes | 524606 | 196 | 151 | 921 | 894 | 0.002129 | 0.001992 |
| SRR1659836 | RNA | CaGg | no | 104278 | 1848 | 1042 | 5 | n/a | 0.01777 | 0.009993 |
| SRR1659837 | RNA | CaGg | no | 44941 | 722 | 399 | 0 | n/a | 0.01607 | 0.008878 |
| SRR1659838 | RNA | CaGg | no | 64218 | 1289 | 727 | 0 | n/a | 0.02007 | 0.01132 |
| SRR1659839 | RNA | CaGg | no | 26630 | 89 | 55 | 0 | n/a | 0.003342 | 0.002065 |
| SRR1659840 | RNA | CaGg | no | 9986 | 2061 | 1364 | 0 | n/a | 0.2064 | 0.1366 |
| SRR1659841 | RNA | CaGg | no | 11821 | 2314 | 1527 | 0 | n/a | 0.1958 | 0.1292 |
| SRR1659842 | RNA | CaGg | no | 23923 | 4761 | 3247 | 0 | n/a | 0.199 | 0.1357 |
| SRR1659846 | RNA | CaGg | no | 12864 | 2582 | 1745 | 1 | n/a | 0.2008 | 0.1356 |
| SRR1659847 | RNA | CaGg | no | 22484 | 4589 | 2992 | 0 | n/a | 0.2041 | 0.1331 |
| SRR1659850 | RNA | CaGg | no | 98 | 3 | 3 | 1 | n/a | 0.04082 | 0.03061 |
| SRR1659851 | RNA | CaGg | no | 125 | 6 | 5 | 0 | n/a | 0.048 | 0.04 |
| SRR15311774 | RNA | CaGg | no | 755 | 77 | 57 | 0 | n/a | 0.102 | 0.0755 |
| SRR15311775 | RNA | CaGg | no | 273 | 36 | 22 | 0 | n/a | 0.1319 | 0.08059 |
| SRR15311776 | RNA | CaGg | no | 2576 | 352 | 224 | 0 | n/a | 0.1366 | 0.08696 |
| SRR6039993 | DNA | CaMg/MRE | yes | 3413552 | 782 | 564 | 4630 | 4399 | 0.001585 | 0.001454 |
| SRR6039994 | DNA | CaMg/MRE | yes | 1853554 | 2888 | 2579 | 7424 | 7140 | 0.005563 | 0.005243 |
| SRR6039995 | DNA | CaMg/MRE | yes | 3132448 | 3675 | 3449 | 24794 | 23814 | 0.009088 | 0.008703 |
| SRR6039996 | DNA | CaMg/MRE | yes | 2995465 | 2821 | 2186 | 4346 | 4205 | 0.002393 | 0.002134 |
| SRR6039997 | DNA | CaMg/MRE | yes | 2182016 | 3117 | 2793 | 8826 | 8555 | 0.005473 | 0.005201 |
| SRR6039998 | DNA | CaMg/MRE | yes | 2591996 | 4112 | 3653 | 10525 | 10144 | 0.005647 | 0.005323 |
| SRR12881897 | DNA | CaMg/MRE | yes | 489398 | 336 | 310 | 759 | 704 | 0.002237 | 0.002072 |
| SRR12881898 | DNA | CaMg/MRE | yes | 503527 | 328 | 304 | 707 | 655 | 0.002056 | 0.001905 |
| SRR12881899 | DNA | CaMg/MRE | yes | 393098 | 264 | 250 | 573 | 534 | 0.002129 | 0.001994 |
| SRR12881900 | DNA | CaMg/MRE | yes | 972422 | 546 | 524 | 1489 | 1413 | 0.002093 | 0.001992 |
| SRR12881901 | DNA | CaMg/MRE | yes | 972034 | 548 | 528 | 1419 | 1339 | 0.002024 | 0.001921 |
| SRR12881902 | DNA | CaMg/MRE | yes | 681857 | 346 | 338 | 981 | 929 | 0.001946 | 0.001858 |
| SRR12881903 | DNA | CaMg/MRE | yes | 1708551 | 858 | 826 | 2533 | 2421 | 0.001985 | 0.0019 |

Of the symbiont SSU reads placed, 3100 had a top hit outside Bacteria and 1211 a top hit in chloroplasts or mitochondria (dropped), and 0 had no hit at 75% or more (kept); `marker_share_entry7.tsv` gives them per run. SingleM is not used in RNA runs (n/a).
RNA runs with more symbiont SSU reads than --rna-sample (500): SRR1659836, SRR1659837, SRR1659838, SRR1659840, SRR1659841, SRR1659842, SRR1659846, SRR1659847. There a random sample was placed (seed 1), and SSU kept applies its kept share to all of the run's symbiont SSU reads. RNA runs are reported only.
