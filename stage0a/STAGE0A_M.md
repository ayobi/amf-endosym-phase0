# Stage 0a, step 2: marker share m

Made by run_stage0a_m.py v0.1 on 2026-11-22 05:15 UTC, with SSUplex 0.2.0 (HMMER 3.4, profiles from Metaxa2 2.2.3), SingleM 0.21.4 with metapackage `S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb`, and minimap2 2.31-r1302 against the full index (refs.fna md5 ecbdc1bd4f9f45a45af3d49dd6585468). Reads: the same as step 1, each mate file on its own.

Rule: PROTOCOL.md 4.1 (marker share) and 4.2 (checks), with section 10, entry 6. For each run, the share of the known symbiont's reads that the route finds: the reads SSUplex calls bacterial SSU, plus the reads in SingleM's OTU table. The found reads were aligned again exactly as in step 1. m is the lowest share among the DNA positives.

**m = 0.001585**: the lowest DNA-positive share, in SRR6039993 (D. epigaea S1, genome reads, mate-pair 3 kb (MP2_3K)).

Checks (section 4.2). A tool finds a symbiont when at least 10 of the reads it found map to that symbiont, in every DNA positive (section 10, entry 6):
- SSUplex finds CaGg and MRE SSU: **passed**. Fewest symbiont reads among the bacterial SSU reads of a DNA positive: 196, in SRR10215375.
- SingleM sees MRE: **yes**. Fewest MRE reads in SingleM's table for a T2 positive: 573, in SRR12881899.

| Run | Test | Kind | Symbiont | Sets m | Symbiont reads | By SSUplex | By SingleM | By either | Share |
|---|---|---|---|---|---|---|---|---|---|
| SRR36098601 | T1 | DNA | CaGg | yes | 2150922 | 811 | 3122 | 3933 | 0.001829 |
| SRR10215375 | T1 | DNA | CaGg | yes | 402536 | 196 | 653 | 849 | 0.002109 |
| SRR10215376 | T1 | DNA | CaGg | yes | 409195 | 205 | 570 | 775 | 0.001894 |
| SRR10215377 | T1 | DNA | CaGg | yes | 524606 | 196 | 921 | 1117 | 0.002129 |
| SRR1659836 | T1 | RNA | CaGg | no | 104278 | 1848 | 5 | 1853 | 0.01777 |
| SRR1659837 | T1 | RNA | CaGg | no | 44941 | 722 | 0 | 722 | 0.01607 |
| SRR1659838 | T1 | RNA | CaGg | no | 64218 | 1289 | 0 | 1289 | 0.02007 |
| SRR1659839 | T1 | RNA | CaGg | no | 26630 | 89 | 0 | 89 | 0.003342 |
| SRR1659840 | T1 | RNA | CaGg | no | 9986 | 2061 | 0 | 2061 | 0.2064 |
| SRR1659841 | T1 | RNA | CaGg | no | 11821 | 2314 | 0 | 2314 | 0.1958 |
| SRR1659842 | T1 | RNA | CaGg | no | 23923 | 4761 | 0 | 4761 | 0.199 |
| SRR1659846 | T1 | RNA | CaGg | no | 12864 | 2582 | 1 | 2583 | 0.2008 |
| SRR1659847 | T1 | RNA | CaGg | no | 22484 | 4589 | 0 | 4589 | 0.2041 |
| SRR1659850 | T1 | RNA | CaGg | no | 98 | 3 | 1 | 4 | 0.04082 |
| SRR1659851 | T1 | RNA | CaGg | no | 125 | 6 | 0 | 6 | 0.048 |
| SRR15311774 | T1 | RNA | CaGg | no | 755 | 77 | 0 | 77 | 0.102 |
| SRR15311775 | T1 | RNA | CaGg | no | 273 | 36 | 0 | 36 | 0.1319 |
| SRR15311776 | T1 | RNA | CaGg | no | 2576 | 352 | 0 | 352 | 0.1366 |
| SRR6039993 | T2 | DNA | CaMg/MRE | yes | 3413552 | 782 | 4630 | 5412 | 0.001585 |
| SRR6039994 | T2 | DNA | CaMg/MRE | yes | 1853554 | 2888 | 7424 | 10312 | 0.005563 |
| SRR6039995 | T2 | DNA | CaMg/MRE | yes | 3132448 | 3675 | 24794 | 28469 | 0.009088 |
| SRR6039996 | T2 | DNA | CaMg/MRE | yes | 2995465 | 2821 | 4346 | 7167 | 0.002393 |
| SRR6039997 | T2 | DNA | CaMg/MRE | yes | 2182016 | 3117 | 8826 | 11943 | 0.005473 |
| SRR6039998 | T2 | DNA | CaMg/MRE | yes | 2591996 | 4112 | 10525 | 14637 | 0.005647 |
| SRR12881897 | T2 | DNA | CaMg/MRE | yes | 489398 | 336 | 759 | 1095 | 0.002237 |
| SRR12881898 | T2 | DNA | CaMg/MRE | yes | 503527 | 328 | 707 | 1035 | 0.002056 |
| SRR12881899 | T2 | DNA | CaMg/MRE | yes | 393098 | 264 | 573 | 837 | 0.002129 |
| SRR12881900 | T2 | DNA | CaMg/MRE | yes | 972422 | 546 | 1489 | 2035 | 0.002093 |
| SRR12881901 | T2 | DNA | CaMg/MRE | yes | 972034 | 548 | 1419 | 1967 | 0.002024 |
| SRR12881902 | T2 | DNA | CaMg/MRE | yes | 681857 | 346 | 981 | 1327 | 0.001946 |
| SRR12881903 | T2 | DNA | CaMg/MRE | yes | 1708551 | 858 | 2533 | 3391 | 0.001985 |

Symbiont reads: those mapping counted for the run's symbiont in step 1. By SSUplex, by SingleM: the symbiont reads among the reads each tool found. By either: a read found by both counts once. All reads found by each tool, symbiont or not: `marker_share.tsv`.

Databases (`get_marker_databases.py`):
- SingleM metapackage `S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb` from doi:10.5281/zenodo.20150069, md5 580f9ae9ae6a0972f15665c7397c6369 as in `references/DATABASES.md`, checked with `singlem data --verify-only`.
- Metaxa2 2.2.3 (https://microbiology.se/sw/Metaxa2_2.2.3.tar.gz, md5 da9c548a8df325d305af30bc65f15ea7), made into SSUplex's profile files: bacteria (B.hmm, 15 profiles, md5 4b48822f26350e545fbdaeb919d5a88e); archaea (A.hmm, 15 profiles, md5 23e94222d87b8eecd34be504bb66b4f8); eukaryota (E.hmm, 15 profiles, md5 2cad87dfdc6d3a8cf9d70a622af05f46); mitochondria (M.hmm + N.hmm, 21 profiles, md5 d1c4f80a60ff853bb7ef4df38fc5ef05); chloroplast (C.hmm, 12 profiles, md5 95d0c15412395f9815fb4712a8427991).
- SSUplex quick start on Metaxa2's test.fasta: 50 of 100 classified, as in the SSUplex README.
