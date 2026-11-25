# Stage 0a, step 1: reference titer t (entry10)

Made by run_stage0a.py v0.2 on 2026-11-24 19:23 UTC, with assign_reads.py v0.2 and minimap2 2.31-r1302, against the full index (refs.fna md5 a2788e50287328b5ff926d49191ebe3e; `references/INDEX.md`). Reads: the first 20,000,000 spots of each run, fetched with fastq-dump 3.4.1. Assignments in the work folder's `stage0a_entry10/assign/`.
Rule: PROTOCOL.md 4.1 (reference titer) with section 10, entry 5: t is the lowest lineage read fraction among the DNA positives, divided by 10. RNA positives are measured and reported; they do not set t.

**t = 0.0009396**: the lowest DNA-positive fraction is 0.009396, in SRR10215376 (G. margarita BEG34 B+, genome reads, B+ line, mate-pair 3 kb (GMAR_GDNA-2-MP-3kb)).
Before this (`fractions.tsv`): t = 0.0009827.

| Run | Test | Kind | Line | Symbiont | Sets t | Reads | Symbiont reads | Fraction | Spread |
|---|---|---|---|---|---|---|---|---|---|
| SRR36098601 | T1 | DNA | G. margarita BEG34 B+ | CaGg | yes | 40000000 | 2079811 | 0.052 | 0.99 |
| SRR10215375 | T1 | DNA | G. margarita BEG34 B+ | CaGg | yes | 40000000 | 379799 | 0.009495 | 0.99 |
| SRR10215376 | T1 | DNA | G. margarita BEG34 B+ | CaGg | yes | 40000000 | 375821 | 0.009396 | 0.99 |
| SRR10215377 | T1 | DNA | G. margarita BEG34 B+ | CaGg | yes | 40000000 | 485853 | 0.01215 | 0.99 |
| SRR1659836 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 32562116 | 65154 | 0.002001 |  |
| SRR1659837 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 31334876 | 22666 | 0.0007233 |  |
| SRR1659838 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 35025904 | 30286 | 0.0008647 |  |
| SRR1659839 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 31780558 | 2980 | 9.377e-05 |  |
| SRR1659840 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 20000000 | 9618 | 0.0004809 |  |
| SRR1659841 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 20000000 | 11505 | 0.0005753 |  |
| SRR1659842 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 20000000 | 23552 | 0.001178 |  |
| SRR1659846 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 16947449 | 12269 | 0.0007239 |  |
| SRR1659847 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 19541952 | 21217 | 0.001086 |  |
| SRR1659850 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 20000000 | 27 | 1.35e-06 |  |
| SRR1659851 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 20000000 | 54 | 2.7e-06 |  |
| SRR15311774 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 40000000 | 322 | 8.05e-06 |  |
| SRR15311775 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 40000000 | 168 | 4.2e-06 |  |
| SRR15311776 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 40000000 | 1935 | 4.837e-05 |  |
| SRR6039993 | T2 | DNA | D. epigaea S1 | CaMg/MRE | yes | 40000000 | 3413545 | 0.08534 | 1.00 |
| SRR6039994 | T2 | DNA | D. epigaea S1 | CaMg/MRE | yes | 40000000 | 1853553 | 0.04634 | 1.00 |
| SRR6039995 | T2 | DNA | D. epigaea S1 | CaMg/MRE | yes | 40000000 | 3132450 | 0.07831 | 1.00 |
| SRR6039996 | T2 | DNA | D. epigaea S1 | CaMg/MRE | yes | 40000000 | 2995463 | 0.07489 | 1.00 |
| SRR6039997 | T2 | DNA | D. epigaea S1 | CaMg/MRE | yes | 40000000 | 2182017 | 0.05455 | 1.00 |
| SRR6039998 | T2 | DNA | D. epigaea S1 | CaMg/MRE | yes | 40000000 | 2591996 | 0.0648 | 1.00 |
| SRR12881897 | T2 | DNA | Rhizophagus sp. LPA8 | CaMg/MRE | yes | 40000000 | 489398 | 0.01223 | 1.00 |
| SRR12881898 | T2 | DNA | Rhizophagus sp. LPA8 | CaMg/MRE | yes | 40000000 | 503527 | 0.01259 | 1.00 |
| SRR12881899 | T2 | DNA | Rhizophagus sp. LPA8 | CaMg/MRE | yes | 40000000 | 393098 | 0.009827 | 1.00 |
| SRR12881900 | T2 | DNA | R. irregularis A2 | CaMg/MRE | yes | 40000000 | 972422 | 0.02431 | 1.00 |
| SRR12881901 | T2 | DNA | R. irregularis A2 | CaMg/MRE | yes | 40000000 | 972034 | 0.0243 | 1.00 |
| SRR12881902 | T2 | DNA | R. irregularis A2 | CaMg/MRE | yes | 40000000 | 681857 | 0.01705 | 1.00 |
| SRR12881903 | T2 | DNA | R. irregularis A2 | CaMg/MRE | yes | 40000000 | 1708551 | 0.04271 | 1.00 |

Reads are counted per read (mates separately). Details per run, with fungal and organelle reads: `fractions_entry10.tsv`. Spread (DNA runs; section 10, entry 10): the highest ratio, among the symbiont's genomes, of the 1-kb windows its mate-1 reads hit to the windows they would hit if spread evenly; a real genome gives about 1, a pile of reads in a few places well under 0.5. Step 2 measures the marker share m (section 10, entry 6).
