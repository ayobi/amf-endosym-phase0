# Stage 0a, step 1: reference titer t

Made by run_stage0a.py v0.1 on 2026-11-21 16:01 UTC, with assign_reads.py v0.1 and minimap2 2.31-r1302, against the full index (refs.fna md5 ecbdc1bd4f9f45a45af3d49dd6585468; `references/INDEX.md`). Reads: the first 20,000,000 spots of each run, fetched with fastq-dump 3.4.1.
Rule: PROTOCOL.md 4.1 (reference titer) with section 10, entry 5: t is the lowest lineage read fraction among the DNA positives, divided by 10. RNA positives are measured and reported; they do not set t.

**t = 0.0009827**: the lowest DNA-positive fraction is 0.009827, in SRR12881899 (Rhizophagus sp. LPA8, plate LPA8_1, metagenome).

| Run | Test | Kind | Line | Symbiont | Sets t | Reads | Symbiont reads | Fraction |
|---|---|---|---|---|---|---|---|---|
| SRR36098601 | T1 | DNA | G. margarita BEG34 B+ | CaGg | yes | 40000000 | 2150922 | 0.05377 |
| SRR10215375 | T1 | DNA | G. margarita BEG34 B+ | CaGg | yes | 40000000 | 402536 | 0.01006 |
| SRR10215376 | T1 | DNA | G. margarita BEG34 B+ | CaGg | yes | 40000000 | 409195 | 0.01023 |
| SRR10215377 | T1 | DNA | G. margarita BEG34 B+ | CaGg | yes | 40000000 | 524606 | 0.01312 |
| SRR1659836 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 32562116 | 104278 | 0.003202 |
| SRR1659837 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 31334876 | 44941 | 0.001434 |
| SRR1659838 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 35025904 | 64218 | 0.001833 |
| SRR1659839 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 31780558 | 26630 | 0.0008379 |
| SRR1659840 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 20000000 | 9986 | 0.0004993 |
| SRR1659841 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 20000000 | 11821 | 0.000591 |
| SRR1659842 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 20000000 | 23923 | 0.001196 |
| SRR1659846 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 16947449 | 12864 | 0.0007591 |
| SRR1659847 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 19541952 | 22484 | 0.001151 |
| SRR1659850 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 20000000 | 98 | 4.9e-06 |
| SRR1659851 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 20000000 | 125 | 6.25e-06 |
| SRR15311774 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 40000000 | 755 | 1.888e-05 |
| SRR15311775 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 40000000 | 273 | 6.825e-06 |
| SRR15311776 | T1 | RNA | G. margarita BEG34 B+ | CaGg | no | 40000000 | 2576 | 6.44e-05 |
| SRR6039993 | T2 | DNA | D. epigaea S1 | CaMg/MRE | yes | 40000000 | 3413552 | 0.08534 |
| SRR6039994 | T2 | DNA | D. epigaea S1 | CaMg/MRE | yes | 40000000 | 1853554 | 0.04634 |
| SRR6039995 | T2 | DNA | D. epigaea S1 | CaMg/MRE | yes | 40000000 | 3132448 | 0.07831 |
| SRR6039996 | T2 | DNA | D. epigaea S1 | CaMg/MRE | yes | 40000000 | 2995465 | 0.07489 |
| SRR6039997 | T2 | DNA | D. epigaea S1 | CaMg/MRE | yes | 40000000 | 2182016 | 0.05455 |
| SRR6039998 | T2 | DNA | D. epigaea S1 | CaMg/MRE | yes | 40000000 | 2591996 | 0.0648 |
| SRR12881897 | T2 | DNA | Rhizophagus sp. LPA8 | CaMg/MRE | yes | 40000000 | 489398 | 0.01223 |
| SRR12881898 | T2 | DNA | Rhizophagus sp. LPA8 | CaMg/MRE | yes | 40000000 | 503527 | 0.01259 |
| SRR12881899 | T2 | DNA | Rhizophagus sp. LPA8 | CaMg/MRE | yes | 40000000 | 393098 | 0.009827 |
| SRR12881900 | T2 | DNA | R. irregularis A2 | CaMg/MRE | yes | 40000000 | 972422 | 0.02431 |
| SRR12881901 | T2 | DNA | R. irregularis A2 | CaMg/MRE | yes | 40000000 | 972034 | 0.0243 |
| SRR12881902 | T2 | DNA | R. irregularis A2 | CaMg/MRE | yes | 40000000 | 681857 | 0.01705 |
| SRR12881903 | T2 | DNA | R. irregularis A2 | CaMg/MRE | yes | 40000000 | 1708551 | 0.04271 |

Reads are counted per read (mates separately). Details per run, with fungal and organelle reads: `fractions.tsv`. Step 2 measures the marker share m (section 10, entry 6).
