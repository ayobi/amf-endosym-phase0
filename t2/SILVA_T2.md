# SILVA for T2

Made by make_t2_silva.py v0.1 on 2026-11-24 17:38 UTC, from `SILVA_144_SSURef_NR99_tax_silva_trunc.fasta.gz` (905,628 sequences; the download's md5 ccead7bba52075f5d8ba2bfab5adc4c1, as SILVA publishes it; the unpacked FASTA's md5 a07fe2ee6ca5e79602221365699d73aa).

Rule: PROTOCOL.md section 5 (T2), as applied in section 10, entry 11. Each header after the sequence id (the seven ranks and the organism name) is searched for the eight taxon names, Moeniiplasma and Mycoplasma-related (followed by any word), ignoring case, and for MRE in capitals, not inside a longer word in capitals. A sequence is removed if any of them is found.

**Removed: 6,604 sequences**; kept: 899,024 (md5 of the T2 FASTA ef134dec00899bb11ac31642543c81cc).

| Term | Sequences it matches |
|---|---|
| Mollicutes | 3,599 |
| Mycoplasmatota | 6,596 |
| Tenericutes | 29 |
| Mycoplasmatales | 405 |
| Entomoplasmatales | 6 |
| Acholeplasmatales | 975 |
| Anaeroplasmatales | 98 |
| Haloplasmatales | 30 |
| Moeniiplasma | 4 |
| Mycoplasma-related | 0 |
| MRE | 6 |

A sequence can match several terms. Every removed sequence, with what matched and its SILVA path: `silva_removed.tsv`. run_t2.py checks that the T2 database was made from this list.

## Check: Mollicutes genera left in

A check of the rule's reach (entry 11); it removes nothing. Kept sequences whose genus rank or organism name holds one of Acholeplasma, Anaeroplasma, Asteroleplasma, Entomoplasma, Haloplasma, Hepatoplasma, Malacoplasma, Mesomycoplasma, Mesoplasma, Metamycoplasma, Mycoplasma, Mycoplasmoides, Mycoplasmopsis, Phytoplasma, Spiroplasma, Ureaplasma: **1,357**.

| Genus named | Sequences |
|---|---|
| Haloplasma | 1,340 |
| Mesomycoplasma | 1 |
| Metamycoplasma | 2 |
| Mycoplasma | 3 |
| Mycoplasmoides | 1 |
| Mycoplasmopsis | 1 |
| Phytoplasma | 9 |

| Sequence | SILVA path |
|---|---|
| FV536182.118321.119816 | Bacteria;Bacillati;Bacillota;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus;Mycoplasmoides pneumoniae |
| KY491109.1.932 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured Methanothermococcus sp. |
| GQ848372.1.914 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| GU270150.1.915 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| GQ927560.1.915 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| JQ815954.1.1399 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| GU363088.1.1325 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| FJ900695.1.923 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| FJ902269.1.1198 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| KC437268.1.1215 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| KC437275.1.1215 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| KY491069.1.932 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured Geoglobus sp. |
| KY491090.1.931 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured Thermogymnomonas sp. |
| KY491093.1.930 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured Methanotorris sp. |
| KY491132.1.926 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured Methanolobus sp. |
| KY491138.1.938 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured Methanocalculus sp. |
| KY491158.1.935 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured Methanotorris sp. |
| KY491159.1.932 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured Thermogymnomonas sp. |
| KY510302.1.919 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured Methanoregula sp. |
| KY510363.1.929 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured Methanoculleus sp. |
| KY510382.1.925 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured Thermococcus sp. |
| KY510384.1.941 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured Methanosarcina sp. |
| KY510388.1.922 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured Methanoculleus sp. |
| KY510397.1.909 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured Methanoculleus sp. |
| KY510420.1.926 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured Thermococcus sp. |
| KY510422.1.931 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured Methanococcus sp. |
| KY510430.1.934 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured Thermococcus sp. |
| FJ604791.1.914 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| GQ927643.1.904 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| AB019737.1.1394 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;unidentified archaeon |
| FJ655609.1.1398 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| FJ655646.1.1399 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| FJ655660.1.1398 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| FJ655721.1.1399 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| HE604515.1.916 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| HE604570.1.914 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| AB801654.1.1398 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| AB801679.1.1398 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| AB803281.1.1398 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| AB803452.1.1398 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| KC682088.1.1397 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| GU553508.1.1398 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| JQ327888.1.919 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| AB355097.1.915 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| JX870241.1.916 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| JX870317.1.915 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| JX870311.1.915 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| JX870313.1.915 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| HF922241.1.1358 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| JN408852.1.1343 | Archaea;Methanobacteriati--other;Thermoplasmatota;Thermoplasmata--other;Incertae Sedis;Incertae Sedis;Haloplasma--other;uncultured archaeon |
| ... | 1,307 more |
