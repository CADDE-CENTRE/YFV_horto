# DATA_DICTIONARY.md

Column summaries for each supplementary CSV.


## SupplementaryData1.csv

- Rows: **818**, Columns: **13**
- Suggested use: Figure 1 maps / site locations (meteorological & geographical context)

| column | dtype | non_null | unique | min | max | example_values |
|---|---|---:|---:|---:|---:|---|
| ID | int64 | 818 | 818 | 1 | 818 | 1, 2, 3, 4, 5 |
| Location | object | 818 | 1 |  |  | Horto, Horto, Horto, Horto, Horto |
| Date | object | 818 | 43 |  |  | 2017-12-22, 2017-12-22, 2017-12-22, 2017-12-22, 2017-12-22 |
| Latitude | float64 | 818 | 33 | -23.46634 | -23.45503 | -23.45906, -23.45906, -23.45906, -23.45906, -23.45906 |
| Longitude | float64 | 818 | 35 | -46.6443 | -46.62975 | -46.63657, -46.63657, -46.63657, -46.63657, -46.63657 |
| Heigh | object | 818 | 2 |  |  | Ground_level, Ground_level, Ground_level, Ground_level, Ground_level |
| Technique | object | 818 | 3 |  |  | Hand_net, Hand_net, Hand_net, Hand_net, Hand_net |
| Pool_origin | object | 818 | 818 |  |  | 3, 4, 11, 5, 13 |
| Genus | object | 818 | 10 |  |  | Culex, Culex, Culex, Aedes, Aedes |
| Species | object | 818 | 30 |  |  | (Cux) sp., (Cux) sp., (Cux) sp., albopictus, albopictus |
| Sex | object | 818 | 3 |  |  | Male, Male, Male, Female_notbloodfed, Female_notbloodfed |
| Quantity | int64 | 818 | 9 | 1 | 10 | 5, 5, 3, 1, 2 |
| PCR_YFV | object | 818 | 3 |  |  | Negative, Negative, Negative, Negative, Negative |

## SupplementaryData2.csv

- Rows: **118**, Columns: **13**
- Suggested use: General analysis inputs (see Data Dictionary)

| column | dtype | non_null | unique | min | max | example_values |
|---|---|---:|---:|---:|---:|---|
| State | object | 118 | 3 |  |  | MG, RJ, RJ, RJ, RJ |
| Municipality | object | 118 | 33 |  |  | Belmiro Braga, Valença, Valença, Nova Iguaçú, Valença |
| Collection date | object | 118 | 44 |  |  | 29/01/2018, 26/01/2018, 24/01/2018, 09/01/2018, 19/01/2018 |
| Code | object | 69 | 69 |  |  | MG 3155, VL 3016, VL 2967, NI 3104, VL 2926 |
| Species | object | 118 | 3 |  |  | janthinomys, janthinomys, janthinomys, janthinomys, janthinomys |
| Species_complete | object | 118 | 14 |  |  | Hg. janthinomys, Hg. janthinomys, Hg. janthinomys, Hg. janthinomys, Hg. janthinomys |
| Number | int64 | 118 | 15 | 1 | 24 | 5, 5, 5, 5, 1 |
| CT_mean | float64 | 118 | 61 | 10.0 | 38.0 | 16.1, 16.5, 16.9, 17.2, 17.5 |
| Study | int64 | 118 | 4 | 1 | 4 | 1, 1, 1, 1, 1 |
| Publication | object | 118 | 4 |  |  | Abreu et al. 2019, Abreu et al. 2019, Abreu et al. 2019, Abreu et al. 2019, Abreu et al. 2019 |
| study_DOI | object | 110 | 3 |  |  | 10.1080/22221751.2019.1568180, 10.1080/22221751.2019.1568180, 10.1080/22221751.2019.1568180, 10.1080/22221751.2019.1568180, 10.1080/22221751.2019.1568180 |
| RT-PCR_protocol | object | 118 | 3 |  |  | Bonaldo et al. 2017, Bonaldo et al. 2017, Bonaldo et al. 2017, Bonaldo et al. 2017, Bonaldo et al. 2017 |
| RT-PCR_protocol_DOI | object | 118 | 3 |  |  | 10.1590/0074-02760170134, 10.1590/0074-02760170134, 10.1590/0074-02760170134, 10.1590/0074-02760170134, 10.1590/0074-02760170134 |

## SupplementaryData3.csv

- Rows: **53**, Columns: **13**
- Suggested use: General analysis inputs (see Data Dictionary)

| column | dtype | non_null | unique | min | max | example_values |
|---|---|---:|---:|---:|---:|---|
| State | object | 53 | 7 |  |  | RJ, ES, MA, ES, MS |
| Species | object | 53 | 3 |  |  | Other, Other, Other, Other, Other |
| Species_complete | object | 53 | 13 |  |  | Sa. Chloropterus, Sa. chloropterus, Sa. Chloropterus, Sa. chloropterus, Sa. Chloropterus |
| Municipality | object | 53 | 22 |  |  | Angra dos Reis_Ilha_Grande, Santa Teresa_a, Pastos Bons, Pancas, Campo Grande |
| include | object | 53 | 2 |  |  | no, no, no, no, no |
| tested_pools | float64 | 46 | 22 | 1.0 | 164.0 | 1.0, 2.0, 3.0, 3.0, 1.0 |
| positive_pools | float64 | 51 | 8 | 0.0 | 35.0 | 1.0, 0.0, 1.0, 1.0, 4.0 |
| total_mosquitoes_captured | float64 | 30 | 24 | 1.0 | 1245.0 | 4.0, 73.0, 12.0, 1.0, 7.0 |
| MIR | float64 | 53 | 35 | 0.0 | 1000.0 | 1000.0, 0.0, 13.699, 83.333, 16.7 |
| MLE | float64 | 13 | 13 | 0.06 | 15.35 | 15.35, 5.55, 2.94, 1.31, 1.09 |
| Collection_start_date | object | 53 | 19 |  |  | 07/02/2018, 08/02/2017, 1994, 08/02/2017, 1992 |
| Collection_end_date | object | 28 | 2 |  |  | 02/03/2017, 02/03/2017, 02/03/2017, 02/03/2017, 02/03/2017 |
| doi | object | 53 | 3 |  |  | 10.3201/eid1612.100608, 10.3390/v14122805, 10.4269/ajtmh.1997.57.132, 10.3390/v14122805, 10.4269/ajtmh.1997.57.132 |

## SupplementaryData4.csv

- Rows: **98**, Columns: **24**
- Suggested use: Figure 2 & sequencing performance (Ct, RPM[log10], N50, coverage); Figure 1 maps / site locations (meteorological & geographical context)

| column | dtype | non_null | unique | min | max | example_values |
|---|---|---:|---:|---:|---:|---|
| id | int64 | 98 | 98 | 1 | 99 | 1, 22, 39, 38, 40 |
| species | object | 98 | 3 |  |  | Alouatta_guariba, Alouatta_guariba, Alouatta_guariba, Alouatta_guariba, Alouatta_guariba |
| host | object | 98 | 2 |  |  | NeotropicalPrimate, NeotropicalPrimate, NeotropicalPrimate, NeotropicalPrimate, NeotropicalPrimate |
| group | object | 98 | 3 |  |  | Group_II, Group_II, Group_III, Group_III, Group_III |
| pool_size | float64 | 9 | 4 | 1.0 | 5.0 | 3.0, 5.0, 1.0, 5.0, 5.0 |
| Pool_origin | float64 | 9 | 9 | 24.0 | 192.0 | 24.0, 66.0, 81.0, 95.0, 120.0 |
| Height_technique | object | 9 | 1 |  |  | Ground_level_Hand_net, Ground_level_Hand_net, Ground_level_Hand_net, Ground_level_Hand_net, Ground_level_Hand_net |
| latitude | float64 | 98 | 54 | -23.4803478 | -23.329005 | -23.464613, -23.459, -23.4189989, -23.4189989, -23.3668848 |
| longitude | float64 | 98 | 57 | -46.7092847 | -46.5776382 | -46.648398, -46.6381, -46.6655516, -46.6655516, -46.6136934 |
| date_identification_collection | object | 98 | 37 |  |  | 09/10/2017, 23/10/2017, 16/11/2017, 16/11/2017, 25/11/2017 |
| decomposition_simple | object | 55 | 3 |  |  | medium_or_advanced, intact_live_2d, intact_live_2d, intact_live_2d, intact_live_2d |
| Result | object | 98 | 1 |  |  | Positive, Positive, Positive, Positive, Positive |
| CT | float64 | 84 | 74 | 5.57 | 31.68 | 24.1, 7.0, 17.5, 18.06, 9.48 |
| amplicon_cov_pct | float64 | 56 | 51 | 37.4804 | 99.3748 | 85.9011, 83.0832, 86.2602, 86.5089, 80.0995 |
| smart9n_Nb_reads  | int64 | 98 | 97 | 15237 | 1164917 | 367336, 142280, 153420, 278911, 136409 |
| N50 | int64 | 98 | 84 | 201 | 945 | 400, 400, 294, 270, 493 |
| smart9n_Nb_of_reads_mapped | int64 | 98 | 97 | 15 | 211290 | 23288, 689, 38, 15073, 48114 |
| RPM | float64 | 98 | 96 | 2.9928 | 6.6813 | 5.8021, 4.6851, 3.3939, 5.7327, 6.5474 |
| smart9n_perc_reads_mapped | float64 | 98 | 96 | 0.01 | 48.01 | 6.34, 0.48, 0.03, 5.4, 35.27 |
| smart9n_average_depth_coverage | float64 | 98 | 97 | 0.69 | 7244.14 | 1175.25, 27.83, 1.2, 819.65, 1522.87 |
| Identity | float64 | 98 | 25 | 92.5 | 95.6 | 95.1, 95.1, 92.9, 94.0, 94.8 |
| smart9n_cov10x_pct | float64 | 98 | 29 | 0.0 | 99.908 | 99.908, 72.355, 0.092, 99.908, 99.908 |
| Final_dataset | object | 97 | 2 |  |  | No, No, Yes, Yes, Yes |
| Accession_number | object | 88 | 88 |  |  | OQ714250.1, OQ714252.1, OQ714294.1, OQ714295.1, OQ714251.1 |