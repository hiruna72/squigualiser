# Profiles

This document lists the different profiles.

## Precomputed kmer lengths and signal moves offsets
The following precomputed kmer length and signal move offset values are available as profiles.
These profiles are used in `reform` tool. More information can be found at [reform](reform.md) and [calculate_offsets](calculate_offsets.md).

| profile                            | kmer length | sig move offset |
|------------------------------------|-------------|-----------------|
| guppy_dna_r9.4.1_450bps_fast_prom  |           3 |               2 |
| guppy_dna_r9.4.1_450bps_hac_prom   |           3 |               2 |
| guppy_dna_r9.4.1_450bps_sup_prom   |           4 |               3 |
| guppy_dna_r10.4.1_e8.2_400bps_fast |           1 |               0 |
| guppy_dna_r10.4.1_e8.2_400bps_hac  |           1 |               0 |
| guppy_dna_r10.4.1_e8.2_400bps_sup  |           1 |               0 |

## Precomputed forward and reverse base shift values
The following precomputed forward and reverse base shift are available as profiles.
These profiles are used in `plot` and `plot_pileup` tool. More information can be found at [pore_model](pore_model.md), [calculate_offsets](calculate_offsets.md), [base_shift_of_reverse_mapped_reads](base_shift_of_reverse_mapped_reads.md), and [base_shift_and_eventalignment](base_shift_and_eventalignment.md).
Note that the default base_shift is 0. Hence, it is not necessary to provide the argument `--profile [name]` where base shift is 0. 

| name                                     | base_shift_forward | base_shift_reverse |
|------------------------------------------|--------------------|--------------------|
| kmer_model_dna_r9.4.1_450bps_5_mer       |         -2         |         -2         |
| kmer_model_dna_r9.4.1_450bps_6_mer       |         -2         |         -3         |
| kmer_model_rna_r9.4.1_70bps_5_mer        |         -1         |         -3         |
| kmer_model_dna_r10.4.1_e8.2_400bps_9_mer |         -6         |         -2         |
| guppy_dna_r9.4.1_450bps_fast_prom        |          0         |          0         |
| guppy_dna_r9.4.1_450bps_hac_prom         |          0         |          0         |
| guppy_dna_r9.4.1_450bps_sup_prom         |          0         |          0         |
| guppy_dna_r10.4.1_e8.2_400bps_fast       |          0         |          0         |
| guppy_dna_r10.4.1_e8.2_400bps_hac        |          0         |          0         |
| guppy_dna_r10.4.1_e8.2_400bps_sup        |          0         |          0         |