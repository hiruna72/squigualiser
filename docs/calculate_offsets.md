# Calculating base shift

As discussed in [pore model](pore_model.md) and [base shift and event alignment](base_shift_and_eventalignment.md) a base shift must be calculated to properly adjust the base color to the raw signal movements.
If the user is planning to use `reform` for a profile that is not listed then he may also use `calculate_offsets.py` to calculate appropriate kmer lenght and signal move offset value to use in the `refom` command.

## Option 1
To calculate the kmer length  and signal move offset for a dataset to use in `reform` use the following command.
```
python calculate_offsets.py -p reform.paf -s reads.blow5 -f reads.fastq
```

## Option 2
To calculate the kmer length  and signal move offset for a specific read and visualise the best base shift use the following command.
```
python calculate_offsets.py -p reform.paf -s reads.blow5 -f reads.fastq -o out.pdf --read_id ${READ_ID} 
```

## Option 3
To calculate and visualise the best base shift for a kmer model use the following command.
```
python calculate_offsets.py --use_model --model pore.model -o out.pdf 
```

The [table](reform.md/#precomputed-kmer-lengths-and-signal-moves-offsets) listing the profiles, recommended kmer lengths and the signal move offsets was calculated using [Option 1](#option-1)

The density plots documented in [pore model](pore_model.md) were calculated using [Option 3](#option-3)

The following density plots were generated using [Option 2](#option-2) for a specific read to find the most significant base offset.

| profile                            | denisty plot file                                                                              |
|------------------------------------|------------------------------------------------------------------------------------------------|
| guppy_dna_r9.4.1_450bps_fast_prom  |  [pdf](density_plots/8e1a33c4-af69-471c-a115-6428c8bf63df_dna_r9.4.1_450bps_fast_prom.cfg.pdf) |
| guppy_dna_r9.4.1_450bps_hac_prom   |   [pdf](density_plots/8e1a33c4-af69-471c-a115-6428c8bf63df_dna_r9.4.1_450bps_hac_prom.cfg.pdf) |
| guppy_dna_r9.4.1_450bps_sup_prom   |   [pdf](density_plots/8e1a33c4-af69-471c-a115-6428c8bf63df_dna_r9.4.1_450bps_sup_prom.cfg.pdf) |
| guppy_dna_r10.4.1_e8.2_400bps_fast | [pdf](density_plots/35142bde-548d-4f55-bf50-21c4cdd254da_dna_r10.4.1_e8.2_400bps_fast.cfg.pdf) |
| guppy_dna_r10.4.1_e8.2_400bps_hac  |  [pdf](density_plots/35142bde-548d-4f55-bf50-21c4cdd254da_dna_r10.4.1_e8.2_400bps_hac.cfg.pdf) |
| guppy_dna_r10.4.1_e8.2_400bps_sup  |  [pdf](density_plots/35142bde-548d-4f55-bf50-21c4cdd254da_dna_r10.4.1_e8.2_400bps_sup.cfg.pdf) |
