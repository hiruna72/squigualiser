# Reform

The format how the moves information stores in the [move table](move_table.md) is not generalised. That is a move that is not a multiple of the stride cannot be stored in that format.
A more generalised format is documented [here](https://hasindu2008.github.io/f5c/docs/output#ss-tag). In addition to the flexibility this format has a greater compression ratio (~1.6) in ASCII format compared to the move table format.
Squigualiser uses this format as the input to generate plots.
Hence, before proceeding to create plots the move table has to be converted to the required format.
This is done using the tool `reform`.
The two important parameters `reform` takes are `kmer_length or k` and `sig_move_offset or m`.
This determines the base colour adjustment to the signal jumps. More information is discussed in [pore_model document](pore_model.md) and [using base shift](base_shift_and_eventalignment.md).

The user can provide `--profile` argument to deteremine the kmer  length and the signal move offset using preset values.
Optionally the user can also provide `-k` and `-m` values.
The argument  `--profile` will override `-k` and/or `-m` values if all are provided.

```
ALIGNMENT=reform_output.paf
PROFILE="guppy_dna_r10.4.1_e8.2_400bps_fast"
squigualiser reform -c --bam out.sam -o ${ALIGNMENT} --profile ${PROFILE} 

KMER_LENGTH=1
SIG_MOVE_OFFSET=0
squigualiser reform --kmer_length ${KMER_LENGTH} --sig_move_offset ${SIG_MOVE_OFFSET} -c --bam out.sam -o ${ALIGNMENT}
```

## Precomputed kmer lengths and signal moves offsets
The following precomputed kmer length and signal move offset values are available as profiles.

| profile                            | kmer length | sig move offset |
|------------------------------------|-------------|-----------------|
| guppy_dna_r9.4.1_450bps_fast_prom  |           3 |               2 |
| guppy_dna_r9.4.1_450bps_hac_prom   |           3 |               2 |
| guppy_dna_r9.4.1_450bps_sup_prom   |           4 |               3 |
| guppy_dna_r10.4.1_e8.2_400bps_fast |           1 |               0 |
| guppy_dna_r10.4.1_e8.2_400bps_hac  |           1 |               0 |
| guppy_dna_r10.4.1_e8.2_400bps_sup  |           1 |               0 |

If `reform` errors out with the message `"Error: specified profile is not found. Please run reform with -k 1 -s 0. Then run calculate_offsets.py and rerun reform with the recommended kmer_length and sig_move_offset."`, then the user is advised to run `calculate_offsets.py` tool.
Information regarding `calcualte_offsets` is documented [calculate_offsets](calculate_offsets.md).

The ASCII BYTE count ratio  between the two formats can be easily calculated using the following bash commands.
```
BASECALLER_MOVES_FILE=
REFORM_FILE=
mv_array_size=$(samtools view ${BASECALLER_MOVES_FILE} | awk '{for(i=5;i<=NF;i++){if($i~/mv/){a=$i}} print a}' | wc -c)
ss_array_size=$(cat ${REFORM_FILE} | awk '{for(i=5;i<=NF;i++){if($i~/ss:Z:/){a=$i}} print a}' | wc -c)
awk -v var1=${mv_array_size} -v var2=${ss_array_size} 'BEGIN { print  ( var1 / var2 ) }'
```
