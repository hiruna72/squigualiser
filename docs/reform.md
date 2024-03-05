# Reform

The [move table](move_table.md) format used by ONT basecallers to store the move information is not generalised. For instance, a move that is not a multiple of the stride, cannot be stored in the move table format. A more generalised format to store the signal to read/reference alignment is documented [here](https://hasindu2008.github.io/f5c/docs/output#ss-tag), which we would refer to as the **ss tag**. In addition to the flexibility, this **ss tag** format has a greater compression ratio (~1.6) in ASCII format compared to the move table format.

Squigualiser uses this format **ss tag** as the input to generate plots. Hence, before proceeding to create plots, the move table has to be converted to the required **ss tag** format. This is done using the subtool `reform` in squigualiser.
The two important parameters `reform` take are `kmer_length (or k)` and `sig_move_offset (or m)`.
This determines the base colour adjustment to the signal jumps. More information is discussed in [pore_model document](pore_model.md) and [using base shift](base_shift_and_eventalignment.md).

The user can provide `--profile` argument to determine the kmer length and the signal move offset using preset values.
Optionally the user can also provide `-k` and `-m` values. The argument  `--profile` will override `-k` and/or `-m` values, if all arguments are provided.

```
ALIGNMENT=reform_output.paf
PROFILE="guppy_dna_r10.4.1_e8.2_400bps_fast"
squigualiser reform -c --bam basecalls.sam -o ${ALIGNMENT} --profile ${PROFILE} 

KMER_LENGTH=1
SIG_MOVE_OFFSET=0
squigualiser reform -k ${KMER_LENGTH} -m ${SIG_MOVE_OFFSET} -c --bam basecalls.sam -o ${ALIGNMENT}
```

## Precomputed kmer lengths and signal move offsets

Precomputed kmer lengths and signal move offsets can be found [here](profiles.md#precomputed-kmer-lengths-and-signal-moves-offsets).

If `reform` errors out with the message `"Error: specified profile is not found. Please run reform with -k 1 -s 0. Then run calculate_offsets.py and rerun reform with the recommended kmer_length and sig_move_offset."`, then the user is advised to run `calculate_offsets` subtool .
Information regarding `calcualte_offsets` is documented at [calculate_offsets](calculate_offsets.md).

## SS tag size vs move table size

The ASCII BYTE count ratio between the two formats can be easily calculated using the following bash commands.
```
BASECALLER_MOVES_FILE=basecalls.sam
REFORM_FILE=reform_output.paf
mv_array_size=$(samtools view ${BASECALLER_MOVES_FILE} | awk '{for(i=5;i<=NF;i++){if($i~/mv/){a=$i}} print a}' | wc -c)
ss_array_size=$(cat ${REFORM_FILE} | awk '{for(i=5;i<=NF;i++){if($i~/ss:Z:/){a=$i}} print a}' | wc -c)
awk -v var1=${mv_array_size} -v var2=${ss_array_size} 'BEGIN { print  ( var1 / var2 ) }'
```

## Reform TSV output 

Reform can also output a custom TSV format for human readability (not supported for plotting):
```
squigualiser reform -k ${KMER_LENGTH} -m ${SIG_MOVE_OFFSET} --bam basecalls.sam -o reform_output.tsv
```
