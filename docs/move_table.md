# Move table

Nanopore basecallers output move arrays in SAM/BAM format. The important fields are listed below.

1. read_id
2. basecalled fastq sequence length
3. basecalled fastq sequence
4. raw signal length in `ns` tag
5. raw signal trim offset in `ts` tag
6. move table in `mv` tag
7. stride used in the neural network (down sampling factor) in `mv` tag

An example move table looks like the following,
```
How the auxiliary field is stored in SAM format -> mv:B:c:5,1,1,0,1,0,0,0,1,0,1,0,1,0,0,0,1,0,1,0,1,1,0,1,0,1,0,1,1,1,1,…
Stride (always the first integer) -> 5
The actual move array (the rest) -> 1,1,0,1,0,0,0,1,0,1,0,1,0,0,0,1,0,1,0,1,1,0,1,0,1,0,1,1,1,1,…
```
The number of ones (1) in the actual move array equals to the fastq sequence length. 
According to the above example the first move corresponds with `1 x stride` signal points. 
The second move corresponds with `2 x stride` signal points. The third with `4 x stride`, the fourth with `2 x stride` and so on (see illustration below).

![image](figures/move_table_annotation.png)

Different models have different strides. Strides from some of the guppy models are listed below.

| Model                            | Stride |
|----------------------------------|--------|
| dna_r10.4.1_e8.2_400bps_fast.cfg | 5      |
| dna_r10.4.1_e8.2_400bps_hac.cfg  | 5      |
| dna_r10.4.1_e8.2_400bps_sup.cfg  | 5      |
| dna_r9.4.1_450bps_fast_prom.cfg  | 5      |
| dna_r9.4.1_450bps_hac_prom.cfg   | 5      |
| dna_r9.4.1_450bps_sup_prom.cfg   | 5      |
| rna_r9.4.1_70bps_hac_prom.cfg    | 10     |
| rna_r9.4.1_70bps_fast_prom.cfg   | 12     |