# ideal-goggles

A simple tool to Visualise nanopore raw signal-base alignment

![image](test/igv.png)
![image](test/pileup_plot.png)

1. In this figure an IGV plot shows the region chr1:6,811,403-6,811,433. 
2. The first read has deletions and the corresponding gaps appear in the signal plot [pileup_plot_0.html](test/pileup_plot_0.html).
3. The last read has an insertion and hence an insertion appears in the last signal plot. 
4. The second read is a reverse mapping and hence its shape is different from the rest.
5. Another example plot [pileup_plot_1.html](test/pileup_plot_1.html)

## INSTALLATION

### using python environment
````
git clone https://github.com/hiruna72/ideal-goggles.git
cd ideal-goggles
python3 -m venv idealg
source idealg/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
````
### using conda environment
````
git clone https://github.com/hiruna72/ideal-goggles.git
cd ideal-goggles
conda create -n idealg python=3.6.10 -y
conda activate idealg
pip install -r requirements.txt
````

## Method 1 - Read to signal visualisation
1. Run basecaller ([slow5-dorado](https://github.com/hiruna72/slow5-dorado), [buttery-eel](https://github.com/Psy-Fer/buttery-eel) or ont-Guppy)
```
# buttery-eel (tested with v0.2.2)
buttery-eel -g [GUPPY exe path] --config [DNA model] -i [INPUT] -o [OUTPUT] --port 5558 --use_tcp -x "cuda:all" --moves_out
e.g buttery-eel -g [GUPPY exe path] --config dna_r10.4.1_e8.2_400bps_sup.cfg -i input_reads.blow5 -o out.sam --port 5558 --use_tcp -x "cuda:all" --moves_out 

# slow5-dorado (tested with v0.2.1)
slow5-dorado basecaller [DNA model] [INPUT] --emit-moves > [OUTPUT]
e.g. slow5-dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.0.0 input_reads.blow5 --emit-moves > out.sam

# ont-guppy (tested with v6.3.7)
guppy_basecaller -c [DNA model] -i [INPUT] --moves_out --bam_out --save_path [OUTPUT]
samtools merge pass/*.bam -o pass_bam.bam # merge passed BAM files to create a single BAM file
```

2. Reformat move table 
```
# PAF output for plotting
REFORMAT_PAF=reform_output.paf
python src/reform.py --sig_move_offset 1 --kmer_length 1 -c --bam out.sam -o ${REFORMAT_PAF}

# For human readability you may prefer the tsv output (not supported for plitting)
python src/reform.py --sig_move_offset 1 --kmer_length 1 --bam out.sam -o reform_output.tsv
```

The output format is explained [here](https://hasindu2008.github.io/f5c/docs/output#resquiggle). `sig_move_offset` in commands above is the number of moves `n` to skip in the signal (`n x stride`) to correct the start of the alignment. This will not skip bases in the fastq sequence.

Pysam does not allow reading SAM/BAM files without a `@SQ` line in the header.
Hence, `reform.py` script might error out with `NotImplementedError: can not iterate over samfile without header`.
Add a fake `@SQ` header line with a zero length reference as follows,

```
echo -e fake_reference'\t'0 > fake_reference.fa.fai
samtools view out.sam -h -t fake_reference.fa.fai -o sq_added_out.sam
```

3. Visualise the signal to sequence alignment
````
FASTA_FILE=read.fasta
SIGNAL_FILE=read.blow5
OUTPUT_HTML=output.html

# use samtools fasta command to create .fasta file from SAM/BAM file
samtools fasta out.sam > ${FASTA_FILE}
# plot it
python src/sqp.py --fasta ${FASTA_FILE} --slow5 ${SIGNAL_FILE} --alignment ${REFORMAT_PAF} --output ${OUTPUT_HTML}
````
## Method 2 - Reference to signal visualisation
The first 3 steps are same as Method 1.
1. Run basecaller
2. Reformat move table
3. Align reads to reference genome
```
REFERENCE=genome.fa
MAPP_SAM=map_output.sam
minimap2 -ax map-ont ${REFERENCE} -t32 --secondary=no pass.fastq -o ${MAP_SAM}

```
5. Realign move array to reference
```
REALIGN_BAM=realign_output.bam
python src/realign.py --bam ${MAPP_SAM} --paf ${REFORMAT_PAF} -o ${REALIGN_BAM}
```

6. Visualise the signal to sequence alignment
````
SIGNAL_FILE=read.slow5
OUTPUT_DIR=output_sqp
REGION=chr1:6811404-6811443

python src/sqp.py --fasta ${REFERENCE} --slow5 ${SIGNAL_FILE} --alignment ${REALIGN_BAM} --output_dir ${OUTPUT_DIR} --tag_name "sqp_fun" --region ${REGION}
````

### Note
1. To get a pileup view, use `scripts/cat_plots.sh` to concatenate multiple `.html` plots in a directory.
2. To skip generating plots for reads mapped in reverse, use `--no_reverse` flag.

## Move table explanation (unconfirmed)
Nanopore basecallers output move arrays in SAM/BAM format. The important fields are listed below.
1. read_id
2. basecalled fastq sequence length
3. basecalled fastq sequence
4. stride used in the neural network (down sampling factor)
5. raw signal length
6. raw signal trim offset
7. move table

An example move array looks like the following,
```
110100010101000101011010101111…
```
The number of ones (1) in the move array equals to the fastq sequence length. 
According to the above example the first move corresponds with `1 x stride` signal points. 
The second move corresponds with `2 x stride` signal points. The third with `4 x stride`, the fourth with `2 x stride` and so on.


## Example
````
EXAMPLE_DIR=test/data/sqp/sigb_formater
FASTA_FILE=${EXAMPLE_DIR}/read.fasta
SIGNAL_FILE=${EXAMPLE_DIR}/read.slow5
ALIGN_FILE=${EXAMPLE_DIR}/r1k1m1.paf
OUTPUT_HTML=output.html

python src/sqp.py --fasta ${FASTA_FILE} --slow5 ${SIGNAL_FILE} --alignment ${ALIGN_FILE} --output ${OUTPUT_HTML}

````
