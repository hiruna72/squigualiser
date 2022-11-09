# ideal-goggles
Visualize nanopore raw signals

![image](test/plot.png)

## INSTALLATION
````
git clone --recursive git@github.com:hiruna72/ideal-goggles.git
cd ideal-goggles
conda create -n idealg python=3.6.10 -y
conda activate idealg
pip install -r requirements.txt
````

## STEPS
1. Basecall a read using ont-guppy 6.3.7 to get the move table
````
CHUNK_SIZE=8000
CHUNKS_PER_RUNNER=100 (decrease this value if you run into CUDA memory problems)
MODEL_R9=dna_r9.4.1_450bps_sup_prom.cfg
MODEL_R10=dna_r10.4.1_e8.2_400bps_sup.cfg
SIGNAL_FILE=signal
OUT_DIR=output

./guppy/basecaller/guppy_basecaller -c ${MODEL_R9} -i ${SIGNAL_FILE} --device cuda:0  --moves_out --bam_out --chunk_size ${CHUNK_SIZE} --chunks_per_runner ${CHUNKS_PER_RUNNER} --save_path ${OUT_DIR} 
````
2. Convert move table to the [f5c resquiggle](https://hasindu2008.github.io/f5c/docs/output) format
````
BAM_FILE=${OUT_DIR}/pass/*.bam
ALIGN_FILE=${OUT_DIR}/pass/alignment.txt

samtools view ${BAM_FILE} | cut -f 1,9,12,19,20 | sed 's/ts:i://' | sed 's/ns:i://' | sed 's/mv:B:c,//g' | sed 's/,/\t/' | sed 's/,//g' | awk '{print $1"\t"$5"\t"$6"\t"$5"\t+\t"$1"\t"$2"\t0\t"$2"\t"$2"\t"$2"\t0\t"$4}' > TEMP0
cut -f 13 TEMP0 | sed 's/1/,1/g' | tr ',' '\n' | tail -n+2 | awk 'BEGIN{STRIDE=5; printf"sc:f:1.000000\tsh:f:0.000000\tss:Z:"} {printf(length($0)*STRIDE",")}' | sed 's/.$//' > TEMP1
cat TEMP0 | cut --complement -f 13 | paste - TEMP1 > ${ALIGN_FILE}
rm TEMP0 TEMP1
````
3. Visualize the signal to sequence alignment
````
FASTQ_FILE=${OUT_DIR}/pass/*.fastq
SIGNAL_FILE
ALIGN_FILE
OUTPUT_HTML=${OUT_DIR}/output.html

python src/sqp.py -f ${FASTQ_FILE} -s ${SIGNAL_FILE} -a ${ALIGN_FILE} -o ${OUTPUT_HTML} 
````
