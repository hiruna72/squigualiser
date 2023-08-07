#!/bin/bash

# steps
# simulate a region using squigulator

RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo ; echo -e "${GREEN}$1${NC}" >&2 ; }
ask() { echo -n "Directory $1 exists. Delete and create again? (y/n)? " ; read answer ; if [ "$answer" != "${answer#[Nn]}" ] ;then exit ; fi ; echo ; }

info "$(date)"

# set -x

RUN_NO="RUN01"

R10_MODEL_FAST="dna_r10.4.1_e8.2_400bps_fast.cfg"
R10_MODEL_HAC="dna_r10.4.1_e8.2_400bps_hac.cfg"
R10_MODEL_SUP="dna_r10.4.1_e8.2_400bps_sup.cfg"

PLOT_TRACK_TOOL="squigualiser plot_tracks"
CALCULATE_OFFSETS_TOOL="squigualiser calculate_offsets"

## variable the user can change start here
SQUIGULATOR=squigulator
SAMTOOLS=samtools
SLOW5TOOLS=slow5tools
SAMTOOLS=samtools
MINIMAP2=minimap2
F5C=f5c

GUPPY=
BUTTERY_EEL_ENV_PATH=
REFERENCE=

MODEL_TO_USE=${R10_MODEL_SUP}
CHUNK_SIZE="--chunk_size 500"

SIMULATING_PROFILE="dna-r10-prom"

PROFILE_TO_DETERMINE_BASE_SHIFT="kmer_model_dna_r10.4.1_e8.2_400bps_9_mer"

SIGNAL_FILE_0="methylated.blow5"
SIGNAL_FILE_1="non_methylated.blow5"

READ_ID="251256d6-1c5d-4756-91bf-8fa5dc6fd1c8"
READ_REGION=${READ_ID}:1-500
REF_REGION="chr21:5,225,759-5,226,051"
SIM_REGION="sim_ref:1-400"
SIG_SCALE="--sig_scale znorm"

## variable the user can change end here

OUTPUT_DIR=${RUN_NO}_${MODEL_TO_USE}
BASECALL_DIR="${OUTPUT_DIR}/basecall_${MODEL_TO_USE}"
SIMULATE_READ_DIR="${OUTPUT_DIR}/simulate_read"
SIMULATE_REF_DIR="${OUTPUT_DIR}/simulate_ref"
SQUIG_PLOT_DIR="${OUTPUT_DIR}/squig_plot_reads"

SEQUENCE_FILE_0="${OUTPUT_DIR}/data/0/fastq/pass.fastq"
SEQUENCE_FILE_1="${OUTPUT_DIR}/data/1/fastq/pass.fastq"

MAPPED_BAM_0="${OUTPUT_DIR}/data/0/bamtx/map.bam"
MAPPED_BAM_1="${OUTPUT_DIR}/data/1/bamtx/map.bam"

EVENTALIGN_BAM_0="${OUTPUT_DIR}/data/0/f5c/eventalign.bam"
EVENTALIGN_BAM_1="${OUTPUT_DIR}/data/1/f5c/eventalign.bam"

EVENTALIGN_TSV_0="${OUTPUT_DIR}/data/0/f5c/eventalign.tsv"
EVENTALIGN_TSV_1="${OUTPUT_DIR}/data/1/f5c/eventalign.tsv"

METH_TSV_0="${OUTPUT_DIR}/data/0/f5c/meth.tsv"
METH_FREQ_TSV_0="${OUTPUT_DIR}/data/0/f5c/meth_freq.tsv"
METH_FREQ_BEDGRAPH_0="${OUTPUT_DIR}/data/0/f5c/meth_freq.bedgraph"
METH_FREQ_BIGWIG_0="${OUTPUT_DIR}/data/0/f5c/meth_freq.bigwig"
METH_BED_0="${OUTPUT_DIR}/data/0/f5c/meth.bed"

METH_TSV_1="${OUTPUT_DIR}/data/1/f5c/meth.tsv"
METH_FREQ_TSV_1="${OUTPUT_DIR}/data/1//f5c/meth_freq.tsv"
METH_FREQ_BEDGRAPH_1="${OUTPUT_DIR}/data/1/f5c/meth_freq.bedgraph"
METH_FREQ_BIGWIG_1="${OUTPUT_DIR}/data/1/f5c/meth_freq.bigwig"
METH_BED_1="${OUTPUT_DIR}/data/1/f5c/meth.bed"

METH_THRESHOLD=0.85

[ "${REFERENCE}" ] || die "edit the reference genome path (variable REF) in the script"

[ "${SQUIGULATOR}" ] || die "edit the executable path of squigulator (variable SQUIGULATOR) in the script"
[ "${SAMTOOLS}" ] || die "edit the executable path of samtools (variable SAMTOOLS) in the script"
[ "${SLOW5TOOLS}" ] || die "edit the executable path of slow5tools (variable SLOW5TOOLS) in the script"
[ "${MINIMAP2}" ] || die "edit the executable path of minimap2 (variable MINIMAP2) in the script"
[ "${F5C}" ] || die "edit the executable path of f5c (variable F5C) in the script"
[ "${GUPPY}" ] || die "edit the path to the dir the of guppy_basecaller (variable GUPPY) in the script"
[ "${BUTTERY_EEL_ENV_PATH}" ] || die "edit the env path of buttery-eel (variable BUTTERY_EEL_ENV_PATH) in the script"

[ ${SIMULATING_PROFILE} ] || die "set the variable SIMULATING_PROFILE"
[ "${READ_ID}" ] || die "READ_ID is not set"

[ "${RUN_NO}" ] || die "RUN_NO is not set"
[ "${OUTPUT_DIR}" ] || die "OUTPUT_DIR is not set"
[ "${BASECALL_DIR}" ] || die "BASECALL_DIR is not set"
[ "${SIMULATE_READ_DIR}" ] || die "SIMULATE_READ_DIR is not set"
[ "${SIMULATE_REF_DIR}" ] || die "SIMULATE_REF_DIR is not set"
[ "${SQUIG_PLOT_DIR}" ] || die "SQUIG_PLOT_DIR is not set"

[ "${PROFILE_TO_DETERMINE_BASE_SHIFT}" ] || die "PROFILE_TO_DETERMINE_BASE_SHIFT is not set"

create_output_dir() {
	test -d "${OUTPUT_DIR}" && ask "${OUTPUT_DIR}" 
	test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
	mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"
	mkdir "$OUTPUT_DIR/data" || die "Failed creating $OUTPUT_DIR"
	mkdir "$OUTPUT_DIR/data/0" || die "Failed creating $OUTPUT_DIR"
	mkdir "$OUTPUT_DIR/data/0/fastq" || die "Failed creating $OUTPUT_DIR"
	mkdir "$OUTPUT_DIR/data/0/bamtx" || die "Failed creating $OUTPUT_DIR"
	mkdir "$OUTPUT_DIR/data/0/f5c" || die "Failed creating $OUTPUT_DIR"
	mkdir "$OUTPUT_DIR/data/1" || die "Failed creating $OUTPUT_DIR"
	mkdir "$OUTPUT_DIR/data/1/fastq" || die "Failed creating $OUTPUT_DIR"
	mkdir "$OUTPUT_DIR/data/1/bamtx" || die "Failed creating $OUTPUT_DIR"
	mkdir "$OUTPUT_DIR/data/1/f5c" || die "Failed creating $OUTPUT_DIR"
}

basecall_fast5() {
	set -x
	info "basecalling using ${MODEL_TO_USE}"
	export CUDA_VISIBLE_DEVICES=0
	test -d "${BASECALL_DIR}" && rm -r "${BASECALL_DIR}"
	mkdir "${BASECALL_DIR}" || die "Failed creating ${BASECALL_DIR}"
	${SLOW5TOOLS} s2f ${SIGNAL_FILE_0} -d "${OUTPUT_DIR}/data/0/fast5"
	${GUPPY}/guppy_basecaller ${CHUNK_SIZE} -x cuda:0 --config ${MODEL_TO_USE} -i ${OUTPUT_DIR}/data/0/fast5 -s ${BASECALL_DIR} || die "guppy_basecaller failed"
	cat ${BASECALL_DIR}/pass/*.fastq > ${SEQUENCE_FILE_0} || die "cat fastq together failed"

	test -d "${BASECALL_DIR}" && rm -r "${BASECALL_DIR}"
	mkdir "${BASECALL_DIR}" || die "Failed creating ${BASECALL_DIR}"
	${SLOW5TOOLS} s2f ${SIGNAL_FILE_1} -d "${OUTPUT_DIR}/data/1/fast5"
	${GUPPY}/guppy_basecaller ${CHUNK_SIZE} -x cuda:0 --config ${MODEL_TO_USE} -i ${OUTPUT_DIR}/data/1/fast5 -s ${BASECALL_DIR} || die "guppy_basecaller failed"
	cat ${BASECALL_DIR}/pass/*.fastq > ${SEQUENCE_FILE_1} || die "cat fastq together failed"
}

minimap2_align() {
	info "aligning reads using minimap2..."
	${MINIMAP2} -ax map-ont ${REFERENCE} -t32 --secondary=no ${SEQUENCE_FILE_0} | ${SAMTOOLS} sort -o ${MAPPED_BAM_0} || die "minimap2_align failed"
	${SAMTOOLS} index ${MAPPED_BAM_0} || die "samtools index failed"

	${MINIMAP2} -ax map-ont ${REFERENCE} -t32 --secondary=no ${SEQUENCE_FILE_1} | ${SAMTOOLS} sort -o ${MAPPED_BAM_1} || die "minimap2_align failed"
	${SAMTOOLS} index ${MAPPED_BAM_1} || die "samtools index failed"
}

f5c_eventalign() {
	info "f5c eventalign..."
	f5c index ${SEQUENCE_FILE_0} --slow5 ${SIGNAL_FILE_0} || die "f5c index failed"
	f5c eventalign -b ${MAPPED_BAM_0} -r ${SEQUENCE_FILE_0} -g ${REFERENCE} --slow5 ${SIGNAL_FILE_0} --sam | head -n -1 | ${SAMTOOLS} sort -o ${EVENTALIGN_BAM_0}
	${SAMTOOLS} index ${EVENTALIGN_BAM_0} || die "samtools index failed"
	f5c eventalign -b ${MAPPED_BAM_0} -r ${SEQUENCE_FILE_0} -g ${REFERENCE} --slow5 ${SIGNAL_FILE_0} -o ${EVENTALIGN_TSV_0} --signal-index --scale-events

	f5c index ${SEQUENCE_FILE_1} --slow5 ${SIGNAL_FILE_1} || die "f5c index failed"
	f5c eventalign -b ${MAPPED_BAM_1} -r ${SEQUENCE_FILE_1} -g ${REFERENCE} --slow5 ${SIGNAL_FILE_1} --sam | head -n -1 | ${SAMTOOLS} sort -o ${EVENTALIGN_BAM_1}
	${SAMTOOLS} index ${EVENTALIGN_BAM_1} || die "samtools index failed"
	f5c eventalign -b ${MAPPED_BAM_1} -r ${SEQUENCE_FILE_1} -g ${REFERENCE} --slow5 ${SIGNAL_FILE_1} -o ${EVENTALIGN_TSV_1} --signal-index --scale-events

}


f5c_methylation() {
	info "f5c methylation..."
	
	${F5C} call-methylation -x hpc-low -b ${MAPPED_BAM_0} -r ${SEQUENCE_FILE_0}  -g ${REFERENCE} --slow5 ${SIGNAL_FILE_0}  > ${METH_TSV_0} || die "f5c methylation failed"
	${F5C} meth-freq -i ${METH_TSV_0} -s > ${METH_FREQ_TSV_0} || die "f5c meth-freq failed"
	tail -n +2 ${METH_FREQ_TSV_0} | awk '{print $1"\t"$2"\t"$3+1"\t"$7}' | sort -k1,1 -k2,2n > ${METH_FREQ_BEDGRAPH_0} || die "meth-freq to bedgraph conversion failed"
	awk -v TH="${METH_THRESHOLD}" '{if($4 > TH) print}' ${METH_FREQ_BEDGRAPH_0} | cut -f 1,2,3  > ${METH_BED_0}

	${F5C} call-methylation -x hpc-low -b ${MAPPED_BAM_1} -r ${SEQUENCE_FILE_1}  -g ${REFERENCE} --slow5 ${SIGNAL_FILE_1}  > ${METH_TSV_1} || die "f5c methylation failed"
	${F5C} meth-freq -i ${METH_TSV_1} -s > ${METH_FREQ_TSV_1} || die "f5c meth-freq failed"
	tail -n +2 ${METH_FREQ_TSV_1} | awk '{print $1"\t"$2"\t"$3+1"\t"$7}' | sort -k1,1 -k2,2n > ${METH_FREQ_BEDGRAPH_1} || die "meth-freq to bedgraph conversion failed"
	awk -v TH="${METH_THRESHOLD}" '{if($4 > TH) print}' ${METH_FREQ_BEDGRAPH_1} | cut -f 1,2,3 > ${METH_BED_1}

}

simulate_ref_signal() {
	info "simulating using ${SIMULATING_PROFILE}"
	
	test -d "$SIMULATE_REF_DIR" && rm -r "$SIMULATE_REF_DIR"
	mkdir "$SIMULATE_REF_DIR" || die "Failed creating $SIMULATE_REF_DIR"
	
	REGION_FILE="${SIMULATE_REF_DIR}/region.file"
	echo ${REF_REGION} > ${REGION_FILE}
	cat ${REGION_FILE}
	${SAMTOOLS} faidx ${REFERENCE} --region-file ${REGION_FILE} -o ${SIMULATE_REF_DIR}/ref.fasta || die "samtools faidx failed"
	sed -i "1s/.*/>sim_ref/" ${SIMULATE_REF_DIR}/ref.fasta

	${SQUIGULATOR} --seed 1 --full-contigs --ideal-time --amp-noise 0 -x ${SIMULATING_PROFILE} ${SIMULATE_REF_DIR}/ref.fasta -o ${SIMULATE_REF_DIR}/sim.slow5 -c ${SIMULATE_REF_DIR}/sim.paf -a ${SIMULATE_REF_DIR}/sim.sam -q ${SIMULATE_REF_DIR}/sim.fasta || die "squigulator failed"
	${SAMTOOLS} sort ${SIMULATE_REF_DIR}/sim.sam -o ${SIMULATE_REF_DIR}/sim.bam || die "samtools sort failed"
	${SAMTOOLS} index ${SIMULATE_REF_DIR}/sim.bam || die "samtools index failed"

}




plot_methylated_vs_non_methylated() {
	mkdir -p "${SQUIG_PLOT_DIR}" || die "Failed creating ${SQUIG_PLOT_DIR}"
	info "plotting eventalign and simulated signals"

	TRACK_COMMAND_FILE="${SQUIG_PLOT_DIR}/track_commands_${FUNCNAME[0]}.txt"
	rm -f ${TRACK_COMMAND_FILE}
	echo "num_commands=2" > ${TRACK_COMMAND_FILE}
	echo "plot_heights=*" >> ${TRACK_COMMAND_FILE}
	echo "squigualiser plot_pileup --bed ${METH_BED_0} --rna -f ${REFERENCE} -s ${SIGNAL_FILE_0} -a ${EVENTALIGN_BAM_0} --region ${REF_REGION} --tag_name methylated --plot_limit 6  --profile ${PROFILE_TO_DETERMINE_BASE_SHIFT} ${SIG_SCALE}" >> ${TRACK_COMMAND_FILE}
	echo "squigualiser plot_pileup --bed ${METH_BED_0} --rna -f ${REFERENCE} -s ${SIGNAL_FILE_1} -a ${EVENTALIGN_BAM_1} --region ${REF_REGION} --tag_name non_methylated --plot_limit 6  --profile ${PROFILE_TO_DETERMINE_BASE_SHIFT} ${SIG_SCALE}" >> ${TRACK_COMMAND_FILE}

	cat ${TRACK_COMMAND_FILE}

	TESTCASE="${MODEL_TO_USE}_methylated_vs_non_methylated"
	OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
	${PLOT_TRACK_TOOL} --shared_x -f ${TRACK_COMMAND_FILE} -o ${SQUIG_PLOT_DIR}/${TESTCASE} --tag_name ${TESTCASE} || die "testcase:$TESTCASE failed"

}

plot_methylated_vs_non_methylated_overlap_only() {
	mkdir -p "${SQUIG_PLOT_DIR}" || die "Failed creating ${SQUIG_PLOT_DIR}"
	info "plotting eventalign and simulated signals"

	TRACK_COMMAND_FILE="${SQUIG_PLOT_DIR}/track_commands_${FUNCNAME[0]}.txt"
	rm -f ${TRACK_COMMAND_FILE}
	echo "num_commands=2" > ${TRACK_COMMAND_FILE}
	echo "plot_heights=*" >> ${TRACK_COMMAND_FILE}
	echo "squigualiser plot_pileup --bed ${METH_BED_0} --rna -f ${REFERENCE} -s ${SIGNAL_FILE_0} -a ${EVENTALIGN_BAM_0} --region ${REF_REGION} --tag_name methylated --plot_limit 20  --profile ${PROFILE_TO_DETERMINE_BASE_SHIFT} ${SIG_SCALE} --overlap_only"  >> ${TRACK_COMMAND_FILE}
	echo "squigualiser plot_pileup --bed ${METH_BED_0} --rna -f ${REFERENCE} -s ${SIGNAL_FILE_1} -a ${EVENTALIGN_BAM_1} --region ${REF_REGION} --tag_name non_methylated --plot_limit 20  --profile ${PROFILE_TO_DETERMINE_BASE_SHIFT} ${SIG_SCALE} --overlap_only"   >> ${TRACK_COMMAND_FILE}

	cat ${TRACK_COMMAND_FILE}

	TESTCASE="${MODEL_TO_USE}_methylated_vs_non_methylated"
	OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
	${PLOT_TRACK_TOOL} --shared_x -f ${TRACK_COMMAND_FILE} -o ${SQUIG_PLOT_DIR}/${TESTCASE} --tag_name ${TESTCASE} || die "testcase:$TESTCASE failed"

}

## stage 1

# create_output_dir
# basecall_fast5
# minimap2_align
# f5c_eventalign
# f5c_methylation
# simulate_ref_signal
# plot_eventalign_and_sim
plot_methylated_vs_non_methylated_overlap_only

info "success"
