#!/bin/bash

# ReadME first!!!
# This is a long pipline and might not work in the first go.
# Go all the way to the bottom and uncomment functions one by one, run and make sure all the commands work before moving to the next function.
# Details of the variables and functions can be found at https://github.com/hiruna72/squigualiser/blob/dev/docs/pipeline_basic.md
# Good luck!

RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo ; echo -e "${GREEN}$1${NC}" >&2 ; }
ask() { echo -n "Directory $1 exists. Delete and create again? (y/n)? " ; read answer ; if [ "$answer" != "${answer#[Nn]}" ] ;then exit ; fi ; echo ; }

info "$(date)"

# set -x

RUN_NO="RUN01"

R9_RNA_MODEL_FAST="rna_r9.4.1_70bps_fast_prom.cfg"
R9_RNA_MODEL_HAC="rna_r9.4.1_70bps_hac_prom.cfg"

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
XPORE_ENV_PATH=

REFERENCE="ref.fa"
MODEL_TO_USE=${R9_RNA_MODEL_HAC}
CHUNK_SIZE="--chunk_size 500"

SIMULATING_PROFILE="rna-r9-prom"

PROFILE_TO_DETERMINE_BASE_SHIFT="kmer_model_rna_r9.4.1_70bps_5_mer"

SIGNAL_FILE_0="HEK293T-METTL3-KO-rep1.blow5"
SIGNAL_FILE_1="HEK293T-WT-rep1.blow5"

CONFIG_YML="Hek293T_config.yml"

READ_ID="251256d6-1c5d-4756-91bf-8fa5dc6fd1c8"
READ_REGION=${READ_ID}:1-500
REF_REGION="ENST00000273480.3:578-687"
SIM_REGION="sim_ref:1-45"
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

DATAPREP_0="${OUTPUT_DIR}/data/0/dataprep"
DATAPREP_1="${OUTPUT_DIR}/data/1/dataprep"

METH_BED="${OUTPUT_DIR}/meth.bed"
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
	${MINIMAP2} -ax splice -uf -k14 ${REFERENCE} -t32 --secondary=no ${SEQUENCE_FILE_0} | ${SAMTOOLS} sort -o ${MAPPED_BAM_0} || die "minimap2_align failed"
	${SAMTOOLS} index ${MAPPED_BAM_0} || die "samtools index failed"

	${MINIMAP2} -ax splice -uf -k14 ${REFERENCE} -t32 --secondary=no ${SEQUENCE_FILE_1} | ${SAMTOOLS} sort -o ${MAPPED_BAM_1} || die "minimap2_align failed"
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

run_xpore() {
	info "running xpore..."
	
	cp ${CONFIG_YML} ${OUTPUT_DIR}/config.yml
	
	source ${XPORE_ENV_PATH}/bin/activate || die "could not activate xpore venv"
	
	xpore dataprep --eventalign ${EVENTALIGN_TSV_0} --transcript_fasta ${REFERENCE} --out_dir ${DATAPREP_0} || die "xpore dataprep failed"
	xpore dataprep --eventalign ${EVENTALIGN_TSV_1} --transcript_fasta ${REFERENCE} --out_dir ${DATAPREP_1} || die "xpore dataprep failed"
	
	cd ${OUTPUT_DIR}
	xpore diffmod --config config.yml || die "xpore diffmod failed"
	cd -
	
	deactivate

	awk -F',' '{print $1"\t"$2"\t"$2+1}' "${OUTPUT_DIR}/xpore_out/diffmod.table" | tail -n +2 > ${METH_BED} || die "creating meth bed file failed"

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




plot_eventalign_and_sim() {
	mkdir -p "${SQUIG_PLOT_DIR}" || die "Failed creating ${SQUIG_PLOT_DIR}"
	info "plotting eventalign and simulated signals"

	TRACK_COMMAND_FILE="${SQUIG_PLOT_DIR}/track_commands_${FUNCNAME[0]}.txt"
	rm -f ${TRACK_COMMAND_FILE}
	echo "num_commands=2" > ${TRACK_COMMAND_FILE}
	echo "plot_heights=*" >> ${TRACK_COMMAND_FILE}
	echo "squigualiser plot_pileup --bed ${METH_BED} --rna -f ${REFERENCE} -s ${SIGNAL_FILE_0} -a ${EVENTALIGN_BAM_0} --region ${REF_REGION} --tag_name HEK293T-METTL3-KO-rep1 --plot_limit 6  --profile ${PROFILE_TO_DETERMINE_BASE_SHIFT} ${SIG_SCALE}" >> ${TRACK_COMMAND_FILE}
	echo "squigualiser plot_pileup --bed ${METH_BED} --rna -f ${REFERENCE} -s ${SIGNAL_FILE_1} -a ${EVENTALIGN_BAM_1} --region ${REF_REGION} --tag_name HEK293T-WT-rep1 --plot_limit 6  --profile ${PROFILE_TO_DETERMINE_BASE_SHIFT} ${SIG_SCALE}" >> ${TRACK_COMMAND_FILE}

	cat ${TRACK_COMMAND_FILE}

	TESTCASE="${MODEL_TO_USE}_evligned_vs_sim"
	OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
	${PLOT_TRACK_TOOL} --shared_x -f ${TRACK_COMMAND_FILE} -o ${SQUIG_PLOT_DIR}/${TESTCASE} --tag_name ${TESTCASE} || die "testcase:$TESTCASE failed"

}

## stage 1

# create_output_dir
# basecall_fast5
# minimap2_align
# f5c_eventalign
# run_xpore
#plot_eventalign_and_sim
# simulate_ref_signal

info "success"
