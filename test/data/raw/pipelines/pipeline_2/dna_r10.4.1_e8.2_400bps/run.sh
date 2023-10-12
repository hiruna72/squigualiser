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

R10_MODEL_FAST="dna_r10.4.1_e8.2_400bps_fast.cfg"
R10_MODEL_HAC="dna_r10.4.1_e8.2_400bps_hac.cfg"
R10_MODEL_SUP="dna_r10.4.1_e8.2_400bps_sup.cfg"

REFORM_TOOL="squigualiser reform"
REALIGN_TOOL="squigualiser realign"
PLOT_TRACK_TOOL="squigualiser plot_tracks"
CALCULATE_OFFSETS_TOOL="squigualiser calculate_offsets"

## variable the user can change start here
SQUIGULATOR=squigulator
SAMTOOLS=samtools
SLOW5TOOLS=slow5tools
SAMTOOLS=samtools
MINIMAP2=minimap2
F5C=f5c
BGZIP=bgzip
TABIX=tabix

GUPPY="none"
BUTTERY_EEL_ENV_PATH="none"
REFERENCE=""

MODEL_TO_USE=${R10_MODEL_SUP}
CHUNK_SIZE="--chunk_size 500"

SIMULATING_PROFILE="dna-r10-prom"

RE_REFORM_K=1
RE_REFORM_M=0
PROFILE_TO_DETERMINE_BASE_SHIFT="kmer_model_dna_r10.4.1_e8.2_400bps_9_mer"

SIGNAL_FILE="reads.blow5"

READ_ID="251256d6-1c5d-4756-91bf-8fa5dc6fd1c8"
READ_REGION=${READ_ID}:1-500
REF_REGION="chr1:92778190-92778210"
SIM_REGION="sim_ref:1-20"
SIG_SCALE="--sig_scale znorm"

## variable the user can change end here

READ_ID_LIST="read_id.list"

OUTPUT_DIR=${RUN_NO}_${MODEL_TO_USE}
BASECALL_DIR="${OUTPUT_DIR}/basecall_${MODEL_TO_USE}"
SIMULATE_READ_DIR="${OUTPUT_DIR}/simulate_read"
SIMULATE_REF_DIR="${OUTPUT_DIR}/simulate_ref"
SQUIG_PLOT_DIR="${OUTPUT_DIR}/squig_plot_reads"

REFORMAT_PAF="${OUTPUT_DIR}/reform.paf"
RE_REFORMAT_PAF="${OUTPUT_DIR}/re_reform.paf"
RE_REFORMAT_PAF_GZ="${OUTPUT_DIR}/re_reform.paf.gz"
MOVES_BAM="${OUTPUT_DIR}/pass_moves.bam"
SEQUENCE_FILE="${OUTPUT_DIR}/pass.fastq"
MAPPED_BAM="${OUTPUT_DIR}/mapped.bam"
REALIGN_BAM="${OUTPUT_DIR}/realigned.bam"
EVENTALIGN_BAM="${OUTPUT_DIR}/eventalign.bam"
EVENTALIGN_PAF="${OUTPUT_DIR}/eventalign.paf"
EVENTALIGN_PAF_FILTERED="${OUTPUT_DIR}/eventalign_filtered.paf"
EVENTALIGN_PAF_FILTERED_GZ="${OUTPUT_DIR}/eventalign_filtered.paf.gz"
SORTED_EVENTALIGN_PAF="${OUTPUT_DIR}/sorted_eventalign.paf"

METH_TSV="${OUTPUT_DIR}/meth.tsv"
METH_FREQ_TSV="${OUTPUT_DIR}/meth-freq.tsv"
METH_FREQ_BEDGRAPH="${OUTPUT_DIR}/meth_freq.bedgraph"
METH_FREQ_BIGWIG="${OUTPUT_DIR}/meth_freq.bigwig"
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
[ "${BGZIP}" ] || die "edit the executable path of bgzip (variable BGZIP) in the script"
[ "${TABIX}" ] || die "edit the executable path of tabix (variable TABIX) in the script"

[ ${SIMULATING_PROFILE} ] || die "set the variable SIMULATING_PROFILE"
[ "${READ_ID}" ] || die "READ_ID is not set"

[ "${RUN_NO}" ] || die "RUN_NO is not set"
[ "${OUTPUT_DIR}" ] || die "OUTPUT_DIR is not set"
[ "${BASECALL_DIR}" ] || die "BASECALL_DIR is not set"
[ "${SIMULATE_READ_DIR}" ] || die "SIMULATE_READ_DIR is not set"
[ "${SIMULATE_REF_DIR}" ] || die "SIMULATE_REF_DIR is not set"
[ "${SQUIG_PLOT_DIR}" ] || die "SQUIG_PLOT_DIR is not set"

[ "${RE_REFORM_K}" ] || die "RE_REFORM_K is not set"
[ "${RE_REFORM_M}" ] || die "RE_REFORM_M is not set"
[ "${PROFILE_TO_DETERMINE_BASE_SHIFT}" ] || die "PROFILE_TO_DETERMINE_BASE_SHIFT is not set"

create_output_dir() {
	test -d "${OUTPUT_DIR}" && ask "${OUTPUT_DIR}" 
	test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
	mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"
}

basecall_eel() {
	info "basecalling using ${MODEL_TO_USE}"
	test -d "${BASECALL_DIR}" && rm -r "${BASECALL_DIR}"
	mkdir "${BASECALL_DIR}" || die "Failed creating ${BASECALL_DIR}"
	export CUDA_VISIBLE_DEVICES=0
	source ${BUTTERY_EEL_ENV_PATH}/bin/activate || die "could not activate buttery-eel venv"
	CUDA_LAUNCH_BLOCKING=1 buttery-eel ${CHUNK_SIZE} -x cuda:0 -g ${GUPPY} --config ${MODEL_TO_USE} -i ${SIGNAL_FILE} -o ${BASECALL_DIR}/moves.sam --log ${OUTPUT_DIR} --moves_out --port 5558 --use_tcp || die "buttery-eel failed"
	${SAMTOOLS} view ${BASECALL_DIR}/moves.sam -o ${MOVES_BAM}  || die "samtools view failed"
	${SAMTOOLS} fastq ${MOVES_BAM} > ${SEQUENCE_FILE}  || die "samtools fastq failed"
	${SAMTOOLS} index ${MOVES_BAM}  || die "samtools index failed"
	deactivate
}

basecall_fast5() {
	set -x
	info "basecalling using ${MODEL_TO_USE}"
	test -d "${BASECALL_DIR}" && rm -r "${BASECALL_DIR}"
	mkdir "${BASECALL_DIR}" || die "Failed creating ${BASECALL_DIR}"
	${SLOW5TOOLS} s2f ${SIGNAL_FILE} -o ${OUTPUT_DIR}/reads.fast5
	export CUDA_VISIBLE_DEVICES=0
	${GUPPY}/guppy_basecaller ${CHUNK_SIZE} -x cuda:0 --config ${MODEL_TO_USE} -i ${OUTPUT_DIR} -s ${BASECALL_DIR} --moves_out --bam_out || die "guppy_basecaller failed"
	cat ${BASECALL_DIR}/pass/*.fastq > ${SEQUENCE_FILE} || die "cat fastq together failed"
	${SAMTOOLS} merge -f ${BASECALL_DIR}/pass/*.bam -o ${MOVES_BAM}  || die "samtools merge failed"
}

reform_move_table() {
	info "running reform..."
	${REFORM_TOOL} -k 1 -m 0 --bam ${MOVES_BAM} -c -o ${REFORMAT_PAF} || die "reform failed"
}

calculate_offsets() {
	info "calculating offsets..."
	${CALCULATE_OFFSETS_TOOL} -p ${REFORMAT_PAF} -f ${SEQUENCE_FILE} -s ${SIGNAL_FILE} || die "calculate_offsets failed"
}

calculate_offsets_read_id() {
	info "calculating offsets for ${READ_ID}..."
	${CALCULATE_OFFSETS_TOOL} --read_id ${READ_ID} -o ${OUTPUT_DIR}/"${READ_ID}_${MODEL_TO_USE}.pdf" -p ${REFORMAT_PAF} -f ${SEQUENCE_FILE} -s ${SIGNAL_FILE} --tag_name ${MODEL_TO_USE} || die "calculate_offsets failed"
}

re_reform() {
	info "running re-reform..."
	${REFORM_TOOL} -k ${RE_REFORM_K} -m ${RE_REFORM_M} --bam ${MOVES_BAM} -c -o ${OUTPUT_DIR}/temp.paf || die "re reform failed"
	sort -k6,6 -k8,8n ${OUTPUT_DIR}/temp.paf -o ${RE_REFORMAT_PAF}
	cp ${RE_REFORMAT_PAF} "${RE_REFORMAT_PAF}.copy" || die "copy failed"
	${BGZIP} ${RE_REFORMAT_PAF} || die "bgzip failed"
	${TABIX} -0 -b 8 -e 9 -s 6 ${RE_REFORMAT_PAF_GZ} || die "tabix failed"
	cp "${RE_REFORMAT_PAF}.copy" ${RE_REFORMAT_PAF} || die "copy failed"
}

simulate_read_signal() {
	info "simulating using ${SIMULATING_PROFILE}"

	test -d "$SIMULATE_READ_DIR" && rm -r "$SIMULATE_READ_DIR"
	mkdir "$SIMULATE_READ_DIR" || die "Failed creating $SIMULATE_READ_DIR"
	
	awk -v var="$READ_ID" 'BEGIN { print ">"var } $0 ~ var{getline; print}' ${SEQUENCE_FILE} > ${SIMULATE_READ_DIR}/ref.fasta
	${SQUIGULATOR} --seed 1 --full-contigs --ideal-time --amp-noise 0 -x ${SIMULATING_PROFILE} ${SIMULATE_READ_DIR}/ref.fasta -o ${SIMULATE_READ_DIR}/sim.slow5 -c ${SIMULATE_READ_DIR}/sim.paf -a ${SIMULATE_READ_DIR}/sim.sam -q ${SIMULATE_READ_DIR}/sim.fasta || die "squigulator failed"
	${SAMTOOLS} sort ${SIMULATE_READ_DIR}/sim.sam -o ${SIMULATE_READ_DIR}/sim.bam || die "samtools sort failed"
	${SAMTOOLS} index ${SIMULATE_READ_DIR}/sim.bam || die "samtools index failed"

}

plot_sim_and_real_signal() {
	info "plotting simulated and real signals"
	
	mkdir -p "${SQUIG_PLOT_DIR}" || die "Failed creating ${SQUIG_PLOT_DIR}"
	
	TRACK_COMMAND_FILE=${SQUIG_PLOT_DIR}/track_commands_${FUNCNAME[0]}.txt
	echo "num_commands=2" > ${TRACK_COMMAND_FILE}
	echo "plot_heights=*" >> ${TRACK_COMMAND_FILE}
	echo "squigualiser plot_pileup -f ${SIMULATE_READ_DIR}/ref.fasta -s ${SIGNAL_FILE} -a ${RE_REFORMAT_PAF_GZ} --region ${READ_REGION} --tag_name ${MODEL_TO_USE}_reformed_read --no_overlap --plot_limit 1  --no_overlap" >> ${TRACK_COMMAND_FILE}
	echo "squigualiser plot_pileup -f ${SIMULATE_READ_DIR}/ref.fasta -s ${SIMULATE_READ_DIR}/sim.slow5 -a ${SIMULATE_READ_DIR}/sim.bam --region ${READ_ID}:1-1000 --tag_name ${MODEL_TO_USE}_sim_read --no_overlap --profile ${PROFILE_TO_DETERMINE_BASE_SHIFT}" >> ${TRACK_COMMAND_FILE}

	cat ${TRACK_COMMAND_FILE}

	TESTCASE="${MODEL_TO_USE}_real_vs_sim"
	OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
	${PLOT_TRACK_TOOL} --shared_x -f ${TRACK_COMMAND_FILE} -o ${SQUIG_PLOT_DIR}/${TESTCASE} --tag_name ${TESTCASE} || die "testcase:$TESTCASE failed"

}

minimap2_align() {
	info "aligning reads using minimap2..."
	${MINIMAP2} -ax map-ont ${REFERENCE} -t32 --secondary=no ${SEQUENCE_FILE} | ${SAMTOOLS} sort -o ${MAPPED_BAM} || die "minimap2_align failed"
	${SAMTOOLS} index ${MAPPED_BAM} || die "samtools index failed"
}

realign() {
	info "re-aligning cigar mapping to moves..."
	${REALIGN_TOOL} --bam ${MAPPED_BAM} --paf ${RE_REFORMAT_PAF} -o ${REALIGN_BAM} || die "realign failed"
	${SAMTOOLS} index ${REALIGN_BAM} || die "samtools index failed"
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

plot_realign_and_sim() {
	mkdir -p "${SQUIG_PLOT_DIR}" || die "Failed creating ${SQUIG_PLOT_DIR}"
	info "plotting realigned and simulated signals"

	TRACK_COMMAND_FILE=${SQUIG_PLOT_DIR}/track_commands_${FUNCNAME[0]}.txt
	rm -f ${TRACK_COMMAND_FILE}
	echo "num_commands=2" > ${TRACK_COMMAND_FILE}
	echo "plot_heights=*" >> ${TRACK_COMMAND_FILE}
	echo "squigualiser plot_pileup -f ${REFERENCE} -s ${SIGNAL_FILE} -a ${REALIGN_BAM} --region ${REF_REGION} --tag_name ${MODEL_TO_USE}_realign --plot_limit 6 ${SIG_SCALE}" >> ${TRACK_COMMAND_FILE}
	echo "squigualiser plot_pileup -f ${SIMULATE_REF_DIR}/ref.fasta -s ${SIMULATE_REF_DIR}/sim.slow5 -a ${SIMULATE_REF_DIR}/sim.bam --region ${SIM_REGION} --tag_name ${MODEL_TO_USE}_sim_read --no_overlap --profile ${PROFILE_TO_DETERMINE_BASE_SHIFT}" >> ${TRACK_COMMAND_FILE}

	cat ${TRACK_COMMAND_FILE}

	TESTCASE="${MODEL_TO_USE}_realigned_vs_sim"
	OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
	${PLOT_TRACK_TOOL} --shared_x -f ${TRACK_COMMAND_FILE} -o ${SQUIG_PLOT_DIR}/${TESTCASE} --tag_name ${TESTCASE} || die "testcase:$TESTCASE failed"

}

f5c_eventalign() {
	info "f5c eventalign..."
	
	f5c index ${SEQUENCE_FILE} --slow5 ${SIGNAL_FILE} || die "f5c index failed"
	# f5c eventalign -b ${MAPPED_BAM} -r ${SEQUENCE_FILE} -g ${REFERENCE} --slow5 ${SIGNAL_FILE} --sam | head -n -1 | ${SAMTOOLS} sort -o ${EVENTALIGN_BAM}
	# ${SAMTOOLS} index ${EVENTALIGN_BAM} || die "samtools index failed"

	f5c eventalign -b ${MAPPED_BAM} -r ${SEQUENCE_FILE} -g ${REFERENCE} --slow5 ${SIGNAL_FILE} -c -o ${EVENTALIGN_PAF}

}

f5c_methylation() {
	info "f5c methylation..."
	
	${F5C} call-methylation -x hpc-low -b ${MAPPED_BAM} -r ${SEQUENCE_FILE}  -g ${REFERENCE} --slow5 ${SIGNAL_FILE}  > ${METH_TSV} || die "f5c methylation failed"
	${F5C} meth-freq -i ${METH_TSV} -s > ${METH_FREQ_TSV} || die "f5c meth-freq failed"

	tail -n +2 ${METH_FREQ_TSV} | awk '{print $1"\t"$2"\t"$3+1"\t"$7}' | sort -k1,1 -k2,2n > ${METH_FREQ_BEDGRAPH} || die "meth-freq to bedgraph conversion failed"
	awk -v TH="${METH_THRESHOLD}" '{if($4 > TH) print}' ${METH_FREQ_BEDGRAPH} > ${METH_BED}

}

plot_eventalign_and_sim() {
	mkdir -p "${SQUIG_PLOT_DIR}" || die "Failed creating ${SQUIG_PLOT_DIR}"
	info "plotting eventalign and simulated signals"

	sort -k6,6 -k8,8n ${EVENTALIGN_PAF} -o ${SORTED_EVENTALIGN_PAF}
	cat ${SORTED_EVENTALIGN_PAF} | grep -f ${READ_ID_LIST} > ${EVENTALIGN_PAF_FILTERED}
	bgzip -f ${EVENTALIGN_PAF_FILTERED}
	tabix -0 -b 8 -e 9 -s 6 ${EVENTALIGN_PAF_FILTERED_GZ}

	TRACK_COMMAND_FILE="${SQUIG_PLOT_DIR}/track_commands_${FUNCNAME[0]}.txt"
	rm -f ${TRACK_COMMAND_FILE}
	echo "num_commands=2" > ${TRACK_COMMAND_FILE}
	echo "plot_heights=*" >> ${TRACK_COMMAND_FILE}
	# echo "squigualiser plot_pileup --bed ${METH_BED} -f ${REFERENCE} -s ${SIGNAL_FILE} -a ${EVENTALIGN_PAF_FILTERED_GZ} --region ${REF_REGION} --tag_name ${MODEL_TO_USE}_eventalign --plot_limit 20  --base_shift -5 ${SIG_SCALE}" >> ${TRACK_COMMAND_FILE}
	echo "squigualiser plot_pileup --bed ${METH_BED} -f ${REFERENCE} -s ${SIGNAL_FILE} -a ${EVENTALIGN_PAF_FILTERED_GZ} --region ${REF_REGION} --tag_name ${MODEL_TO_USE}_eventalign --plot_limit 20  --profile ${PROFILE_TO_DETERMINE_BASE_SHIFT} ${SIG_SCALE}" >> ${TRACK_COMMAND_FILE}
	echo "squigualiser plot_pileup -f ${SIMULATE_REF_DIR}/ref.fasta -s ${SIMULATE_REF_DIR}/sim.slow5 -a ${SIMULATE_REF_DIR}/sim.bam --region ${SIM_REGION} --tag_name ${MODEL_TO_USE}_sim_read --no_overlap --profile ${PROFILE_TO_DETERMINE_BASE_SHIFT}" >> ${TRACK_COMMAND_FILE}

	cat ${TRACK_COMMAND_FILE}

	TESTCASE="${MODEL_TO_USE}_evligned_vs_sim"
	OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
	${PLOT_TRACK_TOOL} --shared_x -f ${TRACK_COMMAND_FILE} -o ${SQUIG_PLOT_DIR}/${TESTCASE} --tag_name ${TESTCASE} || die "testcase:$TESTCASE failed"

}

## stage 1

# create_output_dir
# if there are problems with basecall_eel basecalling use basecall_fast5 basecalling
# basecall_eel
# basecall_fast5
# reform_move_table

# stage 2

# calculate_offsets
# calculate_offsets_read_id
# re_reform

# stage 3

# simulate_read_signal
# plot_sim_and_real_signal
# minimap2_align
# realign
#simulate_ref_signal
# plot_realign_and_sim
# f5c_eventalign
# f5c_methylation
#plot_eventalign_and_sim

info "success"
