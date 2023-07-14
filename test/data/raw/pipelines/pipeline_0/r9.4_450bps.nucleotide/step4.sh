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
[ "${RUN_NO}" ] || die "RUN_NO is not set"

R10_MODEL_FAST="dna_r10.4.1_e8.2_400bps_fast.cfg"
R10_MODEL_HAC="dna_r10.4.1_e8.2_400bps_hac.cfg"
R10_MODEL_SUP="dna_r10.4.1_e8.2_400bps_sup.cfg"

R9_MODEL_FAST="dna_r9.4.1_450bps_fast_prom.cfg"
R9_MODEL_HAC="dna_r9.4.1_450bps_hac_prom.cfg"
R9_MODEL_SUP="dna_r9.4.1_450bps_sup_prom.cfg"

# variable the user can change start here
REF="/media/hiruna/data/basecalling_work/apply_variants_to_genome/genome/hg38noAlt.fa"

MODEL_TO_USE=${R9_MODEL_SUP}
CHUNK_SIZE="--chunk_size 500"

KMER_MODEL="--kmer-model media/hiruna/data/basecalling_work/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model"
SIMULATING_PROFILE="dna-r9-prom"

RE_REFORM_K=4
RE_REFORM_M=3
BASE_SHIFT_FOR_MODEL_FORWARD="-2"
BASE_SHIFT_FOR_MODEL_REVERSE="-3"

REFORM_SCRIPT="/media/hiruna/data/basecalling_work/squigualiser/src/reform.py"
REALIGN_SCRIPT="/media/hiruna/data/basecalling_work/squigualiser/src/realign.py"
TRACK_SCRIPT="/media/hiruna/data/basecalling_work/squigualiser/src/plot_tracks.py"
CALCULATE_OFFSETS_SCRIPT="/media/hiruna/data/basecalling_work/squigualiser/src/calculate_offsets.py"

SQUIGULATOR=squigulator
SAMTOOLS=samtools
MINIMAP2=minimap2
GUPPY="/media/hiruna/data/slow5_work/guppy_integration/slow5-guppy_6.3.7/bin/guppy_basecaller"

SIGNAL_FILE="reads.blow5"

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

READ_ID="8e1a33c4-af69-471c-a115-6428c8bf63df"
READ_REGION=${READ_ID}:1-1000
REF_REGION="chr16:46,389,459-46,390,388"
# REF_REGION="chr1:92,781,684-92,782,120"
SIM_REGION="sim_ref:1-200"
SIG_SCALE="--sig_scale znorm"

[ "${OUTPUT_DIR}" ] || die "OUTPUT_DIR is not set"
# test -d "${OUTPUT_DIR}" && ask "${OUTPUT_DIR}" 
# test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
# mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"

# variable the user can change end here

[ "${SQUIGULATOR}" ] || die "edit the executable path of squigulator in the script"
[ "${SAMTOOLS}" ] || die "edit the executable path of samtools in the script"
[ "${MINIMAP2}" ] || die "edit the executable path of minimap2 in the script"
[ "${GUPPY}" ] || die "edit the executable path of guppy_basecaller in the script"

[ "${KMER_MODEL}" ] && [ ${SIMULATING_PROFILE} ] || die "set one of the variables. KMER_MODEL or SIMULATING_PROFILE"
[ "${READ_ID}" ] || die "READ_ID is not set"

[ "${BASECALL_DIR}" ] || die "BASECALL_DIR is not set"
[ "${SIMULATE_READ_DIR}" ] || die "SIMULATE_READ_DIR is not set"
[ "${SIMULATE_REF_DIR}" ] || die "SIMULATE_REF_DIR is not set"
[ "${SQUIG_PLOT_DIR}" ] || die "SQUIG_PLOT_DIR is not set"

[ "${RE_REFORM_K}" ] || die "RE_REFORM_K is not set"
[ "${RE_REFORM_M}" ] || die "RE_REFORM_M is not set"
[ "${BASE_SHIFT_FOR_MODEL_FORWARD}" ] || die "BASE_SHIFT_FOR_MODEL is not set"
[ "${BASE_SHIFT_FOR_MODEL_REVERSE}" ] || die "BASE_SHIFT_FOR_MODEL is not set"

basecall() {
	info "basecalling using ${MODEL_TO_USE}"
	test -d "$BASECALL_DIR" && rm -r "$BASECALL_DIR"
	mkdir "$BASECALL_DIR" || die "Failed creating $BASECALL_DIR"
	export CUDA_VISIBLE_DEVICES=0
	${GUPPY} -c ${MODEL_TO_USE} -i ${SIGNAL_FILE} --moves_out --bam_out --save_path ${BASECALL_DIR} --device cuda:0 ${CHUNK_SIZE} || die "basecalling failed"
}

samtools_merge() {
	info "samtools merging"
	test -d "$BASECALL_DIR/pass" || die "pass dir not found"
	cat ${BASECALL_DIR}/pass/*.fastq > ${SEQUENCE_FILE}
	${SAMTOOLS} merge -f ${BASECALL_DIR}/pass/*.bam -o ${MOVES_BAM}  || die "samtools merge failed"
}

reform_move_table() {
	info "running reform..."
	python ${REFORM_SCRIPT} -k 1 -m 0 --bam ${MOVES_BAM} -c -o ${REFORMAT_PAF} || die "reform failed"
}

calculate_offsets() {
	info "calculating offsets..."
	python ${CALCULATE_OFFSETS_SCRIPT} -p ${REFORMAT_PAF} -f ${SEQUENCE_FILE} -s ${SIGNAL_FILE} || die "calculate_offsets failed"
}

calculate_offsets_read_id() {
	info "calculating offsets for ${READ_ID}..."
	python ${CALCULATE_OFFSETS_SCRIPT} --read_id ${READ_ID} -o ${OUTPUT_DIR}/"${READ_ID}_${MODEL_TO_USE}.pdf" -p ${REFORMAT_PAF} -f ${SEQUENCE_FILE} -s ${SIGNAL_FILE} --tag_name ${MODEL_TO_USE} || die "calculate_offsets failed"
}

re_reform() {
	info "running re-reform..."
	python ${REFORM_SCRIPT} -k ${RE_REFORM_K} -m ${RE_REFORM_M} --bam ${MOVES_BAM} -c -o ${OUTPUT_DIR}/temp.paf || die "re reform failed"
	sort -k6,6 -k8,8n ${OUTPUT_DIR}/temp.paf -o ${RE_REFORMAT_PAF}
	bgzip -k ${RE_REFORMAT_PAF}
	tabix -0 -b 8 -e 9 -s 6 ${RE_REFORMAT_PAF_GZ}
}

simulate_read_signal() {
	info "simulating using ${KMER_MODEL}"

	test -d "$SIMULATE_READ_DIR" && rm -r "$SIMULATE_READ_DIR"
	mkdir "$SIMULATE_READ_DIR" || die "Failed creating $SIMULATE_READ_DIR"
	
	awk -v var="$READ_ID" 'BEGIN { print ">"var } $0 ~ var{getline; print}' ${SEQUENCE_FILE} > ${SIMULATE_READ_DIR}/ref.fasta
	${SQUIGULATOR} --seed 1 --full-contigs --ideal-time --amp-noise 0 -x ${SIMULATING_PROFILE} ${KMER_MODEL} ${SIMULATE_READ_DIR}/ref.fasta -o ${SIMULATE_READ_DIR}/sim.slow5 -c ${SIMULATE_READ_DIR}/sim.paf -a ${SIMULATE_READ_DIR}/sim.sam -q ${SIMULATE_READ_DIR}/sim.fasta || die "squigulator failed"
	${SAMTOOLS} sort ${SIMULATE_READ_DIR}/sim.sam -o ${SIMULATE_READ_DIR}/sim.bam || die "samtools sort failed"
	${SAMTOOLS} index ${SIMULATE_READ_DIR}/sim.bam || die "samtools index failed"

}

plot_sim_and_real_signal() {
	info "plotting simulated and real signals"
	
	mkdir -p "${SQUIG_PLOT_DIR}" || die "Failed creating ${SQUIG_PLOT_DIR}"
	
	TRACK_COMMAND_FILE=${SQUIG_PLOT_DIR}/track_commands_${FUNCNAME[0]}.txt
	echo "num_commands=2" > ${TRACK_COMMAND_FILE}
	echo "plot_heights=*" >> ${TRACK_COMMAND_FILE}
	echo "python plot_pileup.py -f ${SIMULATE_READ_DIR}/ref.fasta -s ${SIGNAL_FILE} -a ${RE_REFORMAT_PAF_GZ} --region ${READ_REGION} --tag_name ${MODEL_TO_USE}_reformed_read --no_overlap --plot_limit 1  --no_overlap" >> ${TRACK_COMMAND_FILE}
	echo "python plot_pileup.py -f ${SIMULATE_READ_DIR}/ref.fasta -s ${SIMULATE_READ_DIR}/sim.slow5 -a ${SIMULATE_READ_DIR}/sim.bam --region ${READ_ID}:1-1000 --tag_name ${MODEL_TO_USE}_sim_read --no_overlap --base_shift ${BASE_SHIFT_FOR_MODEL_FORWARD}" >> ${TRACK_COMMAND_FILE}

	cat ${TRACK_COMMAND_FILE}

	TESTCASE="${MODEL_TO_USE}_real_vs_sim"
	OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
	python ${TRACK_SCRIPT} -f ${TRACK_COMMAND_FILE} -o ${SQUIG_PLOT_DIR}/${TESTCASE} --tag_name ${TESTCASE} || die "testcase:$TESTCASE failed"

}

minimap2_align() {
	info "aligning reads using minimap2..."
	${MINIMAP2} -ax map-ont ${REF} -t32 --secondary=no ${SEQUENCE_FILE} | ${SAMTOOLS} sort -o ${MAPPED_BAM} || die "minimap2_align failed"
	${SAMTOOLS} index ${MAPPED_BAM} || die "samtools index failed"
}

realign() {
	info "re-aligning cigar mapping to moves..."
	python ${REALIGN_SCRIPT} --bam ${MAPPED_BAM} --paf ${RE_REFORMAT_PAF} -o ${REALIGN_BAM} || die "realign failed"
	${SAMTOOLS} index ${REALIGN_BAM} || die "samtools index failed"
}

simulate_ref_signal() {
	info "simulating using ${KMER_MODEL} ${SIMULATING_PROFILE}"
	
	test -d "$SIMULATE_REF_DIR" && rm -r "$SIMULATE_REF_DIR"
	mkdir "$SIMULATE_REF_DIR" || die "Failed creating $SIMULATE_REF_DIR"
	
	REGION_FILE="${SIMULATE_REF_DIR}/region.file"
	echo ${REF_REGION} > ${REGION_FILE}
	cat ${REGION_FILE}
	${SAMTOOLS} faidx ${REF} --region-file ${REGION_FILE} -o ${SIMULATE_REF_DIR}/ref.fasta || die "samtools faidx failed"
	sed -i "1s/.*/>sim_ref/" ${SIMULATE_REF_DIR}/ref.fasta

	${SQUIGULATOR} --seed 1 --full-contigs --ideal-time --amp-noise 0 -x ${SIMULATING_PROFILE} ${KMER_MODEL} ${SIMULATE_REF_DIR}/ref.fasta -o ${SIMULATE_REF_DIR}/sim.slow5 -c ${SIMULATE_REF_DIR}/sim.paf -a ${SIMULATE_REF_DIR}/sim.sam -q ${SIMULATE_REF_DIR}/sim.fasta || die "squigulator failed"
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
	echo "python plot_pileup.py -f ${REF} -s ${SIGNAL_FILE} -a ${REALIGN_BAM} --region ${REF_REGION} --tag_name ${MODEL_TO_USE}_realign --plot_limit 6 ${SIG_SCALE}" >> ${TRACK_COMMAND_FILE}
	echo "python plot_pileup.py -f ${SIMULATE_REF_DIR}/ref.fasta -s ${SIMULATE_REF_DIR}/sim.slow5 -a ${SIMULATE_REF_DIR}/sim.bam --region ${SIM_REGION} --tag_name ${MODEL_TO_USE}_sim_read --no_overlap --base_shift ${BASE_SHIFT_FOR_MODEL_FORWARD}" >> ${TRACK_COMMAND_FILE}

	cat ${TRACK_COMMAND_FILE}

	TESTCASE="${MODEL_TO_USE}_realigned_vs_sim"
	OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
	python ${TRACK_SCRIPT} --shared_x -f ${TRACK_COMMAND_FILE} -o ${SQUIG_PLOT_DIR}/${TESTCASE} --tag_name ${TESTCASE} || die "testcase:$TESTCASE failed"

}

f5c_eventalign() {
	info "f5c eventalign..."
	
	f5c index ${SEQUENCE_FILE} --slow5 ${SIGNAL_FILE} || die "f5c index failed"
	f5c eventalign -b ${MAPPED_BAM} -r ${SEQUENCE_FILE} -g ${REF} --slow5 ${SIGNAL_FILE} --sam | head -n -1 | ${SAMTOOLS} sort -o ${EVENTALIGN_BAM}
	${SAMTOOLS} index ${EVENTALIGN_BAM} || die "samtools index failed"

}

plot_eventalign_and_sim() {
	mkdir -p "${SQUIG_PLOT_DIR}" || die "Failed creating ${SQUIG_PLOT_DIR}"
	info "plotting eventalign and simulated signals"

	TRACK_COMMAND_FILE="${SQUIG_PLOT_DIR}/track_commands_${FUNCNAME[0]}.txt"
	rm -f ${TRACK_COMMAND_FILE}
	echo "num_commands=2" > ${TRACK_COMMAND_FILE}
	echo "plot_heights=*" >> ${TRACK_COMMAND_FILE}
	echo "python plot_pileup.py -f ${REF} -s ${SIGNAL_FILE} -a ${EVENTALIGN_BAM} --region ${REF_REGION} --tag_name ${MODEL_TO_USE}_eventalign --plot_limit 6  --base_shift ${BASE_SHIFT_FOR_MODEL_FORWARD} ${SIG_SCALE}" >> ${TRACK_COMMAND_FILE}
	echo "python plot_pileup.py -f ${SIMULATE_REF_DIR}/ref.fasta -s ${SIMULATE_REF_DIR}/sim.slow5 -a ${SIMULATE_REF_DIR}/sim.bam --region ${SIM_REGION} --tag_name ${MODEL_TO_USE}_sim_read --no_overlap --base_shift ${BASE_SHIFT_FOR_MODEL_FORWARD}" >> ${TRACK_COMMAND_FILE}

	cat ${TRACK_COMMAND_FILE}

	TESTCASE="${MODEL_TO_USE}_evligned_vs_sim"
	OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
	python ${TRACK_SCRIPT} --shared_x -f ${TRACK_COMMAND_FILE} -o ${SQUIG_PLOT_DIR}/${TESTCASE} --tag_name ${TESTCASE} || die "testcase:$TESTCASE failed"

}

plot_eventalign_forward_reverse() {
	mkdir -p "${SQUIG_PLOT_DIR}" || die "Failed creating ${SQUIG_PLOT_DIR}"
	info "plotting eventalign and simulated signals"

	TRACK_COMMAND_FILE=${SQUIG_PLOT_DIR}/track_commands_${FUNCNAME[0]}.txt
	rm -f ${TRACK_COMMAND_FILE}
	echo "num_commands=2" > ${TRACK_COMMAND_FILE}
	echo "plot_heights=*" >> ${TRACK_COMMAND_FILE}
	echo "python plot_pileup.py -f ${REF} -s ${SIGNAL_FILE} -a ${EVENTALIGN_BAM} --region ${REF_REGION} --tag_name ${MODEL_TO_USE}_eventalign_forward --plot_limit 6  --base_shift ${BASE_SHIFT_FOR_MODEL_FORWARD} ${SIG_SCALE}" >> ${TRACK_COMMAND_FILE}
	echo "python plot_pileup.py -f ${REF} -s ${SIGNAL_FILE} -a ${EVENTALIGN_BAM} --region ${REF_REGION} --tag_name ${MODEL_TO_USE}_eventalign_reverse --plot_limit 6  --base_shift ${BASE_SHIFT_FOR_MODEL_REVERSE} ${SIG_SCALE} --plot_reverse" >> ${TRACK_COMMAND_FILE}

	cat ${TRACK_COMMAND_FILE}

	TESTCASE="${MODEL_TO_USE}_evligned_forward_and_reverse"
	OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
	python ${TRACK_SCRIPT} --shared_x -f ${TRACK_COMMAND_FILE} -o ${SQUIG_PLOT_DIR}/${TESTCASE} --tag_name ${TESTCASE} || die "testcase:$TESTCASE failed"

}

plot_realign_forward_reverse() {
	mkdir -p "${SQUIG_PLOT_DIR}" || die "Failed creating ${SQUIG_PLOT_DIR}"
	info "plotting realigned and simulated signals"

	TRACK_COMMAND_FILE=${SQUIG_PLOT_DIR}/track_commands_${FUNCNAME[0]}.txt
	rm -f ${TRACK_COMMAND_FILE}
	echo "num_commands=2" > ${TRACK_COMMAND_FILE}
	echo "plot_heights=*" >> ${TRACK_COMMAND_FILE}
	echo "python plot_pileup.py -f ${REF} -s ${SIGNAL_FILE} -a ${REALIGN_BAM} --region ${REF_REGION} --tag_name ${MODEL_TO_USE}_realign_forward --plot_limit 6 ${SIG_SCALE}" >> ${TRACK_COMMAND_FILE}
	echo "python plot_pileup.py -f ${REF} -s ${SIGNAL_FILE} -a ${REALIGN_BAM} --region ${REF_REGION} --tag_name ${MODEL_TO_USE}_realign_reverse --plot_limit 6 ${SIG_SCALE} --plot_reverse" >> ${TRACK_COMMAND_FILE}

	cat ${TRACK_COMMAND_FILE}

	TESTCASE="${MODEL_TO_USE}_realigned_forward_and_reverse"
	OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
	python ${TRACK_SCRIPT} --shared_x -f ${TRACK_COMMAND_FILE} -o ${SQUIG_PLOT_DIR}/${TESTCASE} --tag_name ${TESTCASE} || die "testcase:$TESTCASE failed"

}

# basecall
# samtools_merge
# reform_move_table
# calculate_offsets
# calculate_offsets_read_id
# re_reform
# simulate_read_signal
# plot_sim_and_real_signal
# minimap2_align
# realign
# simulate_ref_signal
# plot_realign_and_sim
# f5c_eventalign
# plot_eventalign_and_sim
plot_eventalign_forward_reverse
# plot_realign_forward_reverse

info "success"
