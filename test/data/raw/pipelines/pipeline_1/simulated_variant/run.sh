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

PLOT_PILEUP_TOOL="squigualiser plot_pileup"

## variable the user can change start here
SQUIGULATOR=squigulator
SAMTOOLS=samtools
SLOW5TOOLS=slow5tools

REFERENCE="ref.fasta"

SIMULATING_PROFILE="dna-r10-prom"

PROFILE_TO_DETERMINE_BASE_SHIFT="kmer_model_dna_r10.4.1_e8.2_400bps_9_mer"

SIGNAL_FILE="sim.blow5"

SIM_REGION="sim_ref:300-500"

## variable the user can change end here

OUTPUT_DIR=${RUN_NO}
SIMULATE_READ_DIR="${OUTPUT_DIR}/simulate_read"
SQUIG_PLOT_DIR="${OUTPUT_DIR}/squig_plot_reads"

[ "${REFERENCE}" ] || die "edit the reference genome path (variable REF) in the script"

[ "${SQUIGULATOR}" ] || die "edit the executable path of squigulator (variable SQUIGULATOR) in the script"
[ "${SAMTOOLS}" ] || die "edit the executable path of samtools (variable SAMTOOLS) in the script"
[ "${SLOW5TOOLS}" ] || die "edit the executable path of slow5tools (variable SLOW5TOOLS) in the script"

[ ${SIMULATING_PROFILE} ] || die "set the variable SIMULATING_PROFILE"

[ "${RUN_NO}" ] || die "RUN_NO is not set"
[ "${OUTPUT_DIR}" ] || die "OUTPUT_DIR is not set"
[ "${SIMULATE_READ_DIR}" ] || die "SIMULATE_READ_DIR is not set"
[ "${SQUIG_PLOT_DIR}" ] || die "SQUIG_PLOT_DIR is not set"

[ "${PROFILE_TO_DETERMINE_BASE_SHIFT}" ] || die "PROFILE_TO_DETERMINE_BASE_SHIFT is not set"

create_output_dir() {
	test -d "${OUTPUT_DIR}" && ask "${OUTPUT_DIR}" 
	test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
	mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"
}

simulate_read_signal() {
	info "simulating using ${SIMULATING_PROFILE}"
	
	test -d "$SIMULATE_READ_DIR" && rm -r "$SIMULATE_READ_DIR"
	mkdir "$SIMULATE_READ_DIR" || die "Failed creating $SIMULATE_READ_DIR"
	
	${SQUIGULATOR} --seed 1 -n 100 --ideal-time --amp-noise 0 -x ${SIMULATING_PROFILE} ref.fasta -o ${SIMULATE_READ_DIR}/sim.slow5 -a ${SIMULATE_READ_DIR}/sim.sam || die "squigulator failed"
	${SQUIGULATOR} --seed 499 -n 100 --ideal-time --amp-noise 0 -x ${SIMULATING_PROFILE} ref_variant_added.fasta -o ${SIMULATE_READ_DIR}/sim_variant_added.slow5 -a ${SIMULATE_READ_DIR}/sim_variant_added.sam || die "squigulator failed"
	
	${SAMTOOLS} merge ${SIMULATE_READ_DIR}/sim.sam ${SIMULATE_READ_DIR}/sim_variant_added.sam -o ${SIMULATE_READ_DIR}/sim_merged.bam || die "samtools merege failed"
	${SAMTOOLS} sort ${SIMULATE_READ_DIR}/sim_merged.bam -o ${SIMULATE_READ_DIR}/sim_sorted.bam || die "samtools sort failed"
	${SAMTOOLS} index ${SIMULATE_READ_DIR}/sim_sorted.bam || die "samtools index failed"

	${SLOW5TOOLS} merge ${SIMULATE_READ_DIR}/sim.slow5 ${SIMULATE_READ_DIR}/sim_variant_added.slow5 -o ${SIMULATE_READ_DIR}/sim_merged.slow5 || die "slow5tools merge failed"
	${SLOW5TOOLS} index ${SIMULATE_READ_DIR}/sim_merged.slow5 || die "slow5tools index failed"
}

plot_sim_pileup() {
	set -x
	mkdir -p "${SQUIG_PLOT_DIR}" || die "Failed creating ${SQUIG_PLOT_DIR}"
	info "plotting simulated signal piluep"

	TESTCASE="sim_pileup_forward_mapped"
	OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
	${PLOT_PILEUP_TOOL} -f ref.fasta -s ${SIMULATE_READ_DIR}/sim_merged.slow5 -a ${SIMULATE_READ_DIR}/sim_sorted.bam --region ${SIM_REGION} --tag_name ${TESTCASE} --profile ${PROFILE_TO_DETERMINE_BASE_SHIFT} -o ${OUTPUT} || die "simulated plot pileup failed"
	TESTCASE="sim_pileup_reverse_mapped"
	${PLOT_PILEUP_TOOL} -f ref.fasta -s ${SIMULATE_READ_DIR}/sim_merged.slow5 -a ${SIMULATE_READ_DIR}/sim_sorted.bam --region ${SIM_REGION} --tag_name ${TESTCASE} --profile ${PROFILE_TO_DETERMINE_BASE_SHIFT} -o ${OUTPUT} --plot_reverse || die "simulated plot pileup failed"

}


# create_output_dir
simulate_read_signal
plot_sim_pileup

info "success"
