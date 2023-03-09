#!/bin/bash

RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo ; echo -e "${GREEN}$1${NC}" >&2 ; }

#getting bash arguments
Usage="./simulate_read.sh reference region_file kmer_model output_dir"
NUMBER_ARGS=4
if [[ "$#" -lt ${NUMBER_ARGS} ]]; then
	info "Usage: ${Usage}"
	exit 1
fi
#set -x
REF=${1}
REGION_FILE=${2}
KMER_MODEL=${3}
OUTPUT_DIR=${4}

SAMTOOLS=samtools
SQUIGULATOR=squigulator

[ "${OUTPUT_DIR}" ] || die "OUTPUT_DIR is empty"
test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"

${SAMTOOLS} faidx ${REF} --region-file ${REGION_FILE} -o ${OUTPUT_DIR}/read.fasta || die "samtools faidx failed"

${SQUIGULATOR} --seed 1 --full-contigs --ideal-time --amp-noise 0.4 -x dna-r10-prom --kmer-model ${KMER_MODEL} ${OUTPUT_DIR}/read.fasta -o ${OUTPUT_DIR}/sim.slow5 -q ${OUTPUT_DIR}/sim.fasta -c ${OUTPUT_DIR}/sim.paf || die "squigulator failed"


info "success"
exit 0
