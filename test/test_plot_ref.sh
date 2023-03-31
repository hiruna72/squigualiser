#!/bin/bash

RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo ; echo -e "${GREEN}$1${NC}" >&2 ; }

#redirect
verbose=1
exec 3>&1
exec 4>&2
if ((verbose)); then
  echo "verbose=1"
else
  echo "verbose=0"
  exec 1>/dev/null
  exec 2>/dev/null
fi
#echo "this should be seen if verbose"
#echo "this should always be seen" 1>&3 2>&4

REL_PATH="$(dirname $0)/"
#...directories files tools arguments commands clean
OUTPUT_DIR="${REL_PATH}/data/out/plot"
test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"

GENOME="${REL_PATH}/data/raw/plot/reference_genomes/nCoV-2019.reference.fasta"
RAW_DIR="${REL_PATH}/data/raw/plot"
EXP_DIR="${REL_PATH}/data/exp/plot"

testcase_5s() {
  TESTCASE=5.1
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r2/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r2/sim.sam"
  REGION=""
  PLOT_LIMIT=1
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} || die "testcase:$TESTCASE failed"

  TESTCASE=5.2
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r2/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r2/sim.sam"
  REGION="MN908947.3:1-10000"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" && die "testcase:$TESTCASE failed"

  TESTCASE=5.3
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r2/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r2/sim.bam"
  REGION="MN908947.3:1-10000"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" && die "testcase:$TESTCASE failed"

  TESTCASE=5.4
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r2/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r2/sorted_sim.bam"
  REGION="MN908947.3:14,843-14,914"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup2.html

  TESTCASE=5.5
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r2/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r2/sorted_sim.bam"
  REGION="MN908947.3:14,843-14,914"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" --no_reverse|| die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup2.html

  TESTCASE=5.6
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r2/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r2/sorted_sim.bam"
  REGION="MN908947.3:14,843-14,914"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" --reverse_only|| die "testcase:$TESTCASE failed"

  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup2.html
}
testcase_6s() {
  TESTCASE=6.1
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r3/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r3/sorted_sim.bam"
  REGION="MN908947.3:14,843-14,914"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup2.html

  TESTCASE=6.2
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r3/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r3/sorted_sim.bam"
  REGION="MN908947.3:14,843-14,914"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" --no_reverse|| die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup2.html

  TESTCASE=6.3
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r3/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r3/sorted_sim.bam"
  REGION="MN908947.3:14,843-14,914"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" --reverse_only|| die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup2.html
}

testcase_5s #signal-reference squigulator ideal signals
testcase_6s #signal-reference squigulator

info "all testcases passed"
rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0
