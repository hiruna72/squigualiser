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
OUTPUT_DIR="${REL_PATH}/data/out/plot_annotation"
test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"

RAW_DIR="${REL_PATH}/data/raw/plot"
EXP_DIR="${REL_PATH}/data/exp/plot"

testcase_40s() {
  TESTCASE=40.0
  info "testcase:$TESTCASE - read-signal plot"
  FASTA="-f ${RAW_DIR}/one_read/read.fasta"
  SIGNAL="-s ${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="-a ${RAW_DIR}/one_read/reform/r1k1m0.paf"
  OUTPUT="-o ${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT="--plot_limit 10"
  BED="--bed ${RAW_DIR}/bed_files/sample_1.bed"
  python src/plot.py  ${BED} --no_reverse ${FASTA} ${SIGNAL} ${ALIGNMENT} ${OUTPUT} --tag_name "testcase-${TESTCASE}" ${PLOT_LIMIT}|| die "testcase:$TESTCASE failed"

  TESTCASE=40.1
  info "testcase:$TESTCASE - read-signal plot"
  FASTA="-f ${RAW_DIR}/one_read/read.fasta"
  SIGNAL="-s ${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="-a ${RAW_DIR}/one_read/reform/r1k1m0.paf"
  OUTPUT="-o ${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT="--plot_limit 10"
  BED="--bed ${RAW_DIR}/bed_files/sample_1.bed"
  python src/plot.py  ${BED} --fixed_width --no_reverse ${FASTA} ${SIGNAL} ${ALIGNMENT} ${OUTPUT} --tag_name "testcase-${TESTCASE}" ${PLOT_LIMIT}|| die "testcase:$TESTCASE failed"

  TESTCASE=40.2
  GENOME="/media/hiruna/data/basecalling_work/apply_variants_to_genome/genome/hg38noAlt.fa"
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA="-f ${GENOME}"
  SIGNAL="-s ${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="-a ${RAW_DIR}/f5c_eventalign/sorted_eventalign.paf.gz"
  OUTPUT="-o ${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT="--plot_limit 10"
  BED="--bed ${RAW_DIR}/bed_files/sample_2.bed"
  REGION="--region chr1:6,811,011-6,811,198"
  python src/plot.py  ${BED} ${REGION} --fixed_width --no_reverse ${FASTA} ${SIGNAL} ${ALIGNMENT} ${OUTPUT} --tag_name "testcase-${TESTCASE}" ${PLOT_LIMIT}|| die "testcase:$TESTCASE failed"

}
testcase_40s #annotation

info "all testcases passed"
#rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0
