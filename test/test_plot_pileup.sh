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
OUTPUT_DIR="${REL_PATH}/data/out/plot_sig_ref_pileup"
test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"

RAW_DIR="${REL_PATH}/data/raw/plot"
EXP_DIR="${REL_PATH}/data/exp/plot"

testcase_11s() {
  GENOME="/media/hiruna/data/basecalling_work/apply_variants_to_genome/genome/hg38noAlt.fa"
  TESTCASE=11.0
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="${RAW_DIR}/f5c_eventalign/sorted_eventalign.paf.gz"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  SCALING="znorm"
  python src/plot_pileup.py --no_reverse -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} --sig_scale ${SCALING} && die "testcase:$TESTCASE failed"


  TESTCASE=11.1
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="${RAW_DIR}/f5c_eventalign/sorted_eventalign.paf.gz"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  REGION="chr1:6,811,011-6,811,198"
  SCALING="znorm"
  python src/plot_pileup.py --region ${REGION} --no_reverse -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} --sig_scale ${SCALING} || die "testcase:$TESTCASE failed"

  TESTCASE=11.2
  info "testcase:$TESTCASE - reference-signal plot"
  GENOME="test/data/raw/plot/reference_genomes/chr22_1_5k.fa"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/pileup_view/test_0/na12878_chr22.blow5"
  ALIGNMENT="${RAW_DIR}/pileup_view/test_0/na12878_chr22.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  REGION="chr22_1_5k:4000-4500"
  SCALING="znorm"
  python src/plot_pileup.py --region ${REGION} --no_reverse -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} --sig_scale ${SCALING} || die "testcase:$TESTCASE failed"

}
testcase_11s #pileup

info "all testcases passed"
#rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0
