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

# Relative path to "slow5/tests/"
REL_PATH="$(dirname $0)/"
#...directories files tools arguments commands clean
OUTPUT_DIR="${REL_PATH}/data/out/metric"
test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"
#commands ...

ex() {
  python src/metric.py "$@"
}
RAW_DIR="${REL_PATH}/data/raw/pipelines/pipeline_0/dna_r10.4.1_e8.2_400bps"
#EXP_DIR="${REL_PATH}/data/exp/reform"

[ "${HUMAN_GENOME}" ] || die "Please set the env variable to the human genome path. export HUMAN_GENOME=path/to/file"
GENOME="${HUMAN_GENOME}"
REGION="--region chr1:92,783,745-92,783,946"
PROFILE="--profile guppy_dna_r10.4.1_e8.2_400bps_sup"
TESTCASE=50.1
info "testcase:$TESTCASE - help"
ex && die "testcase:$TESTCASE failed"

TESTCASE=50.2
info "testcase:$TESTCASE"
ex -f ${HUMAN_GENOME} -s ${RAW_DIR}/reads.blow5 ${REGION} ${PROFILE} -a ${RAW_DIR}/nanopolish_signal_projection/realigned.bam -o ${OUTPUT_DIR}/realign_stat.tsv || die "testcase:$TESTCASE failed"

TESTCASE=50.3
info "testcase:$TESTCASE"
ex -f ${HUMAN_GENOME} -s ${RAW_DIR}/reads.blow5 ${REGION} ${PROFILE} -a ${RAW_DIR}/nanopolish_signal_projection/eventalign.bam -o ${OUTPUT_DIR}/eventalign_stat.tsv || die "testcase:$TESTCASE failed"

TESTCASE=50.4
info "testcase:$TESTCASE"
ex -f ${HUMAN_GENOME} -s ${RAW_DIR}/reads.blow5 ${REGION} ${PROFILE} -a ${RAW_DIR}/nanopolish_signal_projection/sorted_projected.paf.gz -o ${OUTPUT_DIR}/naopolish_projection_stat.tsv || die "testcase:$TESTCASE failed"

TESTCASE=50.5
info "testcase:$TESTCASE"
ex --extend_0 -f ${HUMAN_GENOME} -s ${RAW_DIR}/reads.blow5 ${REGION} ${PROFILE} -a ${RAW_DIR}/nanopolish_signal_projection/realigned.bam -o ${OUTPUT_DIR}/realign_stat_ex0.tsv || die "testcase:$TESTCASE failed"

TESTCASE=50.6
info "testcase:$TESTCASE"
ex --extend_0 -f ${HUMAN_GENOME} -s ${RAW_DIR}/reads.blow5 ${REGION} ${PROFILE} -a ${RAW_DIR}/nanopolish_signal_projection/eventalign.bam -o ${OUTPUT_DIR}/eventalign_stat_ex0.tsv || die "testcase:$TESTCASE failed"

TESTCASE=50.7
info "testcase:$TESTCASE"
ex --extend_0 -f ${HUMAN_GENOME} -s ${RAW_DIR}/reads.blow5 ${REGION} ${PROFILE} -a ${RAW_DIR}/nanopolish_signal_projection/sorted_projected.paf.gz -o ${OUTPUT_DIR}/naopolish_projection_stat_ex0.tsv || die "testcase:$TESTCASE failed"

TESTCASE=50.8
info "testcase:$TESTCASE"
ex --extend_0 --plot_reverse -f ${HUMAN_GENOME} -s ${RAW_DIR}/reads.blow5 ${REGION} ${PROFILE} -a ${RAW_DIR}/nanopolish_signal_projection/realigned.bam -o ${OUTPUT_DIR}/realign_stat_rv.tsv || die "testcase:$TESTCASE failed"

TESTCASE=50.9
info "testcase:$TESTCASE"
ex --extend_0 --plot_reverse -f ${HUMAN_GENOME} -s ${RAW_DIR}/reads.blow5 ${REGION} ${PROFILE} -a ${RAW_DIR}/nanopolish_signal_projection/eventalign.bam -o ${OUTPUT_DIR}/eventalign_stat_rv.tsv || die "testcase:$TESTCASE failed"

TESTCASE=50.10
info "testcase:$TESTCASE"
ex --extend_0 --plot_reverse -f ${HUMAN_GENOME} -s ${RAW_DIR}/reads.blow5 ${REGION} ${PROFILE} -a ${RAW_DIR}/nanopolish_signal_projection/sorted_projected.paf.gz -o ${OUTPUT_DIR}/naopolish_projection_stat_rv.tsv || die "testcase:$TESTCASE failed"

TESTCASE=50.11
info "testcase:$TESTCASE"
ex --extend_0 --extend_1 --plot_reverse -f ${HUMAN_GENOME} -s ${RAW_DIR}/reads.blow5 ${REGION} ${PROFILE} -a ${RAW_DIR}/nanopolish_signal_projection/sorted_projected.paf.gz -o ${OUTPUT_DIR}/naopolish_projection_stat_rv.tsv || die "testcase:$TESTCASE failed"

#rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0
