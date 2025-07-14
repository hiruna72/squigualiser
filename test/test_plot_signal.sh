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

SAVE_SVG=""
if [[ "$1" == "--save_svg" ]]; then
  SAVE_SVG="--save_svg"
  shift
fi

REL_PATH="$(dirname $0)/"
#...directories files tools arguments commands clean
OUTPUT_DIR="${REL_PATH}/data/out/plot_signal"
test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"

RAW_DIR="${REL_PATH}/data/raw/plot"

testcase_70s() {
  TESTCASE=70.1
  info "testcase:$TESTCASE - help"
  squigualiser plot_signal && die "testcase:$TESTCASE failed"

  TESTCASE=70.2
  info "testcase:$TESTCASE - help"
  squigualiser plot_signal --help || die "testcase:$TESTCASE failed"

  RAW_DIR="${REL_PATH}/data/raw/plot/f5c_eventalign"
  TESTCASE=70.3
  SIGNAL="${RAW_DIR}/reads.blow5"
  REGION=""
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  squigualiser plot_signal -s ${SIGNAL} -o ${OUTPUT_DIR} &&  die "testcase:$TESTCASE failed"

  RAW_DIR="${REL_PATH}/data/raw/plot/f5c_eventalign"
  info "testcase:$TESTCASE - read-signal plot"
  TESTCASE=70.4
  SIGNAL="${RAW_DIR}/reads.blow5"
  REGION=""
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  READ_ID="64a25d50-50e0-41f8-aed7-2689d566feaa"
  squigualiser plot_signal -s ${SIGNAL} -o ${OUTPUT_DIR} -r ${READ_ID} ||  die "testcase:$TESTCASE failed"

  RAW_DIR="${REL_PATH}/data/raw/plot/f5c_eventalign"
  TESTCASE=70.5
  info "testcase:$TESTCASE - read-signal plot with reverse signal"
  SIGNAL="${RAW_DIR}/reads.blow5"
  REGION=""
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  READ_ID="64a25d50-50e0-41f8-aed7-2689d566feaa"
  squigualiser plot_signal -s ${SIGNAL} -o ${OUTPUT_DIR} -r ${READ_ID} --reverse_signal ||  die "testcase:$TESTCASE failed"

  RAW_DIR="${REL_PATH}/data/raw/plot/f5c_eventalign"
  TESTCASE=70.6
  info "testcase:$TESTCASE - read-signal plot with scaling"
  SIGNAL="${RAW_DIR}/reads.blow5"
  REGION=""
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  READ_ID="64a25d50-50e0-41f8-aed7-2689d566feaa"
  squigualiser plot_signal -s ${SIGNAL} -o ${OUTPUT_DIR} -r ${READ_ID} --sig_scale znorm ||  die "testcase:$TESTCASE failed"

  RAW_DIR="${REL_PATH}/data/raw/plot/f5c_eventalign"
  TESTCASE=70.7
  info "testcase:$TESTCASE - read-signal plot with region"
  SIGNAL="${RAW_DIR}/reads.blow5"
  REGION="300-500"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  READ_ID="64a25d50-50e0-41f8-aed7-2689d566feaa"
  squigualiser plot_signal -s ${SIGNAL} -o ${OUTPUT_DIR} -r ${READ_ID} --region ${REGION} ||  die "testcase:$TESTCASE failed"
  
  RAW_DIR="${REL_PATH}/data/raw/plot"
  RAW_DIR="${RAW_DIR}/rna/t0"
  TESTCASE=70.7
  info "testcase:$TESTCASE - read-signal plot RNA"
  SIGNAL="${RAW_DIR}/sequin_reads.blow5"
  REGION=""
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  READ_ID="00213403-4297-4f03-8412-3cc8b9cb845a"
  squigualiser plot_signal -s ${SIGNAL} -o ${OUTPUT_DIR} -r ${READ_ID} ||  die "testcase:$TESTCASE failed"

}

testcase_70s #basic

info "all testcases passed"
#rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0
