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
OUTPUT_DIR="${REL_PATH}/data/out/plot_tracks"
test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"

RAW_DIR="${REL_PATH}/data/raw/plot_tracks"
EXP_DIR="${REL_PATH}/data/exp/plot_tracks"

testcase_20s() {
  TESTCASE=20.0
  info "testcase:$TESTCASE - plot tracks"
  python src/plot_tracks.py && die "testcase:$TESTCASE failed"

  TESTCASE=20.1
  info "testcase:$TESTCASE - plot tracks"
  COMMAND_FILE="${RAW_DIR}/t_20.1.txt"
  python src/plot_tracks.py -f ${COMMAND_FILE} -o ${OUTPUT_DIR}/${TESTCASE} && die "testcase:$TESTCASE failed"

  TESTCASE=20.2
  info "testcase:$TESTCASE - plot tracks"
  COMMAND_FILE="${RAW_DIR}/t_20.2.txt"
  python src/plot_tracks.py -f ${COMMAND_FILE} -o ${OUTPUT_DIR}/${TESTCASE} && die "testcase:$TESTCASE failed"
}
testcase_20s #pileup

info "all testcases passed"
#rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0