#!/bin/bash

# steps
# cat slow5 files
# cat blow5 files
# convert catenated blow5 to slow5
# compare slow5s against the expected
# additional error catching testcases

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
OUTPUT_DIR="${REL_PATH}/data/out/reform"
test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"
#commands ...

ex() {
  python src/reform.py "$@"
}

RAW_DIR="${REL_PATH}/data/raw/reform"
EXP_DIR="${REL_PATH}/data/exp/reform"

#TESTCASE=0.1
#info "testcase:$TESTCASE - help"
#ex && die "testcase:$TESTCASE failed"
#
#TESTCASE=0.2
#info "testcase:$TESTCASE - read:1,kmer:0,move:1,output:paf"
#ex -k0 -m1 -c "${RAW_DIR}/guppy_one_read.bam" > ${OUTPUT_DIR}/out.paf && die "testcase:$TESTCASE failed"


info "no testcases written here, check pipelines"
rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0
