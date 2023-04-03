#!/bin/bash

# steps

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
OUTPUT_DIR="${REL_PATH}/data/out/realign"
test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"
#commands ...

ex() {
  python src/realign.py "$@"
}

RAW_DIR="${REL_PATH}/data/raw/realign"
EXP_DIR="${REL_PATH}/data/exp/realign"

TESTCASE=0
info "testcase:$TESTCASE - help"
ex && die "testcase:$TESTCASE failed"

TESTCASE=1
info "testcase:$TESTCASE - read:1,forward"
ex --paf ${RAW_DIR}/test_${TESTCASE}/r1k1m1.paf --bam ${RAW_DIR}/test_${TESTCASE}/map.sam --output ${OUTPUT_DIR}/test_${TESTCASE}.sam || die "testcase:${TESTCASE} failed"
#diff "${EXP_DIR}/test_${TESTCASE}.sam" "${OUTPUT_DIR}/test_${TESTCASE}.sam" || die "testcase:${TESTCASE} diff failed"

TESTCASE=2
info "testcase:$TESTCASE - read:1,forward"
ex --paf ${RAW_DIR}/test_${TESTCASE}/p1.paf --bam ${RAW_DIR}/test_${TESTCASE}/m1.sam --output ${OUTPUT_DIR}/test_${TESTCASE}.sam || die "testcase:${TESTCASE} failed"
diff "${EXP_DIR}/test_${TESTCASE}.sam" "${OUTPUT_DIR}/test_${TESTCASE}.sam" || die "testcase:${TESTCASE} diff failed"

TESTCASE=3
info "testcase:$TESTCASE"
ex --paf ${RAW_DIR}/test_${TESTCASE}/move.paf --bam ${RAW_DIR}/test_${TESTCASE}/map.sam --output ${OUTPUT_DIR}/test_${TESTCASE}.sam || die "testcase:${TESTCASE} failed"
diff "${EXP_DIR}/test_${TESTCASE}.sam" "${OUTPUT_DIR}/test_${TESTCASE}.sam" || die "testcase:${TESTCASE} diff failed"

TESTCASE=4
info "testcase:$TESTCASE"
ex --paf ${RAW_DIR}/test_${TESTCASE}/move.paf --bam ${RAW_DIR}/test_${TESTCASE}/map.sam --output ${OUTPUT_DIR}/test_${TESTCASE}.sam
diff "${EXP_DIR}/test_${TESTCASE}.sam" "${OUTPUT_DIR}/test_${TESTCASE}.sam" || die "testcase:${TESTCASE} diff failed"

TESTCASE=5
info "testcase:$TESTCASE - RNA"
ex --paf ${RAW_DIR}/test_${TESTCASE}/move.paf --bam ${RAW_DIR}/test_${TESTCASE}/map.sam --output ${OUTPUT_DIR}/test_${TESTCASE}.sam
diff "${EXP_DIR}/test_${TESTCASE}.sam" "${OUTPUT_DIR}/test_${TESTCASE}.sam" || die "testcase:${TESTCASE} diff failed"

info "all $TESTCASE testcases passed"
rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0
