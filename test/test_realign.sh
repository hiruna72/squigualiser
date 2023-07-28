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

TESTCASE=0.0
info "testcase:$TESTCASE - help"
ex && die "testcase:$TESTCASE failed"

testcase_1s() {
  TESTCASE=1.1
  info "testcase:$TESTCASE - read:1,forward"
  ex --paf ${RAW_DIR}/test_${TESTCASE}/r1k1m0.paf --bam ${RAW_DIR}/test_${TESTCASE}/map.sam --output ${OUTPUT_DIR}/test_${TESTCASE}.sam || die "testcase:${TESTCASE} failed"
  diff "${EXP_DIR}/test_${TESTCASE}.sam" "${OUTPUT_DIR}/test_${TESTCASE}.sam" || die "testcase:${TESTCASE} diff failed"

  TESTCASE=1.2
  info "testcase:$TESTCASE - read:1,forward"
  ex --paf ${RAW_DIR}/test_${TESTCASE}/p1.paf --bam ${RAW_DIR}/test_${TESTCASE}/m1.sam --output ${OUTPUT_DIR}/test_${TESTCASE}.sam || die "testcase:${TESTCASE} failed"
  diff "${EXP_DIR}/test_${TESTCASE}.sam" "${OUTPUT_DIR}/test_${TESTCASE}.sam" || die "testcase:${TESTCASE} diff failed"

  TESTCASE=1.3
  info "testcase:$TESTCASE"
  ex --paf ${RAW_DIR}/test_${TESTCASE}/move.paf --bam ${RAW_DIR}/test_${TESTCASE}/map.sam --output ${OUTPUT_DIR}/test_${TESTCASE}.sam || die "testcase:${TESTCASE} failed"
  diff "${EXP_DIR}/test_${TESTCASE}.sam" "${OUTPUT_DIR}/test_${TESTCASE}.sam" || die "testcase:${TESTCASE} diff failed"

  TESTCASE=1.4
  info "testcase:$TESTCASE"
  ex --paf ${RAW_DIR}/test_${TESTCASE}/move.paf --bam ${RAW_DIR}/test_${TESTCASE}/map.sam --output ${OUTPUT_DIR}/test_${TESTCASE}.sam || die "testcase:${TESTCASE} failed"
  diff "${EXP_DIR}/test_${TESTCASE}.sam" "${OUTPUT_DIR}/test_${TESTCASE}.sam" || die "testcase:${TESTCASE} diff failed"

  TESTCASE=1.5
  info "testcase:$TESTCASE - RNA"
  ex --rna --paf ${RAW_DIR}/test_${TESTCASE}/move.paf --bam ${RAW_DIR}/test_${TESTCASE}/map.bam --output ${OUTPUT_DIR}/test_${TESTCASE}.sam || die "testcase:${TESTCASE} failed"
  diff "${EXP_DIR}/test_${TESTCASE}.sam" "${OUTPUT_DIR}/test_${TESTCASE}.sam" || die "testcase:${TESTCASE} diff failed"

  TESTCASE=1.6
  info "testcase:$TESTCASE - RNA"
  ex --rna --paf "${REL_PATH}/data/raw/plot/rna/t0/reform_rna_2.3.paf" --bam ${RAW_DIR}/test_${TESTCASE}/sorted_map.bam --output ${OUTPUT_DIR}/test_${TESTCASE}.sam || die "testcase:${TESTCASE} failed"
  diff "${EXP_DIR}/test_${TESTCASE}.sam" "${OUTPUT_DIR}/test_${TESTCASE}.sam" || die "testcase:${TESTCASE} diff failed"

  TESTCASE=1.7
  info "testcase:$TESTCASE - DNA kmer size = 1"
  ex --paf "${RAW_DIR}/test_${TESTCASE}/reform.paf" --bam ${RAW_DIR}/test_${TESTCASE}/mapped.bam --output ${OUTPUT_DIR}/test_${TESTCASE}.sam || die "testcase:${TESTCASE} failed"
  diff "${EXP_DIR}/test_${TESTCASE}.sam" "${OUTPUT_DIR}/test_${TESTCASE}.sam" || die "testcase:${TESTCASE} diff failed"

  TESTCASE=1.8
  info "testcase:$TESTCASE - DNA kmer size > 1"
  ex --paf "${RAW_DIR}/test_${TESTCASE}/re_reform.paf" --bam ${RAW_DIR}/test_${TESTCASE}/mapped.bam --output ${OUTPUT_DIR}/test_${TESTCASE}.sam || die "testcase:${TESTCASE} failed"
  diff "${EXP_DIR}/test_${TESTCASE}.sam" "${OUTPUT_DIR}/test_${TESTCASE}.sam" || die "testcase:${TESTCASE} diff failed"

}
testcase_2s() {
  TESTCASE=2.1
  info "testcase:$TESTCASE - read:1,forward"
  ex --paf ${RAW_DIR}/test_${TESTCASE}/r1k1m0.paf --bam ${RAW_DIR}/test_${TESTCASE}/map.sam --output ${OUTPUT_DIR}/test_${TESTCASE}.paf -c || die "testcase:${TESTCASE} failed"
  diff "${EXP_DIR}/test_${TESTCASE}.paf" "${OUTPUT_DIR}/test_${TESTCASE}.paf" || die "testcase:${TESTCASE} diff failed"

  TESTCASE=2.2
  info "testcase:$TESTCASE - read:1,forward"
  ex --paf ${RAW_DIR}/test_${TESTCASE}/p1.paf --bam ${RAW_DIR}/test_${TESTCASE}/m1.sam --output ${OUTPUT_DIR}/test_${TESTCASE}.paf -c || die "testcase:${TESTCASE} failed"
  diff "${EXP_DIR}/test_${TESTCASE}.paf" "${OUTPUT_DIR}/test_${TESTCASE}.paf" || die "testcase:${TESTCASE} diff failed"

  TESTCASE=2.3
  info "testcase:$TESTCASE"
  ex --paf ${RAW_DIR}/test_${TESTCASE}/move.paf --bam ${RAW_DIR}/test_${TESTCASE}/map.sam --output ${OUTPUT_DIR}/test_${TESTCASE}.paf -c || die "testcase:${TESTCASE} failed"
  diff "${EXP_DIR}/test_${TESTCASE}.paf" "${OUTPUT_DIR}/test_${TESTCASE}.paf" || die "testcase:${TESTCASE} diff failed"

  TESTCASE=2.4
  info "testcase:$TESTCASE"
  ex --paf ${RAW_DIR}/test_${TESTCASE}/move.paf --bam ${RAW_DIR}/test_${TESTCASE}/map.sam --output ${OUTPUT_DIR}/test_${TESTCASE}.paf -c || die "testcase:${TESTCASE} failed"
  diff "${EXP_DIR}/test_${TESTCASE}.paf" "${OUTPUT_DIR}/test_${TESTCASE}.paf" || die "testcase:${TESTCASE} diff failed"

  TESTCASE=2.5
  info "testcase:$TESTCASE - RNA"
  ex --rna --paf ${RAW_DIR}/test_${TESTCASE}/move.paf --bam ${RAW_DIR}/test_${TESTCASE}/map.bam --output ${OUTPUT_DIR}/test_${TESTCASE}.paf -c || die "testcase:${TESTCASE} failed"
  diff "${EXP_DIR}/test_${TESTCASE}.paf" "${OUTPUT_DIR}/test_${TESTCASE}.paf" || die "testcase:${TESTCASE} diff failed"

  TESTCASE=2.6
  info "testcase:$TESTCASE - RNA"
  ex --rna --paf "${REL_PATH}/data/raw/plot/rna/t0/reform_rna_2.3.paf" --bam ${RAW_DIR}/test_${TESTCASE}/sorted_map.bam --output ${OUTPUT_DIR}/test_${TESTCASE}.paf -c || die "testcase:${TESTCASE} failed"
  diff "${EXP_DIR}/test_${TESTCASE}.paf" "${OUTPUT_DIR}/test_${TESTCASE}.paf" || die "testcase:${TESTCASE} diff failed"

  TESTCASE=2.7
  info "testcase:$TESTCASE - DNA kmer size = 1"
  ex --paf "${RAW_DIR}/test_1.7/reform.paf" --bam "${RAW_DIR}/test_1.7/mapped.bam" --output ${OUTPUT_DIR}/test_${TESTCASE}.paf -c || die "testcase:${TESTCASE} failed"
  diff "${EXP_DIR}/test_${TESTCASE}.paf" "${OUTPUT_DIR}/test_${TESTCASE}.paf" || die "testcase:${TESTCASE} diff failed"

  TESTCASE=2.8
  info "testcase:$TESTCASE - DNA kmer size > 1"
  ex --paf "${RAW_DIR}/test_1.7/re_reform.paf" --bam "${RAW_DIR}/test_1.7/mapped.bam" --output ${OUTPUT_DIR}/test_${TESTCASE}.paf -c || die "testcase:${TESTCASE} failed"
  diff "${EXP_DIR}/test_${TESTCASE}.paf" "${OUTPUT_DIR}/test_${TESTCASE}.paf" || die "testcase:${TESTCASE} diff failed"

}

testcase_1s #bam output
testcase_2s #paf output

info "all $TESTCASE testcases passed"
rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0
