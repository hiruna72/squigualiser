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

TESTCASE=0.1
info "testcase:$TESTCASE - help"
ex && die "testcase:$TESTCASE failed"

TESTCASE=0.2
info "testcase:$TESTCASE - read:1,kmer:0,move:1,output:paf"
ex -k0 -m1 -c "${RAW_DIR}/guppy_one_read.bam" > ${OUTPUT_DIR}/out.paf && die "testcase:$TESTCASE failed"

TESTCASE=0.3
info "testcase:$TESTCASE - read:1,kmer:9,move:0,output:paf"
ex -k9 -m0 -c "${RAW_DIR}/guppy_one_read.bam" > ${OUTPUT_DIR}/out.paf && die "testcase:$TESTCASE failed"

TESTCASE=0.4
info "testcase:$TESTCASE - read:1,kmer:9,move:10,output:paf"
ex -k9 -m10 -c "${RAW_DIR}/guppy_one_read.bam" > ${OUTPUT_DIR}/out.paf && die "testcase:$TESTCASE failed"

TESTCASE=0.5
info "testcase:$TESTCASE - read:1,kmer:1,move:0,output:paf"
ex -k1 -m0 -c --bam "${RAW_DIR}/guppy_one_read.bam" --output "${OUTPUT_DIR}/r1k1m0.paf" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k1m0.paf" "${OUTPUT_DIR}/r1k1m0.paf" || die "testcase:${TESTCASE} diff failed"

TESTCASE=1.1
info "testcase:$TESTCASE - read:1,kmer:9,move:0,output:paf"
ex -k9 -m0 -c --bam "${RAW_DIR}/guppy_one_read.bam" --output "${OUTPUT_DIR}/r1k9m0.paf" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k9m0.paf" "${OUTPUT_DIR}/r1k9m0.paf" || die "testcase:${TESTCASE} diff failed"

TESTCASE=1.2
info "testcase:$TESTCASE - read:1,kmer:9,move:1,output:paf"
ex -k9 -m1 -c --bam "${RAW_DIR}/guppy_one_read.bam" --output "${OUTPUT_DIR}/r1k9m1.paf" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k9m1.paf" "${OUTPUT_DIR}/r1k9m1.paf" || die "testcase:${TESTCASE} diff failed"

TESTCASE=1.3
info "testcase:$TESTCASE - read:1,kmer:9,move:6,output:paf"
ex -k9 -m6 -c --bam "${RAW_DIR}/guppy_one_read.bam" --output "${OUTPUT_DIR}/r1k9m6.paf" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k9m6.paf" "${OUTPUT_DIR}/r1k9m6.paf" || die "testcase:${TESTCASE} diff failed"

TESTCASE=1.4
info "testcase:$TESTCASE - read:1,kmer:9,move:8,output:paf"
ex -k9 -m8 -c --bam "${RAW_DIR}/guppy_one_read.bam" --output "${OUTPUT_DIR}/r1k9m8.paf" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k9m8.paf" "${OUTPUT_DIR}/r1k9m8.paf" || die "testcase:${TESTCASE} diff failed"

TESTCASE=1.5
info "testcase:$TESTCASE - read:1,kmer:1,move:0,output:tsv"
ex -k1 -m0 --bam "${RAW_DIR}/guppy_one_read.bam" --output "${OUTPUT_DIR}/r1k1m0.tsv" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k1m0.tsv" "${OUTPUT_DIR}/r1k1m0.tsv" || die "testcase:${TESTCASE} diff failed"

TESTCASE=1.6
info "testcase:$TESTCASE - read:1,kmer:9,move:0,output:tsv"
ex -k9 -m0 --bam "${RAW_DIR}/guppy_one_read.bam" --output "${OUTPUT_DIR}/r1k9m0.tsv" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k9m0.tsv" "${OUTPUT_DIR}/r1k9m0.tsv" || die "testcase:${TESTCASE} diff failed"

TESTCASE=1.7
info "testcase:$TESTCASE - read:1,kmer:9,move:1,output:tsv"
ex -k9 -m1 --bam "${RAW_DIR}/guppy_one_read.bam" --output "${OUTPUT_DIR}/r1k9m1.tsv" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k9m1.tsv" "${OUTPUT_DIR}/r1k9m1.tsv" || die "testcase:${TESTCASE} diff failed"

TESTCASE=1.8
info "testcase:$TESTCASE - read:1,kmer:9,move:6,output:tsv"
ex -k9 -m6 --bam "${RAW_DIR}/guppy_one_read.bam" --output "${OUTPUT_DIR}/r1k9m6.tsv" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k9m6.tsv" "${OUTPUT_DIR}/r1k9m6.tsv" || die "testcase:${TESTCASE} diff failed"

TESTCASE=1.9
info "testcase:$TESTCASE - read:1,kmer:9,move:8,output:tsv"
ex -k9 -m8 --bam "${RAW_DIR}/guppy_one_read.bam" --output "${OUTPUT_DIR}/r1k9m8.tsv" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k9m8.tsv" "${OUTPUT_DIR}/r1k9m8.tsv" || die "testcase:${TESTCASE} diff failed"

# dorado output
TESTCASE=1.10
info "testcase:$TESTCASE - read:2,kmer:9,move:8,output:tsv"
ex -k9 -m8 --bam "${RAW_DIR}/slow5-dorado.sam" --output "${OUTPUT_DIR}/dr2k9m8.tsv" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/dr2k9m8.tsv" "${OUTPUT_DIR}/dr2k9m8.tsv" || die "testcase:${TESTCASE} diff failed"

TESTCASE=1.11
info "testcase:$TESTCASE - read:2,kmer:9,move:8,output:paf"
ex -k9 -m8 -c --bam "${RAW_DIR}/slow5-dorado.sam" --output "${OUTPUT_DIR}/dr2k9m8.paf" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/dr2k9m8.paf" "${OUTPUT_DIR}/dr2k9m8.paf" || die "testcase:${TESTCASE} diff failed"

TESTCASE=1.12
info "testcase:$TESTCASE - read:2,kmer:1,move:0,output:tsv"
ex -k1 -m0 --bam "${RAW_DIR}/slow5-dorado.sam" --output "${OUTPUT_DIR}/dr2k1m0.tsv" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/dr2k1m0.tsv" "${OUTPUT_DIR}/dr2k1m0.tsv" || die "testcase:${TESTCASE} diff failed"

TESTCASE=1.13
info "testcase:$TESTCASE - read:2,kmer:1,move:0,output:paf"
ex -k1 -m0 -c --bam "${RAW_DIR}/slow5-dorado.sam" --output "${OUTPUT_DIR}/dr2k1m0.paf" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/dr2k1m0.paf" "${OUTPUT_DIR}/dr2k1m0.paf" || die "testcase:${TESTCASE} diff failed"

RAW_DIR=${RAW_DIR}/rna
EXP_DIR=${EXP_DIR}/rna
TESTCASE=2.1
info "testcase:$TESTCASE - read:2,kmer:1,move:0,output:paf rna"
ex --rna -k1 -m0 -c --bam "${RAW_DIR}/guppy_rna.sam" --output "${OUTPUT_DIR}/rna_2.1.paf" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/rna_2.1.paf" "${OUTPUT_DIR}/rna_2.1.paf" || die "testcase:${TESTCASE} diff failed"

TESTCASE=2.2
info "testcase:$TESTCASE - read:2,kmer:1,move:0,output:tsv rna"
ex --rna -k1 -m0 --bam "${RAW_DIR}/guppy_rna.sam" --output "${OUTPUT_DIR}/rna_2.1.tsv" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/rna_2.1.tsv" "${OUTPUT_DIR}/rna_2.1.tsv" || die "testcase:${TESTCASE} diff failed"

TESTCASE=2.3
info "testcase:$TESTCASE - read:all,kmer:1,move:0,output:paf rna"
ex --rna -k1 -m0 -c --bam "${RAW_DIR}/guppy_rna_all.sam" --output "${OUTPUT_DIR}/rna_2.3.paf" || die "testcase:$TESTCASE failed"

TESTCASE=2.4
info "testcase:$TESTCASE - read:all,kmer:1,move:0,output:tsv rna"
ex --rna -k1 -m0 --bam "${RAW_DIR}/guppy_rna_all.sam" --output "${OUTPUT_DIR}/rna_2.4.tsv" || die "testcase:$TESTCASE failed"

info "all $TESTCASE testcases passed"
rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0
