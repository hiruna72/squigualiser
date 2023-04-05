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

RAW_DIR="${REL_PATH}/data/raw/plot"
EXP_DIR="${REL_PATH}/data/exp/plot"

testcase_0s() {
  TESTCASE=0.1
  info "testcase:$TESTCASE - help"
  python src/plot.py && die "testcase:$TESTCASE failed"

  TESTCASE=0.2
  info "testcase:$TESTCASE - help"
  python src/plot.py --help && die "testcase:$TESTCASE failed"
}
testcase_1s() {
  TESTCASE=1.1
  info "testcase:$TESTCASE - read-signal plot"
  FASTA="${RAW_DIR}/one_read/read.fasta"
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/reform/r1k1m1.paf"
  REGION=""
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=1.2
  info "testcase:$TESTCASE - read-signal plot region specified"
  FASTA="${RAW_DIR}/one_read/read.fasta"
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/reform/r1k1m1.paf"
  REGION="1-20"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=1.3
  info "testcase:$TESTCASE - read-signal plot region specified"
  FASTA="${RAW_DIR}/one_read/read.fasta"
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/reform/r1k1m1.paf"
  REGION="10-40"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=1.4
  info "testcase:$TESTCASE - read-signal plot region specified"
  FASTA="${RAW_DIR}/one_read/read.fasta"
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/reform/r1k9m1.paf"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=1.5
  info "testcase:$TESTCASE - read-signal plot region specified"
  FASTA="${RAW_DIR}/one_read/read.fasta"
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/reform/r1k9m1.paf"
  REGION="10-268"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=1.6
  info "testcase:$TESTCASE - read-signal plot region specified"
  FASTA="${RAW_DIR}/one_read/read.fasta"
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/reform/r1k9m2.paf"
  REGION="10-268"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=1.7
  info "testcase:$TESTCASE - read-signal plot region specified"
  FASTA="${RAW_DIR}/one_read/read.fasta"
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/reform/r1k9m7.paf"
  REGION="10-268"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=1.9
  info "testcase:$TESTCASE - read-signal plot region specified"
  FASTA="${RAW_DIR}/one_read/read.fasta"
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/reform/r1k9m9.paf"
  REGION="10-268"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=1.10
  info "testcase:$TESTCASE - read-signal plot region specified"
  FASTA="${RAW_DIR}/one_read/read.fasta"
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/reform/r1k9m9.paf"
  REGION="10-268"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py --no_pa -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=1.11
  info "testcase:$TESTCASE - read-signal plot region specified"
  FASTA="${RAW_DIR}/simulate_reads/r1/sim.fasta"
  SIGNAL="${RAW_DIR}/simulate_reads/r1/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r1/sim.paf"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py --no_pa -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"
}
testcase_2s() {
  TESTCASE=2.1
  info "testcase:$TESTCASE - resquiggle read-signal plot"
  FASTA="${RAW_DIR}/resquiggle_dna/t1/read.fasta"
  SIGNAL="${RAW_DIR}/resquiggle_dna/t1/read.slow5"
  ALIGNMENT="${RAW_DIR}/resquiggle_dna/t1/resquiggle_move.paf"
  REGION=""
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=2.2
  info "testcase:$TESTCASE - resquiggle read-signal plot"
  FASTA="${RAW_DIR}/resquiggle_dna/t0/read.fastq"
  SIGNAL="${RAW_DIR}/resquiggle_dna/t0/read.blow5"
  ALIGNMENT="${RAW_DIR}/resquiggle_dna/t0/resquiggle_move.paf"
  REGION=""
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=2.3
  info "testcase:$TESTCASE - resquiggle read-signal plot"
  FASTA="${RAW_DIR}/resquiggle_dna/t0/read.fastq"
  SIGNAL="${RAW_DIR}/resquiggle_dna/t0/read.blow5"
  ALIGNMENT="${RAW_DIR}/resquiggle_dna/t0/resquiggle_move.paf"
  REGION="900-1495"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=2.4
  info "testcase:$TESTCASE - resquiggle read-signal plot"
  FASTA="${RAW_DIR}/resquiggle_dna/t0/read.fastq"
  SIGNAL="${RAW_DIR}/resquiggle_dna/t0/read.blow5"
  ALIGNMENT="${RAW_DIR}/resquiggle_dna/t0/resquiggle_move.paf"
  REGION="1483-1495"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"
}
testcase_3s() {
  RAW_DIR="${RAW_DIR}/rna/t0"
  TESTCASE=3.1
  info "testcase:$TESTCASE - resquiggle RNA read-signal plot"
  FASTA="${RAW_DIR}/sequin_reads.fastq"
  SIGNAL="${RAW_DIR}/sequin_reads.blow5"
  ALIGNMENT="${RAW_DIR}/resquiggle_move.paf"
  REGION=""
  READ_ID="00213403-4297-4f03-8412-3cc8b9cb845a"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py --rna -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --read_id ${READ_ID}|| die "testcase:$TESTCASE failed"

  TESTCASE=3.2
  info "testcase:$TESTCASE - resquiggle RNA read-signal plot"
  FASTA="${RAW_DIR}/sequin_reads.fastq"
  SIGNAL="${RAW_DIR}/sequin_reads.blow5"
  ALIGNMENT="${RAW_DIR}/resquiggle_move.paf"
  REGION="1300-1563"
  READ_ID="00213403-4297-4f03-8412-3cc8b9cb845a"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py --rna -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" --read_id ${READ_ID}|| die "testcase:$TESTCASE failed"

  TESTCASE=3.3
  info "testcase:$TESTCASE - resquiggle RNA read-signal plot"
  FASTA="${RAW_DIR}/sequin_reads.fastq"
  SIGNAL="${RAW_DIR}/sequin_reads.blow5"
  ALIGNMENT="${RAW_DIR}/reform_rna_2.3.paf"
  REGION="1300-1563"
  READ_ID="00213403-4297-4f03-8412-3cc8b9cb845a"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py --rna -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" --read_id ${READ_ID}|| die "testcase:$TESTCASE failed"
}

testcase_0s #basic
testcase_1s #signal-read
testcase_2s #resquiggle
testcase_3s #RNA

info "all testcases passed"
rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0
