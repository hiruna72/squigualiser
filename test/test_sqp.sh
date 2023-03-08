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
OUTPUT_DIR="${REL_PATH}/data/out/sqp"
test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"

GENOME="/media/hiruna/data/basecalling_work/apply_variants_to_genome/genome/hg38noAlt.fa"

RAW_DIR="${REL_PATH}/data/raw/sqp"
EXP_DIR="${REL_PATH}/data/exp/sqp"

TESTCASE=01
info "testcase:$TESTCASE - help"
python src/sqp.py && die "testcase:$TESTCASE failed"

TESTCASE=02
info "testcase:$TESTCASE - help"
python src/sqp.py --help || die "testcase:$TESTCASE failed"

TESTCASE=03
info "testcase:$TESTCASE - read-signal plot"
FASTA="${RAW_DIR}/one_read/read.fasta"
SIGNAL="${RAW_DIR}/one_read/read.slow5"
ALIGNMENT="${RAW_DIR}/one_read/reform/r1k1m1.paf"
REGION=""
OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
python src/sqp.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

TESTCASE=04
info "testcase:$TESTCASE - reference-signal plot"
FASTA=${GENOME}
SIGNAL="${RAW_DIR}/one_read/read.slow5"
ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.sam"
REGION=""
OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
python src/sqp.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

TESTCASE=05
info "testcase:$TESTCASE - reference-signal plot"
FASTA=${GENOME}
SIGNAL="${RAW_DIR}/one_read/read.slow5"
ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.sam"
REGION="chr1:92780309-92780570"
OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
python src/sqp.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

TESTCASE=06
info "testcase:$TESTCASE - reference-signal plot"
FASTA=${GENOME}
SIGNAL="${RAW_DIR}/one_read/read.slow5"
ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.sam"
REGION="chr1:92780311-92780570"
OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
python src/sqp.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

TESTCASE=07
info "testcase:$TESTCASE - reference-signal plot"
FASTA=${GENOME}
SIGNAL="${RAW_DIR}/one_read/read.slow5"
ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.sam"
REGION="chr1:92780309-92780569"
OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
python src/sqp.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

TESTCASE=08
info "testcase:$TESTCASE - reference-signal plot"
FASTA=${GENOME}
SIGNAL="${RAW_DIR}/one_read/read.slow5"
ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.sam"
REGION="chr1:92780311-92780569"
OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
python src/sqp.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

TESTCASE=09
info "testcase:$TESTCASE - reference-signal plot"
FASTA=${GENOME}
SIGNAL="${RAW_DIR}/one_read/read.slow5"
ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.sam"
REGION="chr1:92780301-92780308"
OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
python src/sqp.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

TESTCASE=10
info "testcase:$TESTCASE - reference-signal plot"
FASTA=${GENOME}
SIGNAL="${RAW_DIR}/one_read/read.slow5"
ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.sam"
REGION="chr1:92780571-92780579"
OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
python src/sqp.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

TESTCASE=11
info "testcase:$TESTCASE - reference-signal plot"
FASTA=${GENOME}
SIGNAL="${RAW_DIR}/one_read/read.slow5"
ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.sam"
REGION="chr1:92780565-92780569"
OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
python src/sqp.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

TESTCASE=12
info "testcase:$TESTCASE - reference-signal plot"
FASTA=${GENOME}
SIGNAL="${RAW_DIR}/one_read/read.slow5"
ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.sam"
REGION="chr1:92780308-92780570"
OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
python src/sqp.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

TESTCASE=13
info "testcase:$TESTCASE - reference-signal plot"
FASTA=${GENOME}
SIGNAL="${RAW_DIR}/one_read/read.slow5"
ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.sam"
REGION="chr1:92780308-92780569"
OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
python src/sqp.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

TESTCASE=14
info "testcase:$TESTCASE - reference-signal plot"
FASTA=${GENOME}
SIGNAL="${RAW_DIR}/one_read/read.slow5"
ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.sam"
REGION="chr1:92780311-92780571"
OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
python src/sqp.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

TESTCASE=15
info "testcase:$TESTCASE - reference-signal plot"
FASTA=${GENOME}
SIGNAL="${RAW_DIR}/one_read/read.slow5"
ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.sam"
REGION="chr1:92780311-92780569"
OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
python src/sqp.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

cat ${OUTPUT_DIR}/*/*.html > ${OUTPUT_DIR}/pileup.html

TESTCASE=16
info "testcase:$TESTCASE - reference-signal plot"
FASTA=${GENOME}
SIGNAL="${RAW_DIR}/reverse_read/reads.slow5"
ALIGNMENT="${RAW_DIR}/reverse_read/realign.sam"
REGION="chr1:6811011-6811050"
OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
python src/sqp.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

cat ${OUTPUT}/*.html > ${OUTPUT_DIR}/pileup2.html


info "all $TESTCASE testcases passed"
#rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0
# If you want to log to the same file: command1 >> log_file 2>&1
# If you want different files: command1 >> log_file 2>> err_file
# use ANSI syntax format to view stdout/stderr on SublimeText
# use bash -n [script] and shellcheck [script] to check syntax
