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

testcase_5s() {
  GENOME="${REL_PATH}/data/raw/plot/reference_genomes/nCoV-2019.reference.fasta"

  TESTCASE=5.1
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r2/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r2/sim.sam"
  REGION=""
  PLOT_LIMIT=1
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} || die "testcase:$TESTCASE failed"

  TESTCASE=5.2
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r2/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r2/sim.sam"
  REGION="MN908947.3:1-10000"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" && die "testcase:$TESTCASE failed"

  TESTCASE=5.3
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r2/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r2/sim.bam"
  REGION="MN908947.3:1-10000"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" && die "testcase:$TESTCASE failed"

  TESTCASE=5.4
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r2/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r2/sorted_sim.bam"
  REGION="MN908947.3:14,843-14,914"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup2.html

  TESTCASE=5.5
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r2/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r2/sorted_sim.bam"
  REGION="MN908947.3:14,843-14,914"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" --no_reverse|| die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup2.html

  TESTCASE=5.6
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r2/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r2/sorted_sim.bam"
  REGION="MN908947.3:14,843-14,914"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" --reverse_only|| die "testcase:$TESTCASE failed"

  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup2.html

  TESTCASE=5.7
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r5_one/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r5_one/sorted_sim.bam"
  REGION="MN908947.3:3-4"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" --reverse_only || die "testcase:$TESTCASE failed"

  TESTCASE=5.8
  info "testcase:$TESTCASE"
  REGION="MN908947.3:4-4"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" --reverse_only || die "testcase:$TESTCASE failed"

  TESTCASE=5.9
  info "testcase:$TESTCASE"
  REGION="MN908947.3:4-5"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" --reverse_only || die "testcase:$TESTCASE failed"

  TESTCASE=5.10
  info "testcase:$TESTCASE"
  REGION="MN908947.3:5-5"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" --reverse_only || die "testcase:$TESTCASE failed"

  TESTCASE=5.11
  info "testcase:$TESTCASE"
  REGION="MN908947.3:551-552"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" --reverse_only || die "testcase:$TESTCASE failed"

  TESTCASE=5.12
  info "testcase:$TESTCASE"
  REGION="MN908947.3:0-10"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" --reverse_only || die "testcase:$TESTCASE failed"

  TESTCASE=5.13
  info "testcase:$TESTCASE"
  REGION="MN908947.3:552-553"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" --reverse_only || die "testcase:$TESTCASE failed"
}
testcase_6s() {
  GENOME="${REL_PATH}/data/raw/plot/reference_genomes/nCoV-2019.reference.fasta"

  TESTCASE=6.1
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r3/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r3/sorted_sim.bam"
  REGION="MN908947.3:14,843-14,914"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup2.html

  TESTCASE=6.2
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r3/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r3/sorted_sim.bam"
  REGION="MN908947.3:14,843-14,914"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" --no_reverse|| die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup2.html

  TESTCASE=6.3
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r3/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r3/sorted_sim.bam"
  REGION="MN908947.3:14,843-14,914"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}" --reverse_only|| die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup2.html
}
testcase_7s() {
  GENOME="${REL_PATH}/data/raw/plot/reference_genomes/rnasequin_sequences_2.4.fa"

  TESTCASE=7.1
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r4/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r4/sorted_sim.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" && die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup2.html

  TESTCASE=7.2
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/r4/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/r4/sorted_sim.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --rna --plot_limit ${PLOT_LIMIT} || die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup2.html

}
testcase_8s() {
  GENOME="/media/hiruna/data/basecalling_work/apply_variants_to_genome/genome/hg38noAlt.fa"

  TESTCASE=8.1
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/realigned_DNA/reads.blow5"
  ALIGNMENT="${RAW_DIR}/realigned_DNA/realigned.sam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} || die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup2.html

}
testcase_9s() {
  GENOME="${REL_PATH}/data/raw/plot/reference_genomes/rnasequin_sequences_2.4.fa"

  TESTCASE=9.1
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/rna/t0/sequin_reads.blow5"
  ALIGNMENT="${RAW_DIR}/realigned_RNA/one_read/realigned.sam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  python src/plot.py --rna -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} || die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup2.html

}
#testcase_5s #signal-reference squigulator ideal signals
#testcase_6s #signal-reference squigulator
#testcase_7s #signal-reference squigulator RNA
#testcase_8s #signal-reference realigned DNA
testcase_9s #signal-reference realigned RNA

info "all testcases passed"
#rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0
