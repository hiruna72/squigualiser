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

testcase_4s() {
  GENOME="/media/hiruna/data/basecalling_work/apply_variants_to_genome/genome/hg38noAlt.fa"

  TESTCASE=4.1
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.bam"
  REGION=""
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=4.2
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.bam"
  REGION="chr1:92780309-92780570"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=4.3
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.bam"
  REGION="chr1:92780311-92780570"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=4.4
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.bam"
  REGION="chr1:92780309-92780569"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=4.5
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.bam"
  REGION="chr1:92780311-92780569"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=4.6
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.bam"
  REGION="chr1:92780301-92780308"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=4.7
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.bam"
  REGION="chr1:92780571-92780579"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=4.8
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.bam"
  REGION="chr1:92780565-92780569"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=4.9
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.bam"
  REGION="chr1:92780308-92780570"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=4.10
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.bam"
  REGION="chr1:92780308-92780569"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=4.11
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.bam"
  REGION="chr1:92780311-92780571"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  TESTCASE=4.12
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="${RAW_DIR}/one_read/realign/realigned_r1k1m1.bam"
  REGION="chr1:92780311-92780569"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  cat ${OUTPUT_DIR}/*/*.html > ${OUTPUT_DIR}/pileup.html

  TESTCASE=4.13
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA="${RAW_DIR}/simulate_reads/sim.fasta"
  SIGNAL="${RAW_DIR}/simulate_reads/sim.slow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/sim.paf"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  cat ${OUTPUT}/*.html > ${OUTPUT_DIR}/pileup2.html

  TESTCASE=4.14
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/reverse_read/reads.slow5"
  ALIGNMENT="${RAW_DIR}/reverse_read/realign.bam"
  REGION="chr1:6811011-6811050"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --region ${REGION} --tag_name "testcase-${TESTCASE}"|| die "testcase:$TESTCASE failed"

  cat ${OUTPUT}/*.html >> ${OUTPUT_DIR}/pileup2.html
}
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

  TESTCASE=6.4
  info "testcase:$TESTCASE - squigulator sam output"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/reference_paf/one/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/reference_paf/one/sim.sam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" || die "testcase:$TESTCASE failed"

  TESTCASE=6.5
  info "testcase:$TESTCASE - squigulator paf output"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/reference_paf/one/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/reference_paf/one/sim.paf"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" && die "testcase:$TESTCASE failed"

  TESTCASE=6.6
  info "testcase:$TESTCASE - squigulator paf output"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/reference_paf/one/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/reference_paf/one/sim.paf"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --sig_ref && die "testcase:$TESTCASE failed"

  TESTCASE=6.7
  info "testcase:$TESTCASE - squigulator paf output"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/reference_paf/one/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/reference_paf/one/sorted_sim.paf.gz"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --sig_ref || die "testcase:$TESTCASE failed"

  TESTCASE=6.8
  info "testcase:$TESTCASE - squigulator sam output"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/reference_paf/ten/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/reference_paf/ten/sim.sam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=4
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} || die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/${TESTCASE}.html

  TESTCASE=6.9
  info "testcase:$TESTCASE - squigulator paf output"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/simulate_reads/reference_paf/ten/sim.blow5"
  ALIGNMENT="${RAW_DIR}/simulate_reads/reference_paf/ten/sorted_sim.paf.gz"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=4
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --sig_ref --plot_limit ${PLOT_LIMIT} || die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/${TESTCASE}.html

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
  ALIGNMENT="${RAW_DIR}/realigned_DNA/realigned.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  python src/plot.py -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} || die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup_${TESTCASE}.html

  TESTCASE=8.2
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/realigned_DNA/reads.blow5"
  ALIGNMENT="${RAW_DIR}/realigned_DNA/realigned.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  python src/plot.py --fixed_width -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} || die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup_${TESTCASE}.html

  TESTCASE=8.3
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/realigned_DNA/reads.blow5"
  ALIGNMENT="${RAW_DIR}/realigned_DNA/realigned.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  python src/plot.py --pileup -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} && die "testcase:$TESTCASE failed"

  TESTCASE=8.4
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/realigned_DNA/reads.blow5"
  ALIGNMENT="${RAW_DIR}/realigned_DNA/realigned.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  python src/plot.py --pileup --fixed_width -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} && die "testcase:$TESTCASE failed"

  TESTCASE=8.5
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/realigned_DNA/reads.blow5"
  ALIGNMENT="${RAW_DIR}/realigned_DNA/realigned.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  REGION="chr1:6,811,009-6,811,198"
  python src/plot.py --pileup --fixed_width --region ${REGION} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} || die "testcase:$TESTCASE failed"

  TESTCASE=8.6
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/realigned_DNA/reads.blow5"
  ALIGNMENT="${RAW_DIR}/realigned_DNA/realigned.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  REGION="chr1:6,811,009-6,811,198"
  python src/plot.py --pileup --fixed_width --region ${REGION} --no_reverse -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} || die "testcase:$TESTCASE failed"

}
testcase_9s() {
  GENOME="${REL_PATH}/data/raw/plot/reference_genomes/rnasequin_sequences_2.4.fa"

  TESTCASE=9.1
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/rna/t0/sequin_reads.blow5"
  ALIGNMENT="${RAW_DIR}/realigned_RNA/one_read/test_5.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  python src/plot.py --rna -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" || die "testcase:$TESTCASE failed"

  TESTCASE=9.2
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/rna/t0/sequin_reads.blow5"
  ALIGNMENT="${RAW_DIR}/realigned_RNA/all/test_6.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  python src/plot.py --rna -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} || die "testcase:$TESTCASE failed"
  cat ${OUTPUT}/*.html >> ${OUTPUT}/pileup2.html

  TESTCASE=9.3
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/rna/t0/sequin_reads.blow5"
  ALIGNMENT="${RAW_DIR}/realigned_RNA/all/test_6.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  REGION="R2_6_1:472-511"
  READ_ID="dac12fd6-4c9d-4da7-9014-91a2e21e109d"
  python src/plot.py --rna --read_id ${READ_ID} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} --region ${REGION}|| die "testcase:$TESTCASE failed"

  TESTCASE=9.4
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/rna/t0/sequin_reads.blow5"
  ALIGNMENT="${RAW_DIR}/realigned_RNA/all/test_6.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  REGION="R2_6_1:472-511"
  READ_ID="dac12fd6-4c9d-4da7-9014-91a2e21e109d"
  python src/plot.py --fixed_width --rna --read_id ${READ_ID} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} --region ${REGION}|| die "testcase:$TESTCASE failed"

}
testcase_4s #basic
testcase_5s #signal-reference squigulator ideal signals
testcase_6s #signal-reference squigulator
testcase_7s #signal-reference squigulator RNA
testcase_8s #signal-reference realigned DNA
testcase_9s #signal-reference realigned RNA

info "all testcases passed"
#rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0