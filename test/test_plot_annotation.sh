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
OUTPUT_DIR="${REL_PATH}/data/out/plot_annotation"
test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"

RAW_DIR="${REL_PATH}/data/raw/plot"
EXP_DIR="${REL_PATH}/data/exp/plot"

testcase_40s() {
  TESTCASE=40.0
  info "testcase:$TESTCASE - read-signal plot"
  FASTA="-f ${RAW_DIR}/one_read/read.fasta"
  SIGNAL="-s ${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="-a ${RAW_DIR}/one_read/reform/r1k1m0.paf"
  OUTPUT="-o ${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT="--plot_limit 10"
  BED="--bed ${RAW_DIR}/bed_files/DNA_1.bed"
  python src/plot.py  ${BED} ${FASTA} ${SIGNAL} ${ALIGNMENT} ${OUTPUT} --tag_name "testcase-${TESTCASE}" ${PLOT_LIMIT}|| die "testcase:$TESTCASE failed"

  TESTCASE=40.1
  info "testcase:$TESTCASE - read-signal plot"
  FASTA="-f ${RAW_DIR}/one_read/read.fasta"
  SIGNAL="-s ${RAW_DIR}/one_read/read.slow5"
  ALIGNMENT="-a ${RAW_DIR}/one_read/reform/r1k1m0.paf"
  OUTPUT="-o ${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT="--plot_limit 10"
  BED="--bed ${RAW_DIR}/bed_files/DNA_1.bed"
  python src/plot.py  ${BED} --fixed_width ${FASTA} ${SIGNAL} ${ALIGNMENT} ${OUTPUT} --tag_name "testcase-${TESTCASE}" ${PLOT_LIMIT}|| die "testcase:$TESTCASE failed"

  TESTCASE=40.2
  [ "${HUMAN_GENOME}" ] || die "Please set the env variable to the human genome path. export HUMAN_GENOME=path/to/file"
  GENOME="${HUMAN_GENOME}"
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA="-f ${GENOME}"
  SIGNAL="-s ${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="-a ${RAW_DIR}/f5c_eventalign/sorted_eventalign.paf.gz"
  OUTPUT="-o ${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT="--plot_limit 10"
  BED="--bed ${RAW_DIR}/bed_files/DNA_2.bed"
  REGION="--region chr1:6,811,011-6,811,198"
  python src/plot.py  ${BED} ${REGION} --fixed_width ${FASTA} ${SIGNAL} ${ALIGNMENT} ${OUTPUT} --tag_name "testcase-${TESTCASE}" ${PLOT_LIMIT}|| die "testcase:$TESTCASE failed"

  TESTCASE=40.3
  [ "${HUMAN_GENOME}" ] || die "Please set the env variable to the human genome path. export HUMAN_GENOME=path/to/file"
  GENOME="${HUMAN_GENOME}"
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA="-f ${GENOME}"
  SIGNAL="-s ${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="-a ${RAW_DIR}/f5c_eventalign/sorted_eventalign.paf.gz"
  OUTPUT="-o ${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT="--plot_limit 10"
  BED="--bed ${RAW_DIR}/bed_files/DNA_3.bed"
  REGION="--region chr1:6,811,011-6,811,198"
  python src/plot.py  ${BED} ${REGION} --fixed_width ${FASTA} ${SIGNAL} ${ALIGNMENT} ${OUTPUT} --tag_name "testcase-${TESTCASE}" ${PLOT_LIMIT}|| die "testcase:$TESTCASE failed"

  TESTCASE=40.4
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/realigned_DNA/reads.blow5"
  ALIGNMENT="${RAW_DIR}/realigned_DNA/sorted_test_1.4.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  READ_ID="--read_id a2d0e216-8610-40b8-92f0-0a04c4a58e08"
  BED="--bed ${RAW_DIR}/bed_files/DNA_4.bed"
  python src/plot.py ${BED} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" ${READ_ID} || die "testcase:$TESTCASE failed"

}
testcase_41s() {
  RAW_DIR_RNA="${RAW_DIR}/rna/t0"
  TESTCASE=41.0
  info "testcase:$TESTCASE - annotation RNA"
  FASTA="${RAW_DIR_RNA}/sequin_reads.fastq"
  SIGNAL="${RAW_DIR_RNA}/sequin_reads.blow5"
  ALIGNMENT="${RAW_DIR_RNA}/resquiggle_move.paf"
  REGION=""
  READ_ID="00213403-4297-4f03-8412-3cc8b9cb845a"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  BED="--bed ${REL_PATH}/data/raw/plot/bed_files/RNA_1.bed"
  python src/plot.py --rna ${BED} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --read_id ${READ_ID}|| die "testcase:$TESTCASE failed"

  GENOME="${REL_PATH}/data/raw/plot/reference_genomes/rnasequin_sequences_2.4.fa"
  TESTCASE=41.1
  info "testcase:$TESTCASE - annotation RNA"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR_RNA}/sequin_reads.blow5"
  ALIGNMENT="${REL_PATH}/data/raw/plot/realigned_RNA/one_read/sorted_test_1.5.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  BED="--bed ${REL_PATH}/data/raw/plot/bed_files/RNA_2.bed"
  python src/plot.py --rna ${BED} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" || die "testcase:$TESTCASE failed"

}
testcase_42s() {
  TESTCASE=42.0
  [ "${HUMAN_GENOME}" ] || die "Please set the env variable to the human genome path. export HUMAN_GENOME=path/to/file"
  GENOME="${HUMAN_GENOME}"
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA="-f ${GENOME}"
  SIGNAL="-s ${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="-a ${RAW_DIR}/f5c_eventalign/sorted_eventalign.paf.gz"
  OUTPUT="-o ${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT="--plot_limit 10"
  BED="--bed ${RAW_DIR}/bed_files/DNA_2.bed"
  REGION="--region chr1:6,811,011-6,811,198"
  SCALING="--sig_scale znorm"
  python src/plot_pileup.py ${BED} ${REGION} ${FASTA} ${SIGNAL} ${ALIGNMENT} ${OUTPUT} ${PLOT_LIMIT} ${SCALING} --tag_name "testcase-${TESTCASE}" || die "testcase:$TESTCASE failed"

  TESTCASE=42.1
  [ "${HUMAN_GENOME}" ] || die "Please set the env variable to the human genome path. export HUMAN_GENOME=path/to/file"
  GENOME="${HUMAN_GENOME}"
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA="-f ${GENOME}"
  SIGNAL="-s ${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="-a ${RAW_DIR}/f5c_eventalign/sorted_eventalign.paf.gz"
  OUTPUT="-o ${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT="--plot_limit 10"
  BED="--bed ${RAW_DIR}/bed_files/DNA_2.bed"
  REGION="--region chr1:6,811,011-6,811,198"
  SCALING="--sig_scale znorm"
  python src/plot_pileup.py --overlap_bottom ${BED} ${REGION} ${FASTA} ${SIGNAL} ${ALIGNMENT} ${OUTPUT} ${PLOT_LIMIT} ${SCALING} --tag_name "testcase-${TESTCASE}" || die "testcase:$TESTCASE failed"

  TESTCASE=42.2
  [ "${HUMAN_GENOME}" ] || die "Please set the env variable to the human genome path. export HUMAN_GENOME=path/to/file"
  GENOME="${HUMAN_GENOME}"
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA="-f ${GENOME}"
  SIGNAL="-s ${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="-a ${RAW_DIR}/f5c_eventalign/sorted_eventalign.bam"
  OUTPUT="-o ${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT="--plot_limit 10"
  BED="--bed ${RAW_DIR}/bed_files/DNA_2.bed"
  REGION="--region chr1:6,811,011-6,811,198"
  SCALING="--sig_scale znorm"
  python src/plot_pileup.py ${BED} ${REGION} ${FASTA} ${SIGNAL} ${ALIGNMENT} ${OUTPUT} ${PLOT_LIMIT} ${SCALING} --tag_name "testcase-${TESTCASE}" || die "testcase:$TESTCASE failed"

  TESTCASE=42.3
  [ "${HUMAN_GENOME}" ] || die "Please set the env variable to the human genome path. export HUMAN_GENOME=path/to/file"
  GENOME="${HUMAN_GENOME}"
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA="-f ${GENOME}"
  SIGNAL="-s ${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="-a ${RAW_DIR}/f5c_eventalign/sorted_eventalign.bam"
  OUTPUT="-o ${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT="--plot_limit 10"
  BED="--bed ${RAW_DIR}/bed_files/DNA_2.bed"
  REGION="--region chr1:6,811,011-6,811,198"
  SCALING="--sig_scale znorm"
  python src/plot_pileup.py --overlap_bottom ${BED} ${REGION} ${FASTA} ${SIGNAL} ${ALIGNMENT} ${OUTPUT} ${PLOT_LIMIT} ${SCALING} --tag_name "testcase-${TESTCASE}" || die "testcase:$TESTCASE failed"

}
testcase_43s() {
  GENOME="${REL_PATH}/data/raw/plot/reference_genomes/rnasequin_sequences_2.4.fa"

  TESTCASE=43.0
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA="-f ${GENOME}"
  SIGNAL="-s ${RAW_DIR}/f5c_eventalign/rna/reads.blow5"
  ALIGNMENT="-a ${RAW_DIR}/f5c_eventalign/rna/sorted_eventalign.bam"
  OUTPUT="-o ${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT="--plot_limit 10"
  REGION="--region R1_92_1:245-284"
#  SCALING="--sig_scale znorm"
  BED="--bed ${RAW_DIR}/bed_files/RNA_3.bed"
  python src/plot_pileup.py --rna ${BED} ${REGION} ${FASTA} ${SIGNAL} ${ALIGNMENT} ${OUTPUT} ${PLOT_LIMIT} ${SCALING} --tag_name "testcase-${TESTCASE}" || die "testcase:$TESTCASE failed"

  TESTCASE=43.1
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA="-f ${GENOME}"
  SIGNAL="-s ${RAW_DIR}/f5c_eventalign/rna/reads.blow5"
  ALIGNMENT="-a ${RAW_DIR}/f5c_eventalign/rna/sorted_eventalign.bam"
  OUTPUT="-o ${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT="--plot_limit 10"
  REGION="--region R1_92_1:264-500"
  SCALING="--sig_scale znorm"
  BED="--bed ${RAW_DIR}/bed_files/RNA_3.bed"
  python src/plot_pileup.py --rna ${BED} ${REGION} ${FASTA} ${SIGNAL} ${ALIGNMENT} ${OUTPUT} ${PLOT_LIMIT} ${SCALING} --tag_name "testcase-${TESTCASE}" || die "testcase:$TESTCASE failed"
}
testcase_40s #annotation DNA
testcase_41s #annotation RNA
testcase_42s #annotation DNA pileup
testcase_43s #annotation RNA pileup

info "all testcases passed"
#rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0
