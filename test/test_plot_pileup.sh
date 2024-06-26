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
OUTPUT_DIR="${REL_PATH}/data/out/plot_sig_ref_pileup"
test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"

RAW_DIR="${REL_PATH}/data/raw/plot"
EXP_DIR="${REL_PATH}/data/exp/plot"

testcase_20s() {
  [ "${HUMAN_GENOME}" ] || die "Please set the env variable to the human genome path. export HUMAN_GENOME=path/to/file"
  GENOME="${HUMAN_GENOME}"

  TESTCASE=20.0
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="${RAW_DIR}/f5c_eventalign/sorted_eventalign.paf.gz"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  SCALING="znorm"
  squigualiser plot_pileup -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} --sig_scale ${SCALING} && die "testcase:$TESTCASE failed"

  TESTCASE=20.1
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="${RAW_DIR}/f5c_eventalign/sorted_eventalign.paf.gz"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  REGION="chr1:6,811,011-6,811,198"
  SCALING="znorm"
  BASE_SHIFT="--base_shift -6"
  squigualiser plot_pileup ${BASE_SHIFT} --region ${REGION} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} --sig_scale ${SCALING} || die "testcase:$TESTCASE failed"

  TESTCASE=20.2
  info "testcase:$TESTCASE - reference-signal plot"
  GENOME_2="test/data/raw/plot/reference_genomes/chr22_1_5k.fa"
  FASTA=${GENOME_2}
  SIGNAL="${RAW_DIR}/pileup_view/test_0/na12878_chr22.blow5"
  ALIGNMENT="${RAW_DIR}/pileup_view/test_0/na12878_chr22.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  REGION="chr22_1_5k:4000-4500"
  SCALING="znorm"
  BASE_SHIFT="--base_shift -6"
  squigualiser plot_pileup ${BASE_SHIFT} --region ${REGION} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} --sig_scale ${SCALING} || die "testcase:$TESTCASE failed"

  TESTCASE=20.3
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="${RAW_DIR}/f5c_eventalign/sorted_eventalign.paf.gz"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  REGION="chr1:6,811,011-6,811,198"
  SCALING="znorm"
  BASE_SHIFT="--base_shift -6"
  squigualiser plot_pileup --no_overlap ${BASE_SHIFT} --region ${REGION} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} --sig_scale ${SCALING} || die "testcase:$TESTCASE failed"

  TESTCASE=20.4
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="${RAW_DIR}/f5c_eventalign/sorted_eventalign.paf.gz"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  REGION="chr1:6,811,011-6,811,198"
  SCALING="znorm"
  BASE_SHIFT="--base_shift -6"
  squigualiser plot_pileup --overlap_only ${BASE_SHIFT} --region ${REGION} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} --sig_scale ${SCALING} || die "testcase:$TESTCASE failed"

  TESTCASE=20.5
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="${RAW_DIR}/f5c_eventalign/sorted_eventalign.paf.gz"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  REGION="chr1:6,811,011-6,811,198"
  SCALING="znorm"
  BASE_SHIFT="--base_shift -6"
  squigualiser plot_pileup --overlap_only ${BASE_SHIFT} --region ${REGION} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} --sig_scale ${SCALING} && die "testcase:$TESTCASE failed"

  TESTCASE=20.6
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/realigned_DNA/reads.blow5"
  ALIGNMENT="${RAW_DIR}/realigned_DNA/sorted_test_1.4.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  REGION="chr1:6811404-6811443"
  SCALING="znorm"
  BASE_SHIFT="--base_shift 0"
  squigualiser plot_pileup ${BASE_SHIFT} --region ${REGION} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} --sig_scale ${SCALING} || die "testcase:$TESTCASE failed"

  TESTCASE=20.7
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/realigned_DNA/reads.blow5"
  ALIGNMENT="${RAW_DIR}/realigned_DNA/sorted_test_1.4.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  REGION="chr1:6811404-6811443"
  SCALING="znorm"
  BASE_SHIFT="--base_shift 0"
  squigualiser plot_pileup --plot_reverse ${BASE_SHIFT} --region ${REGION} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} --sig_scale ${SCALING} || die "testcase:$TESTCASE failed"

  TESTCASE=20.8
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="${RAW_DIR}/f5c_eventalign/sorted_eventalign.paf.gz"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  REGION="chr1:6,811,011-6,811,198"
  SCALING="znorm"
  BASE_SHIFT="--base_shift 0"
  squigualiser plot_pileup ${BASE_SHIFT} --region ${REGION} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} --sig_scale ${SCALING} || die "testcase:$TESTCASE failed"

  TESTCASE=20.8.1
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="${RAW_DIR}/f5c_eventalign/sorted_eventalign.paf.gz"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  REGION="chr1:6,811,011-6,811,198"
  SCALING="znorm"
  BASE_SHIFT="--base_shift 0"
  squigualiser plot_pileup --auto --region ${REGION} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} --sig_scale ${SCALING} || die "testcase:$TESTCASE failed"


  TESTCASE=20.9
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="${RAW_DIR}/f5c_eventalign/sorted_eventalign.paf.gz"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  REGION="chr1:6,811,011-6,811,198"
  SCALING="znorm"
  PROFILE="--profile kmer_model_dna_r10.4.1_e8.2_400bps_9_mer"
  squigualiser plot_pileup ${PROFILE} --region ${REGION} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} --sig_scale ${SCALING} || die "testcase:$TESTCASE failed"

  TESTCASE=20.10
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/f5c_eventalign/reads.blow5"
  ALIGNMENT="${RAW_DIR}/f5c_eventalign/sorted_eventalign.paf.gz"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT=10
  REGION="chr1:6,811,011-6,811,198"
  SCALING="scaledpA"
  PROFILE="--profile kmer_model_dna_r10.4.1_e8.2_400bps_9_mer"
  squigualiser plot_pileup ${PROFILE} --region ${REGION} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" --plot_limit ${PLOT_LIMIT} --sig_scale ${SCALING} || die "testcase:$TESTCASE failed"


}
testcase_21s() {
  GENOME="${REL_PATH}/data/raw/plot/reference_genomes/rnasequin_sequences_2.4.fa"

  TESTCASE=21.0
  info "testcase:$TESTCASE - RNA"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/f5c_eventalign/rna/reads.blow5"
  ALIGNMENT="${RAW_DIR}/f5c_eventalign/rna/sorted_eventalign.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT="--plot_limit 10"
  SCALING="--sig_scale znorm"
  squigualiser plot_pileup -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" ${PLOT_LIMIT} ${SCALING} && die "testcase:$TESTCASE failed"

  TESTCASE=21.1
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/f5c_eventalign/rna/reads.blow5"
  ALIGNMENT="${RAW_DIR}/f5c_eventalign/rna/sorted_eventalign.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT="--plot_limit 10"
  REGION="R1_92_1:245-284"
  SCALING="--sig_scale znorm"
  BASE_SHIFT="--base_shift -6"
  squigualiser plot_pileup --rna ${BASE_SHIFT} --region ${REGION} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" ${PLOT_LIMIT} ${SCALING} || die "testcase:$TESTCASE failed"

  TESTCASE=21.2
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/f5c_eventalign/rna/reads.blow5"
  ALIGNMENT="${RAW_DIR}/f5c_eventalign/rna/sorted_eventalign.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT="--plot_limit 10"
  REGION="R1_92_1:264-786"
  SCALING="--sig_scale znorm"
  BASE_SHIFT="--base_shift -6"
  squigualiser plot_pileup --rna ${BASE_SHIFT} --region ${REGION} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" ${PLOT_LIMIT} ${SCALING} || die "testcase:$TESTCASE failed"

  TESTCASE=21.3
  info "testcase:$TESTCASE - reference-signal plot"
  FASTA=${GENOME}
  SIGNAL="${RAW_DIR}/f5c_eventalign/rna/reads.blow5"
  ALIGNMENT="${RAW_DIR}/f5c_eventalign/rna/sorted_eventalign.bam"
  OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
  PLOT_LIMIT="--plot_limit 10"
  REGION="R1_92_1:264-786"
  SCALING="--sig_scale znorm"
  BASE_SHIFT="--base_shift -6"
  squigualiser plot_pileup --cprofile --rna ${BASE_SHIFT} --region ${REGION} -f ${FASTA} -s ${SIGNAL} -a ${ALIGNMENT} -o ${OUTPUT} --tag_name "testcase-${TESTCASE}" ${PLOT_LIMIT} ${SCALING} || die "testcase:$TESTCASE failed"

}
testcase_20s #pileup DNA
testcase_21s #pileup RNA

info "all testcases passed"
#rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0
