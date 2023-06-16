#!/bin/bash

# steps
# basecall
# minimap2 align
# f5c index
# f5c eventalign 

RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo ; echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

set -x

REF="../../reference_genomes/rnasequin_sequences_2.4.fa" #reference
# MAP_SAM="mapq60_mapped.bam"
MAP_BAM="mapped.bam"
FASTQ="reads.fastq"
SIGNAL="reads.blow5"
ALIGNMENT_PAF="sorted_eventalign.paf.gz"
ALIGNMENT_BAM="sorted_eventalign.bam"

minimap2 -ax splice -uf -k14  ${REF} -t8 --secondary=no ${FASTQ} -o map.sam
samtools sort map.sam -o ${MAP_BAM}
samtools index ${MAP_BAM}

f5c index ${FASTQ} --slow5 ${SIGNAL}

# BGZIP="/media/hiruna/data/basecalling_work/htslib/bgzip"
# TABIX="/media/hiruna/data/basecalling_work/htslib/tabix"

# f5c eventalign -b ${MAP_BAM} -r ${FASTQ} -g ${REF} --slow5 ${SIGNAL} -c -o full_eventalign.paf --supplementary no
# sort -k6,6 -k9,9n eventalign.paf -o sorted_eventalign.paf
# ${BGZIP} sorted_eventalign.paf --keep -f
# ${TABIX} -0 -b 9 -e 8 -s 6 ${ALIGNMENT} -f

f5c eventalign -b ${MAP_BAM} -r ${FASTQ} -g ${REF} --slow5 ${SIGNAL} --sam -o eventalign.sam || die "f5c evenalign failed"
samtools sort eventalign.sam -o ${ALIGNMENT_BAM}
samtools index ${ALIGNMENT_BAM}


info "success"