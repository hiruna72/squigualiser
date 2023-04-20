#!/bin/bash

# steps
# minimap2 align
# f5c index
# f5c eventalign 

RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo ; echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

set -x

REF=/media/hiruna/data/basecalling_work/apply_variants_to_genome/genome/hg38noAlt.fa #reference
MAP_SAM=mapped.sam
MAP_BAM=mapped.bam
FASTQ=reads.fastq
# samtools fastq basecaller_out.sam > ${FASTQ}
# minimap2 -ax map-ont ${REF} ${FASTQ} -t8 --secondary=no -o ${MAP_SAM}

# samtools sort ${MAP_SAM} -o ${MAP_BAM}
# samtools index ${MAP_BAM}

SIGNAL=reads.blow5
f5c index ${FASTQ} --slow5 ${SIGNAL}

ALIGNMENT=sorted_eventalign.paf.gz
f5c eventalign -b ${MAP_BAM} -r ${FASTQ} -g ${REF} --slow5 ${SIGNAL} -c -o eventalign.paf --supplementary no\
   && sort -k6,6 -nk8,8 eventalign.paf -o sorted_eventalign.paf \
   && bgzip sorted_eventalign.paf \
   && tabix -0 -b 8 -e 9 -s 6 ${ALIGNMENT}


info "success"