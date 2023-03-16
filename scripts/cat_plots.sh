#!/bin/bash

RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo ; echo -e "${GREEN}$1${NC}" >&2 ; }

#getting bash arguments
Usage="./cat_plots.sh INPUT_PLOTS_DIR OUTPUT_HTML"
NUMBER_ARGS=2
if [[ "$#" -lt ${NUMBER_ARGS} ]]; then
	info "Usage: ${Usage}"
	exit 1
fi
#set -x
INPUT_DIR=${1}
OUTPUT_HTML=${2}

counter=0
for file in ${INPUT_DIR}/*.html
do
#  info ${file}
  cat ${file} >> ${OUTPUT_HTML}
  counter=$((counter+1))
done

info "${counter} html files were concatenated successfully."
exit 0
