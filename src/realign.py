import os
import argparse
from readpaf import parse_paf
import pysam

BAM_CMATCH = 0
BAM_CINS = 1
BAM_CDEL = 2
BAM_CREF_SKIP = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD = 6
BAM_CEQUAL = 7
BAM_CDIFF = 8
BAM_CBACK = 9

parser = argparse.ArgumentParser()

parser.add_argument('-p', '--paf', required=True, help="input read-signal alignment .paf file")
parser.add_argument('-b', '--bam', required=True, help="input read-reference alignment SAM/BAM file")
parser.add_argument('-o', '--output', required=True, help="output reference-signal SAM/BAM file")
args = parser.parse_args()

print("input bam: " + args.bam)
print("input paf: " + args.paf)
print("output bam: " + args.output)
print("INFO: supplementary alignments will be skipped.")
print("INFO: move array of a reverse complemented alignment will be treated after reversing.")

samfile = pysam.AlignmentFile(args.bam, mode='r')
fout = ""
if args.output[-4:] == ".bam":
    fout = pysam.AlignmentFile(args.output, "wb", template=samfile)
elif args.output[-4:] == ".sam":
    fout = pysam.AlignmentFile(args.output, "w", template=samfile)
else:
    print("error please provide the output file with correct extension")
    exit()

# inefficient
paf_file = open(args.paf, "r")
paf_dic = {}
for record in parse_paf(paf_file):
    paf_dic[record.query_name] = record

processed_sam_record_count = 0
for sam_record in samfile:
    if sam_record.is_supplementary:
        continue
    # if sam_record.query_name != "63d4a555-4706-4295-b6aa-172b71b9c18f":
    #     continue
    sam_read_id = sam_record.query_name
    if sam_read_id not in paf_dic:
        print("Error associated paf record is missing for the read id: {}".format(sam_read_id))
    paf_read_id = paf_dic[sam_read_id].query_name
    if paf_read_id != sam_read_id:
        print("sam and paf read ids do not match")
        exit(1)
    # print(sam_read_id)
    # print("sam_read.pos: " + str(sam_read.pos+1))
    # print(sam_read.cigarstring)
    cigar_t = sam_record.cigartuples
    # print(cigar_t)

    moves_string = paf_dic[sam_read_id].tags['ss'][2].rstrip(',').split(',')
    if sam_record.is_reverse:
        moves_string.reverse()
    # print(moves_string)

    idx = 0
    ss_string = ""
    count_bases = 0
    for a in cigar_t:
        cig_op = a[0]
        cig_count = a[1]
        if cig_op == BAM_CMATCH:
            for i in range(0, cig_count):
                ss_string = ss_string + moves_string[idx] + ","
                idx = idx + 1
            count_bases += cig_count
        elif cig_op == BAM_CDEL:
            ss_string = ss_string + str(cig_count) + "D"
            count_bases += cig_count
            # print(str(cig_count) + " D " + str(int(sam_read.pos) + 1 + count_bases))
        elif cig_op == BAM_CINS:
            signal_skip = 0
            for i in range(0, cig_count):
                signal_skip = signal_skip + int(moves_string[idx])
                idx = idx + 1
            ss_string = ss_string + str(signal_skip) + "I"
            # print(str(signal_skip) + " I BAM_CINS " + str(int(sam_read.pos) + 1 + count_bases))
        elif cig_op == BAM_CSOFT_CLIP:
            signal_skip = 0
            for i in range(0, cig_count):
                # print(str(idx) + " " + moves_string[idx])
                signal_skip = signal_skip + int(moves_string[idx])
                idx = idx + 1
            ss_string = ss_string + str(signal_skip) + "I"
            # print(str(signal_skip) + " I BAM_CSOFT_CLIP " + str(int(sam_read.pos) + 1 + count_bases))
        elif cig_op == BAM_CHARD_CLIP:
            continue
        else:
            print("error: cigar operation [" + str(cig_op) + "]is not handled yet")

    sam_record.set_tag("ss", ss_string, "Z")
    rq_paf_string = ""
    rq_paf_string = rq_paf_string + paf_dic[paf_read_id].query_name + ","
    rq_paf_string = rq_paf_string + str(paf_dic[paf_read_id].query_length) + ","
    rq_paf_string = rq_paf_string + str(paf_dic[paf_read_id].query_start) + ","
    rq_paf_string = rq_paf_string + str(paf_dic[paf_read_id].query_end) + ","
    rq_paf_string = rq_paf_string + paf_dic[paf_read_id].strand + ","
    rq_paf_string = rq_paf_string + paf_dic[paf_read_id].target_name + ","
    rq_paf_string = rq_paf_string + str(paf_dic[paf_read_id].target_length) + ","
    rq_paf_string = rq_paf_string + str(paf_dic[paf_read_id].target_start) + ","
    rq_paf_string = rq_paf_string + str(paf_dic[paf_read_id].target_end) + ","
    rq_paf_string = rq_paf_string + str(paf_dic[paf_read_id].residue_matches) + ","
    rq_paf_string = rq_paf_string + str(paf_dic[paf_read_id].alignment_block_length) + ","
    rq_paf_string = rq_paf_string + str(paf_dic[paf_read_id].mapping_quality)
    # print(rq_paf_string)
    sam_record.set_tag("rq", rq_paf_string, "Z")
    fout.write(sam_record)
    processed_sam_record_count += 1

print("processed_sam_record_count: " + str(processed_sam_record_count))

paf_file.close()
samfile.close()
fout.close()
