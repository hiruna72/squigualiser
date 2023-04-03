import os
import argparse
from readpaf import parse_paf
import pysam

BAM_CMATCH, BAM_CINS, BAM_CDEL, BAM_CREF_SKIP, BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CPAD, BAM_CEQUAL, BAM_CDIFF, BAM_CBACK = range(10)

def run(args):

    print("input bam: " + args.bam)
    print("input paf: " + args.paf)
    print("output bam: " + args.output)
    print("INFO: supplementary alignments will be skipped.")
    # print("INFO: move array of a reverse complemented alignment will be treated after reversing.")

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

        data_is_rna = False
        if paf_dic[sam_read_id].target_start > paf_dic[sam_read_id].target_end:  # if RNA start_kmer>end_kmer in paf
            data_is_rna = True
            print("Info: data is detected as RNA")
            if not args.rna:
                print("Error: data is not specified as RNA. Please provide the argument --rna ")
                exit(1)

        # print(sam_read_id)
        # print("sam_read.pos: " + str(sam_read.pos+1))
        # print(sam_read.cigarstring)
        cigar_t = sam_record.cigartuples
        if cigar_t is None:
            print("cigartuples for sam record {} is an empty object".format(sam_read_id))
            exit(1)
        # print(cigar_t)

        moves_string = paf_dic[sam_read_id].tags['ss'][2].rstrip(',').split(',')
        raw_start = paf_dic[sam_read_id].query_start
        raw_end = paf_dic[sam_read_id].query_end
        if data_is_rna:
            kmer_end = paf_dic[sam_read_id].target_start
            kmer_start = paf_dic[sam_read_id].target_end
        else:
            kmer_start = paf_dic[sam_read_id].target_start
            kmer_end = paf_dic[sam_read_id].target_end
        # print("{}\t{}\t{}\t{}\t{}\t{}".format(raw_start, raw_end, kmer_start, kmer_end, len(sam_record.seq), abs(kmer_end-kmer_start)))

        strand_dir = "+"
        if sam_record.is_reverse:
            cigar_t.reverse()
            strand_dir = "-"

        idx = 0
        ss_string = ""
        count_bases = 0
        count_bases_seq = 0

        if len(sam_record.query_sequence) != len(moves_string):
            print("Error: the sequence length does not match the number of moves")
            exit(1)

        op_count = 0
        for a in cigar_t:
            cig_op = a[0]
            cig_count = a[1]
            if cig_op == BAM_CMATCH:
                for i in range(0, cig_count):
                    ss_string = ss_string + moves_string[idx] + ","
                    idx = idx + 1
                count_bases += cig_count
                count_bases_seq += cig_count
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
                count_bases_seq += cig_count
                # print(str(signal_skip) + " I BAM_CINS " + str(int(sam_read.pos) + 1 + count_bases))
            elif cig_op == BAM_CSOFT_CLIP:
                signal_skip = 0
                for i in range(0, cig_count):
                    # print(str(idx) + " " + moves_string[idx])
                    signal_skip = signal_skip + int(moves_string[idx])
                    idx = idx + 1
                if sam_record.is_reverse:
                    if op_count == 0:
                        raw_end -= signal_skip
                        kmer_end -= cig_count
                    else:
                        raw_start += signal_skip
                        kmer_start += cig_count
                else:
                    if op_count == 0:
                        raw_start += signal_skip
                        kmer_start += cig_count
                    else:
                        raw_end -= signal_skip
                        kmer_end -= cig_count
                # ss_string = ss_string + str(signal_skip) + "I"
                # print(str(signal_skip) + " I BAM_CSOFT_CLIP " + str(int(sam_read.pos) + 1 + count_bases))
            elif cig_op == BAM_CHARD_CLIP:
                continue
            else:
                print("error: cigar operation [" + str(cig_op) + "]is not handled yet")
            op_count += 1


        sam_record.set_tag("ss", ss_string, "Z")

        #start_raw, end_raw, start_kmer and end_kmer,
        if data_is_rna:
            if not (raw_start < raw_end and kmer_start < kmer_end):
                print("Error: in implementation. Please report on github issues")
            si_string = str(raw_start) + "," + str(raw_end) + "," + str(kmer_end) + "," + str(kmer_start)
        else:
            if not (raw_start < raw_end and kmer_start < kmer_end):
                print("Error: in implementation. Please report on github issues")
            si_string = str(raw_start) + "," + str(raw_end) + "," + str(kmer_start) + "," + str(kmer_end)

        # print("{}\t{}\t{}\t{}\t{}\t{}".format(strand_dir, sam_read_id, si_string, count_bases, count_bases_seq, abs(kmer_end-kmer_start)))
        sam_record.set_tag("si", si_string, "Z")

        fout.write(sam_record)
        processed_sam_record_count += 1

    print("processed_sam_record_count: " + str(processed_sam_record_count))

    paf_file.close()
    samfile.close()
    fout.close()

def argparser():
    # parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    parser.add_argument('-p', '--paf', required=True, help="input read-signal alignment .paf file")
    parser.add_argument('-b', '--bam', required=True, help="input read-reference alignment SAM/BAM file")
    parser.add_argument('-o', '--output', required=True, help="output reference-signal SAM/BAM file")
    parser.add_argument('--rna', required=False, action='store_true', help="specify for RNA reads")

    return parser

if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    run(args)