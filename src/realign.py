import os
import argparse
from readpaf import parse_paf
import pysam

BAM_CMATCH, BAM_CINS, BAM_CDEL, BAM_CREF_SKIP, BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CPAD, BAM_CEQUAL, BAM_CDIFF, BAM_CBACK = range(10)

def run(args):

    print("input bam: " + args.bam)
    print("input paf: " + args.paf)
    print("output file: " + args.output)
    print("INFO: supplementary alignments will be skipped.")
    # print("INFO: move array of a reverse complemented alignment will be treated after reversing.")

    samfile = pysam.AlignmentFile(args.bam, mode='r')
    fout = ""
    if args.output[-4:] == ".bam":
        fout = pysam.AlignmentFile(args.output, "wb", template=samfile)
    elif args.output[-4:] == ".sam":
        fout = pysam.AlignmentFile(args.output, "w", template=samfile)
    elif args.output[-4:] == ".paf":
        if args.c:
            fout = open(args.output, "w")
        else:
            raise Exception("Error: please provide argument '-c' to write in PAF format")
    else:
        raise Exception("Error: please provide the output file with correct extension")

    # inefficient
    paf_file = open(args.paf, "r")
    paf_dic = {}
    for record in parse_paf(paf_file):
        paf_dic[record.query_name] = record

    processed_sam_record_count = 0
    for sam_record in samfile:
        if sam_record.is_unmapped or sam_record.is_supplementary or sam_record.is_secondary:
            continue
        # if sam_record.query_name != "80922061-df02-48ad-bd67-65e5adcc78f5" and sam_record.query_name != "a0d0047d-64d1-4cfd-9792-8c27ad9a8c06":
        #     continue
        sam_read_id = sam_record.query_name
        if sam_read_id not in paf_dic:
            raise Exception("Error: associated paf record is missing for the read id: {}".format(sam_read_id))
        paf_read_id = paf_dic[sam_read_id].query_name
        if paf_read_id != sam_read_id:
            raise Exception("Error: sam and paf read ids do not match")

        data_is_rna = False
        if paf_dic[sam_read_id].target_start > paf_dic[sam_read_id].target_end:  # if RNA start_kmer>end_kmer in paf
            data_is_rna = True
            if not args.rna:
                print("Info: data is detected as RNA")
                raise Exception("Error: data is not specified as RNA. Please provide the argument --rna ")
        if not data_is_rna and args.rna:
            raise Exception("Error: data is not not detected as RNA but the user specified as RNA. Please remove the argument --rna and check dataset")

        # print(sam_read_id)
        # print("sam_read.pos: " + str(sam_read.pos+1))
        # print(sam_read.cigarstring)
        cigar_t = sam_record.cigartuples
        if cigar_t is None:
            raise Exception("Error: cigartuples for sam record {} is an empty object".format(sam_read_id))
        # print(cigar_t)

        moves_string = paf_dic[sam_read_id].tags['ss'][2].rstrip(',').split(',')
        raw_start = paf_dic[sam_read_id].query_start
        raw_end = paf_dic[sam_read_id].query_end
        # print("{}\t{}\t{}\t{}\t{}\t{}".format(raw_start, raw_end, kmer_start, kmer_end, len(sam_record.seq), abs(kmer_end-kmer_start)))
        # print(sam_record.cigarstring)
        # print(sam_record.pos+1)
        strand_dir = "+"
        if sam_record.is_reverse:
            cigar_t.reverse()
            strand_dir = "-"
        if data_is_rna:
            cigar_t.reverse()
            # strand_dir = "RNA"

        idx = 0
        ss_string = ""
        count_bases = 0
        count_bases_seq = 0
        # print(len(moves_string))
        # print(paf_dic[sam_read_id].target_end)
        # print(paf_dic[sam_read_id].target_start)
        # print(len(sam_record.query_sequence))
        # if len(sam_record.query_sequence) != len(moves_string):
        #     raise Exception("Error: the sequence length does not match the number of moves")
        len_moves = len(moves_string)
        op_count = 0
        total_cig_count_seq = 0
        for a in cigar_t:
            cig_op = a[0]
            cig_count = a[1]
            if cig_op == BAM_CMATCH:
                if total_cig_count_seq + cig_count > len_moves:
                    cig_count = len_moves - total_cig_count_seq
                total_cig_count_seq += cig_count
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
                if total_cig_count_seq + cig_count > len_moves:
                    cig_count = len_moves - total_cig_count_seq
                total_cig_count_seq += cig_count
                signal_skip = 0
                for i in range(0, cig_count):
                    signal_skip = signal_skip + int(moves_string[idx])
                    idx = idx + 1
                ss_string = ss_string + str(signal_skip) + "I"
                count_bases_seq += cig_count
                # print(str(signal_skip) + " I BAM_CINS " + str(int(sam_read.pos) + 1 + count_bases))
            elif cig_op == BAM_CSOFT_CLIP:
                if total_cig_count_seq + cig_count > len_moves:
                    cig_count = len_moves - total_cig_count_seq
                total_cig_count_seq += cig_count
                signal_skip = 0
                for i in range(0, cig_count):
                    # print(str(idx) + " " + moves_string[idx])
                    signal_skip = signal_skip + int(moves_string[idx])
                    idx = idx + 1
                if op_count == 0:
                    raw_start += signal_skip
                else:
                    raw_end -= signal_skip
                # ss_string = ss_string + str(signal_skip) + "I"
                # print(str(signal_skip) + " I BAM_CSOFT_CLIP " + str(int(sam_read.pos) + 1 + count_bases))
            elif cig_op == BAM_CHARD_CLIP:
                continue
            elif cig_op == BAM_CREF_SKIP:
                continue
            else:
                raise Exception("Error: cigar operation [" + str(cig_op) + "]is not handled yet found in read : " + sam_read_id)
            op_count += 1

        kmer_start = sam_record.reference_start
        kmer_end = sam_record.reference_start + count_bases
        if args.c:
            if data_is_rna:
                paf_record_str = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tss:Z:{}\n".format(paf_dic[sam_read_id].query_name, paf_dic[sam_read_id].query_length, raw_start, raw_end, strand_dir, sam_record.reference_name, count_bases, kmer_end, kmer_start, count_bases, count_bases, "255", ss_string)
            else:
                paf_record_str = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tss:Z:{}\n".format(paf_dic[sam_read_id].query_name, paf_dic[sam_read_id].query_length, raw_start, raw_end, strand_dir, sam_record.reference_name, count_bases, kmer_start, kmer_end, count_bases, count_bases, "255", ss_string)
            fout.write(paf_record_str)
        else:
            #start_raw, end_raw, start_kmer and end_kmer,
            if data_is_rna:
                si_string = str(raw_start) + "," + str(raw_end) + "," + str(kmer_end) + "," + str(kmer_start)
            else:
                si_string = str(raw_start) + "," + str(raw_end) + "," + str(kmer_start) + "," + str(kmer_end)

            # print("{}\t{}\t{}\t{}\t{}\t{}".format(strand_dir, sam_read_id, si_string, count_bases, count_bases_seq, abs(kmer_end-kmer_start)))
            sam_record.set_tag("ss", ss_string, "Z")
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
    parser.add_argument('-c', action='store_true', help="write move table in paf format")
    parser.add_argument('-o', '--output', required=True, help="output reference-signal SAM/BAM/PAF file")
    parser.add_argument('--rna', required=False, action='store_true', help="specify for RNA reads")

    return parser

if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    try:
        run(args)
    except Exception as e:
        print(str(e))
        exit(1)
