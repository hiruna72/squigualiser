import os
import argparse
from readpaf import parse_paf
import pysam

DEFAULT_KMER_SIZE = 9
DEFAULT_SIG_MOVE_OFFSET = 0
DNA_STRIDE = 5
RNA_STRIDE = 10

profile_dic = {
        "guppy_dna_r9.4.1_450bps_fast_prom": [3, 2],
        "guppy_dna_r9.4.1_450bps_hac_prom": [3, 2],
        "guppy_dna_r9.4.1_450bps_sup_prom": [4, 3],
        "guppy_dna_r10.4.1_e8.2_400bps_fast": [1, 0],
        "guppy_dna_r10.4.1_e8.2_400bps_hac": [1, 0],
        "guppy_dna_r10.4.1_e8.2_400bps_sup": [1, 0]}

def list_profiles():
    # print(profile_dic)
    print("{}\t{}\t{}".format("name", "kmer_length", "sig_move_offset"))
    for profile in profile_dic:
        print("{}\t{}\t{}".format(profile, profile_dic[profile][0], profile_dic[profile][1]))
    print("If the profile you wanted is not listed here, please refer calculate_offsets.md on github to learn how to generate these values for your new data.")

def search_for_profile(profile):
    if profile in profile_dic:
        return profile_dic[profile]
    else:
        raise Exception("Error: specified profile ({}) is not found. Please run reform with -k 1 -s 0. Then run calculate_offsets.py and rerun reform with the recommended kmer_length and sig_move_offset.".format(profile))
def run(args):
    if args.list_profile:
        list_profiles()
        return
    else:
        if args.bam == "":
            raise Exception("Error: the following argument is required: -b/--bam")
        if args.output == "":
            raise Exception("Error: the following argument is required: -o/--output")

    if args.profile == "":
        kmer_length = args.kmer_length
        sig_move_offset = args.sig_move_offset
    else:
        kmer_length, sig_move_offset = search_for_profile(args.profile)
    if kmer_length < 1:
        raise Exception("Error: kmer length must be a positive integer")
    if sig_move_offset < 0:
        raise Exception("Error: signal move offset must not be less than zero")

    if kmer_length <= sig_move_offset:
        raise Exception("Error: signal move offset value must be smaller than the kmer length.")

    if (args.c and args.output[-4:] != ".paf") or (not args.c and args.output[-4:] != ".tsv"):
        raise Exception("Error: please provide the output file with correct extension (.tsv/.paf)")

    print("kmer_length: " + str(kmer_length))
    print("sig_move_offset: " + str(sig_move_offset))
    print("input bam: " + args.bam)
    if args.c:
        print("output format: " + "paf")
    else:
        print("output format: " + "tsv")
    print("output file: " + args.output)

    stride = DNA_STRIDE
    if args.rna:
        stride = RNA_STRIDE
    print("Info: Default stride: {}".format(stride))

    samfile = pysam.AlignmentFile(args.bam, mode='r', check_sq=False)
    fout = open(args.output, "w")
    processed_sam_record_count = 0

    for sam_record in samfile:

        len_seq = len(sam_record.get_forward_sequence()) - kmer_length + 1 # to get the number of kmers

        if not sam_record.has_tag("ns"):
            raise Exception("Error: tag '{}' is not found. Please check your input SAM/BAM file.".format("ns"))
        if not sam_record.has_tag("ts"):
            raise Exception("Error: tag '{}' is not found. Please check your input SAM/BAM file.".format("ts"))
        if not sam_record.has_tag("mv"):
            raise Exception("Error: tag '{}' is not found. Please check your input SAM/BAM file.".format("mv"))

        ns = int(sam_record.get_tag("ns"))
        ts = int(sam_record.get_tag("ts"))
        mv = sam_record.get_tag("mv")

        # print("ns: " + str(ns))
        # print("ts: " + str(ts))
        # print(mv[:5])

        len_mv = len(mv)
        if len_mv == 0:
            raise Exception("Error: mv array length is 0.")

        if mv[0] != stride:
            print("Info: Found stride to be {}, this value will be used.".format(mv[0]))
        stride = int(mv[0])
        if not args.c:
            move_count = 0
            i = 1
            while move_count < sig_move_offset + 1:
                value = mv[i]
                if value == 1:
                    move_count += 1
                i += 1
            end_idx = ts + (i - 1) * stride
            start_idx = end_idx - stride
            kmer_idx = 0
            if args.rna:
                kmer_idx = len_seq - 1

            while i < len_mv:
                value = mv[i]
                if len_seq > 0 and value == 1:
                    fout.write("{}\t".format(sam_record.query_name))
                    fout.write("{}\t".format(kmer_idx))
                    fout.write("{}\t".format(start_idx))
                    fout.write("{}\n".format(end_idx))
                    start_idx = end_idx
                    if args.rna:
                        kmer_idx -= 1
                    else:
                        kmer_idx += 1

                    len_seq -= 1
                if len_seq > 0 and i == len_mv-1:
                    fout.write("{}\t".format(sam_record.query_name))
                    fout.write("{}\t".format(kmer_idx))
                    fout.write("{}\t".format(start_idx))
                    fout.write("{}\n".format(ns))
                    start_idx = ns
                    if args.rna:
                        kmer_idx -= 1
                    else:
                        kmer_idx += 1
                    len_seq -= 1

                end_idx = end_idx + stride
                i += 1
            if len_seq != 0:
                raise Exception("Error: error in the implementation. Please report the command with minimal reproducible data. Read_id: {}".format(sam_record.query_name));

        # write paf format
        else:
            fout.write("{}\t".format(sam_record.query_name)) #1 read_id
            fout.write("{}\t".format(ns)) #2 Raw signal length (number of samples)
            move_count = 0
            i = 1
            start_idx = 0
            kmer_idx = 0

            while move_count < sig_move_offset + 1:
                value = mv[i]
                if value == 1:
                    move_count += 1
                    start_idx = i
                i += 1

            fout.write("{}\t".format(ts + (i - 2) * stride)) #3 Raw signal start index (0-based; BED-like; closed)

            j = 1
            l_end_raw = 0
            len_seq_1 = len_seq + sig_move_offset + 1
            end_idx = j + 1
            while j < len_mv:
                value = mv[j]
                if len_seq_1 > 0 and value == 1:
                    len_seq_1 -= 1
                    end_idx = j
                j += 1
            if len_seq_1 > 0 and j == len_mv:
                l_end_raw = ns
            else:
                l_end_raw = ts + (end_idx-1) * stride

            fout.write("{}\t".format(l_end_raw)) #4 Raw signal end index (0-based; BED-like; open)
            fout.write("+\t") #5 Relative strand: "+" or "-"
            fout.write("{}\t".format(sam_record.query_name)) #6 Same as column 1
            fout.write("{}\t".format(len_seq)) #7 base-called sequence length (no. of k-mers)

            if args.rna:
                fout.write("{}\t".format(len_seq))  # 8 k-mer start index on basecalled sequence (0-based)
                fout.write("{}\t".format(kmer_idx))  # 9 k-mer end index on basecalled sequence (0-based)
            else:
                fout.write("{}\t".format(kmer_idx))  # 8 k-mer start index on basecalled sequence (0-based)
                fout.write("{}\t".format(len_seq))  # 9 k-mer end index on basecalled sequence (0-based)
            fout.write("{}\t".format(len_seq - kmer_idx)) #10 Number of k-mers matched on basecalled sequence
            fout.write("{}\t".format(len_seq)) #11 Same as column 7
            fout.write("{}\t".format("255")) #12 Mapping quality (0-255; 255 for missing)
            fout.write("{}".format("ss:Z:")) #12 Mapping quality (0-255; 255 for missing)

            while i < len_mv:
                value = mv[i]
                # print("{}\t{}".format(i, value))
                if len_seq > 0 and value == 1:
                    fout.write("{},".format((i-start_idx) * stride))  # ss
                    start_idx = i
                    len_seq -= 1
                if len_seq > 0 and i == len_mv-1:
                    if (ns - ((i-1) * stride + ts)) < 0:
                        raise Exception("Error: error in calcuation. (ns - ((i-1)*EXPECTED_STRIDE + ts)) > 0 is not valid")
                    len_seq -= 1
                    l_duration = ((i-start_idx) * stride) + (ns - ((i - 1) * stride + ts))
                    fout.write("{},".format(l_duration))  # ss
                i += 1

            if len_seq != 0:
                raise Exception("Error: error in the implementation. Please report the command with minimal reproducible data. Read_id: {}".format(sam_record.query_name));

            fout.write("{}".format("\n"))  # newline
        processed_sam_record_count += 1

    samfile.close()
    fout.close()
    print("processed_sam_record_count: " + str(processed_sam_record_count))


def argparser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    parser.add_argument('-k', '--kmer_length', required=False, default=DEFAULT_KMER_SIZE, type=int, help="kmer length")
    parser.add_argument('-m', '--sig_move_offset', required=False, default=DEFAULT_SIG_MOVE_OFFSET, type=int, help="signal move offset")
    parser.add_argument('-c', action='store_true', help="write move table in paf format")
    parser.add_argument('-b', '--bam', required=False, default="", type=str, help="input SAM/BAM file produced by the basecaller")
    parser.add_argument('-o', '--output', required=False, default="", type=str, help="output .tsv/.paf file")
    parser.add_argument('--rna', required=False, action='store_true', help="specify for RNA reads")
    parser.add_argument('--profile', required=False, default="", type=str, help="determine -k and -m values using preset values")
    parser.add_argument('--list_profile', action='store_true', help="list the available profiles")
    return parser

if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    try:
        run(args)
    except Exception as e:
        print(str(e))
        exit(1)
