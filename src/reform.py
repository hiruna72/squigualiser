import os
import argparse
from readpaf import parse_paf
import pysam

DEFAULT_KMER_SIZE = 9
DEFAULT_SIG_MOVE_OFFSET = 1
EXPECTED_STRIDE=5

parser = argparse.ArgumentParser()

parser.add_argument('-k', '--kmer_length', required=False, default=DEFAULT_KMER_SIZE, type=int, help="kmer length")
parser.add_argument('-m', '--sig_move_offset', required=False, default=DEFAULT_SIG_MOVE_OFFSET, type=int, help="signal move offset")
parser.add_argument('-c', action='store_true', help="write move table in paf format")
parser.add_argument('-b', '--bam', required=True, help="input SAM/BAM file produced by the basecaller")
parser.add_argument('-o', '--output', required=True, help="output .tsv/.paf file")
args = parser.parse_args()

if args.kmer_length < 1:
    print("kmer length must be a positive integer")
    exit(1)
if args.sig_move_offset < 1:
    print("signal move offset must be a positive integer")
    exit(1)

if args.kmer_length < args.sig_move_offset:
    print("signal move offset value must not be larger than the kmer length.")
    exit(1)

if (args.c and args.output[-4:] != ".paf") or (not args.c and args.output[-4:] != ".tsv"):
    print("error please provide the output file with correct extension (.tsv/.paf)")
    exit(1)

print("kmer_length: " + str(args.kmer_length))
print("sig_move_offset: " + str(args.sig_move_offset))
print("input bam: " + args.bam)
if args.c:
    print("output format: " + "paf")
else:
    print("output format: " + "tsv")
print("output file: " + args.output)


samfile = pysam.AlignmentFile(args.bam, mode='r', check_sq=False)
fout = open(args.output, "w")

for sam_record in samfile:

    len_seq = len(sam_record.get_forward_sequence()) - args.kmer_length + 1 # to get the number of kmers

    if not sam_record.has_tag("ns"):
        print("tag '{}' is not found. Please check your input SAM/BAM file.".format("ns"))
        exit(1)
    if not sam_record.has_tag("ts"):
        print("tag '{}' is not found. Please check your input SAM/BAM file.".format("ts"))
        exit(1)
    if not sam_record.has_tag("mv"):
        print("tag '{}' is not found. Please check your input SAM/BAM file.".format("mv"))
        exit(1)

    ns = int(sam_record.get_tag("ns"))
    ts = int(sam_record.get_tag("ts"))
    mv = sam_record.get_tag("mv")

    # print("ns: " + str(ns))
    # print("ts: " + str(ts))
    # print(mv[:5])

    len_mv = len(mv)
    if len_mv == 0:
        print("mv array length is 0.")
        exit(1)

    if mv[0] != EXPECTED_STRIDE:
        print("expected stride of {} is missing.".format(EXPECTED_STRIDE))
        exit(1)
    if not args.c:
        move_count = 0
        i = 1
        while move_count < args.sig_move_offset:
            value = mv[i]
            if value == 1:
                move_count += 1
            i += 1
        end_idx = ts + (i - 1) * EXPECTED_STRIDE;
        start_idx = end_idx - EXPECTED_STRIDE;
        kmer_idx = 0

        while i < len_mv:
            value = mv[i]
            if len_seq > 0 and value == 1:
                fout.write("{}\t".format(sam_record.query_name))
                fout.write("{}\t".format(kmer_idx))
                fout.write("{}\t".format(start_idx))
                fout.write("{}\n".format(end_idx))
                start_idx = end_idx
                kmer_idx += 1
                len_seq -= 1
            if len_seq > 0 and i == len_mv-1:
                fout.write("{}\t".format(sam_record.query_name))
                fout.write("{}\t".format(kmer_idx))
                fout.write("{}\t".format(start_idx))
                fout.write("{}\n".format(ns))
                start_idx = ns
                kmer_idx += 1
                len_seq -= 1

            end_idx = end_idx + EXPECTED_STRIDE
            i += 1
        if len_seq != 0:
            print("Error in the implementation. Please report the command with minimal reproducible data. Read_id: {}".format(sam_record.query_name));
            exit(1)

    # write paf format
    else:
        fout.write("{}\t".format(sam_record.query_name)) #1 read_id
        fout.write("{}\t".format(ns)) #2 Raw signal length (number of samples)
        move_count = 0
        i = 1
        start_idx = 0
        kmer_idx = 0

        while move_count < args.sig_move_offset:
            value = mv[i]
            if value == 1:
                move_count += 1
                start_idx = i
            i += 1

        fout.write("{}\t".format(ts + (i - 2) * EXPECTED_STRIDE)) #3 Raw signal start index (0-based; BED-like; closed)

        j = 1
        l_end_raw = 0
        len_seq_1 = len_seq + args.sig_move_offset
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
            l_end_raw =  ts+(end_idx-1) * EXPECTED_STRIDE

        fout.write("{}\t".format(l_end_raw)) #4 Raw signal end index (0-based; BED-like; open)
        fout.write("{}\t".format("+")) #5 Relative strand: “+” or “-“
        fout.write("{}\t".format(sam_record.query_name)) #6 Same as column 1
        fout.write("{}\t".format(len_seq + args.kmer_length - 1)) #7 base-called sequence length (no. of k-mers)
        fout.write("{}\t".format(kmer_idx)) #8 k-mer start index on basecalled sequence (0-based)
        fout.write("{}\t".format(len_seq - 1)) #9 k-mer end index on basecalled sequence (0-based)
        fout.write("{}\t".format(len_seq - kmer_idx)) #10 Number of k-mers matched on basecalled sequence
        fout.write("{}\t".format(len_seq + args.kmer_length - 1)) #11 Same as column 7
        fout.write("{}\t".format("255")) #12 Mapping quality (0-255; 255 for missing)
        fout.write("{}".format("ss:Z:")) #12 Mapping quality (0-255; 255 for missing)

        while i < len_mv:
            value = mv[i]
            # print("{}\t{}".format(i, value))
            if len_seq > 0 and value == 1:
                fout.write("{},".format((i-start_idx)*EXPECTED_STRIDE))  # ss
                start_idx = i
                len_seq -= 1
            if len_seq > 0 and i == len_mv-1:
                if (ns - ((i-1)*EXPECTED_STRIDE + ts)) < 0:
                    print("Error in calcuation. (ns - ((i-1)*EXPECTED_STRIDE + ts)) > 0 is not valid")
                    exit(1)
                len_seq -= 1
                l_duration = ((i-start_idx)*EXPECTED_STRIDE) + (ns - ((i-1)*EXPECTED_STRIDE + ts))
                fout.write("{},".format(l_duration))  # ss
            i += 1

        if len_seq != 0:
            print("Error in the implementation. Please report the command with minimal reproducible data. Read_id: {}".format(sam_record.query_name));
            exit(1)

        fout.write("{}".format("\n"))  # newline

samfile.close()
fout.close()