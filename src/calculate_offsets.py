import os
import argparse
from readpaf import parse_paf
import pysam
import pyslow5
from pyfaidx import Fasta
from pyfastx import Fastq
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

DEFAULT_KMER_SIZE = 6
BASE_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
BASE_MAP = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
MAX_DIST_THRESHOLD = 0.0

def create_kmer_model(moves, sequence, raw_signal, kmer_length, sig_move_offset):
    len_seq = len(sequence)
    model = {}
    start_raw = 0
    for j in range(0, sig_move_offset):
        start_raw += int(moves[j])

    for i in range(0, len_seq-kmer_length + 1 - sig_move_offset):
        end_raw = start_raw + int(moves[i + sig_move_offset])
        value = raw_signal[start_raw : end_raw]
        start_raw = end_raw
        key = sequence[i:i+kmer_length]
        if key not in model:
            model[sequence[i:i+kmer_length]] = np.median(value)
        else:
            model[sequence[i:i+kmer_length]] = np.median(value)
    return model
def calculate_offset_values(model, kmer_length):
    test_array = []
    for offset in range(0, kmer_length):
        freq = [[], [], [], []]
        for kmer, value in model.items():
            freq[BASE_INDEX[kmer[offset]]].append(value)

        test_array.append(freq)
    return test_array
def plot_distributions(kmer_length, test_array, output_pdf, plt_title):
    start_offset = 0
    end_offset = kmer_length
    f, axes = plt.subplots(nrows=end_offset-start_offset, ncols=1, figsize=(12, 9))
    for offset in range(start_offset, end_offset):
        i = 0
        for base in test_array[offset]:
            if kmer_length == 1:
                sns.kdeplot(base, label=BASE_MAP[i])
            else:
                sns.kdeplot(base, label=BASE_MAP[i], ax=axes[offset-start_offset])
            i += 1
        if kmer_length == 1:
            axes.set_title('base offset: {}'.format(offset))
        else:
            axes[offset-start_offset].set_title('base offset: {}'.format(offset))
    plt.legend(prop={'size': 10}, title='Base')
    plt.suptitle("{}".format(plt_title), size=16)
    plt.draw()
    plt.savefig(output_pdf, format='pdf')
def calculate_distance(kmer_length, test_array):
    start_offset = 0
    end_offset = kmer_length
    offset_dist = []

    for offset in range(start_offset, end_offset):
        max_mean = -1
        min_mean = 10000
        for base in test_array[offset]:
            median = np.median(base)
            if median < min_mean:
                min_mean = median
            if median > max_mean:
                max_mean = median
        offset_dist.append(max_mean-min_mean)

    # for offset in range(start_offset, end_offset):
    #     std_total = 0
    #     for base in test_array[offset]:
    #         std = np.std(base)
    #         std_total += std
    #     offset_dist.append(1/std_total)

    # offset_dist = []
    # for offset in range(start_offset, end_offset):
    #     base_dist = []
    #     for base in test_array[offset]:
    #         median = np.median(base)
    #         base_dist.append(median)
    #     base_dist.sort()
    #     base_diff = [base_dist[n]-base_dist[n-1] for n in range(1, len(base_dist))]
    #     total_diff = 0
    #     for diff in base_diff:
    #         total_diff += diff
    #     offset_dist.append(total_diff)

    max_dist = -1
    max_offset = 0
    for offset in range(start_offset, end_offset):
        if offset_dist[offset] > max_dist:
            max_offset = offset
            max_dist = offset_dist[offset]
    return max_offset, max_dist
def calculate_offsets(args, sig_move_offset, output_pdf, s5):
    paf_file = args.paf
    sequence_file = args.file
    kmer_length = args.kmer_length

    use_fasta = 0
    if sequence_file[-6:] == ".fastq" or sequence_file[-6:] == ".fq.gz" or sequence_file[-3:] == ".fq":
        use_fasta = 0
    elif sequence_file[-6:] == ".fasta" or sequence_file[-3:] == ".fa":
        use_fasta = 1
    else:
        print("error please provide the sequence file with correct extension")
        return

    max_offset_arr = []
    max_dist_arr = []

    with open(paf_file, "r") as handle:
        if use_fasta:
            sequence_reads = Fasta(sequence_file)
        else:
            sequence_reads = Fastq(sequence_file)
        num_reads = 0
        for paf_record in parse_paf(handle):
            read_id = paf_record.query_name
            if args.read_id != "" and args.read_id != read_id:
                continue
            if use_fasta:
                fasta_seq = sequence_reads.get_seq(name=read_id, start=1, end=int(paf_record.target_length)).seq
            else:
                fasta_seq = sequence_reads[read_id].seq
                if len(fasta_seq) < paf_record.target_length:
                    print("Error: Sequence lengths mismatch. If {} is a multi-line fastq file convert it to a 4-line fastq using seqtk.".format(sequence_file))
                    return
            read = s5.get_read(read_id, pA=True)
            if read is None:
                print("Error: read not found")
                return
            if read is not None:
                start_index = paf_record.query_start
                end_index = read['len_raw_signal']
                y = read['signal'][start_index:end_index]
            moves_string = paf_record.tags['ss'][2]
            moves = re.split(r',+', moves_string)
            if not 'T' in fasta_seq and 'A' in fasta_seq and 'C' in fasta_seq and 'G' in fasta_seq:
                if 'N' in fasta_seq:
                    fasta_seq = fasta_seq.replace('N', 'T')
                elif 'U' in fasta_seq:
                    fasta_seq = fasta_seq.replace('U', 'T')
            if args.rna:
                fasta_seq = fasta_seq[::-1]

            kmer_model = create_kmer_model(moves, fasta_seq, y, kmer_length, sig_move_offset)
            test_array = calculate_offset_values(kmer_model, kmer_length)

            max_offset, max_dist = calculate_distance(kmer_length, test_array)
            # print("{}\tmax_offset:{}\tdist:{}".format(read_id, max_offset, max_dist))
            max_offset_arr.append(max_offset)
            max_dist_arr.append(max_dist)

            if output_pdf is not None:
                plt_title = "{} {}\nkmer_len:{} read_len: {}\nsig_move_offset:{} best_base_offset:{} max_dist:{}".format(args.tag_name, args.read_id, kmer_length, len(fasta_seq), sig_move_offset, max_offset, str(round(max_dist, 4)))
                print(str(round(max_dist, 4)))
                plot_distributions(kmer_length, test_array, output_pdf, plt_title)

            num_reads += 1
            if num_reads == args.read_limit:
                break
    if num_reads == 0:
        raise Exception("Error: no reads were processed. Check the dataset and the read_id if provided.")
    index = np.argsort(max_offset_arr)[len(max_offset_arr)//2]
    return max_offset_arr[index], max_dist_arr[index]
def run(args):
    if args.kmer_length < 1:
        raise Exception("Error: kmer length must be a positive integer")

    if args.use_model:
        if args.model == "":
            raise Exception("Error: please provide a model file")
        print("model file: {}".format(args.model))

        model = {}
        with open(args.model) as f:
            obj = [line.rstrip('\n') for line in f]
            header_line_count = 0
            kmer_length = 0
            while True:
                keyword = obj[header_line_count].split('\t')[0]
                if keyword == '#k':
                    kmer_length = int(obj[header_line_count].split('\t')[1])
                    header_line_count += 1
                elif keyword[0] == '#' or keyword == 'kmer':
                    header_line_count += 1
                else:
                    break
            model_ = obj[header_line_count:]
            for line in model_:
                values_ = line.split('\t')
                model[values_[0]] = float(values_[1])
            print(kmer_length)
            test_array = []
            for base_offset in range(0, kmer_length):
                freq = [[], [], [], []]
                for kmer, value in model.items():
                    freq[BASE_INDEX[kmer[base_offset]]].append(value)
                test_array.append(freq)
            max_offset, max_dist = calculate_distance(kmer_length, test_array)
            print("best_base_offset:{}\tdist:{}".format(max_offset, max_dist))
            if args.output != "":
                output_pdf = PdfPages(args.output)
                print("output file: {}".format(args.output))
                plt_title = "{}\nkmer_len:{}\nbest_base_offset:{} max_dist:{}".format(args.tag_name, kmer_length, max_offset, str(round(max_dist, 4)))
                plot_distributions(kmer_length, test_array, output_pdf, plt_title)
                output_pdf.close()
    else:
        if args.file == "":
            raise Exception("Error: please provide a sequence file")
        if args.paf == "":
            raise Exception("Error: please provide a paf file")
        if args.slow5 == "":
            raise Exception("Error: please provide a signal file")

        print("kmer_length: {}".format(args.kmer_length))
        print("sequence file: {}".format(args.file))
        print("paf file: {}".format(args.paf))
        print("signal file: {}".format(args.slow5))

        output_pdf = None
        if args.read_id != "" and args.output != "":
            output_pdf = PdfPages(args.output)
            print("output file: {}".format(args.output))

        s5 = pyslow5.Open(args.slow5, 'r')

        max_dist = 0
        max_dist_sig_move_offset = 0
        max_dist_base_offset = 0
        base_offset_zero_dist = 0
        base_offset_zero_sig_move_offset = 0
        for sig_move_offset in range(0, args.kmer_length):
            base_offset, dist = calculate_offsets(args, sig_move_offset, output_pdf, s5)
            if dist > max_dist + MAX_DIST_THRESHOLD:
                max_dist = dist
                max_dist_sig_move_offset = sig_move_offset
                max_dist_base_offset = base_offset
            if base_offset == 0:
                base_offset_zero_dist = dist
                base_offset_zero_sig_move_offset = sig_move_offset
        # print(base_offset_zero_sig_move_offset)
        # print(max_dist_sig_move_offset)
        # print(max_dist_base_offset)
        if base_offset_zero_sig_move_offset == max_dist_sig_move_offset - max_dist_base_offset:
            recommended_kmer_length = base_offset_zero_sig_move_offset+1
            recommended_sig_move_offset = base_offset_zero_sig_move_offset
            print("recommended kmer_length:{} recommended sig_move_offset:{}".format(recommended_kmer_length, recommended_sig_move_offset))
            if recommended_kmer_length == 1 and recommended_sig_move_offset == 0:
                print("It is not required to rerun reform.")
            else:
                print("It is recommended to rerun reform with the recommended values. (-k {} -m {})".format(recommended_kmer_length, recommended_sig_move_offset))
        else:
            print("please refer the generated output pdf file to confirm the following results.")
            print("max_dist:{} max_dist_sig_move_offset:{}".format(max_dist, max_dist_sig_move_offset))
            print("recommended kmer_length:{} recommended sig_move_offset:{}".format(max_dist_sig_move_offset+1, max_dist_sig_move_offset))

        if output_pdf is not None:
            output_pdf.close()
        s5.close()


def argparser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    parser.add_argument('-k', '--kmer_length', required=False, default=DEFAULT_KMER_SIZE, type=int, help="kmer length")
    parser.add_argument('-p', '--paf', required=False, type=str, default="", help="input read-signal alignment .paf file")
    parser.add_argument('-f', '--file', required=False, type=str, default="", help="fasta/fa/fastq/fq/fq.gz sequence file")
    parser.add_argument('-s', '--slow5', required=False, type=str, default="", help="slow5 file")
    parser.add_argument('--rna', required=False, action='store_true', help="specify for RNA reads")
    parser.add_argument('--model', required=False, type=str, default="",help="model file")
    parser.add_argument('--read_limit', required=False, type=int, default=100, help="limit the number of reads considered")
    parser.add_argument('-r', '--read_id', required=False, type=str, default="", help="plot the read with read_id")
    parser.add_argument('-o', '--output', required=False, type=str, default="", help="output .pdf file. (works only when a read_id is specified)")
    parser.add_argument('--tag_name', required=False, type=str, default="", help="a tag name to easily identify")
    parser.add_argument('--use_model', required=False, action='store_true', help="calculate offset using the model file. else use the move table.")
    return parser

if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    try:
        run(args)
    except Exception as e:
        print(str(e))
        exit(1)
