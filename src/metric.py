"""
Signal to seQuence alignment Plot - plot
Hiruna Samarakoon - Garvan Medical Institute
hiruna@unsw.edu.au
"""
import numpy as np
import pyslow5
import argparse
import re
from readpaf import parse_paf
from pyfaidx import Fasta
from pyfastx import Fastq
import os
import pysam
from src import plot_utils
import statistics

# ref_start is always 1based closed
# ref_end is always 1based closed
# start_kmer is always 0based closed
# end_kmer is always 0based open

BASE_LIMIT = 1000
SIG_PLOT_LENGTH = 20000
PLOT_BASE_SHIFT = 0
PLOT_LIMIT = 1000
FIXED_BASE_WIDTH = 10

DEFAULT_NUM_BED_COLS = 3
DEFAULT_BED_ANNOTATION_COLOR = (75, 126, 246)
BAM_CMATCH, BAM_CINS, BAM_CDEL, BAM_CREF_SKIP, BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CPAD, BAM_CEQUAL, BAM_CDIFF, BAM_CBACK = range(10)
READ_ID, LEN_RAW_SIGNAL, START_RAW, END_RAW, STRAND, SEQUENCE_ID, LEN_KMER, START_KMER, END_KMER, MATCHES, LEN_KMER, MAPQ = range(12)
SI_START_RAW, SI_END_RAW, SI_START_KMER, SI_END_KMER = range(4)
BED_CHROM, BED_CHROM_START, BED_CHROM_END, BED_NAME, BED_SCORE, BED_STRAND, BED_THICK_START, BED_THICK_END, BED_ITEM_RGB, BED_BLOCK_COUNT, BED_BLOCK_SIZES, BLOCK_STARTS = range(12)

def get_metric(args, fout, metric_record, read_id, signal_tuple, sig_algn_data, fasta_sequence, base_limit, draw_data):
    x = signal_tuple[0]
    x_real = signal_tuple[1]
    y = signal_tuple[2]

    x_coordinate = 0
    initial_x_coordinate = x_coordinate

    # draw moves
    moves = sig_algn_data["ss"]
    base_index = sig_algn_data["start_kmer"]
    total_length_insertions = 0
    total_length_deletions = 0
    flag_base_index_bound = 0
    match_samples = []
    insertion_samples = []
    deletion_bases = []

    for i in moves:
        previous_x_coordinate = x_coordinate
        if 'D' in i:
            i = re.sub('D', '', i)
            n_samples = int(i)
            deletion_bases.append(n_samples)
            prev_x_cord = previous_x_coordinate
            for j in range(0, n_samples):
                base = fasta_sequence[base_index]
                base_index += 1
                total_length_deletions += 1
                if base_index - sig_algn_data["start_kmer"] == base_limit:
                    flag_base_index_bound = 1
                    break
            if flag_base_index_bound == 1:
                break
            x_coordinate = prev_x_cord
        elif 'I' in i:
            i = re.sub('I', '', i)
            n_samples = int(i)
            insertion_samples.append(n_samples)
            total_length_insertions += n_samples
            x_coordinate += n_samples
        else:
            n_samples = int(i)
            match_samples.append(n_samples)
            x_coordinate += n_samples
            base = fasta_sequence[base_index]
            base_index += 1
        if base_index - sig_algn_data["start_kmer"] == base_limit:
            break
        if x_coordinate - initial_x_coordinate > draw_data["sig_plot_limit"]:
            break

    if sig_algn_data["data_is_rna"] == 1:
        region_str = "{}:{}-{}".format(sig_algn_data["ref_name"], sig_algn_data["ref_end"], sig_algn_data["ref_end"] - base_index+1)
        # plot_title = f'{sig_algn_data["tag_name"]}[{sig_algn_data["ref_end"]:,}-{sig_algn_data["ref_end"] - base_index+1:,}]{indt}signal: [{int(x_real[0])}-{int(x_real[x_coordinate - 1])}]{indt}deletions(bases): {total_length_deletions} insertions(samples): {total_length_insertions}{indt}{read_id}{indt}signal dir:{draw_data["sig_dir"]}'
    else:
        region_str = "{}:{}-{}".format(sig_algn_data["ref_name"], sig_algn_data["ref_start"], sig_algn_data["ref_start"] + base_index-1)
        # plot_title = f'{sig_algn_data["tag_name"]}[{sig_algn_data["ref_start"]:,}-{sig_algn_data["ref_start"] + base_index-1:,}]{indt}signal: [{int(x_real[0])}-{int(x_real[x_coordinate - 1])}]{indt}deletions(bases): {total_length_deletions} insertions(samples): {total_length_insertions}{indt}{read_id}{indt}signal dir:{draw_data["sig_dir"]}'

    metric_record['region'] = region_str
    metric_record['total_matches'] = base_index - total_length_deletions
    metric_record['total_deletion_occurrences'] = len(deletion_bases)
    metric_record['total_insertion_occurrences'] = len(insertion_samples)
    metric_record['total_length_deletions'] = total_length_deletions
    metric_record['total_length_insertions'] = total_length_insertions
    if base_index - total_length_deletions > 1:
        metric_record['min_match'] = min(match_samples)
        metric_record['max_match'] = max(match_samples)
        metric_record['mode_match'] = statistics.mode(match_samples)
        metric_record['median_match'] = statistics.median(match_samples)
        metric_record['mean_match'] = statistics.mean(match_samples)
        metric_record['stdev_match'] = statistics.stdev(match_samples)
    else:
        metric_record['min_match'] = "."
        metric_record['max_match'] = "."
        metric_record['mode_match'] = "."
        metric_record['median_match'] = "."
        metric_record['mean_match'] = "."
        metric_record['stdev_match'] = "."
    if len(deletion_bases) > 1:
        metric_record['min_deletion'] = min(deletion_bases)
        metric_record['max_deletion'] = max(deletion_bases)
        metric_record['mode_deletion'] = statistics.mode(deletion_bases)
        metric_record['median_deletion'] = statistics.median(deletion_bases)
        metric_record['mean_deletion'] = statistics.mean(deletion_bases)
        metric_record['stdev_deletion'] = statistics.stdev(deletion_bases)
    else:
        metric_record['min_deletion'] = "."
        metric_record['max_deletion'] = "."
        metric_record['mode_deletion'] = "."
        metric_record['median_deletion'] = "."
        metric_record['mean_deletion'] = "."
        metric_record['stdev_deletion'] = "."
    if len(insertion_samples) > 1:
        metric_record['min_insertion'] = min(insertion_samples)
        metric_record['max_insertion'] = max(insertion_samples)
        metric_record['mode_insertion'] = statistics.mode(insertion_samples)
        metric_record['median_insertion'] = statistics.median(insertion_samples)
        metric_record['mean_insertion'] = statistics.mean(insertion_samples)
        metric_record['stdev_insertion'] = statistics.stdev(insertion_samples)
    else:
        metric_record['min_insertion'] = "."
        metric_record['max_insertion'] = "."
        metric_record['mode_insertion'] = "."
        metric_record['median_insertion'] = "."
        metric_record['mean_insertion'] = "."
        metric_record['stdev_insertion'] = "."

    fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
    metric_record['read_id'],
    metric_record['region'],
    metric_record['total_matches'],
    metric_record['total_deletion_occurrences'],
    metric_record['total_insertion_occurrences'],
    metric_record['total_length_deletions'],
    metric_record['total_length_insertions'],
    metric_record['min_match'],
    metric_record['max_match'],
    metric_record['mode_match'],
    metric_record['median_match'],
    metric_record['mean_match'],
    metric_record['stdev_match'],
    metric_record['min_deletion'],
    metric_record['max_deletion'],
    metric_record['mode_deletion'],
    metric_record['median_deletion'],
    metric_record['mean_deletion'],
    metric_record['stdev_deletion'],
    metric_record['min_insertion'],
    metric_record['max_insertion'],
    metric_record['mode_insertion'],
    metric_record['median_insertion'],
    metric_record['mean_insertion'],
    metric_record['stdev_insertion']))

    if args.extend_0:
        fout.write("\t{}\t{}\t{}".format(
        match_samples,
        deletion_bases,
        insertion_samples))

    fout.write("\n")

def run(args):
    if args.list_profile:
        plot_utils.list_profiles_base_shift()
        return
    else:
        if args.file == "":
            raise Exception("Error: the following argument is required: -f/--file")
        if args.slow5 == "":
            raise Exception("Error: the following argument is required: -s/--slow5")
        if args.alignment == "":
            raise Exception("Error: the following argument is required: -a/--alignment")
        if args.output == "":
            raise Exception("Error:the following argument is required: -o/--output")
    fout = open(args.output, "w")
    if args.read_id != "":
        args.plot_limit = 1

    use_fasta = 0
    if args.file:
        print(f'sequence file: {args.file}')
        if args.file[-6:] == ".fastq" or args.file[-6:] == ".fq.gz" or args.file[-3:] == ".fq":
            use_fasta = 0
        elif args.file[-6:] == ".fasta" or args.file[-3:] == ".fa":
            use_fasta = 1
        else:
            raise Exception("Error: please provide the sequence file with correct extension")

    plot_sig_ref_flag = 0
    if args.sig_ref:
        plot_sig_ref_flag = 1
    use_paf = 0
    if args.alignment:
        print(f'alignment file: {args.alignment}')
        alignment_extension = args.alignment[-4:]
        if alignment_extension == ".bam" or alignment_extension == ".sam":
            use_paf = 0
            plot_sig_ref_flag = 1
        elif alignment_extension == ".paf":
            use_paf = 1
            if args.sig_ref:
                raise Exception("Error: For signal to reference alignment the paf file must be bgzip compressed and tabix indexed")
        elif args.alignment[-7:] == ".paf.gz":
            use_paf = 1
            plot_sig_ref_flag = 1
            index_file = args.alignment + ".tbi"
            if not os.path.exists(index_file):
                raise Exception("Error: please provide a bgzip compressed and tabix indexed paf file to extract the regions")
        else:
            raise Exception("Error: please provide the alignment file with correct extension")

    if use_paf == 0 and use_fasta == 0:
        raise Exception("Error: please provide a .fasta or .fa file when using SAM/BAM")

    if args.base_limit:
        base_limit = args.base_limit
    else:
        base_limit = BASE_LIMIT
    print(f'signal file: {args.slow5}')

    if args.plot_reverse:
        print("Info: reads mapped to the reverse strand will be plotted")

    # open signal file
    s5 = pyslow5.Open(args.slow5, 'r')
    num_plots = 0
    indt = "\t\t\t\t\t\t\t\t"
    draw_data = {}
    draw_data["sig_plot_limit"] = args.sig_plot_limit
    draw_data["fixed_base_width"] = FIXED_BASE_WIDTH
    if args.profile == "":
        draw_data["base_shift"] = args.base_shift
    else:
        if args.plot_reverse:
            draw_data["base_shift"] = plot_utils.search_for_profile_base_shift(args.profile)[1]
        else:
            draw_data["base_shift"] = plot_utils.search_for_profile_base_shift(args.profile)[0]
    kmer_correction = 0
    if args.profile != "":
        kmer_correction = -1*(plot_utils.search_for_profile_base_shift(args.profile)[0] + plot_utils.search_for_profile_base_shift(args.profile)[1])

    metric_header = ("read_id\tregion\t"
                     "total_matches\ttotal_deletion_occurrences\ttotal_insertion_occurrences\t"
                     "total_length_deletions\ttotal_length_insertions\t"
                     "min_match\tmax_match\tmode_match\tmedian_matches\tmean_matches\tstdev_matches\t"
                     "min_deletion\tmax_deletion\tmode_deletion\tmedian_deletion\tmean_deletion\tstdev_deletion\t"
                     "min_insertion\tmax_insertion\tmode_insertion\tmedian_insertion\tmean_insertion\tstdev_insertion")
    if args.extend_0:
        metric_header += "\tmatches\tdeletions\tinsertions"
    fout.write("{}\n".format(metric_header))
    if use_paf == 1 and plot_sig_ref_flag == 0:
        print("Info: Signal to read method using PAF ...")
        with open(args.alignment, "r") as handle:
            if use_fasta:
                sequence_reads = Fasta(args.file)
            else:
                sequence_reads = Fastq(args.file)
            for paf_record in parse_paf(handle):
                metric_record = {}
                if paf_record.query_name != paf_record.target_name:
                    raise Exception("Error: this paf file is a signal to reference mapping. Please provide the argument --sig_ref ")
                read_id = paf_record.query_name
                if args.read_id != "" and read_id != args.read_id:
                    continue
                if read_id not in set(sequence_reads.keys()):
                    raise Exception("Error: read_id {} is not found in {}".format(read_id, args.file))
                if 'ss' not in paf_record.tags:
                    raise Exception("Error: ss string is missing for the read_id {} in {}".format(read_id, args.alignment))

                data_is_rna = 0
                if paf_record.target_start > paf_record.target_end:  # if RNA start_kmer>end_kmer in paf
                    data_is_rna = 1
                    if not args.rna:
                        print("Info: data is detected as RNA")
                        raise Exception("Error: data is not specified as RNA. Please provide the argument --rna ")

                fasta_seq = ""
                if use_fasta:
                    fasta_seq = sequence_reads[read_id][:].seq
                    fasta_seq = fasta_seq.upper()
                else:
                    fasta_seq = sequence_reads[read_id].seq
                    fasta_seq = fasta_seq.upper()
                    if len(fasta_seq) < paf_record.target_length:
                        raise Exception("Error: Sequence lengths mismatch. If {} is a multi-line fastq file convert it to a 4-line fastq using seqtk.".format(args.file))
                if not bool(re.match('^[ACGTUMRWSYKVHDBN]+$', fasta_seq)):
                    raise Exception("Error: base characters other than A,C,G,T/U,M,R,W,S,Y,K,V,H,D,B,N were detected. Please check your sequence files")

                ref_start = -1
                ref_end = -1
                if data_is_rna == 1:
                    fasta_seq = fasta_seq[paf_record.target_end:]
                    ref_start = paf_record.target_end + 1
                else:
                    fasta_seq = fasta_seq[paf_record.target_start:]
                    ref_start = paf_record.target_start + 1

                seq_len = len(fasta_seq)
                if seq_len < base_limit:
                    base_limit = seq_len
                ref_end = base_limit

                if args.region != "":
                    pattern = re.compile("^[0-9]+\-[0-9]+")
                    if not pattern.match(args.region):
                        raise Exception("Error: region provided is not in correct format")
                    ref_start = int(args.region.split("-")[0])
                    ref_end = int(args.region.split("-")[1])

                    if ref_start < 1:
                        raise Exception("Error: region start coordinate ({}) must be positive".format(ref_start))

                    if ref_end < 1:
                        raise Exception("Error: region end coordinate ({}) must be positive".format(ref_end))

                    if data_is_rna == 1 and paf_record.target_start < ref_end:
                        ref_end = paf_record.target_start
                    elif data_is_rna == 0 and paf_record.target_end < ref_end:
                        ref_end = paf_record.target_end

                    # print("ref_start: " + str(ref_start))
                    # print("ref_end: " + str(ref_end))

                    if (ref_end - ref_start + 1) < base_limit:
                        base_limit = ref_end - ref_start + 1

                # print("plot region: {}-{}\tread_id: {}".format(ref_start, ref_end, read_id))
                metric_record['read_id'] = read_id
                metric_record['region'] = "{}:{}-{}".format("", ref_start, ref_end)

                x = []
                x_real = []
                y = []
                read = s5.get_read(read_id, pA=args.no_pa, aux=["read_number", "start_mux"])
                if read is not None:
                    start_index = paf_record.query_start
                    end_index = read['len_raw_signal']
                    x = list(range(1, end_index - start_index + 1))
                    x_real = list(range(start_index + 1, end_index + 1))  # 1based
                    y = read['signal'][start_index:end_index]

                scaling_str = "no scaling"
                if args.sig_scale == "medmad" or args.sig_scale == "znorm" or args.sig_scale == "scaledpA":
                    scaling_str = args.sig_scale
                elif not args.sig_scale == "":
                    raise Exception("Error: given --sig_scale method: {} is not supported".format(args.sig_scale))

                scale_params = {}
                if args.sig_scale == "scaledpA":
                    if not args.no_pa:
                        raise Exception("Error: given --sig_scale method: {} required the signal to be converted to pA levels. Please remove --no_pa argument".format(args.sig_scale))
                    for tag in ["sc", "sh"]:
                        if tag not in paf_record.tags:
                            raise Exception("Error: required tag '{}' for given --sig_scale method: {} is not found in the alignment file".format(tag, args.sig_scale))
                        scale_params[tag] = paf_record.tags[tag][2]

                y = plot_utils.scale_signal(y, args.sig_scale, scale_params)

                strand_dir = "(DNA 5'->3')"
                if data_is_rna == 1:
                    fasta_seq = fasta_seq[:ref_end]
                    fasta_seq = fasta_seq[::-1]
                    strand_dir = "(RNA 3'->5')"
                else:
                    fasta_seq = fasta_seq[ref_start-1:]

                moves_string = paf_record.tags['ss'][2]
                moves_string = re.sub('D', 'D,', moves_string)
                moves_string = re.sub('I', 'I,', moves_string).rstrip(',')
                moves = re.split(r',+', moves_string)

                signal_tuple = (x, x_real, y)
                region_tuple = (ref_start, ref_end, 0, seq_len)

                sig_algn_dic = {}
                sig_algn_dic['start_kmer'] = 0
                sig_algn_dic['ref_start'] = ref_start
                sig_algn_dic['ref_end'] = ref_end
                sig_algn_dic['use_paf'] = use_paf
                sig_algn_dic['plot_sig_ref_flag'] = plot_sig_ref_flag
                sig_algn_dic['data_is_rna'] = data_is_rna
                sig_algn_dic['ss'] = moves
                sig_algn_dic['ref_name'] = ""

                signal_tuple, region_tuple, sig_algn_dic, fasta_seq = plot_utils.adjust_before_plotting(seq_len, signal_tuple, region_tuple, sig_algn_dic, fasta_seq)

                if draw_data["base_shift"] < 0:
                    abs_base_shift = abs(draw_data["base_shift"])
                    x = signal_tuple[0]
                    x_real = signal_tuple[1]
                    y = signal_tuple[2]
                    y_prefix = [np.nan] * abs_base_shift * draw_data["fixed_base_width"]
                    y = np.concatenate((y_prefix, y), axis=0)
                    x_real = np.concatenate(([1] * abs_base_shift * draw_data["fixed_base_width"], x_real), axis=0)
                    x = list(range(1, len(x) + 1 + abs_base_shift * draw_data["fixed_base_width"]))
                    signal_tuple = (x, x_real, y)
                    moves_prefix = [str(draw_data["fixed_base_width"])] * abs_base_shift
                    sig_algn_dic['ss'] = moves_prefix + sig_algn_dic['ss']

                draw_data['y_min'] = np.nanmin(y)
                draw_data['y_max'] = np.nanmax(y)

                get_metric(args=args, fout=fout, metric_record=metric_record, read_id=read_id, signal_tuple=signal_tuple, sig_algn_data=sig_algn_dic, fasta_sequence=fasta_seq, base_limit=base_limit, draw_data=draw_data)

                num_plots += 1
                if num_plots == args.plot_limit:
                    break
    elif use_paf == 0 and plot_sig_ref_flag == 1: # using sam/bam
        print("Info: Signal to reference method using SAM/BAM ...")
        fasta_reads = Fasta(args.file)
        if args.region != "":
            # check if there exists a .bam.bai
            index_file = args.alignment + ".bai"
            if not os.path.exists(index_file):
                raise Exception("Error: please provide a bam file that is sorted and indexed to extract the regions")

            args_region = re.sub(',', '', args.region)
            # print(args_region)
            # pattern = re.compile("^[a-z]+[0-9]+\:[0-9]+\-[0-9]+")
            pattern = re.compile("^.*\:[0-9]+\-[0-9]+")
            if not pattern.match(args_region):
                raise Exception("Error: region provided is not in correct format")
            args_ref_name = args_region.split(":")[0]
            args_ref_start = int(args_region.split(":")[1].split("-")[0])
            args_ref_end = int(args_region.split(":")[1].split("-")[1])

            samfile = pysam.AlignmentFile(args.alignment, mode='rb')
        else:
            args_ref_name = None
            args_ref_start = None
            args_ref_end = None
            samfile = pysam.AlignmentFile(args.alignment, mode='r')

        for sam_record in samfile.fetch(contig=args_ref_name, start=args_ref_start, stop=args_ref_end):
            metric_record = {}
            if args_ref_name is not None and args_ref_name != sam_record.reference_name:
                raise Exception("Error: sam record's reference name [" + sam_record.reference_name + "] and the name specified are different [" + args_ref_name + "]")
            read_id = sam_record.query_name
            if sam_record.is_supplementary or sam_record.is_unmapped or sam_record.is_secondary:
                continue
            if args.plot_reverse is True and sam_record.is_reverse is False:
                continue
            if args.plot_reverse is False and sam_record.is_reverse is True:
                continue
            if args.read_id != "" and read_id != args.read_id:
                continue
            if not sam_record.has_tag("ss"):
                raise Exception("Error: ss string is missing for the read_id {} in {}".format(read_id, args.alignment))
            ref_seq_len = 0
            data_is_rna = 0
            start_index = -1
            if sam_record.has_tag("si"):
                si_tag = sam_record.get_tag("si").split(',')
                start_index = int(si_tag[SI_START_RAW])
                end_index = int(si_tag[SI_END_RAW])
                ref_seq_len = int(si_tag[SI_END_KMER]) - int(si_tag[SI_START_KMER])
                reference_start = int(si_tag[SI_START_KMER])
                if int(si_tag[SI_START_KMER]) > int(si_tag[SI_END_KMER]):  # if RNA start_kmer>end_kmer in paf
                    data_is_rna = 1
                    if not args.rna:
                        print("Info: data is detected as RNA")
                        raise Exception("Error: data is not specified as RNA. Please provide the argument --rna ")
                    ref_seq_len = int(si_tag[SI_START_KMER]) - int(si_tag[SI_END_KMER])
                    reference_start = int(si_tag[SI_END_KMER]) + kmer_correction

            else:
                raise Exception("Error: sam record does not have a 'si' tag.")
            # print("ref_seq_len: " + str(ref_seq_len))

            if ref_seq_len < BASE_LIMIT:
                base_limit = ref_seq_len
            else:
                base_limit = BASE_LIMIT
            sam_record_reference_end = reference_start + ref_seq_len #1based closed
            if not args.loose_bound:
                if data_is_rna == 1:
                    if args_ref_start < reference_start + 1 - kmer_correction:
                        continue
                else:
                    if args_ref_start < reference_start + 1:
                        continue
                if args_ref_end > sam_record_reference_end:
                    continue

            ref_name = sam_record.reference_name
            ref_start = reference_start + 1
            ref_end = ref_start + base_limit - 1 #ref_end is 1based closed
            if args.region != "":
                if args_ref_start > ref_start:
                    ref_start = args_ref_start
                    if (ref_start + base_limit - 1) < sam_record_reference_end:
                        ref_end = ref_start + base_limit - 1
                    else:
                        ref_end = sam_record_reference_end
                if args_ref_end < ref_end:
                    ref_end = args_ref_end
            if ref_end < ref_start:
                print("Warning: a corner case has hit because  the kmer_length used is larger than 1. This alignment will be skipped")
                continue
            base_limit = ref_end - ref_start + 1
            # print("ref_start: {}".format(ref_start))
            # print("ref_end: {}".format(ref_end))
            # print("ref_seq_len: {}".format(ref_seq_len))
            # print("{}".format(sam_record.cigarstring))
            # print("base_limit: {}".format(base_limit))
            if sam_record.is_reverse and data_is_rna == 1:
                raise Exception("Error: A transcript is  always sequenced from 3` to 5`. Squigualiser only supports reads mapped to the transcriptome.")

            if data_is_rna == 1:
                # print("plot (RNA 5'->3') region: {}:{}-{}\tread_id: {}".format(ref_name, ref_end, ref_start, read_id))
                fasta_seq = fasta_reads.get_seq(name=ref_name, start=ref_start, end=ref_end).seq
                fasta_seq = fasta_seq.upper()
            else:
                if sam_record.is_reverse:
                    # print("plot (-) region: {}:{}-{}\tread_id: {}".format(ref_name, ref_start, ref_end, read_id))
                    fasta_seq = fasta_reads.get_seq(name=ref_name, start=ref_start, end=ref_end).seq
                    fasta_seq = fasta_seq.upper()
                    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
                    fasta_seq = "".join(nn[n] for n in fasta_seq)
                else:
                    # print("plot (+) region: {}:{}-{}\tread_id: {}".format(ref_name, ref_start, ref_end, read_id))
                    fasta_seq = fasta_reads.get_seq(name=ref_name, start=ref_start, end=ref_end).seq
                    fasta_seq = fasta_seq.upper()
            if not bool(re.match('^[ACGTUMRWSYKVHDBN]+$', fasta_seq)):
                raise Exception("Error: base characters other than A,C,G,T/U,M,R,W,S,Y,K,V,H,D,B,N were detected. Please check your sequence files")

            metric_record['read_id'] = read_id
            metric_record['region'] = "{}:{}-{}".format(ref_name, ref_start, ref_end)

            x = []
            x_real = []
            y = []

            read = s5.get_read(read_id, pA=args.no_pa, aux=["read_number", "start_mux"])
            if read is not None:
                # print("read_id:", read['read_id'])
                # print("len_raw_signal:", read['len_raw_signal'])
                # end_index = read['len_raw_signal']
                x = list(range(1, end_index - start_index + 1))
                x_real = list(range(start_index+1, end_index+1))             # 1based
                y = read['signal'][start_index:end_index]

            scaling_str = "no scaling"
            if args.sig_scale == "medmad" or args.sig_scale == "znorm" or args.sig_scale == "scaledpA":
                scaling_str = args.sig_scale
            elif not args.sig_scale == "":
                raise Exception("Error: given --sig_scale method: {} is not supported".format(args.sig_scale))

            scale_params = {}
            if args.sig_scale == "scaledpA":
                if not args.no_pa:
                    raise Exception("Error: given --sig_scale method: {} required the signal to be converted to pA levels. Please remove --no_pa argument".format(args.sig_scale))
                for tag in ["sc", "sh"]:
                    if sam_record.has_tag(tag):
                        scale_params[tag] = sam_record.get_tag(tag)
                    else:
                        raise Exception("Error: given --sig_scale method: {} requires {} tag in the alignment file".format(args.sig_scale, tag))

            y = plot_utils.scale_signal(y, args.sig_scale, scale_params)

            moves_string = sam_record.get_tag("ss")
            moves_string = re.sub('D', 'D,', moves_string)
            moves_string = re.sub('I', 'I,', moves_string).rstrip(',')
            moves = re.split(r',+', moves_string)

            if data_is_rna == 0:
                strand_dir = "(DNA +)"
                if sam_record.is_reverse:
                    strand_dir = "(DNA -)"
                    x_real.reverse()
                    y = np.flip(y)
                    moves.reverse()
            if data_is_rna == 1:
                strand_dir = "(RNA 3'->5')"
                fasta_seq = fasta_seq[::-1]
                if sam_record.is_reverse:
                    raise Exception("Error: data is rna and sam record is reverse mapped. This is not implemented yet. Please report")

            signal_tuple = (x, x_real, y)
            region_tuple = (ref_start, ref_end, reference_start, reference_start+ref_seq_len)

            sig_algn_dic = {}
            sig_algn_dic['start_kmer'] = 0
            sig_algn_dic['ref_start'] = ref_start
            sig_algn_dic['ref_end'] = ref_end
            sig_algn_dic['plot_sig_ref_flag'] = plot_sig_ref_flag
            sig_algn_dic['data_is_rna'] = data_is_rna
            sig_algn_dic['ss'] = moves
            sig_algn_dic['ref_name'] = ref_name

            # print(len(moves))
            # print(fasta_seq)
            signal_tuple, region_tuple, sig_algn_dic, fasta_seq = plot_utils.adjust_before_plotting(ref_seq_len, signal_tuple, region_tuple, sig_algn_dic, fasta_seq)

            # if args.auto_base_shift:
            #     draw_data["base_shift"] = plot_utils.calculate_base_shift(signal_tuple[2], fasta_seq, sig_algn_dic['ss'], args)
            #     print("automatically calculated base_shift: {}".format(draw_data["base_shift"]))

            if draw_data["base_shift"] < 0:
                abs_base_shift = abs(draw_data["base_shift"])
                x = signal_tuple[0]
                x_real = signal_tuple[1]
                y = signal_tuple[2]
                y_prefix = [np.nan] * abs_base_shift * draw_data["fixed_base_width"]
                y = np.concatenate((y_prefix, y), axis=0)
                x_real = np.concatenate(([1] * abs_base_shift * draw_data["fixed_base_width"], x_real), axis=0)
                x = list(range(1, len(x) + 1 + abs_base_shift * draw_data["fixed_base_width"]))
                signal_tuple = (x, x_real, y)
                moves_prefix = [str(draw_data["fixed_base_width"])] * abs_base_shift
                sig_algn_dic['ss'] = moves_prefix + sig_algn_dic['ss']

            # print(len(sig_algn_dic['ss']))
            draw_data['y_min'] = np.nanmin(y)
            draw_data['y_max'] = np.nanmax(y)
            get_metric(args=args, fout=fout, metric_record=metric_record, read_id=read_id, signal_tuple=signal_tuple, sig_algn_data=sig_algn_dic, fasta_sequence=fasta_seq, base_limit=base_limit, draw_data=draw_data)

            num_plots += 1
            if num_plots == args.plot_limit:
                break
    elif use_paf == 1 and plot_sig_ref_flag == 1:
        print("Info: Signal to reference method using PAF ...")
        fasta_reads = Fasta(args.file)
        tbxfile = pysam.TabixFile(args.alignment)
        if args.region != "":
            args_region = re.sub(',', '', args.region)
            # print(args_region)
            # pattern = re.compile("^[a-z]+[0-9]+\:[0-9]+\-[0-9]+")
            pattern = re.compile("^.*\:[0-9]+\-[0-9]+")
            if not pattern.match(args_region):
                raise Exception("Error: region provided is not in correct format")
            args_ref_name = args_region.split(":")[0]
            args_ref_start = int(args_region.split(":")[1].split("-")[0])
            args_ref_end = int(args_region.split(":")[1].split("-")[1])
        else:
            args_ref_name = None
            args_ref_start = None
            args_ref_end = None

        for paf_record in tbxfile.fetch(args_ref_name, args_ref_start, args_ref_end, parser=pysam.asTuple()):
            metric_record = {}
            if paf_record[READ_ID] == paf_record[SEQUENCE_ID]:
                raise Exception("Error: this paf file is a signal to read mapping.")
            if args_ref_name is not None and args_ref_name != paf_record[SEQUENCE_ID]:
                raise Exception("Error: sam record's reference name [" + paf_record[SEQUENCE_ID] + "] and the name specified are different [" + args_ref_name + "]")
            read_id = paf_record[READ_ID]
            # if read_id != "285802f0-8f4d-4f03-8d11-ef8a395576e4":
            #     continue
            if args.read_id != "" and read_id != args.read_id:
                continue
            if args.plot_reverse is True and paf_record[STRAND] == "+":
                continue
            if args.plot_reverse is False and paf_record[STRAND] == "-":
                continue
            moves_string = ""
            for i in range(12, len(paf_record)):
                tag = paf_record[i].split(':')[0]
                if tag == "ss":
                    moves_string = paf_record[i].split(':')[2]
            if moves_string == "":
                raise Exception("Error: ss string is missing for the read_id {} in {}".format(read_id, args.alignment))

            data_is_rna = 0
            start_index = int(paf_record[START_RAW])
            end_index = int(paf_record[END_RAW])
            ref_seq_len = int(paf_record[END_KMER]) - int(paf_record[START_KMER])
            reference_start = int(paf_record[START_KMER])
            if int(paf_record[START_KMER]) > int(paf_record[END_KMER]):  # if RNA start_kmer>end_kmer in paf
                data_is_rna = 1
                if not args.rna:
                    print("Info: data is detected as RNA")
                    raise Exception("Error: data is not specified as RNA. Please provide the argument --rna ")
                ref_seq_len = int(paf_record[START_KMER]) - int(paf_record[END_KMER])
                reference_start = int(paf_record[END_KMER]) + kmer_correction
            # print("ref_seq_len: " + str(ref_seq_len))
            if ref_seq_len < BASE_LIMIT:
                base_limit = ref_seq_len
            else:
                base_limit = BASE_LIMIT
            paf_record_reference_end = reference_start + ref_seq_len #1based closed
            if not args.loose_bound:
                if data_is_rna == 1:
                    if args_ref_start < reference_start + 1 - kmer_correction:
                        continue
                else:
                    if args_ref_start < reference_start + 1:
                        continue
                if args_ref_end > paf_record_reference_end:
                    continue
                    
            ref_name = paf_record[SEQUENCE_ID]
            ref_start = reference_start + 1
            ref_end = ref_start + base_limit - 1 #ref_end is 1based closed
            if args.region != "":
                if args_ref_start > ref_start:
                    ref_start = args_ref_start
                    if (ref_start + base_limit - 1) < paf_record_reference_end:
                        ref_end = ref_start + base_limit - 1
                    else:
                        ref_end = paf_record_reference_end
                if args_ref_end < ref_end:
                    ref_end = args_ref_end
            if ref_end < ref_start:
                print("Warning: a corner case has hit because the kmer_length used is larger than 1. This alignment will be skipped")
                continue
            base_limit = ref_end - ref_start + 1
            # print("ref_start: {}".format(ref_start))
            # print("ref_end: {}".format(ref_end))
            # print("ref_seq_len: {}".format(ref_seq_len))
            # print("base_limit: {}".format(base_limit))

            record_is_reverse = 0
            if paf_record[STRAND] == "-":
                record_is_reverse = 1
            if record_is_reverse and data_is_rna == 1:
                raise Exception("Error: A transcript is  always sequenced from 3` to 5`. Squigualiser only supports reads mapped to the transcriptome.")
            if data_is_rna == 1:
                # print("plot (RNA 5'->3') region: {}:{}-{}\tread_id: {}".format(ref_name, ref_end, ref_start, read_id))
                fasta_seq = fasta_reads.get_seq(name=ref_name, start=ref_start, end=ref_end).seq
                fasta_seq = fasta_seq.upper()
            else:
                if record_is_reverse:
                    # print("plot (-) region: {}:{}-{}\tread_id: {}".format(ref_name, ref_start, ref_end, read_id))
                    fasta_seq = fasta_reads.get_seq(name=ref_name, start=ref_start, end=ref_end).seq
                    fasta_seq = fasta_seq.upper()
                    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
                    fasta_seq = "".join(nn[n] for n in fasta_seq)
                else:
                    # print("plot (+) region: {}:{}-{}\tread_id: {}".format(ref_name, ref_start, ref_end, read_id))
                    fasta_seq = fasta_reads.get_seq(name=ref_name, start=ref_start, end=ref_end).seq
                    fasta_seq = fasta_seq.upper()
            if not bool(re.match('^[ACGTUMRWSYKVHDBN]+$', fasta_seq)):
                raise Exception("Error: base characters other than A,C,G,T/U,M,R,W,S,Y,K,V,H,D,B,N were detected. Please check your sequence files")

            metric_record['read_id'] = read_id
            metric_record['region'] = "{}:{}-{}".format(ref_name, ref_start, ref_end)

            x = []
            x_real = []
            y = []

            read = s5.get_read(read_id, pA=args.no_pa, aux=["read_number", "start_mux"])
            if read is not None:
                # print("read_id:", read['read_id'])
                # print("len_raw_signal:", read['len_raw_signal'])
                # end_index = read['len_raw_signal']
                x = list(range(1, end_index - start_index + 1))
                x_real = list(range(start_index + 1, end_index + 1))  # 1based
                y = read['signal'][start_index:end_index]

            scaling_str = "no scaling"
            if args.sig_scale == "medmad" or args.sig_scale == "znorm" or args.sig_scale == "scaledpA":
                scaling_str = args.sig_scale
            elif not args.sig_scale == "":
                raise Exception("Error: given --sig_scale method: {} is not supported".format(args.sig_scale))

            scale_params = {}
            if args.sig_scale == "scaledpA":
                if not args.no_pa:
                    raise Exception("Error: given --sig_scale method: {} required the signal to be converted to pA levels. Please remove --no_pa argument".format(args.sig_scale))
                for tag in ["sc", "sh"]:
                    for i in range(12, len(paf_record)):
                        if tag == paf_record[i].split(':')[0]:
                            scale_params[tag] = float(paf_record[i].split(':')[2])
                    if tag not in scale_params:
                        raise Exception("Error: required tag '{}' for given --sig_scale method: {} is not found in the alignment file".format(tag, args.sig_scale))
            y = plot_utils.scale_signal(y, args.sig_scale, scale_params)

            moves_string = re.sub('D', 'D,', moves_string)
            moves_string = re.sub('I', 'I,', moves_string).rstrip(',')
            moves = re.split(r',+', moves_string)

            if data_is_rna == 0:
                strand_dir = "(DNA +)"
                if record_is_reverse:
                    strand_dir = "(DNA -)"
                    x_real.reverse()
                    y = np.flip(y)
                    moves.reverse()
            if data_is_rna == 1:
                strand_dir = "(RNA 3'->5')"
                fasta_seq = fasta_seq[::-1]

            signal_tuple = (x, x_real, y)
            region_tuple = (ref_start, ref_end, reference_start, reference_start + ref_seq_len)

            sig_algn_dic = {}
            sig_algn_dic['start_kmer'] = 0
            sig_algn_dic['ref_start'] = ref_start
            sig_algn_dic['ref_end'] = ref_end
            sig_algn_dic['pa'] = args.no_pa
            sig_algn_dic['plot_sig_ref_flag'] = plot_sig_ref_flag
            sig_algn_dic['data_is_rna'] = data_is_rna
            sig_algn_dic['ss'] = moves
            sig_algn_dic['ref_name'] = ref_name

            signal_tuple, region_tuple, sig_algn_dic, fasta_seq = plot_utils.adjust_before_plotting(ref_seq_len, signal_tuple, region_tuple, sig_algn_dic, fasta_seq)

            # if args.auto_base_shift:
            #     draw_data["base_shift"] = plot_utils.calculate_base_shift(signal_tuple[2], fasta_seq, sig_algn_dic['ss'], args)
            #     print("automatically calculated base_shift: {}".format(draw_data["base_shift"]))

            if draw_data["base_shift"] < 0:
                abs_base_shift = abs(draw_data["base_shift"])
                x = signal_tuple[0]
                x_real = signal_tuple[1]
                y = signal_tuple[2]
                y_prefix = [np.nan] * abs_base_shift * draw_data["fixed_base_width"]
                y = np.concatenate((y_prefix, y), axis=0)
                x_real = np.concatenate(([1] * abs_base_shift * draw_data["fixed_base_width"], x_real), axis=0)
                x = list(range(1, len(x) + 1 + abs_base_shift * draw_data["fixed_base_width"]))
                signal_tuple = (x, x_real, y)
                moves_prefix = [str(draw_data["fixed_base_width"])] * abs_base_shift
                sig_algn_dic['ss'] = moves_prefix + sig_algn_dic['ss']

            # print(len(moves))
            # print(fasta_seq)
            # print(len(sig_algn_dic['ss']))
            draw_data['y_min'] = np.nanmin(y)
            draw_data['y_max'] = np.nanmax(y)
            get_metric(args=args, fout=fout, metric_record=metric_record, read_id=read_id, signal_tuple=signal_tuple, sig_algn_data=sig_algn_dic, fasta_sequence=fasta_seq, base_limit=base_limit, draw_data=draw_data)

            num_plots += 1
            if num_plots == args.plot_limit:
                break
    else:
        raise Exception("Error: You should not have ended up here. Please check your arguments")

    print("Number of records: {}".format(num_plots))
    if num_plots == 0:
        print("Squigualiser only plots reads that span across the specified region entirely. Reduce the region interval and double check with IGV : {}".format(num_plots))

    s5.close()
def argparser():
    # parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )

    parser.add_argument('-f', '--file', required=False, type=str, default="", help="fasta/fa/fastq/fq/fq.gz sequence file")
    parser.add_argument('-r', '--read_id', required=False, type=str, default="", help="plot the read with read_id")
    parser.add_argument('-s', '--slow5', required=False, type=str, default="", help="slow5 file")
    parser.add_argument('-a', '--alignment', required=False, type=str, default="", help="for read-signal alignment use PAF\nfor reference-signal alignment use SAM/BAM")
    parser.add_argument('--base_limit', required=False, type=int, help="maximum number of bases to plot")
    parser.add_argument('--region', required=False, type=str, default="", help="[start-end] 1-based closed interval region to plot. For SAM/BAM eg: chr1:6811428-6811467 or chr1:6,811,428-6,811,467. For PAF eg:100-200.")
    parser.add_argument('--tag_name', required=False, type=str, default="", help="a tag name to easily identify the plot")
    parser.add_argument('--plot_reverse', required=False, action='store_true', help="plot only the reverse mapped reads.")
    parser.add_argument('--rna', required=False, action='store_true', help="specify for RNA reads")
    parser.add_argument('--sig_ref', required=False, action='store_true', help="plot signal to reference mapping")
    parser.add_argument('--sig_scale', required=False, type=str, default="", help="plot the scaled signal. Supported scalings: [medmad, znorm, scaledpA]")
    # parser.add_argument('--reverse_signal', required=False, action='store_true', help="plot RNA reference/read from 5`-3` and reverse the signal")
    parser.add_argument('--no_pa', required=False, action='store_false', help="skip converting the signal to pA values")
    parser.add_argument('--loose_bound', required=False, action='store_true', help="also plot alignments not completely within the specified region")
    parser.add_argument('--base_shift', required=False, type=int, default=PLOT_BASE_SHIFT, help="the number of bases to shift to align fist signal move")
    parser.add_argument('--profile', required=False, default="", type=str, help="determine base_shift using preset values")
    parser.add_argument('--list_profile', action='store_true', help="list the available profiles")
    parser.add_argument('--plot_limit', required=False, type=int, default=PLOT_LIMIT, help="limit the number of plots generated")
    parser.add_argument('--sig_plot_limit', required=False, type=int, default=SIG_PLOT_LENGTH, help="maximum number of signal samples to plot")
    parser.add_argument('-o', '--output', required=False, type=str, default="", help="output file")
    parser.add_argument('--extend_0', required=False, action='store_true', help="print matches, deletions, insertions arrays")
    return parser

if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    try:
        run(args)
    except Exception as e:
        print(str(e))
        exit(1)


