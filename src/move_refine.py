from src import readpaf_local
import argparse
import numpy as np
import os

fout2 = open("delete", "w")


def get_closest_signal_indices(Q, R):
    signal_dict = {}

    n = len(Q)

    for i in range(n):
        left_bound = Q[i - 1] if i > 0 else float('-inf')  # Previous move or negative infinity
        right_bound = Q[i + 1] if i < n - 1 else float('inf')  # Next move or positive infinity

        # Collect indices of relevant R[j] values for the current Q[i]
        relevant_indices = []
        # for j, r in enumerate(R):
        #     if left_bound <= r <= right_bound:
        #         relevant_indices.append(j)
        
         # If no relevant indices are found, add the nearest left and right segmentations
        if not relevant_indices:
            # Find nearest left segmentation
            left_nearest = max((j for j, r in enumerate(R) if r < Q[i]), default=None)
            # Find nearest right segmentation
            right_nearest = min((j for j, r in enumerate(R) if r > Q[i]), default=None)

            # Find overlap segmentation
            center = min((j for j, r in enumerate(R) if r == Q[i]), default=None)

            if left_nearest is not None:
                relevant_indices.append(left_nearest)
            if right_nearest is not None:
                relevant_indices.append(right_nearest)
            if center is not None:
                relevant_indices.append(center)

        # Add to the dictionary
        signal_dict[i] = relevant_indices
    print(signal_dict)
    return signal_dict

def print_dp_matrix(dp_matrix,i):
    """
    Prints a 2D matrix in a readable format.

    Args:
        dp_matrix (list of lists): The 2D DP matrix to print.
    """
    for row in dp_matrix[i:i+1]:
        print(" ".join(f"{val:6}" for val in row))

def refine_moves(Q, R):

    q_i_dict = get_closest_signal_indices(Q=Q, R=R)
    # print(q_i_dict)
    N, M = len(Q), len(R)

    print("N:{} M:{}".format(N,M))

    # Initialize DP matrix with infinity
    D = np.full((N, M), float('inf'))
    D[0][0] = 0  # start of the signal
    print("101")

    # Fill DP matrix
    for i in range(1, N):
        # for j in range(1, M+1):
        j_indices = q_i_dict[i]
        # print(i)
        # print("{}:{}".format(i,j_indices))
        for j in j_indices:
            if j == 0:
                continue
            cost = abs(Q[i] - R[j])  # Cost function: absolute difference
            D[i][j] = cost + min(
                D[i-1][j],    # Vertical step
                D[i][j-1],    # Horizontal step
                D[i-1][j-1]   # Diagonal step
            )
            fout2.write("D[{}][{}] = {}\n".format(i,j,D[i][j]))
    print("102")
    # print_dp_matrix(D,N+1)
    # print(min(D[N-1]))
    # D[N][M] = min(D[N-1])

    # Backtracking to find the optimal alignment
    i, j = N-1, M-1
    alignment = []
    alignment.append((i,j,Q[i], R[j]))  # Append the aligned pair
    prev_j = M
    while i > 0 and j > 0:
        if i == j:
            i -= 1
            j -= 1
            alignment.append((i,j,Q[i], R[j]))  # Append the aligned pair
            prev_j = j
        elif D[i][j] == D[i-1][j-1] + abs(Q[i] - R[j]):
            i -= 1
            j -= 1
            alignment.append((i,j,Q[i], R[j]))  # Append the aligned pair
            prev_j = j
        elif D[i][j] == D[i-1][j] + abs(Q[i] - R[j]) and prev_j != j:
            i -= 1
            alignment.append((i,j,Q[i], R[j]))  # Append the aligned pair
            prev_j = j
        elif D[i][j] == D[i-1][j] + abs(Q[i] - R[j]) and prev_j == j:
            print("here???")
            i -= 1
            j -= 1
            alignment.append((i,j,Q[i], R[j]))  # Append the aligned pair
            prev_j = j
        else:
            j -= 1
    alignment.reverse()  # Reverse the path to get the correct order

    print(len(alignment))
    print(alignment)

    return [j for _,_,_, j in alignment]

def get_refined_moves_from_alignment(alignment, Q, first_occurrence=True):
    # Dictionary to store the first or last R[j] match for each Q[i]
    matches = {}

    # Iterate over the alignment and update the dictionary based on the flag
    for q, r in alignment:
        if first_occurrence:
            # Store the first occurrence of R[j] for each Q[i]
            if q not in matches:
                matches[q] = r
        else:
            # Store the last occurrence of R[j] for each Q[i]
            matches[q] = r  # This will overwrite with the last match

    # Construct the refined moves array using the selected matches
    refined_moves = [matches[q] for q in Q]

    return refined_moves


def calculate_similarity(Q, refined_moves):
    if len(Q) != len(refined_moves):
        raise ValueError("The lengths of Q and refined_moves must be the same!")

    # Calculate the mean absolute difference
    differences = [abs(q - r) for q, r in zip(Q, refined_moves)]
    mad = sum(differences) / len(differences)

    return mad

def make_ss_string(moves, start_index):
    ss_string = "ss:Z:"
    prev_index = start_index
    for move in moves[1:]:
        ss_string += "{},".format(move-prev_index)
        prev_index = move
    return ss_string

def write_to_file(fout, paf_record, ss_string, end_index):
    fout.write("{}\t".format(paf_record.query_name))
    fout.write("{}\t".format(end_index-paf_record.query_start))
    fout.write("{}\t".format(paf_record.query_start))
    fout.write("{}\t".format(end_index))
    fout.write("{}\t".format(paf_record.strand))
    fout.write("{}\t".format(paf_record.target_name))
    fout.write("{}\t".format(paf_record.target_length))
    fout.write("{}\t".format(paf_record.target_start))
    fout.write("{}\t".format(paf_record.target_end))
    fout.write("{}\t".format(paf_record.residue_matches))
    fout.write("{}\t".format(paf_record.alignment_block_length))
    fout.write("{}\t".format(paf_record.mapping_quality))
    fout.write("{}\t".format(ss_string))

def refine_moves_greedy(Q,R):
    alignment = []
    j = 0  # Pointer for R

    for i, q in enumerate(Q):
        # Advance j to find the two R[j] values surrounding Q[i]
        while j < len(R) - 1 and R[j+1] <= q:
            j += 1
        
        # Determine the closest R[j] or R[j+1] to Q[i]
        if j < len(R) - 1:
            left_dist = abs(q - R[j])
            right_dist = abs(q - R[j + 1])
            if right_dist < left_dist:
                j += 1
        
        # Align Q[i] to the closest R[j]
        alignment.append((i, j, Q[i], R[j]))
        j += 1  # Move to the next R[j] to ensure no reuse

    print(len(alignment))
    print(alignment)

    return [j for _,_,_, j in alignment]


def run(args):
    fout = open(args.output, "w")
    reads = []
    sigann_dict = {}
    if args.sigann:
        print(f'Signal point annotation file: {args.sigann}')
        with open(args.sigann, 'r') as file:
            for line in file:
                line = line.strip()
                if line:  # Skip empty lines
                    read_id, array_str = line.split(' ', 1)  # Split only on the first space
                    array = [int(value) for value in array_str.split(',') if value]  # Parse the array
                    sigann_dict[read_id] = array
    
    moves_dict = {}
    with open(args.alignment, "r") as handle:
        for paf_record in readpaf_local.parse_paf(handle):
            read_id = paf_record.query_name
            start_index = paf_record.query_start
            moves_string = paf_record.tags['ss'][2]
            moves = [int(value) for value in moves_string.split(',') if value]
            moves_array = []
            prev_bound = start_index
            for move in moves:
                prev_bound = prev_bound +  move
                moves_array.append(prev_bound) 
            moves_dict[read_id] = moves_array
            reads.append(read_id)

            print(read_id)
            moves = moves_dict[read_id]
            sigann = sigann_dict[read_id]
            moves = [0] + moves
            sigann = [0] + sigann
            if moves[-1] != sigann[-1]:
                sigann[-1] = moves[-1]
            # print(moves)
            # print(sigann)
            print("len sigann: {} len moves: {}".format(len(sigann), len(moves)))

            refined_moves = refine_moves_greedy(Q=moves, R=sigann)

            # refined_moves = refine_moves(Q=moves, R=sigann)
        
            print(calculate_similarity(Q=moves, refined_moves=refined_moves))

            ss_string = make_ss_string(refined_moves, start_index)

            write_to_file(fout, paf_record, ss_string=ss_string, end_index=refined_moves[-1])

    fout.close()
    fout2.close()


def argparser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )

    parser.add_argument('-r', '--read_id', required=False, type=str, default="", help="plot the read with read_id")
    parser.add_argument('-l', '--read_list', required=False, type=str, default="", help="a file with read_ids to plot")
    parser.add_argument('-a', '--alignment', required=True, type=str, default="", help="read-signal alignment in PAF")
    parser.add_argument('--sigann', required=False, help="file with signal point annotations (0-based)")
    parser.add_argument('-o', '--output', required=True, type=str, default="", help="output file")
    return parser


if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    try:
        run(args)
    except Exception as e:
        print(str(e))
        exit(1)

    