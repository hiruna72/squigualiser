"""
Signal to seQuence alignment Plot - plot
Hiruna Samarakoon - Garvan Medical Institute
hiruna@unsw.edu.au
"""
import numpy as np
import argparse
import re
import os
from src import plot_utils
import statistics
import csv
from sklearn.metrics.pairwise import cosine_similarity
import ast
import seaborn as sns
import matplotlib.pyplot as plt

READ_ID, REF_REGION, SIGNAL_REGION, TOTAL_MATCHES, TOTAL_DELETION_OCCURRENCES, TOTAL_INSERTION_OCCURRENCES, TOTAL_LENGTH_DELETIONS, TOTAL_LENGTH_INSERTIONS, MIN_MATCH, MAX_MATCH, MODE_MATCH, MEDIAN_MATCH, MEAN_MATCH, STDEV_MATCH, MIN_DELETION, MAX_DELETION, MODE_DELETION, MEDIAN_DELETION, MEAN_DELETION, STDEV_DELETION, MIN_INSERTION, MAX_INSERTION, MODE_INSERTION, MEDIAN_INSERTION, MEAN_INSERTION, STDEV_INSERTION,  MATCHES, DELETIONS, INSERTIONS = range(29)

def get_column_from_tsv(file_path, column_index):
    columns = []
    with open(file_path, 'r', newline='', encoding='utf-8') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        # Skip the header (first row)
        next(reader, None)
        for row in reader:
            if row:  # Check if the row is not empty
                columns.append(row[column_index])

    return columns

def calculate_cosine_similarity(list1, list2):
    # Convert lists to NumPy arrays
    array1 = np.array(list1).reshape(1, -1)
    array2 = np.array(list2).reshape(1, -1)

    # Calculate cosine similarity
    similarity_matrix = cosine_similarity(array1, array2)

    # Extract the similarity score
    similarity_score = round(similarity_matrix[0, 0], 4)
    return similarity_score

def plot_density(list1, list2):
    # # Create density plots
    # sns.kdeplot(list1, label='List 1')
    # sns.kdeplot(list2, label='List 2')
    # plt.ylabel('Density')
    # plt.title('Density Plot of List 1 and List 2')
    #
    # Plot histograms
    # plt.hist(list1, bins=20, alpha=0.5, label='List 1', color='blue')
    # plt.hist(list2, bins=20, alpha=0.5, label='List 2', color='orange')
    # plt.ylabel('Frequency')
    # plt.title('Histogram of List 1 and List 2')
    #
    # # Set labels and title
    # plt.xlabel('Values')
    #
    # # Display legend
    # plt.legend()

    # Create a box plot
    plt.boxplot([list1, list2], labels=['List 1', 'List 2'])

    # Set labels and title
    plt.xlabel('Lists')
    plt.ylabel('Values')
    plt.title('Box Plot of List 1 and List 2')

    # Show the plot
    plt.show()

def run(args):
    file_0 = args.f0
    file_1 = args.f1
    readid_0 = get_column_from_tsv(file_0, READ_ID)
    readid_1 = get_column_from_tsv(file_1, READ_ID)
    if readid_0 != readid_1:
        print("read id order is different")

    matches_0 = get_column_from_tsv(file_0, MATCHES)
    matches_1 = get_column_from_tsv(file_1, MATCHES)

    # for i in range(0, len(matches_0)):
    #     # actual_list = ast.literal_eval(string_representation)
    #     cosine_similarity_score = calculate_cosine_similarity(ast.literal_eval(matches_0[i]), ast.literal_eval(matches_1[i]))
    #     print("{}\t{}".format(readid_0[i], cosine_similarity_score))

    matches_0_all = []
    matches_1_all = []
    for i in range(0, len(matches_0)):
        matches_0_all += ast.literal_eval(matches_0[i])
        matches_1_all += ast.literal_eval(matches_1[i])
    plot_density(matches_0_all, matches_1_all)


def argparser():
    # parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )

    parser.add_argument('--f0', required=True, type=str, default="", help="tsv file created using metric")
    parser.add_argument('--f1', required=True, type=str, default="", help="tsv file created using metric")

    return parser

if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    try:
        run(args)
    except Exception as e:
        print(str(e))
        exit(1)


