# Description of each column of the output of metric sub-tool

To run `metric` refer [here](commands.md/#metric).

| Column name                 | Description                                                    |
|-----------------------------|----------------------------------------------------------------|
| read_id                     | read_id of the raw signal                                      |
| ref_region                  | statistics were calculated within this region                  |
| signal_region               | aligned signal region                                          |
| total_matches               | total number of base to base matches                           |
| total_deletion_occurrences  | total number of deletion occurrences                           |
| total_insertion_occurrences | total number of insertion occurrences                          |
| total_length_deletions      | total length of deletions                                      |
| total_length_insertions     | total number of samples counted as insertions to the reference |
| min_match                   | minimum number of samples in a match                           |
| max_match                   | maximum number of samples in a match                           |
| mode_match                  | mode of the number of samples in a match                       |
| median_match                | median of the number of samples in a match                     |
| mean_match                  | mean of the number of samples in a match                       |
| stdev_match                 | standard deviation of the number of samples in a match         |
| min_deletion                | minimum length (in bases) of deletions                         |
| max_deletion                | maximum length (in bases) of deletions                         |
| mode_deletion               | mode of the length (in bases) of deletions                     |
| median_deletion             | median of the length (in bases) of deletions                   |
| mean_deletion               | mean of the length (in bases) of deletions                     |
| stdev_deletion              | standard deviation of the length (in bases) of deletions       |
| min_insertion               | minimum number of samples in an insertion                      |
| max_insertion               | maximum number of samples in an insertion                      |
| mode_insertion              | mode of the number of samples in an insertion                  |
| median_insertion            | median of the number of samples in an insertion                |
| mean_insertion              | mean of the number of samples in an insertion                  |
| stdev_insertion             | standard deviation of the number of samples in an insertion    |
| matches                     | matches array                                                  |
| deletions                   | deletions  array                                               |
| insertions                  | insertions  array                                              |
