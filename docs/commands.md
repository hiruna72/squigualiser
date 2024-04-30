# Commands and Options

* `reform`:<br/>
         Convert basecaller's move table to [ss string](https://hasindu2008.github.io/f5c/docs/output#ss-tag) format.
* `realign`:<br/>
         Realign signal to reference using [cigar string](https://samtools.github.io/hts-specs/SAMv1.pdf) and the move table.
* `plot`:<br/>
         Plot read/reference - signal alignments.
* `plot_pileup`:<br/>
         Plot a reference - signal alignment pileup.
* `plot_tracks`:<br/>
         Plot multiple reference - signal alignment pileup tracks.
* `calculate_offsets`:<br/>
		 A utility program to calculate the most significant base index given a kmer model or a read - signal alignment.
* `metric`:<br/>
		 A utility program to calculate some statistics of the signal alignment.

### reform

```
squigualiser reform [OPTIONS] -c --bam basecall_moves.sam -o reform.paf
squigualiser reform [OPTIONS] --bam basecall_moves.sam -o reform.tsv
```

Convert basecaller's move table to [ss string](https://hasindu2008.github.io/f5c/docs/output#ss-tag) format.
For more information refer [reform](reform.md).

*  `--bam FILE`:<br/>
   The basecaller's move table output in sam or bam format. Guppy outputs multiple sam/bam files. These files have to be merged to create a single sam/bam file.
*  `-o, --output FILE `:<br/>
   Specifies name/location of the output file. A valid relative or absolute path can be provided. Data will be overwritten. The extension should be either `.tsv` or `.paf`.
*  `-k, --kmer_length INT`:<br/>
   (optional) kmer length to consider for the calculation.
*  `-m, --sig_move_offset INT`:<br/>
   (optional) signal move offset to consider for the calculation.
*  `-c`:<br/>
   (optional) write move table in paf format. [default is tsv]
*  `--rna`:<br/>
   (optional) Specify if RNA reads are used. By default, only DNA is accepted [default value: false].
   The tool will not error out if the data is RNA but if the user failed to specify it.
   It will still generate an incorrect output.
*  `--profile STR`:<br/>
   (optional) This is used to determine -k and -m values using preset values. The available profiles can be listed using `--list_profile`  
*  `--list_profile`:<br/>
    (optional) Print the available profiles and exit.

### realign

```
squigualiser realign [OPTIONS] -c --bam minimap2.sam --paf reform_output.paf -o realign_output.paf
squigualiser realign [OPTIONS] -c --bam minimap2.bam --paf reform_output.paf -o realign_output.paf
squigualiser realign [OPTIONS] --bam minimap2.bam --paf reform_output.paf -o realign_output.sam
squigualiser realign [OPTIONS] --bam minimap2.bam --paf reform_output.paf -o realign_output.bam
```

Realign signal to reference using [cigar string](https://samtools.github.io/hts-specs/SAMv1.pdf) and the move table.
For more information refer [realign](realign.md).

* `--paf FILE`:<br/>
   The read-signal alignment file generated using `reform` with `-c` option. 
* `--bam FILE`:<br/>
   The read-reference alignment sam/bam file containing cigar strings (minimap2 sam output). 
*  `-o, --output FILE `:<br/>
   Specifies name/location of the output file. A valid relative or absolute path can be provided. Data will be overwritten. The extension should be `.sam`, `.bam`, or `.paf`.
*  `-c`:<br/>
   (optional) write move table in paf format. [default is sam or bam depending on the extension provided]
*  `--rna`:<br/>
   (optional) Specify if RNA reads are used. By default, only DNA is accepted [default value: false].

### plot

```
squigualiser plot [OPTIONS] -f reads.fasta -s reads.slow5 -a reform.paf -o output_dir
squigualiser plot [OPTIONS] -f reads.fastq -s reads.slow5 -a reform.paf -o output_dir
squigualiser plot [OPTIONS] -f reads.fasta -s reads.slow5 -a realign.sam -o output_dir
squigualiser plot [OPTIONS] -f reads.fasta -s reads.slow5 -a realign.paf.gz -o output_dir
```

Plot read/reference - signal alignments.

* `-f, --file FILE`:<br/>
   The sequence file in `fasta/fa/fastq/fq/fq.gz` format.
* `-r, --read_id STR`
   (optional) Plot only the read with the read id specified.
* `-s, --slow5 FILE`
	Path to slow5 file containing raw signals.
* `-a, --alignment FILE`
	For read-signal alignment plots provide the path to `.paf` file generated using `reform`.
	For reference-signal alignment plots provide the `.sam/.bam` or `.paf.gz` file (generated using `realgin` or `f5c eventalign`).
* `--region STR`
	[start-end] 1-based closed interval region to plot.
	For read-signal alignment eg:`100-200`.
	For reference-signal alignment eg: `chr1:6811428-6811467` or `chr1:6,811,428-6,811,467`.
*  `-o, --output_dir DIR `:<br/>
   Specifies name/location of the output directory. A valid relative or absolute path can be provided. Data will be overwritten but the directory will not be recreated.
* `--tag_name STR`
	(optional) A tag name to easily identify the plot.
* `--plot_reverse`
	(optional) Plot only the reverse mapped reads [default value: false].
* `--rna`
	(optional) Specify if RNA reads are used. By default, only DNA is accepted [default value: false].
* `--base_limit INT`
	(optional) Maximum number of bases to plot.
* `--sig_ref`
	(optional) Plot signal to reference mapping. Can be mostly discarded. Act as a flag to avoid ploting reference-signal alignment using a `.paf` alignment file [default value: false].
* `--fixed_width`
	(optional) Plot with fixed base width. By default, the base widht is stretched to have equal distance between the sample points. By activating this flag, sample points of a base will be squeezed to a fixed width. [default value: false]
* `--sig_scale`
	(optional) Plot the scaled signal.
	By default, the signal is not scaled but converted to pA values.
	Supported scalings are: [medmad, znorm, scaledpA].
	`medmad` is median absolute deviation scaling.  
	`znorm` is zscore normalization scaling.
	`scaledpA` uses `sc` and `sh` tags to scale the raw signal to the pore model.
	The implementation of each method can be found at `src/plot_utils.py/scale_signal()`
* `--no_pa`
	(optional) Do not convert the signal to pA levels. By default, the raw signal is converted to pA levels [default value: false].
* `--loose_bound`
	(optional) Also plot alignments not completely within the specified region but at least part is [default value: false].
* `--point_size INT`
	(optional) Radius of the signal point drawn in the plot [default value: 0.5].
* `--base_width INT`
	(optional) The base width when plotting with fixed base width [default value: 10]. 
* `--base_shift INT`
	(optional) The number of bases to shift to align fist signal move [default value: 0]. More information on this can be found at [here](pore_model.md)
* `--auto`
	(optional) Calculate and automatically set the base shift. [default value: false]
*  `--profile STR`:<br/>
   (optional) This is used to determine base shift using preset values. The available profiles can be listed using `--list_profile`. The precedence is in the order of `--base_shift < --auto < --profile`.
*  `--list_profile`:<br/>
    (optional) Print the available profiles and exit.
* `--plot_limit INT`
	(optional) the number of plots to be generated [default value: 1000]. 
* `--sig_plot_limit INT`
	(optional) The maximum number of signal samples to draw on a plot [default value: 20000].
* `--bed FILE`
	(optional) The bed file with annotations.
* `--no_colours`
	(optional) hide base colours [default value: false].
* `--no_samples`
	(optional) hide sample points [default value: false].
* `--save_svg`
	(optional) save as svg. tweak --region and --xrange to capture the necessary part of the plot [default value: false].
* `--xrange`
	(optional) initial x range [default value: 350].

### plot_pileup

```
squigualiser plot_pileup [OPTIONS] -f genome.fasta -s reads.slow5 -a realign.sam -o output_dir --region chr1:20000-20400
squigualiser plot_pileup [OPTIONS] -f genome.fasta -s reads.slow5 -a realign.paf.gz -o output_dir --region chr1:20000-20400
```

Plot reference - signal alignment pileups.

* `-f, --file FILE`:<br/>
	The sequence file in `fasta/fa/fastq/fq/fq.gz` format.
* `-s, --slow5 FILE`
	Path to slow5 file containing raw signals.
* `-a, --alignment FILE`
	Provide the `.sam/.bam` or `.paf.gz` file (generated using `realgin` or `f5c eventalign`).
* `--region STR`
	[start-end] 1-based closed interval region to plot.
	For read-signal alignment eg:`100-200`.
	For reference-signal alignment eg: `chr1:6811428-6811467` or `chr1:6,811,428-6,811,467`.
*  `-o, --output_dir DIR `:<br/>
	Specifies name/location of the output directory. A valid relative or absolute path can be provided. Data will be overwritten but the directory will not be recreated.
* `-r, --read_id STR`
	(optional) Plot only the read with the read id specified.
* `-l, --read_list STR`
	(optional) Path to a file containing a list of read_ids to plot only the reads listed.
* `--tag_name STR`
	(optional) A tag name to easily identify the plot
* `--plot_reverse`
	(optional) Plot only the reverse mapped reads.
* `--rna`
	(optional) Specify if RNA reads are used. By default, only DNA is accepted.
* `--base_limit INT`
	(optional) Maximum number of bases to plot.
* `--sig_scale`
	(optional) Plot the scaled signal.
	By default, the signal is not scaled but converted to pA values.
	Supported scalings are: [medmad, znorm, scaledpA].
	`medmad` is median absolute deviation scaling.  
	`znorm` is zscore normalization scaling.
	`scaledpA` uses `sc` and `sh` tags to scale the raw signal to the pore model.
	The implementation of each method can be found at `src/plot_utils.py/scale_signal()`
* `--no_pa`
	(optional) Do not convert the signal to pA levels. By default, the raw signal is converted to pA levels [default value: false].
* `--point_size INT`
	(optional) Radius of the signal point drawn in the plot [default value: 0.5].
* `--base_width INT`
	(optional) The base width when plotting with fixed base width [default value: 10]. 
* `--base_shift INT`
	(optional) The number of bases to shift to align fist signal move [default value: 0]. More information on this can be found at [here](pore_model.md)
* `--auto`
	(optional) Calculate and automatically set the base shift. [default value: false]
*  `--profile STR`:<br/>
   (optional) This is used to determine base shift using preset values. The available profiles can be listed using `--list_profile`. The precedence is in the order of `--base_shift < --auto < --profile`.
*  `--list_profile`:<br/>
    (optional) Print the available profiles and exit.
* `--plot_limit INT`
	(optional) the number of plots to be generated [default value: 1000]. 
* `--sig_plot_limit INT`
	(optional) The maximum number of signal samples to draw on a plot [default value: 20000].
* `--bed FILE`
	(optional) The bed file with annotations.
* `--plot_num_samples`
	(optional) Annotate the number of samples for each move. By default, it is not annotated. Annotation can make the plots less responsive.
* `--overlap_bottom`
	(optional) Plot the overlap at the bottom. By default, it is plotted at the top [default value: false].
* `--no_overlap`
	(optional) Do not plot the overlap. By default, it is plotted at the top [default value: false].
* `--overlap_only`
	(optional) Plot only the overlap [default value: false].
* `--cprofile`
	(optional) Create a log file using python cprofile profiler [default value: false].
* `--return_plot`
	(optional) Return the plot object without saving to output. This is used in `plot_tracks` tool [default value: false]. Cannot be used in conjunction with `--output_dir`.
* `--no_colours`
	(optional) hide base colours [default value: false].
* `--save_svg`
	(optional) save as svg. tweak --region and --xrange to capture the necessary part of the plot [default value: false].
* `--xrange`
	(optional) initial x range [default value: 350].

### plot_tracks

```
squigualiser plot_tracks [OPTIONS] -f commands.txt -s reads.slow5 -a reform.paf -o output_dir
```

Plot multiple reference - signal alignment pileup tracks.

* `-f, --file FILE`:<br/>
   The file that contains the `pileup_plot` commands.
*  `-o, --output_dir DIR `:<br/>
   Specifies name/location of the output directory. A valid relative or absolute path can be provided. Data will be overwritten but the directory will not be recreated.
* `--tag_name STR`
	(optional) A tag name to easily identify the plot
* `--shared_x`
	(optional) Share x-axis so that all the plots move together [default value: false].
* `--auto_height`
	(optional) Adjust track height automatically using the number of plots available in each track [default value: false].

### calculate_offsets

```
squigualiser calculate_offsets [OPTIONS] -f reads.fasta -s reads.slow5 -a reform.paf -o output_dir
```

A utility program to calculate the most significant base index given a kmer model or a read - signal alignment.

*  `-k, --kmer_length INT`:<br/>
   (optional) kmer length to consider for the calculation.
* `--paf FILE`:<br/>
   The read-signal alignment file generated using `reform` with `-c` option. 
* `-f, --file FILE`:<br/>
   The sequence file in `fasta/fa/fastq/fq/fq.gz` format.
* `-s, --slow5 FILE`
   Path to slow5 file containing raw signals.
* `--model`
   (optional) Path to the kmer model file if `--use_model` is true.
* `--use_model`
   (optional) Calculate offset using the model file [default value: false]. By default, calculation is done for the reads using move table annotation.
* `--read_limit`
   (optional) Limit the number of reads considered for the calculation [default value: 1000].
* `-r, --read_id STR`
   (optional)Plot only the read with the read id specified.
*  `-o, --output FILE `:<br/>
   Specifies name/location of the output file. A valid relative or absolute path can be provided. Data will be overwritten. The extension should be `.pdf`.
   Output file can be used in `--use_model` mode or when a read id is specified (`--read_id`).
* `--tag_name STR`
	(optional) A tag name to easily identify the plot
*  `--rna`:<br/>
   (optional) Specify if RNA reads are used. By default, only DNA is accepted [default value: false].
   The tool will not error out if the data is RNA but if the user failed to specify it.
   It will still generate an incorrect output.

### metric

```
squigualiser metric [OPTIONS] -f reads.fasta -s reads.slow5 -a reform.paf -o output.tsv
squigualiser metric [OPTIONS] -f reads.fastq -s reads.slow5 -a reform.paf -o output.tsv
squigualiser metric [OPTIONS] -f reads.fasta -s reads.slow5 -a realign.sam -o output.tsv
squigualiser metric [OPTIONS] -f reads.fasta -s reads.slow5 -a realign.paf.gz -o output.tsv
```

Parse the ss string to calculate basic statistics of the read/reference - signal alignments.
The arguments are almost the same as for `plot` and `plot_pileup` tools.
Instead of generating figures `metric` will generate statistics after parsing the ss string.
[These statistics](metric.md) are written to the output file.

* `-f, --file FILE`:<br/>
   The sequence file in `fasta/fa/fastq/fq/fq.gz` format.
* `-r, --read_id STR`
   (optional) Plot only the read with the read id specified.
* `-s, --slow5 FILE`
	Path to slow5 file containing raw signals.
* `-a, --alignment FILE`
	For read-signal alignment plots provide the path to `.paf` file generated using `reform`.
	For reference-signal alignment plots provide the `.sam/.bam` or `.paf.gz` file (generated using `realgin` or `f5c eventalign`).
* `--region STR`
	[start-end] 1-based closed interval region to plot.
	For read-signal alignment eg:`100-200`.
	For reference-signal alignment eg: `chr1:6811428-6811467` or `chr1:6,811,428-6,811,467`.
*  `-o, --output FILE `:<br/>
   Specifies name/location of the output file. A valid relative or absolute path can be provided. Data will be overwritten. [These statistics](metric.md) are written to the output file.
* `--tag_name STR`
	(optional) A tag name to easily identify the plot.
* `--plot_reverse`
	(optional) Plot only the reverse mapped reads [default value: false].
* `--rna`
	(optional) Specify if RNA reads are used. By default, only DNA is accepted [default value: false].
* `--sig_ref`
	(optional) Plot signal to reference mapping. Can be mostly discarded. Act as a flag to avoid ploting reference-signal alignment using a `.paf` alignment file [default value: false].
* `--sig_scale`
	(optional) Plot the scaled signal.
	By default, the signal is not scaled but converted to pA values.
	Supported scalings are: [medmad, znorm, scaledpA].
	`medmad` is median absolute deviation scaling.  
	`znorm` is zscore normalization scaling.
	`scaledpA` uses `sc` and `sh` tags to scale the raw signal to the pore model.
	The implementation of each method can be found at `src/plot_utils.py/scale_signal()`
* `--no_pa`
	(optional) Do not convert the signal to pA levels. By default, the raw signal is converted to pA levels [default value: false].
* `--loose_bound`
	(optional) Also plot alignments not completely within the specified region but at least part is [default value: false].
* `--base_shift INT`
	(optional) The number of bases to shift to align fist signal move [default value: 0]. More information on this can be found at [here](pore_model.md)
*  `--profile STR`:<br/>
   (optional) This is used to determine base shift using preset values. The available profiles can be listed using `--list_profile`  
*  `--list_profile`:<br/>
    (optional) Print the available profiles and exit.
* `--plot_limit INT`
	(optional) the number of plots to be generated [default value: 1000]. 
* `--sig_plot_limit INT`
	(optional) The maximum number of signal samples to draw on a plot [default value: 20000].
