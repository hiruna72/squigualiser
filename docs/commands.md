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
		 Utility program to calculate the most significant base index given a kmer model or a read - signal alignment.

### reform

```
squigualiser reform [OPTIONS] -c --bam basecall_moves.sam -o reform.paf
squigualiser reform [OPTIONS] --bam basecall_moves.sam -o reform.tsv
```

Convert basecaller's move table to [ss string](https://hasindu2008.github.io/f5c/docs/output#ss-tag) format.
For more information refer [reform](reform.md).

*  `--bam FILE`:<br/>
   The basecaller's move table output in sam or bam format. Guppy outputs multiple sam/bam files. These files have to be merged to create a single sam/bam file.
*  `-o, ----output FILE `:<br/>
   Specifies name/location of the output file. A valid relative or absolute path can be provided. Data will be overwritten.
*  `-k, --kmer_length INT`:<br/>
   kmer length to consider for the calculation.
*  `-m, --sig_move_offset INT`:<br/>
   kmer length to consider for the calculation.
*  `-c, --compress compression_type`:<br/>
   Specifies the compression method used for BLOW5 output. `compression_type` can be `none` for uncompressed binary; `zlib` for zlib-based (also known as gzip or DEFLATE) compression; or `zstd` for Z-standard-based compression [default value: zlib]. This option is only valid for BLOW5. `zstd` will only function if slow5tools has been built with zstd support which is turned off by default.
*  `-s, --sig-compress compression_type`:<br/>
   Specifies the raw signal compression method used for BLOW5 output. `compression_type` can be `none` for uncompressed raw signal or `svb-zd` to compress the raw signal using StreamVByte zig-zag delta [default value: svb-zd]. This option is introduced from slow5tools v0.3.0 onwards. Note that record compression (-c option above) is still applied on top of the compressed signal. Signal compression with svb-zd and record compression with zstd is similar to ONT's vbz.  zstd+svb-zd offers slightly smaller file size and slightly better performance compared to the default zlib+svb-zd, however, will be less portable.
*  `-p, --iop INT`:<br/>
    Specifies the number of I/O processes to use during conversion [default value: 8]. Increasing the number of I/O processes makes f2s significantly faster, especially on HPC with RAID systems (multiple disks) where a large value number of processes can be used (e.g., `-p 64`).
*  `--lossless STR`:<br/>
    Retain information in auxiliary fields during FAST5 to SLOW5 conversion. STR can be either `true` or `false`. [default value: true]. This information is generally not required for downstream analysis and can be optionally discarded to reduce filesize. *IMPORTANT: Generated files are only to be used for intermediate analysis and NOT for archiving. You will not be able to convert lossy files back to FAST5*.
* `-a, --allow`:<br/>
   By default f2s will not accept an individual multi-fast5 file or an individual single-fast5 directory containing multiple unique run IDs. When `-a` is specified f2s will allow multiple unique run IDs in an individual multi-fast5 file or single-fast5 directory. In this case, the header of all SLOW5/BLOW5 output files will be determined based on the first occurrence of run ID seen by f2s. This can be used to convert FAST5 files from different samples in a single command if the user does not further require the original run IDs.
*  `--retain`:<br/>
	Retain the same directory structure in the converted output as the input (experimental).
*  `-h, --help`:<br/>
   Prints the help menu.