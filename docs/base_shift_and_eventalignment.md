# Base shift calculation

* User can shift the base sequence to the left by `n` number of bases by providing the argument `--base_shift -n` to `plot` and `plot_pileup` commands.
* This is helpful to correct the signal level to the base. A negative `n` value will shift the base sequence to the left.
* Base shift concept was motivated from the idea presented in [pore_model document](pore_model.md) about the most contributing base to the current level.
* We can programmatically calculate a `base_shift` to nicely align the signal to the base (color). The calculation is implemented [here](https://github.com/hiruna72/squigualiser/blob/2389379fa8898bf78fd695b3bddac982213ea951/src/plot_utils.py#L194).
* However, the user is adviced to use `--profile` (documented [here](docs/profiles.md)) which automatically sets the `--base_shift`.

## Example 1

![image](figures/variants/chr22/fig3.png)
*Figure 1 - variant detection [documented here](pipeline_variant_detection_extended.md)*

* Fig. 1 has a base shift of `-6` respectively.
* At site `A21` there is a heterozygous variant `A/T`.
* If the base shift was some other value instead of `-6` the difference in the signals will not align with base.

## Example 2
* Fig. 2 and 3 have the same `dna_r10.4.1_e8.2_400bps` signal pileup with a base shift of `0` and `-6` respectively.
* The base colors in Fig. 2 do not nicely aligned to the signal. 
* However, in Fig. 3 the signal is moving from low to high when a base `T` is met.
* These signal pileups were generated using f5c eventalign.
* F5c used `dna_r10.4.1_e8.2_400bps` pore model to align these signals. Its most contributing base index is `-6` and hence the appropriate base shift in this scenario is also `-6`.

![image](figures/base_shift/testcase_20.8_base_shift_0.png)

*Figure 2: base shift 0 [link](https://hiruna72.github.io/squigualiser/docs/figures/pileup/pileup_testcase-20.8.html)*

![image](figures/base_shift/testcase_20.1_base_shift_6.png)

*Figure 3: base shift -6 [link](https://hiruna72.github.io/squigualiser/docs/figures/pileup/pileup_testcase-20.1.html)*

## Example 3
* Consider a `dna_r10.4.1_e8.2_400bps` forward and reverse mapped pileups for the same genomic region.
* Fig.4 has `0` base shift for both tracks.
* In Fig. 5, note that the reverse mapped pileup has a `-2` base shift. This is because the signal sequencing direction is from right to left ([more information](base_shift_of_reverse_mapped_reads.md)).
* In Fig. 5, both pileups have the signals going from low to high when a base `T` is met. 

![image](figures/base_shift/plot_tracks_base_shift_0.png)

*Figure 4: base shift 0, 0 [link](https://hiruna72.github.io/squigualiser/docs/figures/plot_tracks/plot_tracks_testcase-30.3.html)*

![image](figures/base_shift/plot_tracks_base_shift_6.png)

*Figure 5: base shift -6, -2 [link](https://hiruna72.github.io/squigualiser/docs/figures/plot_tracks/plot_tracks_testcase-30.6.html)*

## Example 4

Fig. 6 shows the forward and reverse mapped pileups generated using f5c eventalign for `dna_r9.4.1_450bps` data. F5c used the 6mer model ([more information](pore_model.md)).

![image](figures/base_shift/dna_r9.4.1_450bps_fast_forward_reverse.png)
*Figure 6: base shift -2, -3*

## Example 4

* Fig. 7 shows the forward mapped pileups generated using f5c eventalign for `rna_r9.4.1_70bps` data. F5c used the rna 5mer model ([more information](pore_model.md)).
* Note that `squigualiser` always plot the RNA reads in its correct sequencing direction (reverse mapped RNA reads are skipped; reverse mapped RNA reads exist if a genome was used as the reference instead of a transcriptome).

![image](figures/base_shift/rna_r9.4.1_70bps_forward_eventalign.png)
*Figure 7: base shift -3*