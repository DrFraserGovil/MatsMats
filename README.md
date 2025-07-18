# MATSMATS: Matrix-Sequence Matching &amp; Scoring


MATSMATS is part of a project to improve a theoretical model of predicting Transcription Factor Binding Sites (TFBS) by comparing its output to the better-performing models submitted to the [IBIS Challenge](https://ibis.autosome.org/home). Though these models outperformed the theoretical model on some counts, many of them are 'Black Boxes' which do not necessarily provide biological insights; despite performing better on benchmarks. 

## Problem Statement

When computing the 'scores' for a given DNA sequence, we are interested in determining both *which* potential PWM it matches with, and which substring is matching them ost strongly. [Standard tools](https://github.com/jhkorhonen/MOODS/tree/master) for determining this do exist; however for the use case where we want precisely one best-matching score these tools can form a significant bottleneck in the analysis pipeline. 

MATSMATS is therefore a streamlined, highly optimised code which:
* Takes as input:
	1. A set of PWMs in a standard file format
	2. A set of target sequences
* Returns a single file containing for each input sequence
	1. The sequence
	2. The maximal subsequence
	3. A reference to the maximally-activating PWM
	4. The score associated with that matrix


### Licensing

Portions of this code are inspired by the [MOODS repository](https://github.com/jhkorhonen/MOODS/tree/master), released under a GPL-3.0 license. 

Many of the modules in ``src/tools`` are replicated wholesale from [RAMICES](https://github.com/DrFraserGovil/RAMICES_II), a code written by the same authors.

This code is released under a compatible GPL-3.0 license. 