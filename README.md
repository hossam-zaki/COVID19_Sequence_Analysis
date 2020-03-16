# COVID19_Sequence_Analysis

## Overview

This project will study the sequence of the Novel Coronavirus that started in Wuhan, China. We will attempt to compare this sequence to other viruses that have been known to cause other epidemics, such as H1N1, the Spanish Flu, and HIV/AIDS. We aim to learn more which viruses cause epidemics, and to study their comparisons. We also hope to expand this to other viruses, and hopefully gain some more insight to prevent the spread of viruses in the future.

## Usage

This project is in the beginning stages, so we will use standard alignment algorithms.

For global, local, and affine alignment, we will be using the following parameters

```--seqs``` Following this tag should be a file with the two sequences, each on seperate lines. Please see ```data/test_cases/``` for examples

```--sm``` This is a scoring matrix. This will score different nucleotide alignments. please see ``scoringMatrices/``` for examples

```--gapPen``` This is the penalty per gaps.

For affine alignment, there is an added parameter:

```--startPen``` This is the penalty for starting a set of gaps.
