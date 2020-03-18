# COVID19_Sequence_Analysis

## Overview

This project will study the sequence of the Novel Coronavirus that started in Wuhan, China. We will attempt to compare this sequence to other viruses that have been known to cause other epidemics, such as H1N1, the Spanish Flu, and HIV/AIDS. We aim to learn more which viruses cause epidemics, and to study their comparisons. We also hope to expand this to other viruses, and hopefully gain some more insight to prevent the spread of viruses in the future.

## Usage

This project is in the beginning stages, so we will use standard alignment algorithms.

For global, local, and affine alignment, we will be using the following parameters

```--seq1``` Following this tag should be a file with the first sequence. This will be the first sequence in the directory and file name

```--seq2``` Following this tag should be a file with the second sequence

```--sm``` This is a scoring matrix. This will score different nucleotide alignments. please see ``scoringMatrices/``` for examples

```--gapPen``` This is the penalty per gaps.

For affine alignment, there is an added parameter:

```--startPen``` This is the penalty for starting a set of gaps.
