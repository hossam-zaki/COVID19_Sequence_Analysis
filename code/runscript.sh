#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=64G
#SBATCH -t 1:00:00

python3 global_alignment.py --seq1 ../data/sequences/COVID-19_Genome_Sequence.txt  --seq2 ../data/sequences/Lyme_Genome_Sequence.txt --sm ../scoringMatrices/standard.m --gapPen -1
