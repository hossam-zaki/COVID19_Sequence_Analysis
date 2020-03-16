import sys
import argparse
from argparse import ArgumentParser

class affine_align:
    def __init__(self, seqs, scoring_matrix, gap_penalty, start_penalty):
        self.seq_file = seqs
        self.scoring_matrix_file = scoring_matrix
        self.scoring_matrix = None
        self.gap_penalty = int(gap_penalty)
        self.start_penalty = int(start_penalty)
        self.seq1 = None
        self.seq2 = None
        self.scoring_dict = {}
        self.direction_matrix = None
        self.V = None
        self.G = None
        self.E = None
        self.F = None
        self.final_top_seq = None
        self.final_bottom_seq = None
    def retrieve_seqs(self):
        with open(self.seq_file) as f:
            for line in f:
                if self.seq1 == None:
                    self.seq1 = line.strip()
                else:
                    self.seq2 = line.strip()
    def get_scoring_matrix(self):
        with open(self.scoring_matrix_file) as f:
            index = 0
            for line in f:
                if self.scoring_matrix == None:
                    self.scoring_matrix = [[0 for y in range (len(line.split()))] for z in range(len(line.split()))]
                for i in range(0, len(line.split())):
                    self.scoring_matrix[index][i] = line.split()[i]
                index+=1
            for q in range(1, len(self.scoring_matrix)):
                nuc = self.scoring_matrix[q][0].upper()
                for w in range(1, len(self.scoring_matrix)):
                    score = int(self.scoring_matrix[q][w])
                    othernuc = self.scoring_matrix[0][w].upper()
                    self.scoring_dict[(nuc, othernuc)] = score
    def initialize_arrays(self):
        self.V = [[None for y in range (0, len(self.seq2) + 1)] for z in range(0, len(self.seq1) + 1)]
        self.E = [[None for y in range (0, len(self.seq2) + 1)] for z in range(0, len(self.seq1) + 1)]
        self.F = [[None for y in range (0, len(self.seq2) + 1)] for z in range(0, len(self.seq1) + 1)]
        self.G = [[None for y in range (0, len(self.seq2) + 1)] for z in range(0, len(self.seq1) + 1)]
        self.direction_matrix = [[0 for y in range (0, len(self.seq2) + 1)] for z in range(0, len(self.seq1) + 1)]
        self.V[0][0] = 0
        for x in range (0, len(self.seq1)+ 1):
            self.E[x][0] = -float('inf') 
            self.G[x][0] = -float('inf')
            # if x >= 1:
            #     self.F[x][0] = max(self.F[x-1][0] + self.gap_penalty, self.V[x-1][0] + self.gap_penalty + self.start_penalty)
            #     self.V[x][0] = self.F[x][0]
            self.direction_matrix[x][0] = -1
        for y in range (0, len(self.seq2) + 1):
            self.F[0][y] = -float('inf')
            self.G[0][y] = -float('inf')
            # if y >= 1:
            #     self.E[0][y] =  max(self.E[0][y-1] +self.gap_penalty, self.V[0][y-1] + self.gap_penalty + self.start_penalty)
            #     self.V[0][y] = self.E[0][y]
            self.direction_matrix[0][y] = 0
        for a in range(1, len(self.seq1) + 1):
            self.F[a][0] = max(self.F[a-1][0] + self.gap_penalty, self.V[a-1][0] + self.gap_penalty + self.start_penalty)
            self.V[a][0] = self.F[a][0]
        for b in range(1, len(self.seq2) + 1):
            self.E[0][b] =  max(self.E[0][b-1] +self.gap_penalty, self.V[0][b-1] + self.gap_penalty + self.start_penalty)
            self.V[0][b] = self.E[0][b]

    def fill_matrices(self):
        for i in range (1, len(self.V)):
            for j in range (1, len(self.V[0])):
                self.G[i][j] = self.V[i-1][j-1] + self.scoring_dict[(self.seq1[i-1], self.seq2[j-1])]
                self.E[i][j] = max(self.E[i][j-1] + self.gap_penalty, self.V[i][j-1] + self.gap_penalty + self.start_penalty)
                self.F[i][j] = max(self.F[i-1][j] + self.gap_penalty, self.V[i-1][j] + self.gap_penalty + self.start_penalty)
                self.V[i][j] = max(self.F[i][j], self.E[i][j], self.G[i][j])
                if (self.V[i][j] == self.F[i][j]):
                    self.direction_matrix[i][j] = -1
                elif (self.V[i][j] == self.E[i][j]):
                    self.direction_matrix[i][j] = 0
                else:
                    self.direction_matrix[i][j] = 1

    def trace_back(self):
        row_index = len(self.seq1)
        col_index = len(self.seq2)
        top_string = ""
        bottom_string = ""
        while row_index != 0 or col_index != 0:
            value = self.direction_matrix[row_index][col_index]
            if value == -1 :
                bottom_string += "-"
                top_string += self.seq1[row_index -1]
                row_index = row_index - 1
            if value == 0 :
                top_string += "-"
                bottom_string += self.seq2[col_index - 1]
                col_index = col_index - 1
            if value == 1:
                bottom_string += self.seq2[col_index -1]
                top_string += self.seq1[row_index - 1]
                col_index = col_index - 1
                row_index = row_index - 1
        self.final_top_seq=(self.reverse_str(top_string))
        self.final_bottom_seq=(self.reverse_str(bottom_string))

    def affine_aligning(self):
        self.retrieve_seqs()
        self.get_scoring_matrix()
        self.initialize_arrays()
        self.fill_matrices()
        self.trace_back()
        print(self.final_top_seq)
        print(self.final_bottom_seq)
        print(self.V[len(self.seq1)][len(self.seq2)])


    def reverse_str(self, string):
        new_str = ""
        for i in range (len(string) - 1, -1, -1):
            new_str+=string[i]
        return new_str


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--seqs',type=str,required=True)
    parser.add_argument('--sm',type=str,required=True)
    parser.add_argument('--gapPen', type=str, required=True)
    parser.add_argument('--startPen', type=str, required=True)
    config = parser.parse_args()
    align = affine_align(config.seqs, config.sm, config.gapPen, config.startPen)
    align.affine_aligning()