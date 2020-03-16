import sys

class local_aligning:
    def __init__(self, seqs, scoring_matrix, gap_penalty):
        self.seq_file = seqs
        self.scoring_matrix_file = scoring_matrix
        self.scoring_matrix = None
        self.gap_penalty = int(gap_penalty)
        self.seq1 = None
        self.seq2 = None
        self.scoring_dict = {}
        self.alignment_matrix = None
        self.direction_matrix = None
        self.max_Rindex = 0
        self.max_Cindex = 0
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
    def align_locally(self):
        self.alignment_matrix = [[0 for y in range (0, len(self.seq2) + 1)] for z in range(0, len(self.seq1) + 1)]
        self.direction_matrix = [[None for y in range (0, len(self.seq2) + 1)] for z in range(0, len(self.seq1) + 1)]
        for x in range (0, len(self.seq1)+ 1):
            self.alignment_matrix[x][0] = 0
            self.direction_matrix[x][0] = -1
        for y in range (0, len(self.seq2) + 1):
            self.alignment_matrix[0][y] = 0
            self.direction_matrix[0][y] = 0
        max_value = -float('inf')
        for i in range (1, len(self.alignment_matrix)):
            for j in range (1, len(self.alignment_matrix[0])):
                self.alignment_matrix[i][j] = max(0, (self.alignment_matrix[i-1][j] + self.gap_penalty), (self.alignment_matrix[i][j-1] + self.gap_penalty), (self.alignment_matrix[i-1][j-1] + self.scoring_dict[(self.seq1[i-1].upper(), self.seq2[j-1].upper())]))
                if self.alignment_matrix[i][j] > max_value:
                    max_value = self.alignment_matrix[i][j]
                    self.max_Rindex = i
                    self.max_Cindex = j
                if self.alignment_matrix[i][j] == self.alignment_matrix[i-1][j] -1:
                    self.direction_matrix[i][j] = -1
                elif self.alignment_matrix[i][j] == self.alignment_matrix[i][j-1] - 1:
                    self.direction_matrix[i][j] = 0
                else:
                    self.direction_matrix[i][j] = 1
    def trace_back(self):
        row_index = len(self.seq1)
        col_index = len(self.seq2)
        top_string = ""
        bottom_string = ""
        boolean = True
        while boolean:
            value = self.direction_matrix[self.max_Rindex][self.max_Cindex]
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
                col_index= col_index - 1
                row_index = row_index - 1
            if self.alignment_matrix[row_index][col_index] == 0:
                boolean = False
        self.final_top_seq=(self.reverse_str(top_string))
        self.final_bottom_seq=(self.reverse_str(bottom_string))
    def local_align(self):
        self.retrieve_seqs()
        self.get_scoring_matrix()
        self.align_locally()
        self.trace_back()
        print(self.final_top_seq)
        print(self.final_bottom_seq)
        print(self.alignment_matrix[len(self.seq1)][len(self.seq2)])

    def reverse_str(self, string):
        new_str = ""
        for i in range (len(string) - 1, -1, -1):
            new_str+=string[i]
        return new_str
                    
if __name__ == "__main__":
    seqs = sys.argv[1]
    scoring_matrix = sys.argv[2]
    gap_penalty = sys.argv[3]
    align = local_aligning(seqs, scoring_matrix, gap_penalty)
    align.local_align()