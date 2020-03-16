import global_alignment
import sys

def build_files(base, seq1, number):
    f=open(f"{number}seq.txt", "w+")
    base=open(base, "r")
    seq=open(seq1, "r")
    f.write(base.read())
    f.write(seq.read())
    f.close()
    base.close()
    seq.close()
    

if __name__ == "__main__":
    base = sys.argv[1]
    seq1 = sys.argv[2]
    seq2 = sys.argv[3]
    seq3 = sys.argv[4]
    array = [seq1, seq2, seq3]
    number = 1
    for i in array: 
        build_files(base, i, number)
        number+=1

