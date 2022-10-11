"""
BINF3020 tutorial 
Week 5


Q2


Cam McMenamie
"""


import sys

from typing import List, Dict, Tuple


class InvalidPairException(Exception):
  pass

"""
Matrix class used from https://github.com/jwintersinger blosum.py 
"""
class Matrix:
  def __init__(self, matrix_filename: str):
    self._load_matrix(matrix_filename)

  def _load_matrix(self, matrix_filename: str):
    with open(matrix_filename) as matrix_file:
      matrix = matrix_file.read()
    lines = matrix.strip().split('\n')

    header = lines.pop(0)
    columns = header.split()
    matrix = {}

    for row in lines:
      entries = row.split()
      row_name = entries.pop(0)
      matrix[row_name] = {}

      if len(entries) != len(columns):
        raise Exception('Improper entry number in row')
      for column_name in columns:
        matrix[row_name][column_name] = entries.pop(0)

    self._matrix = matrix

  def lookup_score(self, a: str, b: str) -> int:
    a = a.upper()
    b = b.upper()

    if a not in self._matrix or b not in self._matrix[a]:
      raise InvalidPairException('[%s, %s]' % (a, b))
    return self._matrix[a][b]


"""
Generate all k-mers of a sequence.
"""
def generate_kmers(seq: str, k: int) -> List[str]:
    
    kmers = [seq[i:i+k] for i in range(len(seq)-k+1)]
    return kmers

"""
Generate counts for each unique k-mer in a sequence.
"""
def generate_kmer_count(seq: str, k: int) -> Dict[str, int]:

    kmers = generate_kmers(seq, k)

    count = {}
    for k in kmers:
        if k not in count:
            count[k] = 0
        # increment
        count[k] += 1
        
    return count

"""
Score kmers
"""
def match_score(aseq: str, bseq: str, m: Matrix) -> int:

    assert len(aseq) == len(bseq), "Error: k-mers should be of same length"

    score = 0
    for i in range(len(aseq)):
        
        score += int(m.lookup_score(aseq[i], bseq[i]))

    return score

"""
Return k-mer matches above threshold
"""
def threshold_kmer_scores(
    q_kmers: List[str], 
    d_kmers: List[str], 
    matrix: Matrix,
    threshold: int,
) -> List[Tuple]:

    result = []
    for a in q_kmers: 
        for b in d_kmers:
            if match_score(a, b, matrix) > threshold:
                result.append((a, b))
    return result

def main():

    matrix_filename = "blosum62.txt"
    blosum = Matrix(matrix_filename)

    q = "ACEDECADE" # Query sequence
    d = "REDCEDKL"  # Data sequence




    # Get k-mers
    k = 2
    q_2mers = generate_kmers(seq=q, k=k)
    d_2mers = generate_kmers(seq=d, k=k)

    # Print scores 
    
    

    
    

    # Print kmer match scores
    T = 3 # threshold

    print(f"Q: {q}")
    print(f"D: {d}")
    print(f"T = {T}")
    print(f"k = {k}")
    

    if True:
        print("------------------------")
        for i in q_2mers:
            for j in d_2mers:
                score = match_score(i, j, blosum)
                if score > T: print(i, j, score)


    words = threshold_kmer_scores(q_2mers, d_2mers, matrix=blosum, threshold=T)

    print("------------------------")
    print("Total: ", len(words))

if __name__ == "__main__":
    main()