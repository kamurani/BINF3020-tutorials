"""
BINF3020 tutorial 
Week 5


Q2


Cam McMenamie
"""


import sys
import click  

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
            score = match_score(a, b, matrix)
            if score > threshold:
                result.append((a, b, score))
    return result


@click.command()
@click.argument('Q')
@click.argument('D')
@click.option(
    '--threshold', 
    '-t', 
    type=click.INT,
    default=3, 
    help="Threshold score for k-mer pair match.",
)
@click.option(
    '-k',
    type=click.INT,
    default=2,
    help="k-mer length",
)
def main(
    q,
    d,
    threshold, 
    k,
):

    matrix_filename = "blosum62.txt"
    blosum = Matrix(matrix_filename)
    
    if not q:
        q = "ACEDECADE" # Query sequence
    
    if not d:
        d = "REDCEDKL"  # Data sequence

    
    # Get k-mers
    q_2mers = generate_kmers(seq=q, k=k)
    d_2mers = generate_kmers(seq=d, k=k)

    # Get matches 
    words = threshold_kmer_scores(q_2mers, d_2mers, matrix=blosum, threshold=threshold)
    
    
    for w in words:
        (q, d, s) = w
        print(q, d, s) 
        
    
    
    # Optional printing
    if False:
        # Print
        print(f"Q: {q}")
        print(f"D: {d}")
        print(f"T = {threshold}")
        print(f"k = {k}")
        print("------------------------")
        
        for w in words:
            (q, d, s) = w
            print(q, d, s) 

        print("------------------------")
        print("Total: ", len(words))

if __name__ == "__main__":
    main()
