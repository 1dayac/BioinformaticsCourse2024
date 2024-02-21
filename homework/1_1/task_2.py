"""2. Ближайшая подстрока (1 балл)  
Даны 2 строки над алфавитом ATGC, реализуйте алгоритм,
 который позволит находить подстроку в более длинной 
 строке, которая ближе всего по расстоянию Хэминга к 
 боллее короткой. Верните позицию начала этой подстроки, 
 саму подстроку и рассстояние хэмминга.
"""

from task_1 import hamming_distance

def nearest_substring(seq1, seq2):
    if len(seq1) < len(seq2):
        seq1, seq2 = seq2, seq1

    t = len(seq2)
    distance, position = float("inf"), 0
    for i in range(len(seq1)-t+1):
        curr = hamming_distance(seq1[i:i+t], seq2)
        if curr < distance:
            distance, position = curr, i
    return position, seq1[position:position+t], distance


from Bio import SeqIO
f8 = SeqIO.parse("data/f8.fasta", "fasta")

seq = [data.seq for data in f8.records]
print(nearest_substring("FFFFFFFFFFFqwertyyyyyyyyyyyyy", "qqqwerty"))
print(nearest_substring(*seq))
print(nearest_substring(seq[0][:4], seq[1][:4]))
print(nearest_substring(seq[0][:26], seq[1][:4]))
print(nearest_substring(seq[0], seq[1][637:]))

