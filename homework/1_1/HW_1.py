from Bio import SeqIO


def Hamming_distance(seq1, seq2):
    dist = 0
    for a,b in zip(seq1, seq2):
        if a != b:
            dist+=1
    return dist

def Tearest_substring(short, long):
    if len(short) > len(long):
        short, long = long, short

    i = 0
    pos = -1
    min = float('inf')
    while len(short) + i <= len(long):
        dist = Hamming_distance(short,long[i:len(short)+i])
        if min > dist:
            min = dist
            pos = i
        i+=1

    print(f"Позиция начала подстроки: {pos} \n",
          f"Последовательность подстроки: {long[pos:len(short)+pos]}\n",
          f"Расстояние Хэмминга: {min}")
    return

def Levenshtein(seq_1, seq_2):
    n, m = len(seq_1), len(seq_2)
    if n > m:
        seq_1, seq_2 = seq_2, seq_1
        n, m = m, n
    ans = [i for i in range(n+1)]

    for i in range(1, m + 1):
        for j in range(n + 1):
            if j == 0:
                previous = i
            else:
                if seq_1[j-1] == seq_2[i-1]:
                    same = 0
                else:
                    same = 1
                new_previous = min(previous + 1, ans[j] + 1, ans[j-1] + same)
                ans[j-1] = previous
                previous = new_previous
    return previous


# 1.Хэмминг (1 балл)
sequences = list(SeqIO.parse(open(r'data/gattaca.fasta'), "fasta"))
print(Hamming_distance(sequences[0].seq, sequences[1].seq))

# 2. Ближайшая подстрока (1 балл)
sequences = list(SeqIO.parse(open(r'data/f8.fasta'), "fasta"))
Tearest_substring(sequences[0].seq, sequences[1].seq)

# 3. Левенштейн (4 балла)
print(Levenshtein(sequences[0].seq, sequences[1].seq))




