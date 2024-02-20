from Bio import SeqIO, pairwise2
from distance import hamming, levenshtein


def hamming_distance(a, b):
    return sum([a[i] != b[i] for i, _ in enumerate(a)])


def levenshtein_distance(a, b):
    n, m = len(a), len(b)
    d = [[0] * m for _ in range(n)]
    for i in range(n):
        d[i][0] = d[i - 1][0] + 1
    for j in range(m):
        d[0][j] = d[0][j - 1] + 1
    for i in range(n):
        for j in range(m):
            if a[i] != b[j]:
                d[i][j] = min(d[i - 1][j] + 1, d[i][j - 1] + 1,
                              d[i - 1][j - 1] + 1)
            else:
                d[i][j] = d[i - 1][j - 1]
    return d[n - 1][m - 1]


def closest_substring(pattern, s):
    n = len(pattern)
    min_dist = float('inf')
    closest_i = -1
    closest_pattern = ''
    for i in range(len(s) - n + 1):
        cur_pattern = s[i:i + n]
        dist = hamming_distance(pattern, cur_pattern)
        if dist < min_dist:
            min_dist = dist
            closest_i = i
            closest_pattern = cur_pattern
    return closest_i, closest_pattern, min_dist


def __test_task_1(fasta_sequences_1, fasta_sequences_2):
    # Define custom scoring matrix for matches and mismatches
    match_score = 2
    mismatch_penalty = -1
    gap_penalty = -0.5

    # Perform pairwise sequence alignment with custom scoring
    alignments_custom = \
        pairwise2.align.globalms(*fasta_sequences_1, match_score, mismatch_penalty, gap_penalty, gap_penalty,
                                 one_alignment_only=True)[0]

    fasta_sequences_1 = [alignments_custom.seqA, alignments_custom.seqB]
    print('__test_1 (hamming)')
    print('data/f8.fasta aligned results:')
    print('True: ', hamming(*fasta_sequences_1))
    print('Test: ', hamming_distance(*fasta_sequences_1))
    print('data/gattaca.fasta results:')
    print('True: ', hamming(*fasta_sequences_2))
    print('Test: ', hamming_distance(*fasta_sequences_2))


def __test_task_2(fasta_sequences_1, fasta_sequences_2):
    for i, pattern in enumerate(fasta_sequences_2):
        for j, s in enumerate(fasta_sequences_1):
            print(f'__test_2 (gattaca_{i} / f8_{j})')
            ind, found_str, dist = closest_substring(pattern, s)
            print(f'Pattern: {pattern}')
            print(f'Index: {ind}')
            print(f'Found string: {found_str}')
            print(f'Hamming distance: {dist}')
            print()


def __test_task_3(fasta_sequences_1, fasta_sequences_2):
    print('__test_3 (levenshtein)')
    print('data/f8.fasta results:')
    print('True: ', levenshtein(*fasta_sequences_1))
    print('Test: ', levenshtein_distance(*fasta_sequences_1))
    print('data/gattaca.fasta results:')
    print('True: ', levenshtein(*fasta_sequences_2))
    print('Test: ', levenshtein_distance(*fasta_sequences_2))


if __name__ == "__main__":
    fasta_sequence_1 = [record.seq for record in
                        SeqIO.parse("data/f8.fasta", "fasta")]
    fasta_sequence_2 = [record.seq for record in
                        SeqIO.parse("data/gattaca.fasta", "fasta")]
    __test_task_1(fasta_sequence_1, fasta_sequence_2)
    __test_task_2(fasta_sequence_1, fasta_sequence_2)
    __test_task_3(fasta_sequence_1, fasta_sequence_2)
