def hamming_distance(seq1, seq2):
    """Вычисление расстояния Хэмминга"""
    assert len(seq1) == len(seq2)

    count = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            count += 1
    return count

