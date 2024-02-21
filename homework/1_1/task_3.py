"""3. Левенштейн (4 балла)  
Реализуйте алгоритм вычисления расстояния 
редактирования. *O(n•m)* по времени и *O(min(n, m))* по памяти 
(обратите внимание на то, что возвращать способ получения из 
одной строки другую не нужно). Алгоритм возвращаеет одно число - 
расстояние редактирования между входными строками.
"""


import numpy as np

def levenshtein_distance(seq1, seq2):
    n, m = len(seq1), len(seq2)

    dp = np.zeros((n + 1, m + 1), dtype=int)

    for i in range(n + 1):
        dp[i][0] = i
    for j in range(m + 1):
        dp[0][j] = j

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            count = 0 if seq1[i - 1] == seq2[j - 1] else 1
            dp[i][j] = min(dp[i - 1][j] + 1,  # удаление
                          dp[i][j - 1] + 1,  # вставка
                          dp[i - 1][j - 1] + count)  # замена

    return dp[n][m]

print(levenshtein_distance("GATTACA", "AAGAGTAC"))
