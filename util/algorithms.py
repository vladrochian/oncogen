def levenshtein_distance(s1, s2):
    best = [[0 for _ in range(len(s2) + 1)] for _ in range(len(s1) + 1)]
    for i in range(len(s2) + 1):
        best[0][i] = i
    for i in range(1, len(s1) + 1):
        best[i][0] = i
        for j in range(1, len(s2) + 1):
            best[i][j] = min(best[i - 1][j] + 1, best[i][j - 1] + 1,
                             best[i - 1][j - 1] + (0 if s1[i - 1] == s2[j - 1] else 1))
    return best[len(s1)][len(s2)]
