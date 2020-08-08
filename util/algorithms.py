def extract_first_match(s: str, starts_with: str, ends_with: str, search_from=1) -> str:
    """
    Return first substring that starts and ends with given sequences
    :param s: original sequence
    :param starts_with: start of subsequence to match
    :param ends_with: end of subsequence to match
    :param search_from: start position for lookup
    :return: matching subsequence
    """
    suffix = s[s.find(starts_with, search_from - 1):]
    return suffix[:suffix.find(ends_with) + len(ends_with)]


def lev_distance(s1: str, s2: str) -> int:
    """
    Levenshtein distance between two sequences.

    Time complexity: |s1| * |s2|

    Memory complexity: |s2|

    :param s1: first sequence
    :param s2: second sequence
    :return: edit distance
    """
    best = [[0 for _ in range(len(s2) + 1)] for _ in range(2)]
    for i in range(len(s2) + 1):
        best[0][i] = i
    for i in range(1, len(s1) + 1):
        best[i & 1][0] = i
        for j in range(1, len(s2) + 1):
            best[i & 1][j] = min(best[(i & 1) ^ 1][j] + 1, best[i & 1][j - 1] + 1,
                                 best[(i & 1) ^ 1][j - 1] + (0 if s1[i - 1] == s2[j - 1] else 1))
    return best[len(s1) & 1][len(s2)]


def generate_diff(s1: str, s2: str) -> str:
    """
    Generate graphically aligned strings with highlighted differences.

    Time complexity: |s1| * |s2|

    Memory complexity: |s1| * |s2|

    :param s1:
    :param s2:
    :return: printable representation of differences
    """
    best = [[0 for _ in range(len(s2) + 1)] for _ in range(len(s1) + 1)]
    prv = [[(0, 0) for _ in range(len(s2) + 1)] for _ in range(len(s1) + 1)]
    for i in range(len(s2) + 1):
        best[0][i] = i
        prv[0][i] = 0, i - 1
    for i in range(1, len(s1) + 1):
        best[i][0] = i
        prv[i][0] = i - 1, 0
        for j in range(1, len(s2) + 1):
            best[i][j] = best[i - 1][j - 1] + (0 if s1[i - 1] == s2[j - 1] else 1)
            prv[i][j] = i - 1, j - 1
            if best[i - 1][j] + 1 < best[i][j]:
                best[i][j] = best[i - 1][j] + 1
                prv[i][j] = i - 1, j
            if best[i][j - 1] + 1 < best[i][j]:
                best[i][j] = best[i][j - 1] + 1
                prv[i][j] = i, j - 1
    first_line = []
    second_line = []
    third_line = []
    i = len(s1)
    j = len(s2)
    while i > 0 or j > 0:
        prv_i, prv_j = prv[i][j]
        if prv_i == i - 1 and prv_j == j - 1:
            first_line.append(' ' if s1[i - 1] == s2[j - 1] else 'v')
            second_line.append(s1[i - 1])
            third_line.append(s2[j - 1])
        elif prv_i == i - 1:
            first_line.append('*')
            second_line.append(s1[i - 1])
            third_line.append(' ')
        else:
            first_line.append('*')
            second_line.append(' ')
            third_line.append(s2[j - 1])
        i, j = prv_i, prv_j
    first_line.reverse()
    second_line.reverse()
    third_line.reverse()
    return '\n'.join([''.join(first_line), ''.join(second_line), ''.join(third_line)]) + '\n'
