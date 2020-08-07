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
