# This function extracts the longest increasing and decreasing subsequences from
# a sequence of integers.
# The function takes a list of integers and returns two lists of integers 
# containing the longest subsequences. If it gets passed an empty list, None,
# or something other then an integer list, the function returns two times None.
#   
# The subsequence length is not necessarily unique, in this case the elements of
# the longest increasing subsequence should be as small as possible, while the 
# elements of the longest decreasing subsequence should be as large as possible.

import numbers
import math
import array
import random
import time


def is_valid(sequence):
    if (sequence == []) or (sequence is None) or (not isinstance(sequence, list)):
        return False

    for n in sequence:
        if not isinstance(n, numbers.Integral):
            return False

    return True


# def get_lis_or_lds_wiki(sequence, run_lis=True):
#     M = [None] * (len(
#         sequence) + 1)  # at index i, M[i] is the header of the best (in term of longest and smallest) subseq seen so far
#     P = [None] * len(sequence)  # element i of P means at index i, the prev element in subseq is at P[i]
#
#     length = 0
#     for i in range(len(sequence)):
#         lo = 1
#         hi = length
#         while lo <= hi:
#             mid = math.ceil((lo + hi) / 2)
#             is_increasing = (sequence[M[mid]] < sequence[i]) if run_lis else (sequence[M[mid]] > sequence[i])
#             if is_increasing:
#                 lo = mid + 1
#             else:
#                 hi = mid - 1
#
#         new_lo = lo
#         P[i] = M[new_lo - 1]
#         M[new_lo] = i
#
#         if new_lo > length:
#             length = new_lo
#
#     # increasing_subseq = [-1] * len(sequence)
#     longest_subseq = []
#     k = M[length]
#     for i in range(length, 0, -1):
#         # increasing_subseq[i] = sequence[k]
#         longest_subseq.insert(0, sequence[k])
#         k = P[k]
#     return longest_subseq


def bin_search(arr, x, rev=False):
    l = 0
    r = len(arr)
    while l < r:
        mid = (l + r) // 2
        # if arr[mid] == x:
        #     return mid + 1
        pred = arr[mid] < x if not rev else arr[mid] > x
        if pred:
            l = mid + 1
        else:
            r = mid
    return l


def get_lis_or_lds(sequence, lis=True):
    best_subseqs_by_len = [[]] * (len(sequence) + 1)  # best_subseqs_by_len[i] stores the best subseq of length i
    longest = 0
    for i, num in enumerate(sequence):
        if i == 0:
            longest = 1
            best_subseqs_by_len[longest] = [num]
        else:
            for l in range(1, longest + 1):
                curr_best = best_subseqs_by_len[l]
                pos = bin_search(curr_best, num, rev=(not lis))

                if pos == len(curr_best):
                    candidate = best_subseqs_by_len[l] + [num]
                    if best_subseqs_by_len[l + 1]:
                        next_len_best = best_subseqs_by_len[l + 1]
                        assert len(candidate) == len(next_len_best)
                        is_increasing = num > next_len_best[len(next_len_best) - 1] if lis else num < next_len_best[
                            len(next_len_best) - 1]
                        if is_increasing:
                            continue
                        for j in range(len(candidate)):
                            is_increasing = next_len_best[j] > candidate[j] if lis else next_len_best[j] < candidate[j]
                            if is_increasing:
                                best_subseqs_by_len[l + 1] = candidate
                                longest = max(len(candidate), longest)
                                break
                    else:
                        best_subseqs_by_len[l + 1] = candidate
                        longest = max(len(candidate), longest)

                elif pos == len(curr_best) - 1:
                    curr_best[pos] = num

    return best_subseqs_by_len[longest]


def get_longest_subsequences(sequence):
    if not is_valid(sequence):
        return None, None
    # return get_lis_or_lds(sequence), get_lis_or_lds(sequence, run_lis=False)
    return get_lis_or_lds(sequence), get_lis_or_lds(sequence, False)


def print_all(seq):
    if len(seq) <= 1000:
        print("Self ", get_lis_or_lds(seq), get_lis_or_lds(seq, lis=False))
        # print("Wiki ", get_lis_or_lds_wiki(seq), get_lis_or_lds_wiki(seq, False))
        print("")
    else:
        start = time.time()
        res_my_lis = get_lis_or_lds(seq)
        my_lis_runtime = time.time() - start

        start = time.time()
        # res_wiki_lis = get_lis_or_lds_wiki(seq)
        wiki_lis_runtime = time.time() - start
        # print(f"LIS func performs correctly: {res_my_lis == res_wiki_lis}")
        print(f"Runtime (my_lis, wiki_lis): {my_lis_runtime} s, {wiki_lis_runtime} s")

        start = time.time()
        res_my_lds = get_lis_or_lds(seq, False)
        my_lds_runtime = time.time() - start
        start = time.time()
        # res_wiki_lds = get_lis_or_lds_wiki(seq, False)
        wiki_lds_runtime = time.time() - start
        # print(f"LDS func performs correctly: {res_my_lds == res_wiki_lds}")
        print(f"Runtime (my_lds, wiki_lds): {my_lds_runtime} s, {wiki_lds_runtime} s")

        # Last run with 100_000 numbers:
        # Runtime (my_lis, wiki_lis): 219.30497789382935 s, 0.7740004062652588 s
        # Runtime (my_lis, wiki_lis): 261.25200057029724 s, 0.7490057945251465 s


if __name__ == '__main__':
    # seq = [0, 8, 4, 12, 2, 10]
    print_all([0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15])
    # seq = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11]
    print_all([-2202, -9644, 6681, -374, 4020, 9549, 3069, 9205, 2988, 4206])

    print_all([-1599, 4755, -9887, 9436, -5354, -9141, 659, 4247, -5844, -1521])

    print_all([5649, 6081, 2164, 3617, 1207, -4687, 402, 7447, -6240, 4155])
    print_all([15, 27, 14, 38, 63, 55, 46, 65, 85])
    print_all([50, 3, 10, 7, 40, 80])

    print_all([49, 40, 93, 24, 34, 83, 74, 12, 89, 69, 100, 25, 75, 67, 45, 88])

    seq = [49, 40, 93, 24, 34, 83, 74, 12, 89, 69, 100, 25, 75, 67, 45, 88, 36, 7, 50, 14, 59, 70, 68, 13, 99, 15, 80,
           97, 39, 90, 82, 43, 48, 71, 47, 16, 96, 23, 9, 17, 29, 1, 22, 5, 77, 86, 63, 42, 20, 10, 95, 53, 35, 57, 55,
           33, 18, 4, 62, 8, 73, 54, 87, 37, 76, 64, 2, 72, 19, 79, 27, 38, 46, 6, 11, 26, 3, 32, 30, 66, 58, 31, 98,
           52, 44, 56, 91, 94, 41, 84, 85, 78, 65, 51, 21, 81, 61, 28, 92, 60]
    print_all(seq)

    # seq = random.sample(range(0, 100000), 100000)
    # print_all(seq)
