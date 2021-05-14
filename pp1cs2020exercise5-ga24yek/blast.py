from collections import defaultdict, OrderedDict, namedtuple
from dataclasses import dataclass
from itertools import chain
from itertools import product
import itertools
from operator import itemgetter

import numpy as np

from pathlib import Path
import pickle
import json

"""
ATTENTION: Use the following dictionaries to get the correct index for each
           amino acid when accessing any type of matrix (PSSM or substitution
           matrix) parameters. Failure to do so will most likely result in not
           passing the tests.
"""
ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'
AA_TO_INT = {aa: index for index, aa in enumerate(ALPHABET)}
INT_TO_AA = {index: aa for index, aa in enumerate(ALPHABET)}


@dataclass
class HSP:
    # HSP = namedtuple('HSP', ['seq', 'query', 'score', 'q_pos', 'target', 't_pos', 'mat', 'X'])
    seq: str
    score: float
    query: str
    q_pos: int
    target: str
    t_pos: int
    mat: np.ndarray
    X: int


class BlastDB:

    def __init__(self):
        """
        Initialize the BlastDB class.
        """
        self.seq_indexer = defaultdict(set)
        # self.db = OrderedDict()
        self.seq_store = list()
        self.seq2words = defaultdict(set)

    def add_sequence(self, sequence):
        """
        Add a sequence to the database.

        :param sequence: a protein sequence (string).
        """
        # self.db[sequence] = None
        self.seq_store.append(sequence)
        seq_index = len(self.seq_store) - 1
        for i in range(len(sequence) - 2):
            word = sequence[i:i + 3]
            self.seq_indexer[word].add(seq_index)
            self.seq2words[sequence].add(word)

    def get_sequences(self, word):
        """
        Return all sequences in the database containing a given word.

        :param word: a word (string).

        :return: List with sequences.
        """
        ret_seqs = []
        if isinstance(word, str) and len(word) == 3:
            indexes = self.seq_indexer[word]
            if indexes:
                # _db = list(list(zip(*self.db.items()))[0])
                _db = self.seq_store
                ret_seqs = list(map(list(_db).__getitem__, indexes))
                # ret_seqs = list(itemgetter(*indexes)(_db))
        return ret_seqs

    def get_db_stats(self):
        """
        Return some database statistics:
            - Number of sequences in database
            - Number of different words in database
            - Average number of words per sequence (rounded to nearest int)
            - Average number of sequences per word (rounded to nearest int)

        :return: Tuple with four integer numbers corrsponding to the mentioned
                 statistics (in order of listing above).
        """
        seq_count = len(self.seq_store)
        word_count = len(self.seq_indexer.keys())
        index_count = len(list(chain(*self.seq_indexer.values())))

        avg_words_per_seq = int(round(float(index_count) / seq_count))
        avg_seqs_per_word = int(round(float(index_count) / word_count))

        return seq_count, word_count, avg_words_per_seq, avg_seqs_per_word





class Blast:

    def __init__(self, substitution_matrix):
        """
        Initialize the Blast class with the given substitution_matrix.

        :param substitution_matrix: 20x20 amino acid substitution score matrix.
        """
        self.submat = substitution_matrix
        self.db = BlastDB()
        self.words = defaultdict(set)
        self.p_ws = defaultdict(set)

    @staticmethod
    def __get_qualified_words(rows, T, pos):
        row0 = rows[0]
        row1 = rows[1]
        row2 = rows[2]
        words = defaultdict(set)
        for i, score_i in enumerate(row0):
            for j, score_j in enumerate(row1):
                for k, score_k in enumerate(row2):
                    score = score_i + score_j + score_k
                    if score >= T:
                        word = f'{INT_TO_AA[i]}{INT_TO_AA[j]}{INT_TO_AA[k]}'
                        words[word].add((pos, score))
        return words

    def get_words(self, *, sequence=None, pssm=None, T=11):
        """
        Return all words with score >= T for given protein sequence or PSSM.
        Only a sequence or PSSM will be provided, not both at the same time.
        A word may only appear once in the list.

        :param sequence: a protein sequence (string).
        :param pssm: a PSSM (Lx20 matrix, where L is length of sequence).
        :param T: score threshold T for the words.

        :return: List of unique words.
        """
        self.words = defaultdict(set)
        self.p_ws = defaultdict(set)
        seen_subseq = set()
        if sequence and not pssm:
            for a in range(len(sequence) - 2):
                a0 = sequence[a]
                a1 = sequence[a + 1]
                a2 = sequence[a + 2]
                subseq = sequence[a: a + 3]
                # if seen_subseq.intersection({''.join(per) for per in set(itertools.permutations(list(subseq)))}):
                #     continue

                # if subseq in seen_subseq:
                #     continue

                seen_subseq.add(subseq)

                rows = self.submat[[AA_TO_INT[a0], AA_TO_INT[a1], AA_TO_INT[a2]], :]
                # self.words.update(self.__get_qualified_words(rows, T, a))
                new_words = self.__get_qualified_words(rows, T, a)
                for word, pos_score_set in new_words.items():
                    self.words[word].update(pos_score_set)
                    for pos, score in pos_score_set:
                        self.p_ws[pos].add((word, score))

        elif not sequence and pssm is not None:
            row_count = pssm.shape[0]
            for a in range(row_count - 2):
                rows = pssm[a:a + 3, :]
                new_words = self.__get_qualified_words(rows, T, a)
                for word, pos_score_set in new_words.items():
                    self.words[word].update(pos_score_set)
                    for pos, score in pos_score_set:
                        self.p_ws[pos].add((word, score))

        return list(self.words.keys())

    @staticmethod
    def __find_all(string, substr, step=1):
        start = 0
        while True:
            start = string.find(substr, start)
            if start == -1:
                return
            yield start
            start += step

    def search_one_hit(self, blast_db, *, query=None, pssm=None, T=13, X=5, S=30):
        """
        Search a database for target sequences with a given query sequence or
        PSSM. Return a dictionary where the keys are the target sequences for
        which HSPs have been found and the corresponding values are lists of
        tuples. Each tuple is a HSP with the following elements (and order):
            - Start position of HSP in query sequence
            - Start position of HSP in target sequence
            - Length of the HSP
            - Total score of the HSP
        The same HSP may not appear twice in the list (remove duplictes).
        Only a sequence or PSSM will be provided, not both at the same time.

        :param blast_db: BlastDB class object with protein sequences.
        :param query: query protein sequence.
        :param pssm: query PSSM (Lx20 matrix, where L is length of sequence).
        :param T: score threshold T for the words.
        :param X: drop-off threshold X during extension.
        :param S: score threshold S for the HSP.

        :return: dictionary of target sequences and list of HSP tuples.
        """

        HSPs = defaultdict(set)
        query_words = self.get_words(sequence=query, T=T) if query else self.get_words(pssm=pssm, T=T)
        mat = self.submat if query else pssm

        for word in query_words:
            targets = blast_db.get_sequences(word)
            for q_pos, score in self.words[word]:
                for target in targets:
                    target_positions = list(self.__find_all(target, word, 1))
                    for t_pos in target_positions:
                        hsp = HSP(seq=word, score=score, query=query, q_pos=q_pos, target=target, t_pos=t_pos,
                                  mat=mat, X=X)
                        hsp = self.__rexpand(hsp)
                        hsp = self.__lexpand(hsp)
                        if hsp.score >= S:
                            HSPs[hsp.target].add((hsp.q_pos, hsp.t_pos, len(hsp.seq), hsp.score))
        return HSPs

    def __lexpand(self, p):
        i = -1
        best_score, best_seq = p.score, p.seq
        best_q_pos, best_t_pos = -1, -1
        while p.q_pos + i >= 0 and p.t_pos + i >= 0:
            best_seq_copied = best_seq
            p, best_score, best_seq = self.__update_cur_best(p, best_score, best_seq, i, left=True)
            if best_seq_copied != best_seq:
                best_q_pos = p.q_pos + i
                best_t_pos = p.t_pos + i

            if p.score <= best_score - p.X:
                break
            i -= 1

        p.seq = best_seq
        p.score = best_score
        if best_q_pos != -1 and best_t_pos != -1:
            p.q_pos, p.t_pos = best_q_pos, best_t_pos

        return p

    def __rexpand(self, p):
        best_score, best_seq = p.score, p.seq
        i = 3
        right_end = len(p.query) if p.query else p.mat.shape[0]
        while p.q_pos + i < right_end and p.t_pos + i < len(p.target):
            p, best_score, best_seq = self.__update_cur_best(p, best_score, best_seq, i)

            i += 1
            if p.score <= best_score - p.X:
                break

        p.score = best_score
        p.seq = best_seq
        return p

    @staticmethod
    def __update_cur_best(p, best_score, best_seq, i, left=False):
        _best_score = best_score
        _best_seq = best_seq

        row_index = AA_TO_INT[p.query[p.q_pos + i]] if p.query else p.q_pos + i
        cur_t_residue = p.target[p.t_pos + i]

        p.score += p.mat[row_index][AA_TO_INT[cur_t_residue]]
        p.seq = str(p.seq) + str(cur_t_residue) if not left else str(cur_t_residue) + str(p.seq)

        if (p.score > best_score) or (p.score == best_score and len(p.seq) < len(best_seq)):
            _best_score, _best_seq = p.score, p.seq

        return p, _best_score, _best_seq

    def search_two_hit(self, blast_db, *, query=None, pssm=None, T=11, X=5, S=30, A=40):
        """
        Search a database for target sequences with a given query sequence or
        PSSM. Return a dictionary where the keys are the target sequences for
        which HSPs have been found and the corresponding values are lists of
        tuples. Each tuple is a HSP with the following elements (and order):
            - Start position of HSP in query sequence
            - Start position of HSP in target sequence
            - Length of the HSP
            - Total score of the HSP
        The same HSP may not appear twice in the list (remove duplictes).
        Only a sequence or PSSM will be provided, not both at the same time.

        :param blast_db: BlastDB class object with protein sequences.
        :param query: query protein sequence.
        :param pssm: query PSSM (Lx20 matrix, where L is length of sequence).
        :param T: score threshold T for the words.
        :param X: drop-off threshold X during extension.
        :param S: score threshold S for the HSP.
        :param A: max distance A between two hits for the two-hit method.

        :return: dictionary of target sequences and list of HSP tuples.
        """
        import bisect
        import timeit

        HSPs = defaultdict(set)
        right_end = len(query) if query else pssm.shape[0]
        query_words = set(self.get_words(sequence=query, T=T) if query else self.get_words(pssm=pssm, T=T))
        mat = self.submat if query else pssm

        start = timeit.default_timer()

        all_target_indices = set()
        for word in query_words:
            all_target_indices.update(blast_db.seq_indexer[word])

        # hit_db = defaultdict(lambda : defaultdict(list))
        hsp_tracker = defaultdict(lambda: defaultdict(list))
        for i in all_target_indices:
            target = blast_db.seq_store[i]
            diags_hits_tracker = defaultdict(list)
            hit_words = query_words.intersection(blast_db.seq2words[target])
            if len(hit_words) >= 2:
                psw_tups = list()
                for word in hit_words:
                    for q_pos, score in self.words[word]:
                        bisect.insort(psw_tups, (q_pos, score, word))

                for q_pos, score, word in psw_tups:
                    hit_poses = self.__find_all(target, word)
                    for t_pos in hit_poses:
                        diag = t_pos - q_pos
                        rhit = (q_pos, t_pos, word, score)

                        # used_hits = defaultdict(list)
                        hsp_added = False
                        lh_to_remove = -1
                        if diags_hits_tracker[diag]:
                            for lh_i, lhit in enumerate(diags_hits_tracker[diag]):  #CAN DO BIN SEARCH TO HALVE THE LIST
                                if hsp_added:
                                    break

                                if rhit[0] - lhit[0] >= 3 and rhit[1] - lhit[1] >= 3:  # L < R
                                    distance = rhit[0] - lhit[0]
                                    if abs(distance) <= A and distance == rhit[1] - lhit[1]:
                                        is_included = False
                                        for hsp in hsp_tracker[target][diag]:
                                            if hsp[1] <= lhit[1] < rhit[1] < hsp[1] + hsp[2]:
                                                is_included = True
                                                break
                                        if not is_included:
                                            hsp = HSP(seq=word, score=score, query=query, q_pos=q_pos, target=target,
                                                      t_pos=t_pos, mat=mat, X=X)
                                            hsp = self.__lexpand(hsp)
                                            if hsp.t_pos <= lhit[1] + 2:  # check if part of lhit included on target seq
                                                hsp_q_pos = hsp.q_pos
                                                hsp_t_pos = hsp.t_pos

                                                hsp.q_pos = rhit[0]
                                                hsp.t_pos = rhit[1]

                                                hsp = self.__rexpand(hsp)

                                                hsp_tracker[target][diag].append((hsp_q_pos, hsp_t_pos, len(hsp.seq)))
                                                if hsp.score >= S:
                                                    HSPs[target].add((hsp_q_pos, hsp_t_pos, len(hsp.seq), hsp.score))
                                                    lh_to_remove = lh_i
                                                    hsp_added = True
                        if not hsp_added:
                            diags_hits_tracker[diag].append(rhit)
                        else:
                            del diags_hits_tracker[diag][lh_to_remove]

        # for target, hits in hit_db.items():

        stop = timeit.default_timer()
        duration = stop - start
        print('Time: ', duration)
        print()

        # seen_target_indexes = set()
        # for word in query_words:
        #     targets = blast_db.get_sequences(word)
        #     for target in targets:
        #         self.__find_all(target, word)
        #         for q_pos, score in self.words[word]:
        #             hits[target].append((word, q_pos, score))


        diagonals = defaultdict(dict)

        # for l_pos, l_word_score_tups in self.p_ws.items():
        #     for r_pos in range(l_pos + 3, min(l_pos + A, right_end - 2)):
        #         distance = r_pos - l_pos
        #         r_word_score_tups = self.p_ws.get(r_pos)
        #         if not r_word_score_tups:
        #             continue
        #
        #         for l_word, l_score in l_word_score_tups:
        #             for r_word, r_score in r_word_score_tups:
        #                 l_targets = set(blast_db.get_sequences(l_word))
        #                 r_targets = set(blast_db.get_sequences(r_word))
        #
        #                 targets = l_targets.intersection(r_targets)
        #                 valid_targets = defaultdict(set)
        #                 for target in targets:
        #                     l_pos_t = set(self.__find_all(target, l_word, 1))
        #                     r_pos_t = set(self.__find_all(target, r_word, 1))
        #
        #                     for l in l_pos_t:
        #                         r = distance + l
        #                         if r in r_pos_t:
        #                             is_included = False
        #                             diag_dict = diagonals.get(target)
        #                             if diag_dict:
        #                                 if (r_pos - r) in set(diag_dict.keys()): # same diagonal
        #                                     hsp = diag_dict[r_pos - r]
        #                                     if hsp[1] <= l < hsp[1] + hsp[2] or hsp[1] <= r < hsp[1] + hsp[2]:
        #                                         is_included = True
        #
        #                             if not is_included:
        #                                 valid_targets[target].add((l, r))
        #                 ## DEBUG
        #                 # if len(valid_targets) != 0:
        #                 #     print(f'valid_targets: {len(valid_targets)}')
        #
        #                 for target, pos_set in valid_targets.items():
        #                     for l, r in pos_set:
        #                         q_pos = r_pos
        #                         t_pos = r
        #                         hsp = HSP(seq=r_word, score=r_score, query=query, q_pos=q_pos, target=target,
        #                                   t_pos=t_pos, mat=mat, X=X)
        #                         hsp = self.__lexpand(hsp)
        #                         if hsp.t_pos <= l + 2:
        #                             hsp_q_pos = hsp.q_pos
        #                             hsp_t_pos = hsp.t_pos
        #
        #                             hsp.q_pos = r_pos
        #                             hsp.t_pos = r
        #
        #                             hsp = self.__rexpand(hsp)
        #                             if hsp.score >= S:
        #                                 tup = (hsp_q_pos, hsp_t_pos, len(hsp.seq))
        #                                 diagonals[hsp.target][r_pos - r] = tup
        #                                 HSPs[hsp.target].add((hsp_q_pos, hsp_t_pos, len(hsp.seq), hsp.score))
        return HSPs


if __name__ == '__main__':
    _blast_db = BlastDB()

    # test_json = 'blast_test_debug.json'
    test_json = 'blast_test_full.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    for s in json_data['db_sequences']:
        _blast_db.add_sequence(s)

    sub_matrix = np.array(json_data['sub_matrix'], dtype=np.int64)
    _pssm = np.array(json_data['query_pssm'], dtype=np.int64)
    _blast = Blast(sub_matrix)
    _query_seq = json_data['query_seq']

    load = False
    if load:
        with open('blast_db.pkl', 'rb') as file:
            _blast_db = pickle.load(file)

    results = None
    blast_results = None

    # test_search_one_hit = False
    # if test_search_one_hit:
    #     blast_results = json_data['blast_hsp_one_hit']
    #     for key, value in blast_results.items():
    #         blast_results[key] = [tuple(x) for x in value]
    #
    #     results = blast.search_one_hit(blast_db, query=_query_seq, T=13, X=5, S=30)
    #################################################################################################
    # test_search_one_hit_pssm = True
    # if test_search_one_hit_pssm:
    #     blast_results = json_data['blast_hsp_one_hit_pssm']
    #     for key, value in blast_results.items():
    #         blast_results[key] = [tuple(x) for x in value]
    #
    #     results = _blast.search_one_hit(_blast_db, pssm=_pssm, T=13, X=5, S=30)

    test_search_2_hits = False
    if test_search_2_hits:
        from time import gmtime, strftime
        blast_results = json_data['blast_hsp_two_hit']
        for key, value in blast_results.items():
            blast_results[key] = [tuple(x) for x in value]
        results = _blast.search_two_hit(_blast_db, query=_query_seq, T=11, X=5, S=30, A=40)

        if 'full' in test_json:
            with open(f'results_{strftime("%Y_%m_%d_%H_%M_%S", gmtime())}.pkl', 'ab') as file:
                pickle.dump(results, file)

    else:
        blast_results = json_data['blast_hsp_two_hit']
        for key, value in blast_results.items():
            blast_results[key] = [tuple(x) for x in value]

        with open(f'results_2020_06_24_23_38_43.pkl', 'rb') as file:
            results = pickle.load(file)

        faulty = {}
        for target in results.keys():
            bres = blast_results.get(target)
            if not bres:
                faulty[target] = None
                continue
            gt_set = set(bres)
            my_set = set(results[target])
            if gt_set != my_set:
                faulty[target] = gt_set.difference(my_set)
        print(faulty)

    assert results is not None
    assert (len(blast_results) == len(results))
    assert (set(blast_results) == set(results))
    print('OK')
