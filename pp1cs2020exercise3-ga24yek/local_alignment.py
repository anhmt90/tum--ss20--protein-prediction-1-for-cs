import numpy as np


class LocalAlignment:
    def __init__(self, string1, string2, gap_penalty, matrix):
        """
        :param string1: first string to be aligned, string
        :param string2: second string to be aligned, string
        :param gap_penalty: gap penalty, integer
        :param matrix: substitution matrix containing scores for amino acid
                       matches and mismatches, dict

        Attention! string1 is used to index columns, string2 is used to index rows
        """
        self.gap_penalty = gap_penalty
        self.substitution_matrix = matrix
        self.score_matrix = np.zeros((len(string2) + 1, len(string1) + 1), dtype=np.int)

        self.__alignment = ('', '')
        self.__s1h = ['-'] + list(string1)
        self.__s2v = ['-'] + list(string2)

        self.align()

    def align(self):
        """
        Align given strings using the Smith-Waterman algorithm.
        NB: score matrix and the substitution matrix are different matrices!
        """
        num_rows, num_cols = self.score_matrix.shape

        for m in range(1, len(self.__s2v)):
            for n in range(1, len(self.__s1h)):
                d_m_n = self.substitution_matrix[self.__s2v[m]][self.__s1h[n]]
                score_m1_n1 = self.score_matrix[m - 1, n - 1]
                score_m1_n = self.score_matrix[m - 1, n]
                score_m_n1 = self.score_matrix[m, n - 1]

                from_topleft = score_m1_n1 + d_m_n
                from_top = score_m1_n + self.gap_penalty
                from_left = score_m_n1 + self.gap_penalty

                score_m_n = max(0, max(from_topleft, max(from_left, from_top)))
                self.score_matrix[m, n] = score_m_n

        backtrace_start_cell = np.unravel_index(self.score_matrix.argmax(), self.score_matrix.shape)

        self._backtrace(cur_cell_idx=backtrace_start_cell, str1_aligned='', str2_aligned='')

    def has_alignment(self):
        """
        :return: True if a local alignment has been found, False otherwise
        """
        return self.__alignment != ('', '')

    def get_alignment(self):
        """
        :return: alignment represented as a tuple of aligned strings
        """
        return self.__alignment

    def is_residue_aligned(self, string_number, residue_index):
        """
        :param string_number: number of the string (1 for string1, 2 for string2) to check
        :param residue_index: index of the residue to check
        :return: True if the residue with a given index in a given string has been aligned
                 False otherwise
        """
        seq = self.__s1h[1:] if string_number == 1 else self.__s2v[1:]
        seq_aligned = self.__alignment[string_number - 1]

        subseq_aligned = seq_aligned.replace('-', '')
        if not seq_aligned:
            return False
        seq_parts = ''.join(seq).partition(subseq_aligned)
        return len(seq_parts[0]) - 1 < residue_index < len(seq_parts[0]) + len(seq_parts[1])

    def _backtrace(self, cur_cell_idx, str1_aligned, str2_aligned):
        m = cur_cell_idx[0]
        n = cur_cell_idx[1]
        if self.score_matrix[m, n] == 0:
            self.__alignment = (str1_aligned, str2_aligned)
        else:
            score_m_n = self.score_matrix[m, n]
            score_m1_n1 = self.score_matrix[m - 1, n - 1]
            score_m1_n = self.score_matrix[m - 1, n]
            score_m_n1 = self.score_matrix[m, n - 1]

            d_m_n = self.substitution_matrix[self.__s2v[m]][self.__s1h[n]]

            from_topleft = (score_m1_n1 + d_m_n) == score_m_n
            from_top = (score_m1_n + self.gap_penalty) == score_m_n
            from_left = (score_m_n1 + self.gap_penalty) == score_m_n

            if from_topleft:
                __str1_aligned = self.__s1h[n] + str1_aligned
                __str2_aligned = self.__s2v[m] + str2_aligned
                self._backtrace((m - 1, n - 1), str1_aligned=__str1_aligned, str2_aligned=__str2_aligned)
            if from_top:
                __str1_aligned = '-' + str1_aligned
                __str2_aligned = self.__s2v[m] + str2_aligned
                self._backtrace((m - 1, n), str1_aligned=__str1_aligned, str2_aligned=__str2_aligned)
            if from_left:
                __str1_aligned = self.__s1h[n] + str1_aligned
                __str2_aligned = '-' + str2_aligned
                self._backtrace((m, n - 1), str1_aligned=__str1_aligned, str2_aligned=__str2_aligned)


if __name__ == '__main__':
    import tests.matrices as mat

    la = LocalAlignment("ARNDCEQGHI", "DDCEQHG", -6, mat.MATRICES['blosum'])
    print(la.get_alignment())
    print(la.is_residue_aligned(string_number=2, residue_index=6))
