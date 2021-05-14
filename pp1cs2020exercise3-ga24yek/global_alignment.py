import numpy as np


class GlobalAlignment:
    def __init__(self, string1, string2, gap_penalty, matrix):
        """
        :param string1: first string to be aligned, string
        :param string2: second string to be aligned, string
        :param gap_penalty: gap penalty, integer
        :param matrix: substitution matrix containing scores for amino acid
                       matches and mismatches, dict

        Attention! string1 is used to index columns, string2 is used to index rows
        """
        # self.string1 = string1
        # self.string2 = string2
        self.gap_penalty = gap_penalty
        self.substitution_matrix = matrix
        self.score_matrix = np.zeros((len(string2) + 1, len(string1) + 1), dtype=np.int)

        self.__alignments = set()
        self.__s1h = ['-'] + list(string1)
        self.__s2v = ['-'] + list(string2)

        self.align()

    def align(self):
        """
        Align given strings using the Needleman-Wunsch algorithm,
        store the alignments and the score matrix used to compute those alignments.
        NB: score matrix and the substitution matrix are different matrices!
        """
        num_rows, num_cols = self.score_matrix.shape
        self.score_matrix[0, :] = np.arange(start=0, stop=num_cols) * self.gap_penalty
        self.score_matrix[:, 0] = np.arange(start=0, stop=num_rows) * self.gap_penalty

        for m in range(1, len(self.__s2v)):
            for n in range(1, len(self.__s1h)):
                d_m_n = self.substitution_matrix[self.__s2v[m]][self.__s1h[n]]
                score_m1_n1 = self.score_matrix[m - 1, n - 1]
                score_m1_n = self.score_matrix[m - 1, n]
                score_m_n1 = self.score_matrix[m, n - 1]

                from_topleft = score_m1_n1 + d_m_n
                from_top = score_m1_n + self.gap_penalty
                from_left = score_m_n1 + self.gap_penalty

                score_m_n = max(from_topleft, max(from_left, from_top))
                self.score_matrix[m, n] = score_m_n

        backtrace_start_cell = (num_rows - 1, num_cols - 1)
        self._backtrace(cur_cell_idx=backtrace_start_cell, str1_aligned='', str2_aligned='')

    def get_best_score(self):
        """
        :return: the highest score for the aligned strings, int
        """
        num_rows, num_cols = self.score_matrix.shape
        return self.score_matrix[num_rows - 1, num_cols - 1]

    def get_number_of_alignments(self):
        """
        :return: number of found alignments with the best score
        """
        return len(self.__alignments)

    def get_alignments(self):
        """
        :return: list of alignments, where each alignment is represented
                 as a tuple of aligned strings
        """
        return list(self.__alignments)

    def get_score_matrix(self):
        """
        :return: matrix built during the alignment process as a list of lists
        """
        return self.score_matrix.tolist()

    def _backtrace(self, cur_cell_idx, str1_aligned, str2_aligned):
        if cur_cell_idx == (0, 0):
            self.__alignments.add((str1_aligned, str2_aligned))
        else:
            m = cur_cell_idx[0]
            n = cur_cell_idx[1]
            if m != 0 and n != 0:
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

    ga = GlobalAlignment("AVNCCEGQHI", "ARNDEQ", -1, mat.MATRICES['identity'])
    print(ga.get_best_score())
    print(ga.get_alignments())
    print(ga.get_score_matrix())
