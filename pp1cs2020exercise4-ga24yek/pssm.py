import numpy as np

"""
ATTENTION: Use the following dictionaries to get the correct index for each
           amino acid when accessing any type of matrix or array provided as
           parameters. Further, use those indices when generating or returning
           any matrices or arrays. Failure to do so will most likely result in
           not passing the tests.
EXAMPLE: To access the substitution frequency from alanine 'A' to proline 'P'
         in the bg_matrix use bg_matrix[AA_TO_INT['A'], AA_TO_INT['P']].
"""
ALPHABET = 'ACDEFGHIKLMNPQRSTVWY-'
AA_TO_INT = {aa: index for index, aa in enumerate(ALPHABET)}
INT_TO_AA = {index: aa for index, aa in enumerate(ALPHABET)}
GAP_INDEX = AA_TO_INT['-']


def _check_empty(sequences):
    if isinstance(sequences, list) and len(sequences) == 0:
        raise TypeError('No sequence provided')


def _check_same_length(sequences):
    seq_iter = iter(sequences)
    _len = len(next(seq_iter))
    if not all(len(seq) == _len for seq in seq_iter):
        raise TypeError(f'Sequences are not the same length.')


def _check_valid(sequences):
    all_seqs = ''.join(sequences)
    for letter in all_seqs:
        if letter not in ALPHABET:
            raise TypeError(f'Invalid letter in sequence ')


class MSA:
    def __init__(self, sequences):
        """
        Initialize the MSA class with the provided list of sequences. Check the
        sequences for correctness. Pre-calculate any statistics you seem fit.

        :param sequences: List containing the MSA sequences.
        """
        _check_empty(sequences)
        _check_same_length(sequences)
        _check_valid(sequences)

        self.seqs = sequences
        self.MSA_len = len(self.seqs[0])
        self.MSA_mat = np.array([list(seq) for seq in sequences])

        self.weights = self.get_sequence_weights()

    def get_pssm(self, *, bg_matrix=None, beta=10, use_sequence_weights=False,
                 redistribute_gaps=False, add_pseudocounts=False):
        """
        Return a PSSM for the underlying MSA. Use the appropriate refinements
        according to the parameters. If no bg_matrix is specified, use uniform
        background and pair frequencies.
        Every row in the resulting PSSM corresponds to a non-gap position in
        the primary sequence of the MSA (i.e. the first one).
        Every column in the PSSM corresponds to one of the 20 amino acids.
        Values that would be -inf must be replaced by -20 in the final PSSM.
        Before casting to dtype=numpy.int64, round all values to the nearest
        integer (do not just FLOOR all values).

        :param bg_matrix: Amino acid pair frequencies as numpy array (20, 20).
                          Access the matrix using the indices from AA_TO_INT.
        :param beta: Beta value (float) used to weight the pseudocounts
                     against the observed amino acids in the MSA.
        :param use_sequence_weights: Calculate and apply sequence weights.
        :param redistribute_gaps: Redistribute the gaps according to the
                                  background frequencies.
        :param add_pseudocounts: Calculate and add pseudocounts according
                                 to the background frequencies.

        :return: PSSM as numpy array of shape (L x 20, dtype=numpy.int64).
                 L = ungapped length of the primary sequence.
        """
        # pssm = np.zeros((self.MSA_len, 20))
        pssm = np.zeros((self.MSA_len, 21), dtype=np.float64)
        if bg_matrix is None:
            bg_matrix = np.empty((20, 20))
            bg_matrix.fill(1.0 / (20 * 20))

        pssm = self.__count_observations(pssm, use_sequence_weights)

        if redistribute_gaps:
            pssm = self.__redistribute_gaps(pssm, bg_matrix)

        pssm = pssm[:, :20]
        if add_pseudocounts:
            pssm = self.__compute_adjusted_freqs(pssm, bg_matrix, beta)

        pssm = self.__normalize(pssm)
        pssm = self.__divide_by_bg_freq(pssm, bg_matrix)
        pssm = self.__log2(pssm)
        pssm = self.__remove_gap_rows(pssm)
        pssm[pssm < -20] = -20

        return np.rint(pssm).astype(np.int64)

    def __remove_gap_rows(self, pssm):
        primary_seq = self.MSA_mat[0]
        target_row_idxes = np.where(primary_seq == '-')[0]
        return np.delete(pssm, target_row_idxes, axis=0)

    @staticmethod
    def __log2(pssm):
        return 2 * np.log2(pssm)

    @staticmethod
    def __divide_by_bg_freq(pssm, bg_mat):
        p = np.sum(bg_mat, axis=0)
        return pssm[:, :20] / p.reshape((1, -1))

    @staticmethod
    def __normalize(pssm):
        row_sums = pssm.sum(axis=1)
        return pssm / row_sums[:, np.newaxis]

    def __compute_adjusted_freqs(self, pssm, bg_mat, beta):
        pseudocount_mat = self.__compute_pseudocounts(pssm, bg_mat)
        alpha = self.get_number_of_observations() - 1
        return (alpha * pssm + float(beta) * pseudocount_mat) / (alpha + beta)

    def __compute_pseudocounts(self, pssm, bg_mat):
        p = np.sum(bg_mat, axis=0)
        pseudocount_mat = np.zeros(pssm.shape)
        for i in range(pssm.shape[0]):
            for j in range(pssm.shape[1]):
                if pssm[i, j] != 0.0:
                    for a in range(20):
                        pseudocount_mat[i, a] += bg_mat[j, a] * (pssm[i, j] / p[j])
        return pseudocount_mat

    @staticmethod
    def __redistribute_gaps(pssm, bg_mat):
        expected_bg_freq = np.sum(bg_mat, axis=0)
        pssm[:, :20] = pssm[:, :20] + pssm[:, 20].reshape(-1, 1) * expected_bg_freq.reshape(1, -1)
        return pssm

    def __count_observations(self, pssm, use_sequence_weights):
        if not use_sequence_weights:
            for i in range(self.MSA_len):
                unique_AAs, counts = np.unique(self.MSA_mat[:, i], return_counts=True)
                for j, aa in enumerate(unique_AAs):
                    pssm[i, AA_TO_INT[aa]] = counts[j]
        else:
            seq_weights = self.get_sequence_weights()
            for i in range(self.MSA_len):
                unique_AAs, seq_indices = np.unique(self.MSA_mat[:, i], return_inverse=True)
                for j, aa in enumerate(unique_AAs):
                    weight_indices = np.where(seq_indices == j)[0]
                    pssm[i, AA_TO_INT[aa]] = np.sum(seq_weights[weight_indices])
        return pssm

    def get_size(self):
        """
        Return the number of sequences in the MSA and the MSA length, i.e.
        the number of columns in the MSA. This includes gaps.

        :return: Tuple of two integers. First element is the number of
                 sequences in the MSA, second element is the MSA length.
        """
        return len(self.seqs), self.MSA_len

    def get_primary_sequence(self):
        """
        Return the primary sequence of the MSA. In this exercise, the primary
        sequence is always the first sequence of the MSA. The returned
        sequence must NOT include gap characters.

        :return: String containing the ungapped primary sequence.
        """
        primary_seq = self.seqs[0].replace('-', '')
        return primary_seq

    def get_sequence_weights(self):
        """
        Return the calculated sequence weights for all sequences in the MSA.
        The order of weights in the array must be equal to the order of the
        sequences in the MSA.

        :return: Numpy array (dtype=numpy.float64) containing the weights for
                 all sequences in the MSA.
        """

        weights = np.sum(self.__calc_weight_matrix(), axis=1)
        return weights.astype(np.float64)

    def get_number_of_observations(self):
        """
        Return the estimated number of independent observations in the MSA.

        :return: Estimate of independent observation (dtype=numpy.float64).
        """
        L = self.MSA_len
        r = np.zeros(L, dtype=np.float64)
        for i in range(L):
            r[i] = len(np.unique(self.MSA_mat[:, i]))

        num_obs = np.sum(r) / L
        return num_obs.astype(np.float64)

    def __calc_weight_matrix(self):
        weight_mat = np.zeros(self.get_size(), dtype=np.float64)
        r = np.zeros(self.MSA_len)
        for i in range(self.MSA_len):
            unique_AAs, seq_indices = np.unique(self.MSA_mat[:, i], return_inverse=True)
            r[i] = len(unique_AAs)
            if r[i] > 1:
                for u in range(len(unique_AAs)):
                    indices = np.where(seq_indices == u)
                    s = len(indices[0])
                    weight_mat[indices, i] = 1.0 / (r[i] * s)

        return weight_mat


if __name__ == '__main__':
    import tests.pssm_test as dat

    my_seqs = ['SE-AN', 'SE-ES', 'SEVEN', 'SE-AS']
    my_bg_matrix = np.array(dat.get_bg())
    msa = MSA(sequences=my_seqs)
    my_pssm = msa.get_pssm(bg_matrix=my_bg_matrix, redistribute_gaps=True, add_pseudocounts=True,
                           use_sequence_weights=True)
    print(my_pssm)
