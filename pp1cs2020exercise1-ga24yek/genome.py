from itertools import permutations
import itertools
from collections import defaultdict

inversed_genetic_code_table = {
    'F': ['TTT', 'TTC'],
    'L': ['TTA', 'TTG'] + ['CTT', 'CTA', 'CTC', 'CTG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG'] + ['AGT', 'AGC'],
    'Y': ['TAT', 'TAC'],
    'C': ['TGT', 'TGC'],
    'W': ['TGG'],
    'P': ['CCT', 'CCA', 'CCC', 'CCG'],
    'H': ['CAT', 'CAC'],
    'Q': ['CAA', 'CAG'],
    'R': ['CGA', 'CGT', 'CGC', 'CGG'] + ['AGG', 'AGA'],
    'I': ['ATT', 'ATC', 'ATA'],
    'M': ['ATG'],
    'T': ['ACA', 'ACT', 'ACC', 'ACG'],
    'N': ['AAC', 'AAT'],
    'K': ['AAG', 'AAA'],
    'V': ['GTA', 'GTT', 'GTC', 'GTG'],
    'A': ['GCA', 'GCT', 'GCC', 'GCG'],
    'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'],
    'G': ['GGA', 'GGT', 'GGC', 'GGG'],
}

base_complements = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}


class Genome:

    def __init__(self, genome):
        """
        Initialize the Genome class with the provided genome sequence.

        :param genome: String with the genome sequence.
        """
        self.genome = genome.upper()

        g_len = len(self.genome)
        self.base_content = {
            'A': self.genome.count('A') / g_len,
            'T': self.genome.count('T') / g_len,
            'G': self.genome.count('G') / g_len,
            'C': self.genome.count('C') / g_len
        }

    def get_gc_content(self):
        """
        Return the GC content of the genome sequence, i.e. the combined
        fraction of 'G' and 'C' in the entire genome sequence.

        :return: GC content (float, rounded to 6 digits)
        """
        return round(self.base_content['G'] + self.base_content['C'], 6)

    def get_codon_dist(self):
        """
        Return the expected codon distribution (fractions) based on the
        distribution (fractions) of the four different nucleotides (ATGC).

        :return: Tree-like structure made out of nested dictionaries. The nodes
                 represent the different nucleotides and the path down the tree
                 forms the corresponding codons. The leafs contain the expected
                 codon frequencies (rounded to 6 digits).
        """
        all_codons = dict(
            zip([''.join(codon) for codon in list(itertools.product(list('ATCG'), repeat=3))], [0] * (4 * 4 * 4)))

        callable_dict = lambda: defaultdict(callable_dict)
        codon_dist = callable_dict()
        for codon, freq in all_codons.items():
            freq = round(self.base_content[codon[0]] * self.base_content[codon[1]] * self.base_content[codon[2]], 6)
            codon_dist[codon[0]][codon[1]][codon[2]] = freq

        return convert_from_defaultdict(codon_dist)

    def get_amino_acid_dist(self):
        """
        Return the expected amino acid distribution (fractions) based on the
        expected distribution (fractions) of the different codons.

        :return: Dictionary that contains the expected amino acid distribution.
                 The keys are the 20 different amino acids, the values are the
                 corresponding frequencies (rounded to 6 digits).
        """
        all_codons = dict(
            zip([''.join(codon) for codon in list(itertools.product(list('ATCG'), repeat=3))], [0] * (4 * 4 * 4)))

        for i in range(len(self.genome) - 2):
            all_codons[''.join([self.genome[i], self.genome[i + 1], self.genome[i + 2]])] += 1

        condon_dist = self.get_codon_dist()
        aa_dist = {}
        for aa, codons in inversed_genetic_code_table.items():
            aa_freq = 0
            for c in codons:
                cdist = condon_dist[c[0]][c[1]][c[2]]
                ccount = all_codons[f'{c[0]}{c[1]}{c[2]}']
                aa_freq = aa_freq + cdist
            aa_dist[aa] = aa_freq

        s = sum(aa_dist.values())
        for k, v in aa_dist.items():
            aa_dist[k] = round(v / s, 6)

        return aa_dist


def convert_from_defaultdict(d):
    if isinstance(d, defaultdict):
        d = {k: convert_from_defaultdict(v) for k, v in d.items()}
    return d


if __name__ == "__main__":
    seq = 'TTATTTAAGAGCAATGGCCAACAAGTAAAAACGGTTAGCAGGGTTAGGGATATGTTTGTTGACTCTAAAGAAGAGTACGCAAAGCACTATCAAGAAAAGTATTACAATGAGTATTGTCCGTTTTACAGAGATCTCCCGCCTGCATGGGTAGCCATTGAGTTAATGACTTTCGGCAACGTAGTGAAGTTAATTCAAAACATCAGTGATGATAAAATTCAATCACTTAAGATGGATAGATTTTCTAAGAAGTTCAATATTCAGAAATTTCAGACATTAATTAGTTGGATGAATGTGCTGCACCAGATGAGAAACTACTGTGGGCATCATAACCGACTGTTTAATCGAAACTTCCCTGCTCCAACAGCGATTAAAAAGAGCTTGTCTGATGCAATTCCTCTTGTCAGGACCAAACCAAATCCAGATAAGCGTGAAGAGGATCAGTTAAACCGACTTTATACAGCTCTTGCTGCATTACAATGTATATATTCAGGGCTTGGTTTCGATGAAAAAATAGGACCAAAAATCTCTGATTTATTTGATAATTATACAACAACACAGCGGTTTAGCTTATCAATGGGTTTCCCTAATGGTTGGAAAGAAGAGCCGCTTTTTTTTGATTTATAATCGCTATACTTAACATAAAACCCGTTCTACGAATCGTAAAAGGTCGCCTATTTAGGTGGCCTTTTTTTATGAAAACTGTTTTTTAAGTGATAAAAATGTCGTTTGTCCTAACGAAATAGTATTTTTAACAGGCTTATTCAGGCGTTTTGATTTTTAGCTATCTGGTAATATTTAATTATTTTTTATTGGAGTTTGTATATGCCGCCACCTAAAGGTCGAAGTTTTCCATTTGCACCGCGTCACTCTGCTGATTGGTTAGTGAGTCATGTAACGTATGACCAAGCCGTTGATATGTTTTTTAATCAGACCGCAACACAACAGAATTTGGGGCATGATCCTTTAGTTTATTCAAAAGTATTTAGAGGTGTGACATACGCAACGCTCCAAGAAGCTCAACAAGTGTTTACTGAAACGATGAATGCAGAGTATGAAGTGAGAGAGCAAAGGGATTTAGCAGATGAAAATCGAGGGTGTAAGGCTTCTAAAATTCTCAATGATAATATTCGGAATCGTATTGTACCAACAGAGGATGGTTTAGAGTCATTAACTGCATATCGAAAGCCTTTTGGTGAAAAGTGTGTATCTCCTCTTTCTTTGTTTACCAAATCGTTAAATGGCGGTTCAAACAGTCATATTCAAGCAAACGAAACTCAAGCTCAATTAGTTACGCAACAAGCTTATGATTTCCCACTCGTCACTAAAGGTTCTAAAGCGATTCAGGCATTAGATGATTCTTGGGGTTTACAGTGCTCTTTAGAGGATGAATTAAACCAAAATCACCAAGATGTTGTTATTTTAGATTTGGAAGACTCACATAAATTATTAAATAATTTTTGGGTGGATTTAGGCAAAGATACAGCGTTAATAGCAACAAGTTTAAATACGGCTAACGGTTGGTTTCAACATACAACACCGATATTTGATGCTAAAGGTTTAGTTAAGCAATTCGGCGATATAAATATTAAAGCGGATATTGTAGAGAGTAAGGGGAAACAGTTTATTGCTTTCTCTGGTAAGAAGAACGGCAAAGAAATTTAACACGCACTTGTTAATGGCACTCGAATTAACATGAACGGTAAGAAGTACCCAATTAACAGCCCTAAAGTTCAACAAGTCGGTTTGAGTCCTAAAGCACGAGCAAATGGATTTAAAGGCGCAGGTGTATTGACGTTTGTGGTATCGGCAGCGATAGCAACGACAGATCTTGTTTTTAAAGACGATTATCACTTAGTTGATTGGTTCGGTAATGTTGGGGCAGATATGTTTAAGGCATTGCTACAGTTTGGTGCAGGAGAGGCTATATTGTTTGGAATTATTGCGATGACTGGTTACGTTACTCTGGGCTTGATTGCTGTGTTTTTTGTCTATGTGTCTATTGAGTGGATATGGAGTGAATACAAGGTAAATGATGAGGTAGTTAAGGGATTAGAAAGTGTTATCAGTTAAAATGAGGTTTTTAGGAATAGTATTTATACCTTTATTGATTTTACTTTGGTGGATACCTACAAATGGAGTCTTGGGGGATTACCAAGATTTATTAAATCAGACAGATGAAGTTAGATTATCACTGATAACATTAATCGTCCCCATTGGTGGTTTTATCCCTTTACTTATTTCCGTTGCTTTAATTTCTATTGCTATTTATTCGGGGAAACAAGCTCGTTTGGTCATTGGTGAAAAATGGACTAGTACTATCAATAGAACGTGTATTTATAGCATGATATTAGGTGTAGTATTTGCGGTTGTATTTGCACTATATTGTATAAAATTGCTTGATGAAAACGGTTATGAATACAGTTATAACCTAACCCAAATAACCCCGACAGGTATTCATTTAATGTATGTAAAATCACACAAAGAATAGAAACAAAAAAGCCACTCAAAAGAGTGGCTTTTTTGTCTATGAACCAGAGGATCATAAAATTGATCATTGAGATCTTGAGCATTTAGCCGTAACCTTATGCCCAGATGTGAAAAAACCAAGCCTGCAAGCTTGGTTTTTTCGAAATACTTAACTACGTTTAGACGGCGTTGTAAGTCCTTTCATGATAGATCATATCACATAAGAACTGCAAGTATTTCCAAGCAATCAGTGGTTGGGTTTAGAGTCACTAGGGTTGTGGGACAGCATCATGCGGTTGGTTTTCATACCATATTTCGCCCAGATACTGGACGAGAGACTAGCGTAACATCCGATAACGCACACCAAATGAGTGGACTGAAGCAAGAACTAGCGTAAGCAGCTAGGCAAAGAGGCTTGTAGATTGGGGAAACCCTTGCAAGGTACTGATGCAGAGGTCGGCAGAGAACGGTGAGGATGGGGGCAATACAGCCCGTAATCTCTGCACTTGAATAGTGTATTAGAATTCAAGCGAATGATGGCGATGAACGGCTCCAGTCGGCAAGCAATCAGGGTTCACACTAGAGTGATAAACTTAGCGAAAGCGGTTATCCTCTAGGGTGAGTATTACCGAAATCTAGCTCAAAACTCACCTATCAGCATACATGCCGACTTCTTTTGCTCTTCAATTCAAAAAGAATTAAAAAGAAAATATTAAGAGCTGCACGATTAGAAATGAAAGTTGCCGCAACGAGCATAAGCGAGTTAAGCCCACTTTCTTTTCTAATGTGCGAACAAGCGAAGCGCGTTAGTTGGGTCGATTATCTTCTTTTATGGATCGATGCTGACCAACTGAATGTGCTTTTTTCTTCTTTCGGTAACTTTTTCAATTCCGCTCTCGCTGCCCTTGCTATGGTGAGTCGTTTTAATTTTTCTCTGTTCCTTTTGAAAATCGTTTCAAAGGCTTGGTCTCTATAATCTTCACCGTTGCTGCTGCAGAGTGAATCTATAAGAGCAAAGTCTCTTGGCTTGAGTCCGTATACTTTTTCATAAATTGCTTTGGCCGCTTCTTTTCGTGCCACTCCTTTTTCTATTTGTAGTTCAATAAGAGTGCGGTCTTCTTGAGTAAATCGTCTTTTCTTCATTCTCTTTTTTGCTCTTAGGTTTCAAAGGACATCCTTTGGCGCGGCATGTTTTTCCAAATACGTATGTAGAGGGTGAAAACGTGCCCGTCGTGCCTATCTAAAAAATAGATATCACTGAGTAATCAATAATAAAATGCTAACTCCTTGTATTAATACGTTAATTTATGGGTAACGTTATTTTTTAAAGAAACTATCCAAAAGCGTCTTTTTTCCTTCTTTTTTCGTTGGTGTTAATACCGGCTCTTTATGCTCTCTTATCATGGCTTCTAGTTGGTCAAAACGAGTGGCGATGTAATGTAAGTCATCATCAAGCATTTCTCTTTCGTAGGCGTAATACGATTTATCAATTGAGTGATGCACTCTATCGCCATGAGTTTGATACAAATCCCTTTTATAGGCGAACAAACCTCTTACTAGCAACGCCCTTAACAGTTCTGCTTTATCTTTATGCGTGTCGTTTAACAGCTCTTCAATTTCAGCATTCACTCGACCATCCACTCGAACGGACACCAAATGAGTCTTGGATTTTGGTTTATTCATTTTCTAGTTTTCTCTTTTGCCTAATTCGATAAAGAACCGTCTTCAAATAGTCTCTCGATACCGTGATATTTTTTGTCTTTAAAAACGCAACGATGTGCTCTAGTGAGATACCTTCACTCATCAAGATGGTTATCACGTCATATAACTCCGTAATGACGGCGGTTTTAGATTTCGTTTGATGCTTTTTCGTGTACGAAAGAAGCTCTTTTTCTAACGGTGTTTTCTTTTTCACTGAATACTTTATTGATTACTGATTGATACTTTAGAGACTAGCAATAAAAAAGAAGAACGAGAATGGTTATTTCAGTTGGGCTGAAATTGGCATGATTTATCTTAGTAACCTACCTATAAACTCATTGATTCAATGGGTTAGCTTTTTAATGTTAGTTACGTTGTGATTACCCGTTGATTTCTATTTTCAAGATAGGCACGACGGGGACGTTTTCACCCCCTGCTACGCATTGGATAAAAACATCCCGTGCCAAAGGTGAAACCTTTGAAACCCTGACGCACTTCGTTTGTCTCCAGTACGCCCAAAAAAGCCAGCTACGCTTCGGGCTTTTTCGAGCGCACCGGCTAATTGTTGAAAATAGACGTTACAGGGTGATTATTGATTGAAACGGATTGTGATTTTATTGCCCTGTAAATCGCCTAGAAAGCGTCATAAACTTTCAAACTGGTTTTGGGGTATGTTGGTGTGGTTTATATCTTAAAAACGCCTCTATGGGCTTCTGAGGCTGCGAAGGCATGGTTTTTATAGGCTAACTTACGGATTTAACATAAGGGGCATAATGAACACCAATTTTGTGCGTACCTAGATTCAACACAAGACAACCTACTTCACCCGCGGGGATCTTAGGTTCTGCGTTACGGAGTCACGTGATTGAAAATAACTGGACTCTGTATATGTAACATGTATAATTAACTTAACAATGACTTCCTCCTGAACGCAAGTTCGTAAGAGTGGAAGTTTTTTTTGCCTAGAGAAAATGGATATGACTACGGCTCCCCTCGTTCCTTATAAAAAACCTTACCTATCTAGTTCTCAGTTATGTAAAAAACTCATAGACCAAGGGCTAATCATTGATGATGAAGATTTTGCTGAAAAAGTCTTAAATCGTTGTAGTTACTATCGATTTAAAGCTTATCTATCACCCTTTAAAGATAAAGGGACTAAGAAATTTTCAGAACATACTACGTTTCATAATGGTTATGAACTTTATATGTTTGATAGTGAGTTGCGAAGTTATATCTTTGATATAATTGAAAAAGTTGAGATAGGCGTTCGTTCGGCGCTAGATCAGTGGATTACGAAACAAACAGATAATCCATTCTGGTATTTAGATGCATCT'
    genome = Genome(seq)

    import pprint

    pp = pprint.PrettyPrinter()
    pp.pprint(genome.get_codon_dist())
    pp.pprint(genome.get_amino_acid_dist())
    print(sum(genome.get_amino_acid_dist().values()))
