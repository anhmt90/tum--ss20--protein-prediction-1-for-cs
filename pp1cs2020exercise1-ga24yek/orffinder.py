from genome import base_complements, inversed_genetic_code_table
from collections import deque


def is_start_codon(codon):
    return genetic_code_table.get(codon) == 'M'


def is_stop_codon(codon):
    return genetic_code_table.get(codon) is None


def get_genetic_code_table():
    gc_table = {}
    for key, val in inversed_genetic_code_table.items():
        for v in val:
            gc_table[f'{v}'] = key

    assert (gc_table.get('TAA') is None) and (gc_table.get('TAG') is None) and (gc_table.get('TGA') is None)
    return gc_table


def is_valid(seq, min_num_aa):
    valid = isinstance(seq, str) and isinstance(min_num_aa, int) and (len(seq) >= 3)
    if valid:
        for base in seq:
            valid = valid and (base_complements.get(base) is not None)
    return valid


def get_reverse_complement_in_reading_direction(dna):
    complement = []
    for base in dna:
        complement.append(base_complements.get(base))
    return ''.join(complement[::-1])


def get_reading_frame(seq, start, rev=False):
    """
    This function also considers circular DNA
    :param seq:
    :param start:
    :param rev:
    :return:
    """
    assert 0 <= start <= 2
    i = start
    while i <= (len(seq) // 3):
        codon = f'{seq[i]}{seq[i + 1]}{seq[i + 2]}'
        i += 3


def get_polypeptide(orf):
    assert len(orf) % 3 == 0
    polypeptide = []
    i = 0
    for i in range(0, len(orf), 3):
        codon = f'{orf[i]}{orf[i + 1]}{orf[i + 2]}'
        aa = genetic_code_table.get(codon)
        if aa is not None:
            polypeptide.append(aa)
            # assert i == (len(orf) - 3)
            # break
        i += 3
    return ''.join(polypeptide)


def get_orf_len(start, stop, seq_len):
    num_nt = ((stop - start) if start <= stop else (stop + seq_len - start)) + 1
    num_nt -= 3
    assert num_nt % 3 == 0
    return num_nt // 3


def is_overlapped(start1, stop1, start2, stop2):
    return max(start1, start2) < min(stop1, stop2)


def is_longer_frame(frame1, frame2, seq_len):
    return get_orf_len(frame1[0], frame1[1], seq_len) > get_orf_len(frame2[0], frame2[1], seq_len)


def find_orfs(seq, min_num_aa):
    # orfs_in_3_rfs = set()
    orfs_in_3_rfs = {}
    seq_len = len(seq)
    for f in range(3):
        orfs_ = {}
        start = -1
        search_seq = seq_len
        circular = False
        c = f
        # circular_seq = seq[f:seq_len] + seq[:f] + seq[f:seq_len]
        while c <= search_seq:
            c0 = c % seq_len
            c1 = (c + 1) % seq_len
            c2 = (c + 2) % seq_len
            codon = f'{seq[c0]}{seq[c1]}{seq[c2]}'

            if ((c + 2) >= seq_len) and (start != -1 or is_start_codon(codon)) and not circular:
                search_seq += search_seq
                circular = True

            if is_start_codon(codon) and start == -1:
                start = c0
            elif is_stop_codon(codon) and start != -1:
                stop = c2
                num_aa = get_orf_len(start, stop, seq_len)

                if num_aa >= min_num_aa:
                    if orfs_.get(stop) is None:
                        orfs_[stop] = start
                    else:
                        current_num_aa = get_orf_len(orfs_.get(stop), stop, seq_len)
                        if num_aa > current_num_aa:
                            orfs_[stop] = start
                    if circular:
                        break
                start = -1
            c += 3

        for stop, start in orfs_.items():
            if orfs_in_3_rfs.get(stop) is None:
                orfs_in_3_rfs[stop] = start
            else:
                curr_start = orfs_in_3_rfs[stop]
                if is_longer_frame((start, stop), (curr_start, stop), seq_len):
                    orfs_in_3_rfs[stop] = start

    return orfs_in_3_rfs


def get_orfs(seq, min_num_aa):
    """
    Find and return all ORFs within the provided genome sequence.

    :param seq: String with the genome sequence.
    :param min_num_aa: Minimum number of amino acids in the translated protein
                       sequence.

    :return: List of tuples representing the ORFs (details in worksheet).
    """

    if not is_valid(seq, min_num_aa):
        raise TypeError("Invalid genome input!")

    leading_orfs = find_orfs(seq, min_num_aa)
    seq_rc_rd = get_reverse_complement_in_reading_direction(seq)
    lagging_orfs = find_orfs(seq_rc_rd, min_num_aa)
    # lagging_orfs = []

    ORFs = []
    for stop, start in leading_orfs.items():
        aa_orf = seq[start:(stop + 1)] if start <= stop else (seq[start:] + seq[:(stop + 1)])
        ORFs.append((start, stop, get_polypeptide(aa_orf), False))

    for stop, start in lagging_orfs.items():
        aa_orf = seq_rc_rd[start:(stop + 1)] if start <= stop else (seq_rc_rd[start:] + seq_rc_rd[:(stop + 1)])
        ORFs.append((len(seq) - start - 1, len(seq) - stop - 1, get_polypeptide(aa_orf), True))

    return ORFs


genetic_code_table = get_genetic_code_table()

if __name__ == "__main__":
    sequence = 'TTATTTAAGAGCAATGGCCAACAAGTAAAAACGGTTAGCAGGGTTAGGGATATGTTTGTTGACTCTAAAGAAGAGTACGCAAAGCACTATCAAGAAAAGTATTACAATGAGTATTGTCCGTTTTACAGAGATCTCCCGCCTGCATGGGTAGCCATTGAGTTAATGACTTTCGGCAACGTAGTGAAGTTAATTCAAAACATCAGTGATGATAAAATTCAATCACTTAAGATGGATAGATTTTCTAAGAAGTTCAATATTCAGAAATTTCAGACATTAATTAGTTGGATGAATGTGCTGCACCAGATGAGAAACTACTGTGGGCATCATAACCGACTGTTTAATCGAAACTTCCCTGCTCCAACAGCGATTAAAAAGAGCTTGTCTGATGCAATTCCTCTTGTCAGGACCAAACCAAATCCAGATAAGCGTGAAGAGGATCAGTTAAACCGACTTTATACAGCTCTTGCTGCATTACAATGTATATATTCAGGGCTTGGTTTCGATGAAAAAATAGGACCAAAAATCTCTGATTTATTTGATAATTATACAACAACACAGCGGTTTAGCTTATCAATGGGTTTCCCTAATGGTTGGAAAGAAGAGCCGCTTTTTTTTGATTTATAATCGCTATACTTAACATAAAACCCGTTCTACGAATCGTAAAAGGTCGCCTATTTAGGTGGCCTTTTTTTATGAAAACTGTTTTTTAAGTGATAAAAATGTCGTTTGTCCTAACGAAATAGTATTTTTAACAGGCTTATTCAGGCGTTTTGATTTTTAGCTATCTGGTAATATTTAATTATTTTTTATTGGAGTTTGTATATGCCGCCACCTAAAGGTCGAAGTTTTCCATTTGCACCGCGTCACTCTGCTGATTGGTTAGTGAGTCATGTAACGTATGACCAAGCCGTTGATATGTTTTTTAATCAGACCGCAACACAACAGAATTTGGGGCATGATCCTTTAGTTTATTCAAAAGTATTTAGAGGTGTGACATACGCAACGCTCCAAGAAGCTCAACAAGTGTTTACTGAAACGATGAATGCAGAGTATGAAGTGAGAGAGCAAAGGGATTTAGCAGATGAAAATCGAGGGTGTAAGGCTTCTAAAATTCTCAATGATAATATTCGGAATCGTATTGTACCAACAGAGGATGGTTTAGAGTCATTAACTGCATATCGAAAGCCTTTTGGTGAAAAGTGTGTATCTCCTCTTTCTTTGTTTACCAAATCGTTAAATGGCGGTTCAAACAGTCATATTCAAGCAAACGAAACTCAAGCTCAATTAGTTACGCAACAAGCTTATGATTTCCCACTCGTCACTAAAGGTTCTAAAGCGATTCAGGCATTAGATGATTCTTGGGGTTTACAGTGCTCTTTAGAGGATGAATTAAACCAAAATCACCAAGATGTTGTTATTTTAGATTTGGAAGACTCACATAAATTATTAAATAATTTTTGGGTGGATTTAGGCAAAGATACAGCGTTAATAGCAACAAGTTTAAATACGGCTAACGGTTGGTTTCAACATACAACACCGATATTTGATGCTAAAGGTTTAGTTAAGCAATTCGGCGATATAAATATTAAAGCGGATATTGTAGAGAGTAAGGGGAAACAGTTTATTGCTTTCTCTGGTAAGAAGAACGGCAAAGAAATTTAACACGCACTTGTTAATGGCACTCGAATTAACATGAACGGTAAGAAGTACCCAATTAACAGCCCTAAAGTTCAACAAGTCGGTTTGAGTCCTAAAGCACGAGCAAATGGATTTAAAGGCGCAGGTGTATTGACGTTTGTGGTATCGGCAGCGATAGCAACGACAGATCTTGTTTTTAAAGACGATTATCACTTAGTTGATTGGTTCGGTAATGTTGGGGCAGATATGTTTAAGGCATTGCTACAGTTTGGTGCAGGAGAGGCTATATTGTTTGGAATTATTGCGATGACTGGTTACGTTACTCTGGGCTTGATTGCTGTGTTTTTTGTCTATGTGTCTATTGAGTGGATATGGAGTGAATACAAGGTAAATGATGAGGTAGTTAAGGGATTAGAAAGTGTTATCAGTTAAAATGAGGTTTTTAGGAATAGTATTTATACCTTTATTGATTTTACTTTGGTGGATACCTACAAATGGAGTCTTGGGGGATTACCAAGATTTATTAAATCAGACAGATGAAGTTAGATTATCACTGATAACATTAATCGTCCCCATTGGTGGTTTTATCCCTTTACTTATTTCCGTTGCTTTAATTTCTATTGCTATTTATTCGGGGAAACAAGCTCGTTTGGTCATTGGTGAAAAATGGACTAGTACTATCAATAGAACGTGTATTTATAGCATGATATTAGGTGTAGTATTTGCGGTTGTATTTGCACTATATTGTATAAAATTGCTTGATGAAAACGGTTATGAATACAGTTATAACCTAACCCAAATAACCCCGACAGGTATTCATTTAATGTATGTAAAATCACACAAAGAATAGAAACAAAAAAGCCACTCAAAAGAGTGGCTTTTTTGTCTATGAACCAGAGGATCATAAAATTGATCATTGAGATCTTGAGCATTTAGCCGTAACCTTATGCCCAGATGTGAAAAAACCAAGCCTGCAAGCTTGGTTTTTTCGAAATACTTAACTACGTTTAGACGGCGTTGTAAGTCCTTTCATGATAGATCATATCACATAAGAACTGCAAGTATTTCCAAGCAATCAGTGGTTGGGTTTAGAGTCACTAGGGTTGTGGGACAGCATCATGCGGTTGGTTTTCATACCATATTTCGCCCAGATACTGGACGAGAGACTAGCGTAACATCCGATAACGCACACCAAATGAGTGGACTGAAGCAAGAACTAGCGTAAGCAGCTAGGCAAAGAGGCTTGTAGATTGGGGAAACCCTTGCAAGGTACTGATGCAGAGGTCGGCAGAGAACGGTGAGGATGGGGGCAATACAGCCCGTAATCTCTGCACTTGAATAGTGTATTAGAATTCAAGCGAATGATGGCGATGAACGGCTCCAGTCGGCAAGCAATCAGGGTTCACACTAGAGTGATAAACTTAGCGAAAGCGGTTATCCTCTAGGGTGAGTATTACCGAAATCTAGCTCAAAACTCACCTATCAGCATACATGCCGACTTCTTTTGCTCTTCAATTCAAAAAGAATTAAAAAGAAAATATTAAGAGCTGCACGATTAGAAATGAAAGTTGCCGCAACGAGCATAAGCGAGTTAAGCCCACTTTCTTTTCTAATGTGCGAACAAGCGAAGCGCGTTAGTTGGGTCGATTATCTTCTTTTATGGATCGATGCTGACCAACTGAATGTGCTTTTTTCTTCTTTCGGTAACTTTTTCAATTCCGCTCTCGCTGCCCTTGCTATGGTGAGTCGTTTTAATTTTTCTCTGTTCCTTTTGAAAATCGTTTCAAAGGCTTGGTCTCTATAATCTTCACCGTTGCTGCTGCAGAGTGAATCTATAAGAGCAAAGTCTCTTGGCTTGAGTCCGTATACTTTTTCATAAATTGCTTTGGCCGCTTCTTTTCGTGCCACTCCTTTTTCTATTTGTAGTTCAATAAGAGTGCGGTCTTCTTGAGTAAATCGTCTTTTCTTCATTCTCTTTTTTGCTCTTAGGTTTCAAAGGACATCCTTTGGCGCGGCATGTTTTTCCAAATACGTATGTAGAGGGTGAAAACGTGCCCGTCGTGCCTATCTAAAAAATAGATATCACTGAGTAATCAATAATAAAATGCTAACTCCTTGTATTAATACGTTAATTTATGGGTAACGTTATTTTTTAAAGAAACTATCCAAAAGCGTCTTTTTTCCTTCTTTTTTCGTTGGTGTTAATACCGGCTCTTTATGCTCTCTTATCATGGCTTCTAGTTGGTCAAAACGAGTGGCGATGTAATGTAAGTCATCATCAAGCATTTCTCTTTCGTAGGCGTAATACGATTTATCAATTGAGTGATGCACTCTATCGCCATGAGTTTGATACAAATCCCTTTTATAGGCGAACAAACCTCTTACTAGCAACGCCCTTAACAGTTCTGCTTTATCTTTATGCGTGTCGTTTAACAGCTCTTCAATTTCAGCATTCACTCGACCATCCACTCGAACGGACACCAAATGAGTCTTGGATTTTGGTTTATTCATTTTCTAGTTTTCTCTTTTGCCTAATTCGATAAAGAACCGTCTTCAAATAGTCTCTCGATACCGTGATATTTTTTGTCTTTAAAAACGCAACGATGTGCTCTAGTGAGATACCTTCACTCATCAAGATGGTTATCACGTCATATAACTCCGTAATGACGGCGGTTTTAGATTTCGTTTGATGCTTTTTCGTGTACGAAAGAAGCTCTTTTTCTAACGGTGTTTTCTTTTTCACTGAATACTTTATTGATTACTGATTGATACTTTAGAGACTAGCAATAAAAAAGAAGAACGAGAATGGTTATTTCAGTTGGGCTGAAATTGGCATGATTTATCTTAGTAACCTACCTATAAACTCATTGATTCAATGGGTTAGCTTTTTAATGTTAGTTACGTTGTGATTACCCGTTGATTTCTATTTTCAAGATAGGCACGACGGGGACGTTTTCACCCCCTGCTACGCATTGGATAAAAACATCCCGTGCCAAAGGTGAAACCTTTGAAACCCTGACGCACTTCGTTTGTCTCCAGTACGCCCAAAAAAGCCAGCTACGCTTCGGGCTTTTTCGAGCGCACCGGCTAATTGTTGAAAATAGACGTTACAGGGTGATTATTGATTGAAACGGATTGTGATTTTATTGCCCTGTAAATCGCCTAGAAAGCGTCATAAACTTTCAAACTGGTTTTGGGGTATGTTGGTGTGGTTTATATCTTAAAAACGCCTCTATGGGCTTCTGAGGCTGCGAAGGCATGGTTTTTATAGGCTAACTTACGGATTTAACATAAGGGGCATAATGAACACCAATTTTGTGCGTACCTAGATTCAACACAAGACAACCTACTTCACCCGCGGGGATCTTAGGTTCTGCGTTACGGAGTCACGTGATTGAAAATAACTGGACTCTGTATATGTAACATGTATAATTAACTTAACAATGACTTCCTCCTGAACGCAAGTTCGTAAGAGTGGAAGTTTTTTTTGCCTAGAGAAAATGGATATGACTACGGCTCCCCTCGTTCCTTATAAAAAACCTTACCTATCTAGTTCTCAGTTATGTAAAAAACTCATAGACCAAGGGCTAATCATTGATGATGAAGATTTTGCTGAAAAAGTCTTAAATCGTTGTAGTTACTATCGATTTAAAGCTTATCTATCACCCTTTAAAGATAAAGGGACTAAGAAATTTTCAGAACATACTACGTTTCATAATGGTTATGAACTTTATATGTTTGATAGTGAGTTGCGAAGTTATATCTTTGATATAATTGAAAAAGTTGAGATAGGCGTTCGTTCGGCGCTAGATCAGTGGATTACGAAACAAACAGATAATCCATTCTGGTATTTAGATGCATCT'
    orfs = get_orfs(sequence, 34)

    import pprint

    pp = pprint.PrettyPrinter()
    pp.pprint(orfs)
    print("Num:", len(orfs))
