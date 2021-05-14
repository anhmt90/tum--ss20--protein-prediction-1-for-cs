# The following method should take a strand as an argument and
# return a complementary strand


def complementary(strand):
    strand = strand.upper()
    complementary_strand = []
    complement = ''
    for nucleobase in strand:
        if nucleobase == 'A':
            complement = 'T'
        elif nucleobase == 'T':
            complement = 'A'
        elif nucleobase == 'C':
            complement = 'G'
        elif nucleobase == 'G':
            complement = 'C'
        else:
            raise InvalidNucleobaseException(f"{nucleobase} is an invalid nucleobase")
        complementary_strand.append(complement)

    return "".join(complementary_strand)


class InvalidNucleobaseException(Exception):
    pass
