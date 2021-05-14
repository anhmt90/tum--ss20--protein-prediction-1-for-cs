# The following function should take a string and a sub-string and return 
# a list of starting positions for the substring


def get_sub_position(search_string, needle):
    sub_positions = set()
    for i, _ in enumerate(search_string):
        if (i + len(needle)) > len(search_string):
            break
        pos = search_string.find(needle, i)
        if pos > -1 and pos not in sub_positions:
            sub_positions.add(pos)
    return sorted(sub_positions)


if __name__ == '__main__':
    print(get_sub_position('AAAA', 'AAA'))
    print(get_sub_position("ATGAACCAACCTTATGAAGGAATTAACCTG", "CCC"))
    print(get_sub_position("ATGAACCAACCTTATGAAGGAATTAACCTG", "TGA"))
