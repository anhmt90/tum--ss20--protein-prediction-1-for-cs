# -*- coding: utf-8 -*-
from Bio import SeqIO  # Tip: This module might be useful for parsing...
from Bio.SeqIO.UniprotIO import Parser

############ Exercise 3: SwissProt ##########
class SwissProt_Parser:
    PARSER = SeqIO

    def __init__(self, path, frmt='uniprot-xml'):
        """
            Initialize every SwissProt_Parser with a path to a XML-formatted UniProt file.
            An example file is included in the repository (P09616.xml).
            Tip: Store the parsed XML entry in an object variable instead of parsing it
            again & again ...
        """

        self.sp_record = next(SeqIO.parse(path, frmt))

        # parser = Parser(path)
        # self.sp_anno = parser.parse()

    # 2.2 SwissProt Identifiers
    def get_sp_identifier(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                Unique SwissProt identifier for the given xml file
        """

        return self.sp_record.id

    # 2.3 SwissProt Sequence length
    def get_sp_sequence_length(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                Return sequence length of the UniProt entry as an integer.
        """

        return self.sp_record.annotations['sequence_length']

    # 2.4 Organism 
    def get_organism(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                Return the name of the organsim as stated in the corresponding field
                of the XML data. Return value has to be a string.
        """

        return self.sp_record.annotations['organism']

    # 2.5 Localizations
    def get_localization(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                Return the name of the subcellular localization as stated in the
                corresponding field.
                Return value has to be a list of strings.
        """

        return self.sp_record.annotations['comment_subcellularlocation_location']

    # 2.6 Cross-references to PDB
    def get_pdb_support(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                Returns a list of all PDB IDs which support the annotation of the
                given SwissProt XML file. Return the PDB IDs as list.
        """
        pdb_ids = []
        for ref in self.sp_record.dbxrefs:
            if 'PDB:' in ref:
                pdb_ids.append(ref.partition(':')[2])

        return pdb_ids


def main():
    print('SwissProt XML Parser class')
    sp = SwissProt_Parser('./tests/P09616.xml')
    print(sp.sp_record)
    print("\n\n\n\n\n")
    print(sp.get_sp_identifier())
    print(sp.get_sp_sequence_length())
    print(sp.get_localization())
    return None


if __name__ == '__main__':
    main()
