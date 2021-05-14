# -*- coding: utf-8 -*-

from Bio.PDB.MMCIFParser import MMCIFParser  # Tip: This module might be useful for parsing...
import numpy as np

aa_dict = {
    'ALA': 'A',
    'CYS': 'C',
    'ASP': 'D',
    'GLU': 'E',
    'PHE': 'F',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LYS': 'K',
    'LEU': 'L',
    'MET': 'M',
    'ASN': 'N',
    'PRO': 'P',
    'GLN': 'Q',
    'ARG': 'R',
    'SER': 'S',
    'THR': 'T',
    'VAL': 'V',
    'TRP': 'W',
    'TYR': 'Y',
}


############# Exercise 2: Protein Data Bank #############
# General remark: In our exercise every structure will have EXACTLY ONE model.
# This is true for nearly all X-Ray structures. NMR structures have several models.
class PDB_Parser:
    CIF_PARSER = MMCIFParser(QUIET=True)  # parser object for reading in structure in CIF format

    def __init__(self, path):
        """
            Initialize every PDB_Parser with a path to a structure-file in CIF format.
            An example file is included in the repository (7ahl.cif).
            Tip: Store the parsed structure in an object variable instead of parsing it
            again & again ...
        """
        self.structure = self.CIF_PARSER.get_structure('7AHL',
                                                       path)  # Parse the structure once and re-use it in the functions below

    # 2.8 Chains    
    def get_number_of_chains(self):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
            Return:
                Number of chains in this structure as integer.
        """
        n_chains = list(self.structure.get_chains())
        return len(n_chains)

    # 2.9 Sequence  
    def get_sequence(self, chain_id):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the amino acid sequence (single-letter alphabet!) of a given chain (chain_id)
                in a Biopython.PDB structure as a string.
        """

        sequence = []
        chain = next(ch for ch in self.structure.get_chains() if ch.id == chain_id)
        for res in chain.child_list:
            mol = aa_dict.get(res.resname)
            if mol is not None:
                sequence.append(mol)
        return ''.join(sequence)

    # 2.10 Water molecules
    def get_number_of_water_molecules(self, chain_id):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the number of water molecules of a given chain (chain_id)
                in a Biopython.PDB structure as an integer.
        """
        n_waters = 0
        chain = next(ch for ch in self.structure.get_chains() if ch.id == chain_id)
        for res in chain.child_list:
            if res.resname == 'HOH':
                n_waters += 1
        return n_waters

    # 2.11 C-Alpha distance    
    def get_ca_distance(self, chain_id_1, index_1, chain_id_2, index_2):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id_1 : String (usually in ['A','B', 'C' ...]. The number of chains
                                depends on the specific protein and the resulting structure)
                index_1    : index of a residue in a given chain in a Biopython.PDB structure
                chain_id_2 : String (usually in ['A','B', 'C' ...]. The number of chains
                            depends on the specific protein and the resulting structure)
                index_2    : index of a residue in a given chain in a Biopython.PDB structure

                chain_id_1 and index_1 describe precisely one residue in a PDB structure,
                chain_id_2 and index_2 describe the second residue.

            Return:
                Return the C-alpha (!) distance between the two residues, described by
                chain_id_1/index_1 and chain_id_2/index_2. Round the returned value via int().

            The reason for using two different chains as an input is that also the distance
            between residues of different chains can be interesting.
            Different chains in a PDB structure can either occur between two different proteins
            (Heterodimers) or between different copies of the same protein (Homodimers).
        """
        ca_1 = next(ch for ch in self.structure.get_chains() if ch.id == chain_id_1)[index_1]['CA']
        ca_2 = next(ch for ch in self.structure.get_chains() if ch.id == chain_id_2)[index_2]['CA']

        ca_distance = np.linalg.norm(ca_1.coord - ca_2.coord)
        return int(ca_distance)

    def get_residue_sequence(self, chain_id):
        chain = next(ch for ch in self.structure.get_chains() if ch.id == chain_id)
        seq = []
        for res in chain.child_list:
            if aa_dict.get(res.resname) is not None:
                seq.append(res)
        return seq

    # 2.12 Contact Map  
    def get_contact_map(self, chain_id):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return a complete contact map (see description in exercise sheet)
                for a given chain in a Biopython.PDB structure as numpy array.
                The values in the matrix describe the c-alpha distance between all residues
                in a chain of a Biopython.PDB structure.
                Only integer values of the distance have to be given (see below).
        """

        seq = self.get_residue_sequence(chain_id)
        length = len(seq)
        contact_map = np.zeros((length, length), dtype=np.float32)

        for i in range(length):
            for j in range(length):
                contact_map[i][j] = np.linalg.norm(seq[i].child_dict.get('CA').coord - seq[j].child_dict.get('CA').coord)

        return contact_map.astype(np.int64)  # return rounded (integer) values

    # 2.13 B-Factors    
    def get_bfactors(self, chain_id):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the B-Factors for all residues in a chain of a Biopython.PDB structure.
                The B-Factors describe the mobility of an atom or a residue.
                In a Biopython.PDB structure B-Factors are given for each atom in a residue.
                Calculate the mean B-Factor for a residue by averaging over the B-Factor
                of all atoms in a residue.
                Sometimes B-Factors are not available for a certain residue;
                (e.g. the residue was not resolved); insert np.nan for those cases.

                Finally normalize your B-Factors using Standard scores (zero mean, unit variance).
                You have to use np.nanmean, np.nanvar etc. if you have nan values in your array.
                The returned data structure has to be a numpy array rounded again to integer.
        """
        seq = self.get_residue_sequence(chain_id)
        b_factors = np.zeros(len(seq), dtype=np.float32)
        for r, res in enumerate(seq):
            b_factors[r] = np.mean([res.child_list[atom_idx].bfactor for atom_idx in range(len(res.child_list))])
        #
        # for r, res in enumerate(seq):
        #     atom_bfs = np.zeros(len(res.child_list), dtype=np.float32)
        #     for a, atom in enumerate(res.child_list):
        #         atom_bfs[a] = atom.bfactor
        #     b_factors[r] = np.mean(atom_bfs)

        b_factors = (b_factors - np.nanmean(b_factors)) / np.nanstd(b_factors)
        return b_factors.astype(np.int64)  # return rounded (integer) values


def main():
    print('PDB parser class.')
    parser = PDB_Parser('./tests/7ahl.cif')
    # print(parser.get_number_of_chains())
    # print(parser.get_sequence('C'))
    # print(parser.get_number_of_water_molecules('D'))
    # print(parser.get_ca_distance('A', 20, 'A', 55))
    # print(parser.get_ca_distance('A', 121, 'E', 120))

    # import pprint
    # pp = pprint.PrettyPrinter()
    with np.printoptions(threshold=np.inf):
        # test_cm = np.load('./tests/contact_map_1.npy')
        # print('JSON shape:', test_cm.shape)
        # matching_arr = parser.get_contact_map('A') == test_cm
        # print(f"contact maps matched: {matching_arr.all()}")

        test_bf = np.load('./tests/bfactors_5.npy')
        # print(test_bf)
        bfs = parser.get_bfactors('E')
        # print(bfs)
        matching_arr = bfs == test_bf
        print(f"B-Factors matched: {matching_arr.all()}")
    return None


if __name__ == '__main__':
    main()
