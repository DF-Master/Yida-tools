from Bio import PDB

pdb_ub = 'G:/MS_data/220125-UBBSA/mono_ub.pdb'


def cal_distance(animo_pos1,
                 animo_pos2,
                 pdb,
                 name='default',
                 animo_core1='CA',
                 animo_core2='CA',
                 in_one_domain=True,
                 autoprint=True):
    parser = PDB.PDBParser()
    structure = parser.get_structure(name, pdb)
    model = structure[0]
    chain = model['A']
    residue1 = chain[animo_pos1]  # Start from 1 not 0
    residue2 = chain[animo_pos2]
    atom1 = residue1[animo_core1]
    atom2 = residue2[animo_core2]
    distance = atom1 - atom2
    if autoprint == True:
        print('Distance between', residue1.get_resname(), animo_core1, 'and',
              residue2.get_resname(), animo_core2, 'is ', distance)
    return distance


if __name__ == '__main__':
    cal_distance(1, 2, pdb_ub)
