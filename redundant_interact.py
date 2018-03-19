from Bio.PDB import *
import sys
import os
import copy
from Bio import pairwise2

if len(sys.argv)<2:
    inter = None
else:
    inter = sys.argv[1]

if os.path.isdir(inter):
    files_list = []
    for files in os.listdir(inter):
        if files.endswith(".pdb"):
            files_list.append(os.path.join(inter, files))

elif inter is None:
    path = "./"
    files_list = []
    for files in os.listdir(path):
        if files.endswith(".pdb"):
            files_list.append(files)

def GetStructures(pdbfile):
    """
    Given a pdbfile, gets it's structure.
    :param pdbfile: pdb file to open
    :return: structure object
    """
    parser = PDBParser()
    structure = parser.get_structure(pdbfile[0:-4], pdbfile)
    return structure

def Alignsequence(structure1, structure2):
    """
    Given 2 structures, gets the sequence and aligns it
    :param structure1: structure 1 to align
    :param structure2: structure 2 to align
    :return: Alignment
    """
    ppb = PPBuilder()
    for pp in ppb.build_peptides(structure1):
        sequence1 = pp.get_sequence()
    for pp in ppb.build_peptides(structure2):
        sequence2 = pp.get_sequence()

    alignment = pairwise2.align.globalxx(sequence1, sequence2)
    return alignment

class Superimpose(object):

    def __init__(self, structure1, structure2):
        self.structure1 = structure1
        self.structure2 = structure2
        self.si = Superimposer()

    def SuperimposeStructures(self):
        """
        :param structure1: first structure to superimpose
        :param structure2: second structure to superimpose
        :return: rmsd object
            """
        atoms_a = list(self.structure1.get_atoms())
        atoms_b = list(self.structure2.get_atoms())
        if len(atoms_a) > len(atoms_b):
            atoms_a = atoms_a[:len(atoms_b)]
        else:
            atoms_b = atoms_b[:len(atoms_a)]

        self.si.set_atoms(atoms_a, atoms_b)

        return self.si

    def getRMSD(self):

        superimpose = self.SuperimposeStructures()
        rmsd = superimpose.rms

        return rmsd

    def Rotate(self, atoms_b):

        superimpose = self.SuperimposeStructures()
        return superimpose.apply(atoms_b)

def get_interactions(list_atoms1, list_atoms2):
    """given 2 lists of atoms corresponding to 2 different chains,
    returns a tuple with 3 elements:
        1. tuple of chains that interact, i.e. ("A","B")
        2. tuple of residue number that interact, i.e. (125, 543)
        3. distance
    """
    beta_carbons1 = list(filter(lambda x: x.get_id() == "CB", list_atoms1))
    beta_carbons2 = list(filter(lambda x: x.get_id() == "CB", list_atoms2))
    ns = NeighborSearch(beta_carbons1)
    interactions = []

    for atom in beta_carbons2:
        interact = ns.search(atom.get_coord(), 7.0)
        interactions.extend([tuple(sorted([str(atom.get_parent().resname), str(x.get_parent().resname)])) for x in interact])        
    return interactions


structures = {}

for file in files_list:
    structures[file] = GetStructures(file)


# print(structures)

structures2 = copy.deepcopy(structures)
scores = {}
all_files = set()
discarded_files = set()
superimpositions = []


for file, structure in structures2.items():
    chains1 = list(structure[0].get_chains())
    all_files.add(file)
    for file2, structure2 in structures2.items():
        chains2 = list(structure2[0].get_chains())
        if file2 not in all_files:
            scores[file + file2] = []
            if file != file2:

                for chain in chains1:
                #     #chains_copy.remove(chain)
                     for chain2 in chains2:
                         Alignment = Alignsequence(chain, chain2)
                         scores[file + file2].append(Alignment[0][2]/len(list(chain.get_residues())))
               
                same_chains = list(filter(lambda x: x[1] > 0.95, enumerate(scores[file + file2])))
                
                # print(chains1[0].id + chains1[1].id, chains2[0].id + chains2[1].id)
                if len(same_chains) == 2 or len(same_chains) == 4:
                    distances1 = set(get_interactions(list(chains1[0].get_atoms()), list(chains1[1].get_atoms())))
                    distances2 = set(get_interactions(list(chains2[0].get_atoms()), list(chains2[1].get_atoms())))
                    distances_union = distances1.union(distances2)
                    distances_intersect = distances1.intersection(distances2)
                    # distances_diff = distances_union.difference(distances_intersect)
                    percent_similarity = len(distances_intersect) / len(distances_union) * 100
                    if percent_similarity >= 80:
                        discarded_files.add(file2)



selected_files = list(all_files.difference(discarded_files))
print(selected_files)

io = PDBIO()
i = 0
for i in range(len(selected_files)):
    with open("./selected_pdbs/dimer" + str(i+1) + ".pdb", "w") as fp:
        io.set_structure(structures[selected_files[i]])
        io.save(fp, write_end = 1)

 


                    # if len(distances_diff) / len(dista)
                    # if distancpees:
                    #     discarded_files.add(file2)



# selected_files = list(all_files.difference(discarded_files))
# print(selected_files)
    #     io.set_structure(structure2[0])
                #     io.save(fp, write_end = 1)


                #final_sup = Superimpose(structure2[0], structure[0])
                
                #print(rotated_prot)
                # # superimposed_chains_str.apply(structure2[0])
                # superimpositions.append(superimposed_chains_str.apply(structure2[0]))









#                 chains1_interact = get_interactions(list(chains1[0].get_atoms()), list(chains1[1].get_atoms()))
#                 chains2_interact = get_interactions(list(chains2[0].get_atoms()), list(chains2[1].get_atoms()))
#                 print(chains1_interact)
#                 print(chains2_interact)
# #                 if set(chains1_interact) != set(chains2_interact):
# #                     desired_interactions.append((file, file2))

# sup = Superimposer()
# # Specify the atom lists
# # 'fixed' and 'moving' are lists of Atom objects
# # The moving atoms will be put on the fixed atoms
# sup.set_atoms(fixed, moving)
# # Print rotation/translation/rmsd
# print sup.rotran
# print sup.rms
# # Apply rotation/translation to the moving atoms
# sup.apply(moving)
# # print(desired_interactions)
                


# print(scores)











