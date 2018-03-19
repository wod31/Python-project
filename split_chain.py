from Bio.PDB import *
import sys


parser = PDBParser()
structure = parser.get_structure("2f1d", "/Users/titu/OneDrive/master_bioinformatics/2nd term/PYT/project/2f1d.pdb")
chains = list(structure[0].get_chains())

def get_interactions(list_atoms1, list_atoms2):
	"""given 2 lists of atoms corresponding to 2 different chains,
	returns a tuple with 3 elements:
		1. tuple of chains that interact, i.e. ("A","B")
		2. tuple of residue number that interact, i.e. (125, 543)
		3. distance
	"""
	ns = NeighborSearch(list_atoms1)
	interactions = []

	for atom in list_atoms2:
		interact = ns.search(atom.get_coord(), 8.0)
		interactions.extend([(str(atom.get_parent().get_parent().get_id())+ str(atom.get_parent().id[1]), str(x.get_parent().get_parent().get_id()) + str(x.get_parent().id[1])) for x in interact])
	
	return interactions

between_chains = []

for index1, chain1 in enumerate(chains):
	atoms_chain1 = Selection.unfold_entities(chain1, 'A')
	beta_carbons_chain1 = list(filter(lambda x: x.get_id() == "CB", atoms_chain1))

	for chain2 in chains[index1+1:]:
		atoms_chain2 = Selection.unfold_entities(chain2, 'A')
		beta_carbons_chain2 = filter(lambda x: x.get_id() == "CB", atoms_chain2)
		provi = get_interactions(list(beta_carbons_chain1), list(beta_carbons_chain2))
		if len(provi) > 0:
			between_chains.extend(provi)

print(between_chains)
set_dimers = set([x[0] + y[0] for x, y in between_chains])


for dimer in set_dimers:
	pdb_chain_file = './split_chains/chain_{}_{}.pdb'.format(dimer[0], dimer[1])
	chain1 = list(filter(lambda x: x.id == dimer[0], chains))[0]
	chain2 = list(filter(lambda x: x.id == dimer[1], chains))[0]
	with open(pdb_chain_file, "w") as fp:
		io = PDBIO()
		io.set_structure(chain1)
		io.save(fp, write_end = 0)
		io.set_structure(chain2)
		io.save(fp, write_end = 1)

