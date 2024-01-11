import Bio.PDB
import Bio.SeqRecord

import h5py
from sklearn.manifold import TSNE


# Wrapper class for single atoms
class Atom:
    def __init__(self, bio_atom, residue, index):
        self.bio_atom = bio_atom
        self.residue = residue
        self.index = index


# Wrapper class for collections of atoms
class Residue:
    def __init__(self, atoms, bio_residue, index):
        self.atoms = atoms
        self.bio_residue = bio_residue
        self.index = index
        self.color = [0, 0, 0]


class Protein:
    def __init__(self, structure_path, sequence_path, embedding_path):
        parser = Bio.PDB.PDBParser(QUIET=True)
        bio_structure = parser.get_structure("id", structure_path)
        self.residues = []
        self.atoms = []

        # Initialize data from PDB file
        atom_index = 0
        for i, bio_residue in enumerate(bio_structure.get_residues()):
            residue_atoms = []
            residue = Residue(residue_atoms, bio_residue, i)
            for bio_atom in bio_residue.get_atoms():
                residue_atoms.append(Atom(bio_atom, residue, atom_index))
                atom_index += 1
            self.residues.append(residue)
            self.atoms.extend(residue_atoms)

        embedding_file = h5py.File(embedding_path, 'r')
        sequence_file = open(sequence_path, 'r')

        title = sequence_file.readline()[1:].strip()
        self.sequence = sequence_file.readline().strip()

        embeddings = embedding_file[title][()]

        # Use the t-SNE algorithm to transform the embeddings into 2D vectors
        transform = TSNE(n_components=2, learning_rate='auto', init='random', perplexity=3).fit_transform(embeddings)
        self.embedding_points = transform.flatten()
