import Bio.PDB
import Bio.SeqRecord
import colour

import h5py
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN


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
        self.highlighted = False


class Protein:
    RESIDUE_INDEX = 1
    CLUSTER_INDEX = 2
    ATOM_TYPE = 3

    RAINBOW = 8
    POISSON = 9

    poisson_palette = [(187, 176, 148), (128, 118, 101), (89, 82, 70), (51, 51, 51), (25, 31, 34), (47, 68, 67),
                       (59, 94, 88), (90, 140, 108), (139, 180, 141), (192, 208, 165), (247, 239, 199),
                       (161, 205, 176), (112, 147, 149), (74, 120, 123), (56, 49, 64), (115, 77, 92),
                       (167, 103, 114), (204, 134, 125), (224, 186, 139), (195, 130, 82), (161, 86, 60),
                       (111, 52, 45), (68, 39, 31)]

    cpk_colors = {"C": (64, 58, 64), "O": (219, 73, 70), "N": (70, 110, 219), "S": (235, 208, 56),
                  "P": (235, 145, 56),
                  "_": (255, 255, 255)}

    def __init__(self, structure_path, embedding_path,
                 color_mode=CLUSTER_INDEX, color_palette=RAINBOW, cluster_max_distance=6.0):
        self.color_mode = color_mode
        self.color_palette = color_palette

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
        embeddings = embedding_file[list(embedding_file.keys())[0]][()]

        # Use the t-SNE algorithm to transform the embeddings into 2D vectors
        transform = TSNE(n_components=2, learning_rate='auto', init='random', perplexity=3).fit_transform(embeddings)
        self.embedding_points = transform.flatten()

        # Calculate residue clusters
        db = DBSCAN(eps=cluster_max_distance).fit(transform)

        self.cluster_index = db.labels_
        self.cluster_count = len(set(self.cluster_index)) - (1 if -1 in self.cluster_index else 0)

        self.update_colors()

    def update_colors(self, new_color_mode=None, new_color_palette=None):
        if new_color_mode:
            self.color_mode = new_color_mode
        if new_color_palette:
            self.color_palette = new_color_palette
        for residue in self.residues:
            self.color_residue(residue)

    def color_residue(self, residue):
        match self.color_mode:
            case self.RESIDUE_INDEX:
                color = self.get_color(residue.index / len(self.residues), residue.highlighted)

                for atom in residue.atoms:
                    atom.color = color
                residue.color = color
            case self.CLUSTER_INDEX:
                color = self.get_color(self.cluster_index[residue.index] / self.cluster_count, residue.highlighted)
                for atom in residue.atoms:
                    atom.color = color
                residue.color = color
            case self.ATOM_TYPE:
                residue.color = self.get_color(self.cluster_index[residue.index] / self.cluster_count, residue.highlighted)
                for atom in residue.atoms:
                    atom.color = self.cpk_colors[atom.bio_atom.get_id()[0]]

    def get_color(self, x, highlight=False):
        def to_rgb(c):
            return [int(b * 255) for b in c.rgb]

        def get_color_from_palette(palette):
            step = 1.0 / len(palette)
            i = int(x * len(palette))
            x_i, x_j = step * i, step * (i + 1)
            z = (x - x_i) / (x_j - x_i)
            a, b = palette[i], palette[(i + 1) % len(palette)]
            return [(b[j] * z + a[j] * (1.0 - z)) / 255.0 for j in range(3)]

        color = colour.Color(hue=0, saturation=0.0, luminance=0.0)

        if x >= 0:
            match self.color_palette:
                case self.RAINBOW:
                    color = colour.Color(hue=x * 0.8, saturation=0.65, luminance=0.5)
                case self.POISSON:
                    color.rgb = get_color_from_palette(self.poisson_palette)

        if highlight:
            color.set_luminance(color.get_luminance() + 0.25)
        return to_rgb(color)
