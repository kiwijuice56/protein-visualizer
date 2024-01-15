import os
import subprocess
import warnings

import Bio.PDB
import Bio.SeqRecord
import numpy as np
from Bio import SeqIO

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
        self.color = [0, 0, 0]


# Wrapper class for collections of atoms
class Residue:
    def __init__(self, atoms, bio_residue, index, seq_index):
        self.atoms = atoms
        self.bio_residue = bio_residue
        self.index = index
        self.seq_index = seq_index
        self.color = [0, 0, 0]
        self.highlighted = False


# Contains information about a protein, such as its 3D structure and node embeddings
class Protein:
    # Color modes
    RESIDUE_INDEX = 1  # Color residues by their index in the amino acid sequence
    CLUSTER_INDEX = 2  # Color residues by their associated cluster in the embedding space
    ATOM_TYPE = 3  # Color residues by their atoms using CPK coloring

    # Color palettes (for color modes 1 and 2)
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

    def __init__(self, pdb_path, chain_id=None, color_mode=CLUSTER_INDEX, color_palette=RAINBOW, cluster_distance=6.0):
        """
        @param pdb_path: Filepath to .pdb file containing the protein structure
        @param chain_id: Which chain to load from the .pdb file. Defaults to first chain
        @param color_mode: Options: Protein.RESIDUE_INDEX, Protein.CLUSTER_INDEX, Protein.ATOM_TYPE
        @param color_palette: Options: Protein.RAINBOW, Protein.POISSON
        @param cluster_distance: Maximum distance between nodes of a single cluster within the embedding space
        """
        self.color_mode = color_mode
        self.color_palette = color_palette

        # Get the full 3D structure from the PDB file
        warnings.filterwarnings("ignore", "'HEADER' line not found; can't determine PDB ID.")
        bio_structure = Bio.PDB.PDBParser(QUIET=True).get_structure("struct", pdb_path)

        # Choose default chain if none specified
        chains = [c for c in bio_structure.get_chains()]
        if chain_id is None:
            chain_id = chains[0].get_id()

        # Retrieve the requested chain from the PDB file
        chain = None
        for other_chain in chains:
            if other_chain.get_id() == chain_id:
                chain = other_chain
                break

        if chain is None:
            raise Exception(
                f"Chain with ID '{chain_id}' could not be found within .pdb file."
                f"Available chains: {[chain.get_id() for chain in chains]}.")

        # There are two methods of attaining an amino acid sequence from a PDB file: `atom` and `seqres`
        # `atom` calculates the sequence from the available residues found in the 3D space
        # `seqres` retrieves the sequence from the PDB file metadata

        # We will use the full, accurate sequence from the `seqres` record in order to calculate embeddings,
        # but we also need the `atom` sequence in order to match the 3D data with those same embeddings
        physical_record, full_record = None, None
        for other_record in [r for r in SeqIO.parse(pdb_path, "pdb-atom")]:
            if other_record.annotations["chain"] == chain_id:
                physical_record = other_record
        for other_record in [r for r in SeqIO.parse(pdb_path, "pdb-seqres")]:
            if other_record.annotations["chain"] == chain_id:
                full_record = other_record
        if full_record is None:
            warnings.warn("Could not read SEQRES information; Continuing with sequence as determined from the 3D "
                          "structure.", stacklevel=2)
            full_record = physical_record

        prot_name = "".join(x for x in full_record.name if x.isalnum())

        # ProSE accepts sequences as .fa files only
        SeqIO.write(full_record, f"data/{prot_name}_{chain_id}.fa", "fasta")

        # Parse the PDB structure data
        self.residues = []
        self.atoms = []
        atom_index = 0
        for i, bio_residue in enumerate(chain.get_residues()):
            seq_index = bio_residue.get_id()[1]
            # Remove any residues not in range (often, these are water molecules)
            if seq_index < physical_record.annotations["start"]:
                continue
            if seq_index >= physical_record.annotations["end"]:
                continue
            residue_atoms = []
            residue = Residue(residue_atoms, bio_residue, i, seq_index)
            for bio_atom in bio_residue.get_atoms():
                residue_atoms.append(Atom(bio_atom, residue, atom_index))
                atom_index += 1
            self.residues.append(residue)
            self.atoms.extend(residue_atoms)

        # Calculate embeddings with ProSE
        if not os.path.isfile(f"data/{prot_name}_{chain_id}_embeddings.h5"):
            print("Embeddings not found in `data` directory. Generating new embeddings.")
            subprocess.call(["python", "prose/embed_sequences.py", '-o', f"data/{prot_name}_{chain_id}_embeddings.h5",
                             f"data/{prot_name}_{chain_id}.fa"], shell=True)
        else:
            print("Embeddings from previous session found in `data` directory.")

        # Parse the embeddings file for the residues we will render
        embedding_file = h5py.File(f"data/{prot_name}_{chain_id}_embeddings.h5", 'r')
        embeddings = embedding_file[list(embedding_file.keys())[0]][()]
        realized_embeddings = []
        for residue in self.residues:
            realized_embeddings.append(embeddings[residue.seq_index])

        # Use the t-SNE algorithm to transform the embeddings into 2D vectors
        transform = TSNE(n_components=2, perplexity=3).fit_transform(np.array(realized_embeddings))
        self.embedding_points = transform.flatten()

        # Calculate residue clusters
        db = DBSCAN(eps=cluster_distance).fit(transform)
        self.cluster_index = db.labels_
        self.cluster_count = len(set(self.cluster_index)) - (1 if -1 in self.cluster_index else 0)

        self.residues[0].highlighted = True
        self.update_colors()

    def update_colors(self, new_color_mode=None, new_color_palette=None):
        """
        Updates the color of the residues within this protein. Slow performance, do not call regularly
        @param new_color_mode: Options: Protein.RESIDUE_INDEX, Protein.CLUSTER_INDEX, Protein.ATOM_TYPE
        @param new_color_palette: Options: Protein.RAINBOW, Protein.POISSON
        """
        if new_color_mode:
            self.color_mode = new_color_mode
        if new_color_palette:
            self.color_palette = new_color_palette
        for residue in self.residues:
            self.color_residue(residue)

    def color_residue(self, residue):
        """
        Colors the given residue (of this protein) using the color settings
        @param residue: A residue of this protein.
        """
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
                residue.color = self.get_color(self.cluster_index[residue.index] / self.cluster_count,
                                               residue.highlighted)
                for atom in residue.atoms:
                    atom.color = self.cpk_colors[atom.bio_atom.get_id()[0]]
                    if residue.highlighted:
                        atom.color = [min(255, c + 128) for c in atom.color]

    def get_color(self, x, highlight=False):
        """
        Internal function to retrieve a color from the current palette
        @param x: [0, 1], different context depending on the current color mode
        @param highlight: Whether this residue is currently highlighted
        @return: An RGB array.
        """
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
                    color = colour.Color(hue=x * 0.8, saturation=0.75, luminance=0.5)
                case self.POISSON:
                    color.rgb = get_color_from_palette(self.poisson_palette)

        if highlight:
            color.set_luminance(min(1.0, color.get_luminance() + 0.25))
        return to_rgb(color)
