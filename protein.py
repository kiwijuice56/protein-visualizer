from os.path import isfile
from warnings import filterwarnings
from json import load, dump

from colour import Color
import numpy as np


# Wrapper class for single atoms
class Atom:
    def __init__(self, bio_atom, residue, index):
        """
        @param bio_atom: A Bio.PDB.Atom object
        @param residue: The Residue object that owns this atom
        @param index: The index of this atom out of all atoms in the ordered protein sequence
        """
        self.bio_atom = bio_atom
        self.residue = residue
        self.index = index

        self.color = [0, 0, 0]
        self.outline_color = [32, 32, 32]


# Wrapper class for collections of atoms
class Residue:
    def __init__(self, atoms, bio_residue, index):
        """
        @param atoms: A list of Atom objects that make up this residue
        @param bio_residue: A Bio.PDB.Residue object
        @param index: The index of this residue out of all residues in the ordered protein sequence
        """
        self.atoms = atoms
        self.bio_residue = bio_residue
        self.index = index

        self.color = [0, 0, 0]
        self.outline_color = [32, 32, 32]
        self.highlighted = False

        # A string [GO id] : float [0.0, 1.0] pair of how strongly this residue contributed to a certain GO prediction
        self.go_map = {}


# Contains information about a protein, such as its 3D structure and node embeddings
class Protein:
    # Color modes
    CLUSTER_INDEX = 1  # Color residues by their associated cluster in the embedding space
    RESIDUE_INDEX = 2  # Color residues by their index in the amino acid sequence
    ATOM_TYPE = 3  # Color residues by their atoms using CPK coloring

    # Color palettes (for color modes 1 and 2)
    RAINBOW = 6
    MONOCOLOR = 7
    POISSON = 8
    GRAPE = 9
    PENGUIN = 10
    LEMON = 11

    poisson_palette = [(187, 176, 148), (128, 118, 101), (89, 82, 70), (51, 51, 51), (25, 31, 34), (47, 68, 67),
                       (59, 94, 88), (90, 140, 108), (139, 180, 141), (192, 208, 165), (247, 239, 199),
                       (161, 205, 176), (112, 147, 149), (74, 120, 123), (56, 49, 64), (115, 77, 92),
                       (167, 103, 114), (204, 134, 125), (224, 186, 139), (195, 130, 82), (161, 86, 60),
                       (111, 52, 45), (68, 39, 31)]

    grape_palette = [(3, 6, 55), (60, 7, 83), (114, 4, 85), (145, 10, 103), (194, 48, 131)]

    monocolor_palette = [(59, 212, 59)]

    penguin_palette = [(43, 48, 58), (146, 220, 229), (238, 229, 233), (124, 124, 124), (214, 73, 51)]

    lemon_palette = [(191, 174, 72), (95, 173, 65), (45, 147, 108), (57, 20, 99), (58, 8, 66)]

    cpk_colors = {"C": (64, 58, 64), "O": (219, 73, 70), "N": (70, 110, 219), "S": (235, 208, 56),
                  "P": (235, 145, 56), "_": (255, 255, 255)}

    # Codes to color the terminal text
    output_color = "\033[96m"
    output_color2 = "\u001b[33m"

    def __init__(self, pdb_path, chain_id=None, verbose=False):
        """
        @param pdb_path: Filepath to .pdb file containing the protein structure
        @param chain_id: Which chain to load from the .pdb file. Defaults to first chain
        @param verbose: Whether extraneous output such as internal warnings are printed
        """
        self.color_mode = self.CLUSTER_INDEX
        self.color_palette = self.RAINBOW

        if not verbose:
            filterwarnings("ignore")

        print(self.output_color + "Loading protein structure.")

        # Get the full 3D structure from the PDB file
        from Bio.PDB import PDBParser
        bio_structure = PDBParser(QUIET=not verbose).get_structure("struct", pdb_path)

        # Parse the protein name from the file name
        protein_name = "".join(x for x in pdb_path[pdb_path.rfind("/"):-4] if x.isalnum() or x in "_-").lower()

        # Choose default chain if none specified
        chains = [c for c in bio_structure.get_chains()]
        if chain_id is None:
            chain_id = chains[0].get_id()

        # Retrieve the requested chain from the .pdb file
        chain = None
        for other_chain in chains:
            if other_chain.get_id() == chain_id:
                chain = other_chain
                break

        print(self.output_color + f"Rendering chain with ID '{chain_id}' from {[chain.get_id() for chain in chains]}")

        # Parse the .pdb file for the protein sequence
        physical_record = None
        from Bio.SeqIO import parse
        for other_record in [r for r in parse(pdb_path, "pdb-atom")]:
            if other_record.annotations["chain"] == chain_id:
                physical_record = other_record
        self.sequence = ''.join([r for r in str(physical_record.seq) if r not in "-*X"])

        # Parse the .pdb file for structure data
        self.residues = []
        self.atoms = []

        atom_index = 0
        for i, bio_residue in enumerate(chain.get_residues()):
            true_index = bio_residue.get_id()[1]
            # Remove any residues not in range (often, these are water molecules)
            if true_index < physical_record.annotations["start"]:
                continue
            if true_index > physical_record.annotations["end"]:
                continue
            residue_atoms = []
            residue = Residue(residue_atoms, bio_residue, i)
            for bio_atom in bio_residue.get_atoms():
                residue_atoms.append(Atom(bio_atom, residue, atom_index))
                atom_index += 1
            self.residues.append(residue)
            self.atoms.extend(residue_atoms)

        if not isfile(f"data/{protein_name}_data.json"):
            print(self.output_color2 + "Cache for this protein not found in the 'data' directory.")
            print(self.output_color2 + "Generating embeddings.")

            data = self.generate_embeddings(self.sequence)

            print(self.output_color2 + "Generating contact map.")
            contact_map = self.generate_contact_map(self.residues)
            print(self.output_color2 + "Generating GO annotations.")

            data.update(self.generate_go_annotations(self.sequence, contact_map)["query_prot"])

            with open(f"data/{protein_name}_data.json", 'w') as f:
                dump(data, f, indent=1)
        else:
            print(self.output_color + "Cache for this protein was found in the 'data' directory.")

        with open(f"data/{protein_name}_data.json", 'r') as f:
            data = load(f)
            self.embedding_points = data["embedding_points"]
            self.cluster_index = data["cluster_indices"]
            self.cluster_count = len(set(self.cluster_index)) - (1 if -1 in self.cluster_index else 0)
            self.go_ids = data["GO_ids"]
            self.go_names = data["GO_names"]
            self.scores = data["confidence"]
            self.current_go_id = self.go_ids[0]

            # Assign saliency to each residue in the protein sequence
            for i, annotation in enumerate(self.go_ids):
                for residue in self.residues:
                    residue.go_map[annotation] = data["saliency_maps"][i][residue.index]

        print(self.output_color + "All calculations complete. Preparing to render.")
        self.update_colors()

    def update_colors(self, new_color_mode=None, new_color_palette=None):
        """
        Updates the color of the residues within this protein. Slow performance, do not call regularly
        @param new_color_mode: Options: Any of the Protein color modes
        @param new_color_palette: Options: Any of the Protein color palettes
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

        def get_color(x):
            """
            Internal function to retrieve a color from the current palette
            @param x: [0, 1], different context depending on the current color mode
            @return: An RGB array.
            """

            def get_color_from_palette(palette):
                step = 1.0 / len(palette)
                i = int(x * len(palette))
                x_i, x_j = step * i, step * (i + 1)
                z = (x - x_i) / (x_j - x_i)
                a, b = palette[i % len(palette)], palette[(i + 1) % len(palette)]
                return [(b[j] * z + a[j] * (1.0 - z)) / 255.0 for j in range(3)]

            new_color = Color()
            if x >= 0:
                match self.color_palette:
                    case self.RAINBOW:
                        new_color = Color(hue=x * 0.7, saturation=0.6, luminance=0.5)
                    case self.POISSON:
                        new_color.rgb = get_color_from_palette(self.poisson_palette)
                    case self.GRAPE:
                        new_color.rgb = get_color_from_palette(self.grape_palette)
                    case self.MONOCOLOR:
                        new_color.rgb = get_color_from_palette(self.monocolor_palette)
                    case self.PENGUIN:
                        new_color.rgb = get_color_from_palette(self.penguin_palette)
                    case self.LEMON:
                        new_color.rgb = get_color_from_palette(self.lemon_palette)
                new_color.set_luminance(new_color.get_luminance() * luminance)
            if residue.highlighted:
                new_color.set_luminance(min(1.0, new_color.get_luminance() + 0.35))
            return [int(b * 255) for b in new_color.rgb]

        luminance = residue.go_map[self.current_go_id]
        luminance *= luminance

        residue.outline_color = [int(luminance * 86)] * 3

        for atom in residue.atoms:
            atom.outline_color = residue.outline_color

        match self.color_mode:
            case self.RESIDUE_INDEX:
                color = get_color(residue.index / len(self.residues))

                for atom in residue.atoms:
                    atom.color = color
                residue.color = color
            case self.CLUSTER_INDEX:
                color = get_color(self.cluster_index[residue.index] / self.cluster_count)
                for atom in residue.atoms:
                    atom.color = color
                residue.color = color
            case self.ATOM_TYPE:
                residue.color = get_color(self.cluster_index[residue.index] / self.cluster_count)

                for atom in residue.atoms:
                    color = Color()
                    color.rgb = [c / 255.0 for c in self.cpk_colors[atom.bio_atom.get_id()[0]]]
                    color.set_luminance(color.get_luminance() * luminance)
                    if residue.highlighted:
                        color.set_luminance(min(1.0, color.get_luminance() + 0.35))
                    atom.color = [int(c * 255) for c in color.rgb]

    # Taken from https://github.com/tbepler/prose
    @staticmethod
    def generate_embeddings(sequence):
        """
        Use the ProSE model to generate embeddings
        @param sequence: The amino acid sequence (with FASTA amino acid names) as a string
        @return A dictionary with the embedding points and cluster indices
        """
        from torch import from_numpy, no_grad
        from sklearn.cluster import DBSCAN
        from sklearn.manifold import TSNE

        from deep_learning.prose.alphabets import Uniprot21
        from deep_learning.prose.models.multitask import ProSEMT

        def embed_sequence(x):
            if len(x) == 0:
                n = model.embedding.proj.weight.size(1)
                y = np.zeros((1, n), dtype=np.float32)
                return y

            alphabet = Uniprot21()
            x = x.upper()

            # Convert to alphabet index
            x = alphabet.encode(x)
            x = from_numpy(x)

            # Embed the sequence
            with no_grad():
                x = x.long().unsqueeze(0)
                y = model.transform(x)
                y = y.squeeze()
                y = y.cpu().numpy()

            return y

        model = ProSEMT.load_pretrained()
        model.eval()

        generated_embeddings = embed_sequence(sequence.encode("utf-8"))

        # Use the t-SNE algorithm to transform the embeddings into 2D vectors
        transform = TSNE(n_components=2, perplexity=30).fit_transform(np.array(generated_embeddings))
        embedding_points = transform.flatten()

        # Use the DBSCAN algorithm to locate clusters in the embedding space
        db = DBSCAN(eps=2.8).fit(transform)
        cluster_index = db.labels_

        return {"embedding_points": [float(x) for x in embedding_points],
                "cluster_indices": [int(x) for x in cluster_index]}

    @staticmethod
    def generate_go_annotations(sequence, contact_map):
        """
        Use the DeepFRI model to generate GO annotations
        @param contact_map: The NxN numpy matrix containing the distance between each residue
        @param sequence: The amino acid sequence (with FASTA amino acid names) as a string
        @return A dictionary with the predicted GO annotations, IDs, confidence levels, and saliency
        """
        from deep_learning.deepfrier.predictor import Predictor

        with open("saved_models/model_config.json") as f:
            params = load(f)

        params = params["gcn"]
        gcn = params["gcn"]
        layer_name = params["layer_name"]
        models = params["models"]

        predictor = Predictor(models['mf'], gcn=gcn)
        predictor.predict(contact_map, sequence)

        predictor.compute_GradCAM(layer_name=layer_name, use_guided_grads=False)
        return predictor.pdb2cam

    # https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/
    @staticmethod
    def generate_contact_map(residues):
        """
        Generate a C-Alpha contact map given a list of residues
        @param residues: The list to create the contact map for
        @return: An NxN numpy matrix of the distance between each residue
        """
        contact_map = np.full((len(residues), len(residues)), 0, float)
        for row in range(len(residues)):
            for col in range(row, len(residues)):
                bio_residue1, bio_residue2 = residues[row].bio_residue, residues[col].bio_residue
                if bio_residue1['CA'] is None or bio_residue2['CA'] is None:
                    diff_vector = bio_residue1.get_atom().coord - bio_residue2.get_atom().coord
                else:
                    diff_vector = bio_residue1['CA'].coord - bio_residue2['CA'].coord
                distance = np.sqrt(np.sum(diff_vector * diff_vector))
                contact_map[row, col] = contact_map[col, row] = distance
        return contact_map
