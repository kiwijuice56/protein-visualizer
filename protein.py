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
        self.go_map = {"GO:-------": 1.0}


# Contains information about a protein, such as its 3D structure and node embeddings
class Protein:
    # Color modes
    CLUSTER_INDEX = 1  # Color residues by their associated cluster in the embedding space
    RESIDUE_INDEX = 2  # Color residues by their index in the amino acid sequence
    ATOM_TYPE = 3  # Color residues by their atoms using CPK coloring
    RESIDUE_TYPE = 4 # Colors by identity of each amino acid

    # Color palettes (for color modes 1 and 2)
    RAINBOW = 6
    MONOCOLOR = 7
    POISSON = 8
    PENGUIN = 9
    GRAPE = 10
    LEMON = 11
    MULBERRY = 12

    POISSON_PALETTE = [(187, 176, 148), (128, 118, 101), (89, 82, 70), (51, 51, 51), (25, 31, 34), (47, 68, 67),
                       (59, 94, 88), (90, 140, 108), (139, 180, 141), (192, 208, 165), (247, 239, 199),
                       (161, 205, 176), (112, 147, 149), (74, 120, 123), (56, 49, 64), (115, 77, 92),
                       (167, 103, 114), (204, 134, 125), (224, 186, 139), (195, 130, 82), (161, 86, 60),
                       (111, 52, 45), (68, 39, 31)]

    GRAPE_PALETTE = [(3, 6, 55), (60, 7, 83), (114, 4, 85), (145, 10, 103), (194, 48, 131)]

    MULBERRY_PALETTE = [(197, 235, 195), (183, 200, 181), (167, 144, 165), (135, 92, 116), (84, 65, 78), (177, 133, 167), (195, 162, 158), (232, 218, 197), (255, 244, 233)]

    MONOCOLOR_PALETTE = [(59, 212, 59)]

    PENGUIN_PALETTE = [(43, 48, 58), (146, 220, 229), (238, 229, 233), (124, 124, 124), (214, 73, 51)]

    LEMON_PALETTE = [(191, 174, 72), (95, 173, 65), (45, 147, 108), (57, 20, 99), (58, 8, 66)]

    CPK_COLORS = {"C": (64, 58, 64), "O": (219, 73, 70), "N": (70, 110, 219), "S": (235, 208, 56),
                  "P": (235, 145, 56), "_": (255, 255, 255)}

    # Codes to color the terminal text
    OUTPUT_COLOR_GOOD = "\033[96m"
    OUTPUT_COLOR_WORKING = "\u001b[33m"

    MAX_OUTLINE_BRIGHTNESS = 90  # 0 to 255
    HIGHLIGHT_LUMINANCE = 0.35  # 0 to 1

    def __init__(self, protein_path, prompt_for_chain=True, chain_id=None, verbose=False):
        """
        @param protein_path: Filepath to .pdb or .cif file containing the protein structure
        @param chain_id: The specific protein chain to draw as a string
        @param prompt_for_chain: Whether the program should ask for chain input or simply default to the first
        @param verbose: Whether extraneous output such as internal warnings are printed
        """
        self.color_mode = self.CLUSTER_INDEX
        self.color_palette = self.RAINBOW

        if not verbose:
            filterwarnings("ignore")

        # Get the full 3D structure from the protein file
        if str(protein_path).endswith(".cif"):
            from Bio.PDB import MMCIFParser
            bio_structure = MMCIFParser(QUIET=not verbose).get_structure("struct", protein_path)
        else:
            from Bio.PDB import PDBParser
            bio_structure = PDBParser(QUIET=not verbose).get_structure("struct", protein_path)

        # Parse the protein name from the file name
        protein_name = "".join(x for x in protein_path[protein_path.rfind("/"):-4] if x.isalnum() or x in "_-").lower()

        print(self.OUTPUT_COLOR_GOOD + f"Loading {protein_name} structure.")

        from Bio.SeqIO import parse
        records = {r.annotations["chain"]: r for r in parse(protein_path, "cif-atom" if str(protein_path).endswith(".cif") else "pdb-atom")}
        if chain_id is None:
            if len(records) > 1 and prompt_for_chain:
                print(f"Multiple chains found: {list(records.keys())}. Please select one by typing its name.")
                chain_id = input()
            else:
                chain_id = list(records.keys())[0]
        protein_name += '_' + chain_id

        # Retrieve the requested chain from the .pdb file
        chain = None
        for bio_chain in [c for c in bio_structure.get_chains()]:
            if bio_chain.get_id() == chain_id:
                chain = bio_chain
                break

        physical_record = records[chain_id]

        print(self.OUTPUT_COLOR_GOOD + f"Loading chain with ID '{chain_id}'.")

        for other_record in [r for r in parse(protein_path, "cif-atom" if str(protein_path).endswith(".cif") else "pdb-atom")]:
            if other_record.annotations["chain"] == chain_id:
                physical_record = other_record
        self.sequence = ''.join([r for r in str(physical_record.seq) if r not in "-*X"])

        # Parse the .pdb file for structure data
        self.residues = []
        self.atoms = []
        self.residue_map = {}

        atom_index = 0
        res_index = 0
        for bio_residue in chain.get_residues():
            true_index = bio_residue.get_id()[1]
            # Remove any residues not in range (often, these are water molecules)
            if true_index < physical_record.annotations["start"]:
                continue
            if true_index > physical_record.annotations["end"]:
                continue
            residue_atoms = []
            residue = Residue(residue_atoms, bio_residue, res_index)
            for bio_atom in bio_residue.get_atoms():
                residue_atoms.append(Atom(bio_atom, residue, atom_index))
                atom_index += 1
            self.residues.append(residue)
            self.atoms.extend(residue_atoms)
            if not bio_residue.get_resname() in self.residue_map:
                self.residue_map[bio_residue.get_resname()] = len(self.residue_map)
            res_index += 1

        if not isfile(f"data/{protein_name}_data.json"):
            print(self.OUTPUT_COLOR_WORKING + "Cache for this protein not found in the 'data' directory.")
            print(self.OUTPUT_COLOR_WORKING + "Generating embeddings.")

            data = self.generate_embeddings(self.sequence)

            print(self.OUTPUT_COLOR_WORKING + "Generating contact map.")
            contact_map = self.generate_contact_map(self.residues)
            print(self.OUTPUT_COLOR_WORKING + "Generating GO annotations.")

            annotations = self.generate_go_annotations(self.sequence, contact_map)
            if "query_prot" in annotations:
                data.update(annotations["query_prot"])
            else:
                data.update({"GO_ids": [], "GO_names": [], "confidence": []})

            with open(f"data/{protein_name}_data.json", 'w') as f:
                dump(data, f, indent=1)
        else:
            print(self.OUTPUT_COLOR_GOOD + "Cache for this protein was found in the 'data' directory.")

        with open(f"data/{protein_name}_data.json", 'r') as f:
            data = load(f)
            self.embedding_points = data["embedding_points"]
            self.cluster_index = data["cluster_indices"]
            self.cluster_count = len(set(self.cluster_index)) - (1 if -1 in self.cluster_index else 0)
            self.go_ids = data["GO_ids"]
            self.go_names = data["GO_names"]
            self.scores = data["confidence"]

            self.go_ids.insert(0, "GO:-------")
            self.go_names.insert(0, "")
            self.scores.insert(0, 0.0)

            self.current_go_id = self.go_ids[0]

            # Assign saliency to each residue in the protein sequence
            for i, annotation in enumerate(self.go_ids):
                if i == 0:
                    continue
                for residue in self.residues:
                    residue.go_map[annotation] = data["saliency_maps"][i - 1][residue.index]

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
                        new_color.rgb = get_color_from_palette(self.POISSON_PALETTE)
                    case self.GRAPE:
                        new_color.rgb = get_color_from_palette(self.GRAPE_PALETTE)
                    case self.MONOCOLOR:
                        new_color.rgb = get_color_from_palette(self.MONOCOLOR_PALETTE)
                    case self.PENGUIN:
                        new_color.rgb = get_color_from_palette(self.PENGUIN_PALETTE)
                    case self.LEMON:
                        new_color.rgb = get_color_from_palette(self.LEMON_PALETTE)
                    case self.MULBERRY:
                        new_color.rgb = get_color_from_palette(self.MULBERRY_PALETTE)
                new_color.set_luminance(new_color.get_luminance() * luminance)
            if residue.highlighted:
                new_color.set_luminance(min(1.0, new_color.get_luminance() + self.HIGHLIGHT_LUMINANCE))
            return [int(b * 255) for b in new_color.rgb]

        luminance = residue.go_map[self.current_go_id]
        luminance *= luminance
        luminance = min(1.0, luminance + 0.1)

        residue.outline_color = [int(luminance * self.MAX_OUTLINE_BRIGHTNESS)] * 3

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
                    color.rgb = [c / 255.0 for c in self.CPK_COLORS[atom.bio_atom.get_id()[0]]]
                    color.set_luminance(color.get_luminance() * luminance)
                    if residue.highlighted:
                        color.set_luminance(min(1.0, color.get_luminance() + self.HIGHLIGHT_LUMINANCE))
                    atom.color = [int(c * 255) for c in color.rgb]
            case self.RESIDUE_TYPE:
                color = get_color(self.residue_map[residue.bio_residue.get_resname()] / len(self.residue_map))
                for atom in residue.atoms:
                    atom.color = color
                residue.color = color

    # Adapted from https://github.com/tbepler/prose
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
