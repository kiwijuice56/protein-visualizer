import numpy as np
import h5py

import torch

from prose import fasta
from prose.alphabets import Uniprot21


# Taken from https://github.com/tbepler/prose
def embed_sequence(model, x):
    if len(x) == 0:
        n = model.embedding.proj.weight.size(1)
        z = np.zeros((1, n), dtype=np.float32)
        return z

    alphabet = Uniprot21()
    x = x.upper()

    # Convert to alphabet index
    x = alphabet.encode(x)
    x = torch.from_numpy(x)

    # Embed the sequence
    with torch.no_grad():
        x = x.long().unsqueeze(0)
        z = model.transform(x)
        z = z.squeeze()
        z = z.cpu().numpy()

    return z


def generate_embeddings(output_path, input_path):
    from prose.models.multitask import ProSEMT
    model = ProSEMT.load_pretrained()
    model.eval()

    # Parse the sequences and embed them
    # Write them to hdf5 file
    h5 = h5py.File(output_path, 'w')

    with open(input_path, 'rb') as f:
        for name, sequence in fasta.parse_stream(f):
            pid = name.decode('utf-8')
            z = embed_sequence(model, sequence)
            # Write as hdf5 dataset
            h5.create_dataset(pid, data=z)
