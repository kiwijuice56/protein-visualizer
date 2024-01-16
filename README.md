# protein-visualizer
A protein visualization program that uses deep learning to map 3D structure to amino acid function.
See [website](https://kiwijuice56.github.io/protein-visualizer/) for installation and usage details.

## Libraries
- `pyglet 1.5.28` as an OpenGL interface
- `bio 1.6.2` to parse PDB files
- `scikit-learn 1.3.2` for the t-SNE algorithm
- `h5py 3.10.0` to parse h5 database files
- `colour 0.1.5` for convenient color arithmetic
- `numpy 1.26.2`

## Attribution
- Bepler, T., Berger, B. Learning the protein language: evolution, structure, and function. Cell Systems 12, 6 (2021). https://doi.org/10.1016/j.cels.2021.05.017
- Bepler, T., Berger, B. Learning protein sequence embeddings using information from structure. International Conference on Learning Representations (2019). https://openreview.net/pdf?id=SygLehCqtm