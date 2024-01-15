# protein-visualizer
A protein visualization program that uses deep learning to map 3D structure to amino acid function.

## Usage
1. Download this repository
2. Download and configure the pretrained models in `prose/` according to the 
instructions on https://github.com/tbepler/prose.
3. Navigate to the main directory and run: `python main.py [path to .pdb file]`

## Controls
- `Left Mouse Drag`: Move the camera in both 3D and 2D space
- `Right Mouse Drag`: Rotate the camera in 3D space
- `Left/Right Arrow`: Shift which residue is currently highlighted
- `Up/Down Arrow`: Increase and decrease atom point size
- `1/2/3`: Changes color mode
- `8/9` Changes color palette
- `O`: Toggles atom outline
- `Esc`: Close the program

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