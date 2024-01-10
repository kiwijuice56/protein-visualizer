# protein-visualizer
A program to visualize how amino acid residues contribute to protein structure and function.
Work in progress.

![Protein rendered using this program](picture_demonstration.png)

## Controls
- `WASD`: Translate the camera
- `Shift/Space`: Translate the camera vertically
- `Mouse`: Rotate the camera
- `Left/Right Arrow`: Shift which residue is currently highlighted
- `Up/Down Arrow`: Increase and decrease atom point size
- `O`: Toggles atom outline
- `Esc`: Close the program

## Libraries
- `pyglet 1.5.28` as an OpenGL interface
- `bio 1.6.2` to parse PDB files
- `scikit-learn 1.3.2` for the t-SNE algorithm
- `h5py` to parse h5 database files
- `colour` for convenient color arithmetic
- `camera.py` script from https://gist.github.com/mr-linch/f6dacd2a069887a47fbc