# PDB-Renderer
A simple program to display Protein Data Bank (PDB) protein files in 3D. Run `pdb-renderer.py` after installing the 
dependencies. Sample proteins in `proteins` are sourced from RCSB PDB (except the generation, created by myself
using AlphaFold 2)

![Protein rendered using this program](picture_demonstration.png)

## Controls
- `WASD`: Translate the camera
- `Shift/Space`: Translate the camera vertically
- `Mouse`: Rotate the camera
- `1/2`: Shift which residue is currently highlighted
- `Esc`: Close the program

## Libraries
- `pyglet 1.5.28` as an OpenGL interface
- `bio 1.6.2` to parse PDB files
- `camera.py` script from https://gist.github.com/mr-linch/f6dacd2a069887a47fbc