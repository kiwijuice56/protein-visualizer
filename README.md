# PDB-Renderer
A simple program to display Protein Data Bank (PDB) protein files in 3D. See `main.py` for example usage. Sample proteins in `proteins` are sourced from RCSB PDB (except the generation, created by myself
using AlphaFold 2)

![Protein rendered using this program](picture_demonstration.png)

## Controls
- `WASD`: Translate the camera
- `Shift/Space`: Translate the camera vertically
- `Mouse`: Rotate the camera
- `Left/Right Arrow`: Shift which residue is currently highlighted
- `Up/Down Arrow`: Increase and decrease atom point size
- `1`: Set color mode to `atomic` in which atoms are colored according to CPK guidelines
- `2`: Set color mode to `residue` in which atoms are colored according to which residue they belong to
- `3`: Set color mode to `contrast` in which highlighted atoms are much brighter than the surrounding atoms
- `Esc`: Close the program

## Libraries
- `pyglet 1.5.28` as an OpenGL interface
- `bio 1.6.2` to parse PDB files
- `camera.py` script from https://gist.github.com/mr-linch/f6dacd2a069887a47fbc