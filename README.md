# PDB-Renderer
A simple program to display Protein Data Bank (PDB) protein files in 3D.

![Protein rendered using this program](picture_demonstration.png)

## Controls
- WASD: Translate the camera
- Shift/Space: Translate the camera vertically
- Mouse: Rotate the camera
- 1/2: Shift which residue is currently highlighted

## Libraries
- `pyglet 1.5.28` as an OpenGL interface
- `bio 1.6.2` to parse PDB files
- `camera.py` script from https://gist.github.com/mr-linch/f6dacd2a069887a47fbc