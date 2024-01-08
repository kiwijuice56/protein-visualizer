import Bio.PDB
import Bio.SeqRecord

# Pyglet version 1.5.28 (https://pyglet.readthedocs.io/en/pyglet-1.5-maintenance/)
import pyglet
from pyglet.window import key
from pyglet.gl import *

import math

# Camera code from https://gist.github.com/mr-linch/f6dacd2a069887a47fbc
from camera import FirstPersonCamera

pdbparser = Bio.PDB.PDBParser(QUIET=True)
prot = pdbparser.get_structure("example protein", "alphafold-generation.pdb")

atoms = [atom for atom in prot.get_atoms()]
residues = [residue for residue in prot.get_residues()]
residue_index_range = {}

# Determine which atoms belong to each residue
atom_idx = 0
for i, residue in enumerate(prot.get_residues()):
    residue_atoms = [atom for atom in residue.get_atoms()]
    residue_index_range[i] = [atom_idx, atom_idx + len(residue_atoms)]
    atom_idx += len(residue_atoms)

# Initialize the window and camera
win = pyglet.window.Window()
win.set_exclusive_mouse(True)
win.set_fullscreen(True)
cam = FirstPersonCamera(win, movement_speed=16)


def color_atom(atom, highlight=False):
    # Based on CPK coloring (https://en.wikipedia.org/wiki/CPK_coloring)
    colors = {"C": (64, 58, 64), "O": (219, 73, 70), "N": (70, 110, 219), "S": (235, 208, 56), "P": (235, 145, 56),
              "_": (255, 255, 255)}
    if highlight:
        return [min(255, c + 128) for c in colors[atom.get_id()[0]]]
    elif atom.get_id()[0] in colors:
        return colors[atom.get_id()[0]]
    else:
        return colors['_']


# Load data into a format that OpenGL can parse
point_coord = []
draw_color = []
residue_idx = 0

for atom in prot.get_atoms():
    point_coord.extend(atom.get_coord())
    draw_color.extend(color_atom(atom, atom in atoms[residue_index_range[0][0]: residue_index_range[0][1]]))

vertex_list = pyglet.graphics.vertex_list(len(atoms), ('v3f', point_coord), ('c3B', draw_color))


# Add extra inputs for cycling through residues
def on_key_press(symbol, modifiers):
    global residue_idx
    old_idx = residue_idx

    if symbol == key._1:
        residue_idx -= 1
    if symbol == key._2:
        residue_idx += 1

    if residue_idx >= len(residue_index_range):
        residue_idx = 0
    elif residue_idx < 0:
        residue_idx = len(residue_index_range) - 1

    if not old_idx == residue_idx:
        for i in range(residue_index_range[old_idx][0], residue_index_range[old_idx][1]):
            vertex_list.colors[i * 3: i * 3 + 3] = color_atom(atoms[i])
        for i in range(residue_index_range[residue_idx][0], residue_index_range[residue_idx][1]):
            vertex_list.colors[i * 3: i * 3 + 3] = color_atom(atoms[i], True)


win.push_handlers(on_key_press)


@win.event
def on_draw():
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    # Initialize OpenGL context
    # glClearColor(1.0, 1.0, 1.0, 1.0)
    glLoadIdentity()
    glEnable(GL_POINT_SMOOTH)
    glEnable(GL_DEPTH_TEST)
    glMatrixMode(GL_PROJECTION)
    gluPerspective(65, win.width / float(win.height), 0.01, 512)
    glPointSize(8)

    # Draw the protein
    cam.draw()
    vertex_list.draw(pyglet.gl.GL_POINTS)

    # Reset projection to 2D for UI
    glLoadIdentity()
    glOrtho(0, win.width, -win.height, 0, 0, 1000)

    document = pyglet.text.document.FormattedDocument()
    document.insert_text(0, f"{prot.get_id()}\n", {"font_name": "Consolas", "font_size": 16, "color": (255, 255, 255, 255)})

    # Draw residue label
    for i in range(-2, 4):
        adj = residue_idx + i
        if adj < 0:
            adj = len(residues) + adj
        if adj >= len(residues):
            adj -= len(residues)
        snippet = f"{f'%0{int(math.log10(len(residues))) + 1}d' % adj}:{'%3s' % residues[adj].get_resname()} "
        color = tuple([255] * 4 if i == 0 else [64 - abs(i) * 12] * 3 + [255])
        document.insert_text(len(document.text), snippet, {"color": color})

    layout = pyglet.text.layout.TextLayout(document, multiline=True, width=win.width, height=win.height)
    layout.anchor_x = "left"
    layout.anchor_y = "top"

    layout.draw()

    return pyglet.event.EVENT_HANDLED


# Main update loop
def on_update(delta_time):
    cam.update(delta_time)


if __name__ == "__main__":
    pyglet.clock.schedule(on_update)
    pyglet.app.run()
