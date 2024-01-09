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
prot = pdbparser.get_structure("example protein", "proteins/alphafold-generation.pdb")

atoms = [atom for atom in prot.get_atoms()]
residues = [residue for residue in prot.get_residues()]
residue_index_range = {}
atom_residue_map = {}

# Determine which atoms belong to each residue
atom_idx = 0
for i, residue in enumerate(prot.get_residues()):
    for atom in residue.get_atoms():
        atom_residue_map[atom] = i
    residue_atoms = [atom for atom in residue.get_atoms()]
    residue_index_range[i] = [atom_idx, atom_idx + len(residue_atoms)]
    atom_idx += len(residue_atoms)

# Initialize the window and camera
win = pyglet.window.Window()
win.set_exclusive_mouse(True)
win.set_fullscreen(True)
cam = FirstPersonCamera(win, movement_speed=16)

residue_idx = 0
point_size = 8
color_mode = "atomic"


def color_atom(atom, highlight=False, mode="atomic"):
    # Based on CPK coloring (https://en.wikipedia.org/wiki/CPK_coloring)
    atomic_colors = {"C": (64, 58, 64), "O": (219, 73, 70), "N": (70, 110, 219), "S": (235, 208, 56),
                     "P": (235, 145, 56),
                     "_": (255, 255, 255)}
    residue_colors = [(191, 53, 15), (204, 135, 24), (212, 209, 25), (110, 204, 33), (33, 204, 130), (33, 178, 204),
                      (36, 51, 212), (133, 36, 212), (212, 36, 145)]
    base_color = (128, 128, 128)

    if mode == "atomic":
        base_color = atomic_colors[atom.get_id()[0]] if atom.get_id()[0] in atomic_colors else atomic_colors['_']
    elif mode == "residue":
        base_color = residue_colors[atom_residue_map[atom] % len(residue_colors)]
    elif mode == "contrast":
        base_color = (255, 213, 25) if highlight else (46, 29, 115)

    if highlight:
        return tuple([min(255, c + 128) for c in base_color])
    else:
        return base_color


# Load data into a format that OpenGL can parse
point_coord = []
draw_color = []

first_residue = atoms[residue_index_range[0][0]: residue_index_range[0][1]]
for atom in prot.get_atoms():
    point_coord.extend(atom.get_coord())
    draw_color.extend(color_atom(atom, atom in first_residue, mode=color_mode))

vertex_list = pyglet.graphics.vertex_list(len(atoms), ('v3f', point_coord), ('c3B', draw_color))


# Implement inputs for controlling the UI
def on_key_press(symbol, modifiers):
    global color_mode
    global residue_idx
    global point_size

    old_idx = residue_idx

    if symbol == key.LEFT:
        residue_idx -= 1
    if symbol == key.RIGHT:
        residue_idx += 1

    if symbol == key.UP:
        point_size += 1
    if symbol == key.DOWN:
        point_size -= 1
    point_size = max(1, min(10, point_size))

    old_mode = color_mode
    if symbol == key._1:
        color_mode = "atomic"
    if symbol == key._2:
        color_mode = "residue"
    if symbol == key._3:
        color_mode = "contrast"
    if not color_mode == old_mode:
        for i in range(0, len(atoms)):
            draw_color[i * 3: i * 3 + 3] = color_atom(atoms[i], mode=color_mode)
        for i in range(residue_index_range[residue_idx][0], residue_index_range[residue_idx][1]):
            draw_color[i * 3: i * 3 + 3] = color_atom(atoms[i], highlight=True, mode=color_mode)

    if residue_idx >= len(residue_index_range):
        residue_idx = 0
    elif residue_idx < 0:
        residue_idx = len(residue_index_range) - 1

    if not old_idx == residue_idx:
        for i in range(residue_index_range[old_idx][0], residue_index_range[old_idx][1]):
            draw_color[i * 3: i * 3 + 3] = color_atom(atoms[i], mode=color_mode)
        for i in range(residue_index_range[residue_idx][0], residue_index_range[residue_idx][1]):
            draw_color[i * 3: i * 3 + 3] = color_atom(atoms[i], highlight=True, mode=color_mode)


win.push_handlers(on_key_press)


@win.event
def on_draw():
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    # Initialize OpenGL context
    glClearColor(1.0, 1.0, 1.0, 1.0)
    glLoadIdentity()
    glEnable(GL_POINT_SMOOTH)
    glDisable(GL_DEPTH_TEST)
    glMatrixMode(GL_PROJECTION)
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE)
    gluPerspective(65, win.width / float(win.height), 0.5, 512)

    # Draw the protein
    cam.draw()

    glPointSize(point_size * 2)
    vertex_list.colors = [32] * len(vertex_list.colors)
    vertex_list.draw(pyglet.gl.GL_POINTS)

    glPointSize(point_size)
    glEnable(GL_DEPTH_TEST)
    vertex_list.colors = draw_color
    vertex_list.draw(pyglet.gl.GL_POINTS)

    # Reset projection to 2D for UI
    glLoadIdentity()
    glDisable(GL_DEPTH_TEST)
    glEnable(GL_BLEND)
    glOrtho(0, win.width, -win.height, 0, 0, 1000)

    document = pyglet.text.document.FormattedDocument()
    document.insert_text(0, f"{prot.get_id()}\n",
                         {"font_name": "Consolas", "font_size": 16, "background_color": (0, 0, 0, 180),
                          "color": (255, 255, 255, 255)})

    # Draw residue label
    for i in range(-2, 4):
        adj = residue_idx + i
        if adj < 0:
            adj = len(residues) + adj
        if adj >= len(residues):
            adj -= len(residues)
        snippet = f"{f'%0{int(math.log10(len(residues))) + 1}d' % adj}:{'%3s' % residues[adj].get_resname()} "
        color = tuple((255, 230, 0, 255) if i == 0 else [200 - abs(i) * 24] * 3 + [255])
        document.insert_text(len(document.text), snippet, {"color": color})
    document.insert_text(len(document.text), "\nhttps://github.com/kiwijuice56/pdb-renderer", {"color": (255, 255, 255, 255)})

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
