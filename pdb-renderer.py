import Bio.PDB
import Bio.SeqRecord

# Pyglet version 1.5.x (https://pyglet.readthedocs.io/en/pyglet-1.5-maintenance/)
import pyglet
from pyglet.window import key
from pyglet.gl import *

# Camera code from https://gist.github.com/mr-linch/f6dacd2a069887a47fbc
from camera import FirstPersonCamera

pdbparser = Bio.PDB.PDBParser(QUIET=True)
prot = pdbparser.get_structure("alfaro_worm", "alfaro_worm.pdb")
atoms = [atom for atom in prot.get_atoms()]
amino_acids = {}

atom_idx = 0
for i, amino_acid in enumerate(prot.get_residues()):
    residue_atoms = [atom for atom in amino_acid.get_atoms()]
    amino_acids[i] = [atom_idx, atom_idx + len(residue_atoms)]
    atom_idx += len(residue_atoms)


win = pyglet.window.Window()
win.set_exclusive_mouse(True)
cam = FirstPersonCamera(win, movement_speed=16)

# Based on CPK coloring (https://en.wikipedia.org/wiki/CPK_coloring)
colors = {"C": (45, 45, 45), "O": (219, 73, 70), "N": (70, 110, 219), "S": (235, 208, 56),
          "P": (235, 145, 56),
          "_": (100, 100, 100), "*": (246, 77, 255)}

# Load atom data from the pdb file (coordinates and color based on atom type)
points = []
draw_color = []
residue_idx = 0
vertex_list = None  # Final data to be sent to OpenGL Renderer


def color_atom(atom):
    if atom.get_id()[0] in colors:
        return colors[atom.get_id()[0]]
    else:
        draw_color.extend(colors["_"])
        return colors["_"]


for atom in prot.get_atoms():
    points.extend(atom.get_coord())
    draw_color.extend(color_atom(atom))

vertex_list = pyglet.graphics.vertex_list(len(points) // 3, ('v3f', points), ('c3B', draw_color))


def on_key_press(symbol, modifiers):
    global residue_idx
    old_idx = residue_idx

    if symbol == key._1:
        residue_idx += 1
    if symbol == key._2:
        residue_idx -= 1

    if not old_idx == residue_idx:
        for i in range(amino_acids[old_idx][0], amino_acids[old_idx][1]):
            color = color_atom(atoms[i])
            vertex_list.colors[i * 3] = color[0]
            vertex_list.colors[i * 3 + 1] = color[1]
            vertex_list.colors[i * 3 + 2] = color[2]
        for i in range(amino_acids[residue_idx][0], amino_acids[residue_idx][1]):
            color = colors['*']
            vertex_list.colors[i * 3] = color[0]
            vertex_list.colors[i * 3 + 1] = color[1]
            vertex_list.colors[i * 3 + 2] = color[2]


win.push_handlers(on_key_press)


@win.event
def on_draw():
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    # Initialize OpenGL context
    glClearColor(1.0, 1.0, 1.0, 1.0)
    glLoadIdentity()
    glEnable(GL_POINT_SMOOTH)
    glEnable(GL_DEPTH_TEST)
    glMatrixMode(GL_PROJECTION)
    gluPerspective(65, win.width / float(win.height), 0.1, 1000)
    glPointSize(5)

    # Draw the protein
    cam.draw()

    vertex_list.draw(pyglet.gl.GL_POINTS)

    return pyglet.event.EVENT_HANDLED


# Main update loop
def on_update(delta_time):
    cam.update(delta_time)


if __name__ == "__main__":
    pyglet.clock.schedule(on_update)
    pyglet.app.run()
