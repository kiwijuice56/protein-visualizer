import Bio.PDB
import Bio.SeqRecord

# Pyglet version 1.5.x (https://pyglet.readthedocs.io/en/pyglet-1.5-maintenance/)
import pyglet
from pyglet.gl import *

# Camera code from https://gist.github.com/mr-linch/f6dacd2a069887a47fbc
from camera import FirstPersonCamera

pdbparser = Bio.PDB.PDBParser(QUIET=True)
prot = pdbparser.get_structure("alfaro_worm", "alfaro_worm.pdb")

win = pyglet.window.Window()
win.set_exclusive_mouse(True)
cam = FirstPersonCamera(win, movement_speed=16)

# Based on CPK coloring (https://en.wikipedia.org/wiki/CPK_coloring)
colors = {"C": (45, 45, 45), "O": (219, 73, 70), "N": (70, 110, 219), "S": (235, 208, 56),
          "P": (235, 145, 56), "_": (100, 100, 100)}

# Load atom data from the pdb file (coordinates and color based on atom type)
points = []
draw_color = []
for atom in prot.get_atoms():
    points.extend(atom.get_coord())
    if atom.get_id()[0] in colors:
        draw_color.extend(colors[atom.get_id()[0]])
    else:
        draw_color.extend(colors["_"])

vertex_list = pyglet.graphics.vertex_list(len(points) // 3, ('v3f', points), ('c3B', draw_color))


@win.event
def on_draw():
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    # Initialize OpenGL context
    glClearColor(1.0, 1.0, 1.0, 1.0)
    glLoadIdentity()
    glEnable(GL_POINT_SMOOTH)
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
