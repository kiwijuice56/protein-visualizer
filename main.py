import pyglet
from pyglet.window import key
from pyglet.gl import *

import math

# Camera code from https://gist.github.com/mr-linch/f6dacd2a069887a47fbc
from camera import FirstPersonCamera

from pdb_renderer import PDBRenderer


# Initialize the window and camera
window = pyglet.window.Window()
window.set_exclusive_mouse(True)
window.set_fullscreen(True)
cam = FirstPersonCamera(window, movement_speed=16)
pdb_renderer = PDBRenderer("proteins/alphafold-generation.pdb", window)


# Add inputs for the renderer parameters
def on_key_press(symbol, modifiers):
    if symbol == key.LEFT:
        pdb_renderer.set_highlighted_index(pdb_renderer.highlighted_index - 1)
    if symbol == key.RIGHT:
        pdb_renderer.set_highlighted_index(pdb_renderer.highlighted_index + 1)

    if symbol == key.UP:
        pdb_renderer.set_point_size(pdb_renderer.point_size + 1)
    if symbol == key.DOWN:
        pdb_renderer.set_point_size(pdb_renderer.point_size - 1)

    if symbol == key.O:
        pdb_renderer.outline = not pdb_renderer.outline

    if symbol == key._1:
        pdb_renderer.set_color_mode(pdb_renderer.ColorMode.CPK)
    if symbol == key._2:
        pdb_renderer.set_color_mode(pdb_renderer.ColorMode.CHAINBOW)
    if symbol == key._3:
        pdb_renderer.set_color_mode(pdb_renderer.ColorMode.CONTRAST)
    if symbol == key._4:
        pdb_renderer.set_color_mode(pdb_renderer.ColorMode.RESIDUE)


window.push_handlers(on_key_press)


@window.event
def on_draw():
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glLoadIdentity()
    glMatrixMode(GL_PROJECTION)
    gluPerspective(65, window.width / float(window.height), 0.8, 512)

    # Draw the 3D protein
    cam.draw()
    pdb_renderer.draw()

    # Reset projection to 2D for UI
    glLoadIdentity()
    glDisable(GL_DEPTH_TEST)
    glEnable(GL_BLEND)
    glOrtho(0, window.width, -window.height, 0, 0, 1000)

    document = pyglet.text.document.FormattedDocument()

    # Draw the title label
    document.insert_text(0, "example protein\n",
                         {"font_name": "Consolas",
                          "font_size": 16,
                          "background_color": (0, 0, 0, 180),
                          "color": (255, 255, 255, 255)})

    # Draw the adjacent residue label
    for i in range(-2, 4):
        adj = pdb_renderer.highlighted_index + i
        if adj < 0:
            adj = len(pdb_renderer.residues) + adj
        if adj >= len(pdb_renderer.residues):
            adj -= len(pdb_renderer.residues)
        digit_length = int(math.log10(len(pdb_renderer.residues))) + 1
        snippet = f"{f'%0{digit_length}d' % adj}:{'%3s' % pdb_renderer.residues[adj].bio_residue.get_resname()} "
        color = tuple((255, 230, 0, 255) if i == 0 else [200 - abs(i) * 24] * 3 + [255])
        document.insert_text(len(document.text), snippet, {"color": color})

    layout = pyglet.text.layout.TextLayout(document, multiline=True, width=window.width, height=window.height)
    layout.anchor_x = "left"
    layout.anchor_y = "top"

    # Padding
    layout.x += 4
    layout.y -= 4

    layout.draw()

    return pyglet.event.EVENT_HANDLED


def on_update(delta_time):
    cam.update(delta_time)


if __name__ == "__main__":
    pyglet.clock.schedule(on_update)
    pyglet.app.run()
