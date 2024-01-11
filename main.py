from pdb_renderer import PDBRenderer
from embedding_renderer import EmbeddingRenderer
from protein import Protein

import pyglet
from pyglet.window import key
from pyglet.gl import *

import math

import colour

# Camera code from https://gist.github.com/mr-linch/f6dacd2a069887a47fbc
from camera import FirstPersonCamera


# Initialize window
window = pyglet.window.Window(fullscreen=True)
window.set_exclusive_mouse(True)

# Initialize camera and protein
cam = FirstPersonCamera(window, movement_speed=16)
protein = Protein("proteins/alphafold_generation.pdb",
                  "proteins/alphafold_generation_sequence.fa",
                  "proteins/alphafold_generation_embeddings.h5")

initial_color = colour.Color(hue=0, saturation=0.8, luminance=0.5)
for i, residue in enumerate(protein.residues):
    initial_color.hue = (i / len(protein.residues)) * 0.8
    residue.color = [int(b * 255) for b in initial_color.rgb]

pdb_renderer = PDBRenderer(protein, window)
embedding_renderer = EmbeddingRenderer(protein, window)

highlight_index = 0


# Add inputs for the renderer parameters
def on_key_press(symbol, modifiers):
    global highlight_index
    old_index = highlight_index
    if symbol == key.LEFT:
        highlight_index -= 1
    if symbol == key.RIGHT:
        highlight_index += 1
    if highlight_index < 0:
        highlight_index = len(protein.residues) + highlight_index
    elif highlight_index >= len(protein.residues):
        highlight_index -= len(protein.residues)

    new_color = colour.Color(hue=(float(old_index) / len(protein.residues)) * 0.8, saturation=0.8, luminance=0.5)
    protein.residues[old_index].color = [int(b * 255) for b in new_color.rgb]

    new_color = colour.Color(hue=(float(highlight_index) / len(protein.residues)) * 0.8, saturation=0.8, luminance=0.5)
    protein.residues[highlight_index].color = [min(255, int(b * 255) + 128) for b in new_color.rgb]

    if old_index != highlight_index:
        pdb_renderer.update_colors(protein.residues[old_index].atoms[0].index,
                                   protein.residues[old_index].atoms[-1].index + 1)
        pdb_renderer.update_colors(protein.residues[highlight_index].atoms[0].index,
                                   protein.residues[highlight_index].atoms[-1].index + 1)

        embedding_renderer.update_colors(old_index, old_index + 1)
        embedding_renderer.update_colors(highlight_index, highlight_index + 1)

    if symbol == key.UP:
        pdb_renderer.set_point_size(pdb_renderer.point_size + 1)
    if symbol == key.DOWN:
        pdb_renderer.set_point_size(pdb_renderer.point_size - 1)

    if symbol == key.O:
        pdb_renderer.outline = not pdb_renderer.outline


window.push_handlers(on_key_press)


@window.event
def on_draw():
    global highlight_index
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glLoadIdentity()
    glMatrixMode(GL_PROJECTION)
    gluPerspective(65, window.width / float(window.height), 0.8, 512)

    # Draw the 3D protein
    cam.draw()

    pdb_renderer.draw()

    # Draw the 2D embeddings
    embedding_renderer.set_bounding_box([int(window.width * 0.6) - 16, -window.height + 16, int(window.width * 0.4), window.height])
    embedding_renderer.draw()

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
        adj = highlight_index + i
        if adj < 0:
            adj = len(protein.residues) + adj
        if adj >= len(protein.residues):
            adj -= len(protein.residues)
        digit_length = int(math.log10(len(protein.residues))) + 1
        snippet = f"{f'%0{digit_length}d' % adj}:{'%3s' % protein.residues[adj].bio_residue.get_resname()} "
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
