from pdb_renderer import PDBRenderer
from embedding_renderer import EmbeddingRenderer
from protein import Protein

import pyglet
from pyglet.window import key
from pyglet.gl import *

import math


# Initialize window
window = pyglet.window.Window(resizable=True)
window.set_exclusive_mouse(False)

# Initialize protein
protein = Protein("proteins/alphafold_generation.pdb",
                  "proteins/alphafold_generation_sequence.fa",
                  "proteins/alphafold_generation_embeddings.h5")

# Initialize rendering windows
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

    protein.residues[old_index].highlighted = False
    protein.color_residue(protein.residues[old_index])

    protein.residues[highlight_index].highlighted = True
    protein.color_residue(protein.residues[highlight_index])

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

    if symbol == key._1:
        protein.update_colors(protein.CLUSTER_INDEX)
        pdb_renderer.update_colors()
        pdb_renderer.update_colors()
    if symbol == key._2:
        protein.update_colors(protein.RESIDUE_INDEX)
        pdb_renderer.update_colors()
        pdb_renderer.update_colors()

    if symbol == key._8:
        protein.update_colors(new_color_palette=protein.RAINBOW)
        pdb_renderer.update_colors()
        pdb_renderer.update_colors()
    if symbol == key._9:
        protein.update_colors(new_color_palette=protein.POISSON)
        pdb_renderer.update_colors()
        pdb_renderer.update_colors()

    if symbol == key.O:
        pdb_renderer.outline = not pdb_renderer.outline


window.push_handlers(on_key_press)


@window.event
def on_draw():
    global highlight_index

    # Draw the 3D protein
    pdb_renderer.set_bounding_box([0, 0, window.width, window.height])
    pdb_renderer.draw()

    # Draw the 2D embeddings
    embedding_renderer.set_bounding_box([int(window.width * 0.7) - 16, 16, int(window.width * 0.3), int(window.width * 0.3)])
    embedding_renderer.draw()

    # Reset projection to 2D for UI
    glLoadIdentity()
    glDisable(GL_DEPTH_TEST)
    glEnable(GL_BLEND)
    glOrtho(0, window.width, -window.height, 0, 0, 1000)

    document = pyglet.text.document.FormattedDocument()

    text_style = {"font_name": "Consolas",
                  "font_size": 16,
                  "background_color": (0, 0, 0, 180),
                  "color": (255, 255, 255, 255)}

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
        text_style["color"] = color
        document.insert_text(len(document.text), snippet, text_style)

    layout = pyglet.text.layout.TextLayout(document, multiline=True, width=window.width, height=window.height)
    layout.anchor_x = "left"
    layout.anchor_y = "top"

    # Padding
    layout.x += 4
    layout.y -= 4

    layout.draw()

    return pyglet.event.EVENT_HANDLED


def on_update(delta_time):
    embedding_renderer.camera.update()
    pdb_renderer.camera.update()


if __name__ == "__main__":
    pyglet.clock.schedule(on_update)
    pyglet.app.run()
