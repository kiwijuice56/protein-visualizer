import sys

from pdb_renderer import PDBRenderer
from embedding_renderer import EmbeddingRenderer
from protein import Protein

import pyglet
from pyglet.window import key
from pyglet.gl import *

import math

# Initialize window
window = pyglet.window.Window(resizable=True)
window.set_caption("protein-visualizer")
fps_display = pyglet.window.FPSDisplay(window=window)

# Initialize protein
default_path = "data/alphafold_generation.pdb"
protein = Protein(sys.argv[1] if len(sys.argv) > 1 else default_path)

# Initialize rendering windows
pdb_renderer = PDBRenderer(protein, window)
embedding_renderer = EmbeddingRenderer(protein, window)


# Initialize UI
residue_doc = pyglet.text.document.FormattedDocument()
highlight_index = -1
remembered_index = 0


def update_residue_label():
    residue_doc.text = ""

    text_style = {"font_name": "Consolas",
                  "font_size": 16,
                  "background_color": (0, 0, 0, 140),
                  "color": (255, 255, 255, 255),
                  "wrap": False}

    # Rough estimate of how many amino acid labels can fit on screen
    labels_fit = math.ceil(window.width / 32.0)

    # Draw the residue label indices
    for i in range(-2, labels_fit + 1):
        adj = remembered_index + i
        if adj < 0:
            adj = len(protein.residues) + adj
        if adj >= len(protein.residues):
            adj -= len(protein.residues)
        snippet = f"{'%-8d' % protein.residues[adj].seq_index}"
        color = tuple((255, 230, 0, 255) if i == 0 and not highlight_index == -1 else [255, 255, 255, 255])
        text_style["color"] = color
        text_style["font_size"] = 8
        residue_doc.insert_text(len(residue_doc.text), snippet, text_style)
    residue_doc.insert_text(len(residue_doc.text), '\n')

    # Draw the residue label names
    for i in range(-2, labels_fit + 1):
        adj = remembered_index + i
        if adj < 0:
            adj = len(protein.residues) + adj
        if adj >= len(protein.residues):
            adj -= len(protein.residues)
        snippet = f"{'%3s' % protein.residues[adj].bio_residue.get_resname()} "
        color = tuple((255, 230, 0, 255) if i == 0 and not highlight_index == -1 else [255, 255, 255, 255])
        text_style["color"] = color
        text_style["font_size"] = 16
        residue_doc.insert_text(len(residue_doc.text), snippet, text_style)


# Add inputs for the renderer parameters
def on_key_press(symbol, modifiers):
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
    if symbol == key._3:
        protein.update_colors(protein.ATOM_TYPE)
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


def on_mouse_motion(x, y, dx, dy):
    global highlight_index
    global remembered_index
    old_index = highlight_index
    highlight_index = pdb_renderer.hovered_residue

    if not old_index == -1:
        protein.residues[old_index].highlighted = False
        protein.color_residue(protein.residues[old_index])

    if not highlight_index == -1:
        remembered_index = highlight_index
        protein.residues[highlight_index].highlighted = True
        protein.color_residue(protein.residues[highlight_index])

    if old_index != highlight_index:
        update_residue_label()
        pdb_renderer.update_colors(protein.residues[old_index].atoms[0].index,
                                   protein.residues[old_index].atoms[-1].index + 1)
        pdb_renderer.update_colors(protein.residues[highlight_index].atoms[0].index,
                                   protein.residues[highlight_index].atoms[-1].index + 1)

        embedding_renderer.update_colors(old_index, old_index + 1)
        embedding_renderer.update_colors(highlight_index, highlight_index + 1)


window.push_handlers(on_key_press)
window.push_handlers(on_mouse_motion)


@window.event
def on_draw():
    global highlight_index

    # Draw the 3D protein
    pdb_renderer.set_bounding_box([0, 0, window.width, window.height])
    pdb_renderer.draw()

    # Draw the 2D embeddings
    embedding_renderer.set_bounding_box(
        [int(window.width * 0.7) - 8, 47, int(window.width * 0.3), int(window.width * 0.3)])
    embedding_renderer.draw()

    # Reset projection to 2D for UI
    glLoadIdentity()
    glDisable(GL_DEPTH_TEST)
    glEnable(GL_BLEND)
    glOrtho(0, window.width, -window.height, 0, 0, 1000)

    fps_doc = pyglet.text.Label(text=f"fps: {round(pyglet.clock.get_fps())}",
                                font_name="Consolas",
                                font_size=16,
                                x=16, y=-16,
                                anchor_x="left", anchor_y="top")

    residue_layout = pyglet.text.layout.TextLayout(residue_doc, multiline=True, width=window.width)
    residue_layout.anchor_x = "left"
    residue_layout.anchor_y = "bottom"
    residue_layout.y = -window.height

    residue_layout.draw()
    fps_doc.draw()

    return pyglet.event.EVENT_HANDLED


def on_update(delta_time):
    embedding_renderer.camera.update()
    pdb_renderer.camera.update()


if __name__ == "__main__":
    update_residue_label()
    pyglet.clock.schedule(on_update)
    pyglet.app.run()
