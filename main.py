from pdb_renderer import PDBRenderer
from embedding_renderer import EmbeddingRenderer
from protein import Protein
from user_interface import UserInterface

import sys

import pyglet
from pyglet.gl import *

# Initialize window
window = pyglet.window.Window(resizable=True)
window.set_caption("protein-visualizer")

# Initialize protein
default_path = "data/alphafold_generation.pdb"
protein = Protein(sys.argv[1] if len(sys.argv) > 1 else default_path)

# Initialize renderers
pdb_renderer = PDBRenderer(protein, window)
embedding_renderer = EmbeddingRenderer(protein, window)

# Initialize UI
ui = UserInterface(protein, window, pdb_renderer, embedding_renderer)


@window.event
def on_draw():
    # Draw the 3D protein
    pdb_renderer.set_bounding_box([0, 0, window.width, window.height])
    pdb_renderer.draw()

    # Draw the 2D embeddings
    embedding_renderer.set_bounding_box([int(window.width * 0.7) - 8, 47, int(window.width * 0.3), int(window.width * 0.3)])
    embedding_renderer.draw()

    # Draw the user interface
    ui.draw()

    return pyglet.event.EVENT_HANDLED


@window.event
def on_resize(width, height):
    if ui.res_layout:
        ui.update_residue_label()


def on_update(delta_time):
    embedding_renderer.camera.update()
    pdb_renderer.camera.update()


if __name__ == "__main__":
    pyglet.clock.schedule(on_update)
    pyglet.app.run()
