import os

from renderers.pdb_renderer import PDBRenderer
from renderers.embedding_renderer import EmbeddingRenderer
from protein import Protein
from gui.user_interface import UserInterface

import pyglet
from pyglet.gl import *

import tkinter as tk
from tkinter import filedialog

from ctypes import windll

# Enable terminal color on Windows
os.system("color")

print("protein-visualizer version 1.0, documentation: \u001b[35mhttps://kiwijuice56.github.io/protein-visualizer/")
print("\u001b[37m-" * 8)
print("Select a .pdb file to render.")

# Fix screen resolution
windll.shcore.SetProcessDpiAwareness(1)

root = tk.Tk()
root.iconbitmap("img/icon.ico")
root.withdraw()

# Initialize protein
protein = Protein(filedialog.askopenfilename(title="Select a protein file", filetypes=[('Protein Data Bank', '.pdb')]))

# Initialize window
window = pyglet.window.Window(resizable=True, vsync=0)
window.set_caption("protein-visualizer")
icon_ico = pyglet.image.load("img/icon.ico")
icon_png = pyglet.image.load("img/icon.png")
window.set_icon(icon_ico, icon_png)

root.destroy()

# Initialize renderers
pdb_renderer = PDBRenderer(protein, window)
embedding_renderer = EmbeddingRenderer(protein, window)
embedding_renderer.set_bounding_box([int(window.width * 0.7) - 8, 47, int(window.width * 0.3), int(window.width * 0.3)])

# Initialize UI
ui = UserInterface(protein, window, pdb_renderer, embedding_renderer)


@window.event
def on_draw():
    embedding_renderer.detect_mouse()
    pdb_renderer.detect_mouse()

    # Draw the 3D protein
    pdb_renderer.draw()

    # Draw the 2D embeddings
    embedding_renderer.draw()

    # Draw the user interface
    ui.draw()

    return pyglet.event.EVENT_HANDLED


@window.event
def on_resize(width, height):
    if ui.res_layout:
        ui.update_residue_label()
    pdb_renderer.set_bounding_box([0, 0, window.width, window.height])
    embedding_renderer.set_bounding_box([int(window.width * 0.7) - 8, 47, int(window.width * 0.3),
                                         int(window.width * 0.3)])


def on_update(delta_time):
    embedding_renderer.camera.update()
    pdb_renderer.camera.update()


if __name__ == "__main__":
    pyglet.clock.schedule(on_update)
    pyglet.app.run()
