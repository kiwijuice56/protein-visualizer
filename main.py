import os
import sys

from renderers.pdb_renderer import PDBRenderer
from renderers.embedding_renderer import EmbeddingRenderer
from protein import Protein
from gui.user_interface import UserInterface

import pyglet
from pyglet.gl import *

import tkinter as tk
from tkinter import filedialog

from ctypes import windll


VERSION = "1.0"
URL_COLOR = "\u001b[35m"
DEFAULT_COLOR = "\u001b[37m"
EMBEDDING_SPACE_WIDTH = 0.3

# Fix screen resolution issue on high DPI screens
windll.shcore.SetProcessDpiAwareness(1)

# Enable terminal color on Windows
os.system("color")

print(f"protein-visualizer version {VERSION}, documentation: {URL_COLOR}https://kiwijuice56.github.io/protein"
      f"-visualizer/")
print(f"{DEFAULT_COLOR}-" * 8)


# Accept file input using either Tkinter GUI or command line arguments
if len(sys.argv) > 1:
    file = sys.argv[1]
    chain_id = sys.argv[2] if len(sys.argv) > 2 else None
    visible_window = False
else:
    print("Select a .pdb file to render.")

    root = tk.Tk()
    root.iconbitmap("img/icon.ico")
    root.withdraw()

    file = filedialog.askopenfilename(title="Select a protein file", filetypes=[('Protein Data Bank', '.pdb')])
    root.destroy()

    chain_id = None
    visible_window = True

# Initialize the protein
protein = Protein(file, chain_id=chain_id, prompt_for_chain=visible_window)


# If running the program from a command line, quit and do not show the window
if not visible_window:
    exit()

# Initialize window
window = pyglet.window.Window(resizable=True, vsync=0)
window.set_caption("protein-visualizer")
icon_ico = pyglet.image.load("img/icon.ico")
icon_png = pyglet.image.load("img/icon.png")
window.set_icon(icon_ico, icon_png)

# Initialize renderers
pdb_renderer = PDBRenderer(protein, window)
embedding_renderer = EmbeddingRenderer(protein, window)
embedding_renderer.set_bounding_box([int(window.width * (1.0 - EMBEDDING_SPACE_WIDTH)) - 8,
                                     48, int(window.width * EMBEDDING_SPACE_WIDTH),
                                     int(window.width * EMBEDDING_SPACE_WIDTH)])

# Initialize UI
ui = UserInterface(protein, window, pdb_renderer, embedding_renderer)


@window.event
def on_draw():
    # Draw a layer to detect the residues the mouse hovers over
    embedding_renderer.detect_mouse()
    pdb_renderer.detect_mouse()

    pdb_renderer.draw()
    embedding_renderer.draw()

    ui.draw()

    return pyglet.event.EVENT_HANDLED


@window.event
def on_resize(_width, _height):
    if ui.res_layout:
        ui.update_residue_label()
    pdb_renderer.set_bounding_box([0, 0, window.width, window.height])
    embedding_renderer.set_bounding_box([int(window.width * (1.0 - EMBEDDING_SPACE_WIDTH)) - 8,
                                         48, int(window.width * EMBEDDING_SPACE_WIDTH),
                                         int(window.width * EMBEDDING_SPACE_WIDTH)])


def on_update(_delta_time):
    embedding_renderer.camera.update()
    pdb_renderer.camera.update()


if __name__ == "__main__":
    pyglet.clock.schedule(on_update)
    pyglet.app.run()
