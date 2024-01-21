import math

import pyglet
from pyglet.gl import *


class UserInterface:
    def __init__(self, protein, window, pdb_renderer, embedding_renderer):
        self.protein = protein
        self.window = window
        self.pdb_renderer = pdb_renderer
        self.embedding_renderer = embedding_renderer

        window.push_handlers(self.on_key_press)
        window.push_handlers(self.on_mouse_motion)
        window.push_handlers(self.on_mouse_scroll)

        self.hl_idx = -1
        self.stored_hl_idx = -1

        self.go_idx = 0

        self.res_doc = pyglet.text.document.FormattedDocument()
        self.res_layout = None

        self.update_residue_label()

    def on_key_press(self, symbol, modifiers):
        if symbol == pyglet.window.key.UP:
            self.pdb_renderer.set_point_size(self.pdb_renderer.point_size + 1)
        if symbol == pyglet.window.key.DOWN:
            self.pdb_renderer.set_point_size(self.pdb_renderer.point_size - 1)

        if symbol == pyglet.window.key.RIGHT:
            self.go_idx += 1
            self.go_idx = self.go_idx % len(self.protein.go_ids)
            self.protein.current_go_id = self.protein.go_ids[self.go_idx]

            self.protein.update_colors(self.protein.GO_ANNOTATION)
            self.pdb_renderer.update_colors()
            self.pdb_renderer.update_colors()

            print(f"Viewing {self.protein.go_names[self.go_idx]}.")

        if symbol == pyglet.window.key.LEFT:
            self.go_idx -= 1
            self.go_idx = self.go_idx % len(self.protein.go_ids)
            self.protein.current_go_id = self.protein.go_ids[self.go_idx]

            self.protein.update_colors(self.protein.GO_ANNOTATION)
            self.pdb_renderer.update_colors()
            self.pdb_renderer.update_colors()

            print(f"Viewing {self.protein.go_names[self.go_idx]}.")

        if symbol == pyglet.window.key._1:
            self.protein.update_colors(self.protein.CLUSTER_INDEX)
            self.pdb_renderer.update_colors()
            self.pdb_renderer.update_colors()
        if symbol == pyglet.window.key._2:
            self.protein.update_colors(self.protein.RESIDUE_INDEX)
            self.pdb_renderer.update_colors()
            self.pdb_renderer.update_colors()
        if symbol == pyglet.window.key._3:
            self.protein.update_colors(self.protein.ATOM_TYPE)
            self.pdb_renderer.update_colors()
            self.pdb_renderer.update_colors()
        if symbol == pyglet.window.key._4:
            self.protein.update_colors(self.protein.GO_ANNOTATION)
            self.pdb_renderer.update_colors()
            self.pdb_renderer.update_colors()

        if symbol == pyglet.window.key._8:
            self.protein.update_colors(new_color_palette=self.protein.RAINBOW)
            self.pdb_renderer.update_colors()
            self.pdb_renderer.update_colors()
        if symbol == pyglet.window.key._9:
            self.protein.update_colors(new_color_palette=self.protein.POISSON)
            self.pdb_renderer.update_colors()
            self.pdb_renderer.update_colors()

        if symbol == pyglet.window.key.O:
            self.pdb_renderer.outline = not self.pdb_renderer.outline

    def on_mouse_scroll(self, x, y, scroll_x, scroll_y):
        self.on_mouse_motion(x, y, 0, 0)

    def on_mouse_motion(self, x, y, dx, dy):
        prev_hl_idx = self.hl_idx
        self.hl_idx = self.pdb_renderer.hovered_residue

        if not prev_hl_idx == -1:
            self.protein.residues[prev_hl_idx].highlighted = False
            self.protein.color_residue(self.protein.residues[prev_hl_idx])

        if not self.hl_idx == -1:
            self.stored_hl_idx = self.hl_idx
            self.protein.residues[self.hl_idx].highlighted = True
            self.protein.color_residue(self.protein.residues[self.hl_idx])

        if prev_hl_idx != self.hl_idx:
            self.update_residue_label()
            self.pdb_renderer.update_colors(self.protein.residues[prev_hl_idx].atoms[0].index,
                                            self.protein.residues[prev_hl_idx].atoms[-1].index + 1)
            self.pdb_renderer.update_colors(self.protein.residues[self.hl_idx].atoms[0].index,
                                            self.protein.residues[self.hl_idx].atoms[-1].index + 1)

            self.embedding_renderer.update_colors(prev_hl_idx, prev_hl_idx + 1)
            self.embedding_renderer.update_colors(self.hl_idx, self.hl_idx + 1)

    def update_residue_label(self):
        if self.res_layout:
            self.res_layout.begin_update()
        self.res_doc.text = ""

        text_style = {"font_name": "Consolas",
                      "font_size": 16,
                      "background_color": (0, 0, 0, 140),
                      "color": (255, 255, 255, 255),
                      "wrap": False}

        # Rough estimate of how many amino acid labels can fit on screen
        labels_fit = math.ceil(self.window.width / 32.0)

        # Draw the residue label indices
        for i in range(-2, labels_fit + 1):
            adj = self.stored_hl_idx + i
            if adj < 0:
                adj = len(self.protein.residues) + adj
            if adj >= len(self.protein.residues):
                adj -= len(self.protein.residues)
            snippet = f"{'%-8d' % self.protein.residues[adj].seq_index}"
            color = tuple((255, 230, 0, 255) if i == 0 and not self.hl_idx == -1 else [255, 255, 255, 255])
            text_style["color"] = color
            text_style["font_size"] = 8
            self.res_doc.insert_text(len(self.res_doc.text), snippet, text_style)
        self.res_doc.insert_text(len(self.res_doc.text), '\n')

        # Draw the residue label names
        for i in range(-2, labels_fit + 1):
            adj = self.stored_hl_idx + i
            if adj < 0:
                adj = len(self.protein.residues) + adj
            if adj >= len(self.protein.residues):
                adj -= len(self.protein.residues)
            snippet = f"{'%3s' % self.protein.residues[adj].bio_residue.get_resname()} "
            color = tuple((255, 230, 0, 255) if i == 0 and not self.hl_idx == -1 else [255, 255, 255, 255])
            text_style["color"] = color
            text_style["font_size"] = 16
            self.res_doc.insert_text(len(self.res_doc.text), snippet, text_style)
        self.res_layout = pyglet.text.layout.TextLayout(self.res_doc, multiline=True, width=self.window.width)
        self.res_layout.anchor_x = "left"
        self.res_layout.anchor_y = "bottom"
        self.res_layout.y = -self.window.height

    def draw(self):
        glLoadIdentity()
        glDisable(GL_DEPTH_TEST)
        glEnable(GL_BLEND)
        glOrtho(0, self.window.width, -self.window.height, 0, 0, 1000)

        fps_doc = pyglet.text.Label(text=f"fps: {round(pyglet.clock.get_fps())}",
                                    font_name="Consolas",
                                    font_size=16,
                                    x=16, y=-16,
                                    anchor_x="left", anchor_y="top")
        if self.res_layout:
            self.res_layout.draw()
        fps_doc.draw()
