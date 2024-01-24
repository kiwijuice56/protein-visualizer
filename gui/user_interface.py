import math

import pyglet
from pyglet.gl import *


class DropDown(pyglet.gui.WidgetBase):
    BACKGROUND_COLOR = (0, 0, 0, 140)

    def __init__(self, x, y, width, height, title, options, window, batch):
        pyglet.gui.WidgetBase.__init__(self, x, y, width, height)
        self.batch = batch
        self.window = window
        self.window.push_handlers(self.on_mouse_press)

        self.bounding_box = [x, y, width, height]

        self.title = title
        self.options = options
        self.selected_option = options[0]
        self.open = False

        self.batch = batch
        self.buttons = []
        for i, option in enumerate(self.options):
            button = Button(self.x, self.y - 32 * (i + 2), self.width, 32, self.options[i], i, self.window, self.batch)
            button.hide()
            self.buttons.append(button)

        self.title = pyglet.text.Label(
            text=self.title, font_name="Consolas", multiline=True,
            font_size=16, x=self.x, y=-20, width=self.width,
            anchor_x="left", anchor_y="top", batch=self.batch)

        self.current_option = pyglet.text.Label(
            text='%-21s ▼' % self.selected_option, font_name="Consolas", multiline=True,
            font_size=16, x=self.x, y=-52, width=self.width,
            anchor_x="left", anchor_y="top", batch=self.batch)

        self.background = pyglet.shapes.Rectangle(*self.bounding_box, color=self.BACKGROUND_COLOR[0:3],
                                                  batch=self.batch)
        self.background.y -= 32
        self.background.opacity = self.BACKGROUND_COLOR[3]

    def on_mouse_press(self, x, y, button, modifiers):
        y = y - self.window.height
        box = [c for c in self.bounding_box]
        box[1] -= 32
        if not (x > box[0] and x - box[0] < box[2] and
                y > box[1] and y - box[1] < box[3]) or self.open:
            for button in self.buttons:
                self.open = False
                button.hide()
        else:
            for button in self.buttons:
                self.open = True
                button.show()
            return pyglet.event.EVENT_HANDLED


class Button(pyglet.gui.WidgetBase):
    BACKGROUND_COLOR = (0, 0, 0, 140)

    def __init__(self, x, y, width, height, text, index, window, batch):
        pyglet.gui.WidgetBase.__init__(self, x, y, width, height)
        self.batch = batch
        self.window = window
        window.push_handlers(self.on_mouse_press)

        self.bounding_box = [x, y, width, height]
        self.text = text
        self.index = index

        self.background = pyglet.shapes.Rectangle(*self.bounding_box, color=self.BACKGROUND_COLOR[0:3],
                                                  batch=self.batch)
        self.background.opacity = self.BACKGROUND_COLOR[3]
        self.label = pyglet.text.Label(
            text=self.text, font_name="Consolas", multiline=True,
            font_size=16, x=self.x, y=self.y + 28, width=self.width,
            anchor_x="left", anchor_y="top", batch=self.batch)

        self.visible = False
        self.pressed = False

    def show(self):
        self.label.batch = self.batch
        self.background.batch = self.batch
        self.background.opacity = self.BACKGROUND_COLOR[3]
        self.visible = True

    def hide(self):
        self.label.batch = None
        self.background.opacity = 0
        self.visible = False

    def on_mouse_press(self, x, y, button, modifiers):
        if not self.visible:
            return

        y = y - self.window.height
        box = [c for c in self.bounding_box]
        if not (x > box[0] and x - box[0] < box[2] and
                y > box[1] and y - box[1] < box[3]):
            return
        else:
            self.pressed = True


class UserInterface:
    def __init__(self, protein, window, pdb_renderer, embedding_renderer):
        self.protein = protein
        self.window = window
        self.pdb_renderer = pdb_renderer
        self.embedding_renderer = embedding_renderer

        self.batch = pyglet.graphics.Batch()

        window.push_handlers(self.on_key_press)
        window.push_handlers(self.on_mouse_motion)
        window.push_handlers(self.on_mouse_scroll)

        self.hl_idx = -1
        self.stored_hl_idx = -1

        self.go_idx = 0

        self.res_doc = pyglet.text.document.FormattedDocument()
        self.res_layout = None

        self.color_mode = DropDown(x=16, y=-48, width=292, height=32, title="Coloring Mode",
                                   options=["Functional Similarity", "Amino Acid Order", "Atom Type"],
                                   window=window, batch=self.batch)
        self.color_palette = DropDown(x=324, y=-48, width=292, height=32, title="Color Palette",
                                      options=["Rainbow", "Monocolor", "Poisson", "Grape"],
                                      window=window, batch=self.batch)

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

            self.protein.update_colors()
            self.pdb_renderer.update_colors()
            self.pdb_renderer.update_colors()

        if symbol == pyglet.window.key.LEFT:
            self.go_idx -= 1
            self.go_idx = self.go_idx % len(self.protein.go_ids)
            self.protein.current_go_id = self.protein.go_ids[self.go_idx]

            self.protein.update_colors()
            self.pdb_renderer.update_colors()
            self.embedding_renderer.update_colors()

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
            self.res_layout.delete()
            self.res_layout.begin_update()
        self.res_doc.text = ""

        text_style = {"font_name": "Consolas", "font_size": 16, "wrap": False,
                      "background_color": (0, 0, 0, 140), "color": (255, 255, 255, 255)}

        # Rough estimate of how many amino acid labels can fit on screen
        labels_fit = math.ceil(self.window.width / 52.0)

        # Draw the residue label indices
        for i in range(-2, labels_fit + 1):
            adj = self.stored_hl_idx + i
            if adj < 0:
                adj = len(self.protein.residues) + adj
            if adj >= len(self.protein.residues):
                adj -= len(self.protein.residues)
            snippet = '%-8d ' % self.protein.residues[adj].index
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
            snippet = '%3s ' % self.protein.residues[adj].bio_residue.get_resname()
            color = tuple((255, 230, 0, 255) if i == 0 and not self.hl_idx == -1 else [255, 255, 255, 255])
            text_style["color"] = color
            text_style["font_size"] = 16
            self.res_doc.insert_text(len(self.res_doc.text), snippet, text_style)
        self.res_layout = pyglet.text.layout.TextLayout(self.res_doc, multiline=True, width=self.window.width,
                                                        batch=self.batch)
        self.res_layout.anchor_x = "left"
        self.res_layout.anchor_y = "bottom"
        self.res_layout.y = -self.window.height

    def draw(self):
        glLoadIdentity()
        glDisable(GL_DEPTH_TEST)
        glEnable(GL_BLEND)
        glOrtho(0, self.window.width, -self.window.height, 0, 0, 1000)

        for button in self.color_mode.buttons:
            if button.pressed:
                button.pressed = False
                self.color_mode.current_option.text = '%-21s ▼' % button.label.text
                self.protein.update_colors(new_color_mode=button.index + 1)
                self.pdb_renderer.update_colors()
                self.embedding_renderer.update_colors()
        for button in self.color_palette.buttons:
            if button.pressed:
                button.pressed = False
                self.color_palette.current_option.text = '%-21s ▼' % button.label.text
                self.protein.update_colors(new_color_palette=button.index + 6)
                self.pdb_renderer.update_colors()
                self.embedding_renderer.update_colors()

        go_doc = pyglet.text.Label(
            text=f"{self.protein.go_ids[self.go_idx]}: {self.protein.go_names[self.go_idx]}\n{'%.2f' % (self.protein.scores[self.go_idx] * 100)}% confidence",
            font_name="Consolas", multiline=True,
            font_size=16, x=632, y=-16, width=self.window.width,
            anchor_x="left", anchor_y="top")
        go_doc.draw()

        self.batch.draw()
