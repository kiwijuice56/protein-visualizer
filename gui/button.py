import pyglet


class Button(pyglet.gui.WidgetBase):
    BACKGROUND_COLOR = (0, 0, 0, 140)
    HOVER_COLOR = (64, 64, 64, 140)

    def __init__(self, bounding_box, text, index, window, batch, bg_batch, text_width=21):
        pyglet.gui.WidgetBase.__init__(self, *bounding_box)
        self.batch = batch
        self.window = window
        window.push_handlers(self.on_mouse_press, self.on_mouse_motion)

        self.bounding_box = bounding_box

        self.text = text
        self.index = index

        self.background = pyglet.shapes.Rectangle(*self.bounding_box, color=self.BACKGROUND_COLOR[0:3], batch=bg_batch)
        self.background.opacity = self.BACKGROUND_COLOR[3]
        self.label = pyglet.text.Label(
            text=self.text, font_name="Consolas", multiline=True,
            font_size=16, x=self.x, y=self.y + 28, width=self.width,
            anchor_x="left", anchor_y="top", batch=self.batch)
        self.visible = False
        self.pressed = False
        self.hovered = False

    def show(self):
        self.background.opacity = self.BACKGROUND_COLOR[3]
        self.label.batch = self.batch
        self.visible = True

    def hide(self):
        self.background.color = self.BACKGROUND_COLOR[0: 3]
        self.background.opacity = 0
        self.label.batch = None
        self.visible = False
        self.hovered = False

    def on_mouse_press(self, x, y, button, modifiers):
        if not self.visible:
            return
        b = [c for c in self.bounding_box]
        y -= self.window.height
        if x > b[0] and x - b[0] < b[2] and y > b[1] and y - b[1] < b[3]:
            self.pressed = True
            return pyglet.event.EVENT_HANDLED

    def on_mouse_motion(self, x, y, dx, dy):
        if not self.visible:
            return
        b = [c for c in self.bounding_box]
        y -= self.window.height
        if x > b[0] and x - b[0] < b[2] and y > b[1] and y - b[1] < b[3]:
            self.background.color = self.HOVER_COLOR[0: 3]
            self.hovered = True
        elif self.hovered:
            self.background.color = self.BACKGROUND_COLOR[0: 3]
            self.hovered = False


class DropDown(Button):
    def __init__(self, bounding_box, title, options, window, batch, bg_batch, text_width=21):
        Button.__init__(self, bounding_box, f'%-{text_width}s ▾' % options[0], -1, window, batch, bg_batch)
        self.batch = batch
        self.window = window
        self.window.push_handlers(self.on_mouse_press, self.on_mouse_motion)

        self.options = options
        self.buttons = []
        for i, option in enumerate(self.options):
            button = Button([self.x, self.y - 32 * (i + 1), self.width, 32], self.options[i], i, self.window, self.batch, bg_batch)
            self.buttons.append(button)

        self.title = pyglet.text.Label(
            text=title, font_name="Consolas", multiline=True,
            font_size=16, x=self.x, y=self.bounding_box[1] + 64, width=self.width,
            anchor_x="left", anchor_y="top", batch=self.batch)

        self.is_open = False
        self.close()
        self.show()

    def open(self):
        self.is_open = True
        for button in self.buttons:
            button.show()

    def close(self):
        self.is_open = False
        for button in self.buttons:
            button.hide()
