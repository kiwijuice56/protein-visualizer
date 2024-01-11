import pyglet


class InputHandler(object):
    def __init__(self):
        self.dx_left = 0
        self.dy_left = 0
        self.dx_middle = 0
        self.dy_middle = 0
        self.scroll_y = 0
        self.gesture_pos = [0, 0]
        self.bounding_box = [0, 0, 0, 0]

    def on_mouse_press(self, x, y, button, modifiers):
        self.gesture_pos = [x, y]

    def on_mouse_scroll(self, x, y, scroll_x, scroll_y):
        self.scroll_y = scroll_y

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
        if self.gesture_pos[0] < self.bounding_box[0] or self.gesture_pos[0] - self.bounding_box[0] > self.bounding_box[2]:
            return
        if self.gesture_pos[1] < self.bounding_box[1] or self.gesture_pos[1] - self.bounding_box[1] > self.bounding_box[3]:
            return
        if buttons & pyglet.window.mouse.LEFT:
            self.dx_left = dx
            self.dy_left = dy
        if buttons & (pyglet.window.mouse.MIDDLE | pyglet.window.mouse.RIGHT):
            self.dx_middle = dx
            self.dy_middle = dy
        return pyglet.event.EVENT_HANDLED
