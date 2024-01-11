import math

import pyglet
from pyglet.gl import *


class Camera3D(object):
    class InputHandler(object):
        def __init__(self):
            self.dx = 0
            self.dy = 0
            self.dx2 = 0
            self.dy2 = 0
            self.scroll_y = 0

        def on_mouse_scroll(self, x, y, scroll_x, scroll_y):
            self.scroll_y = scroll_y

        def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
            if buttons & pyglet.window.mouse.LEFT:
                self.dx = dx
                self.dy = dy
            if buttons & pyglet.window.mouse.MIDDLE:
                self.dx2 = dx
                self.dy2 = dy

    def __init__(self, window, movement_speed=16, mouse_sensitivity=0.01):
        self.pivot_pos = [0, 0, -32]
        self.camera_pos = [12, 0, math.pi/2]
        self.forward_dir = [0, 0, 1]
        self.global_up_dir = [0, -1, 0]
        self.local_up_dir = [0, -1, 0]
        self.right_dir = [0, 0, 0]

        self.input_handler = Camera3D.InputHandler()

        window.push_handlers(self.input_handler)

        self.movement_speed = movement_speed
        self.mouse_sensitivity = mouse_sensitivity

    def move_horizontal(self, distance):
        self.pivot_pos = [self.right_dir[i] * distance + self.pivot_pos[i] for i in range(3)]

    def move_vertical(self, distance):
        self.pivot_pos = [self.local_up_dir[i] * -distance + self.pivot_pos[i] for i in range(3)]

    def update(self, delta_time):
        self.move_horizontal(self.input_handler.dx * self.movement_speed)
        self.input_handler.dx = 0

        self.move_vertical(self.input_handler.dy * self.movement_speed)
        self.input_handler.dy = 0

        self.camera_pos[1] -= self.input_handler.dx2 * self.mouse_sensitivity
        self.input_handler.dx2 = 0

        self.camera_pos[2] += self.input_handler.dy2 * self.mouse_sensitivity
        self.camera_pos[2] = max(0.1, min(math.pi - 0.1, self.camera_pos[2]))
        self.input_handler.dy2 = 0

        self.camera_pos[0] -= self.input_handler.scroll_y
        self.camera_pos[0] = max(0.05, self.camera_pos[0])
        self.input_handler.scroll_y = 0

    @staticmethod
    def cross(a, b):
        c = [a[1] * b[2] - a[2] * b[1],
             a[2] * b[0] - a[0] * b[2],
             a[0] * b[1] - a[1] * b[0]]

        return c

    def draw(self):
        camera_trans = [0, 0, 0]
        camera_trans[0] = self.camera_pos[0] * math.sin(self.camera_pos[2]) * math.cos(self.camera_pos[1])
        camera_trans[2] = self.camera_pos[0] * math.sin(self.camera_pos[2]) * math.sin(self.camera_pos[1])
        camera_trans[1] = self.camera_pos[0] * math.cos(self.camera_pos[2])

        self.forward_dir = [component / self.camera_pos[0] for component in camera_trans]
        self.right_dir = self.cross(self.forward_dir, self.global_up_dir)
        self.local_up_dir = self.cross(self.forward_dir, self.right_dir)

        gluLookAt(0, 0, 0, *self.forward_dir, *self.global_up_dir)
        glTranslatef(*camera_trans)
        glTranslatef(*self.pivot_pos)
