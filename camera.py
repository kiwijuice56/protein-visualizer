import numpy as np

import pyglet
from pyglet.gl import *


class Camera3D(object):
    class InputHandler(object):
        def __init__(self):
            self.dx_left = 0
            self.dy_left = 0
            self.dx_middle = 0
            self.dy_middle = 0
            self.scroll_y = 0

        def on_mouse_scroll(self, x, y, scroll_x, scroll_y):
            self.scroll_y = scroll_y

        def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
            if buttons & pyglet.window.mouse.LEFT:
                self.dx_left = dx
                self.dy_left = dy
            if buttons & pyglet.window.mouse.RIGHT:
                self.dx_middle = dx
                self.dy_middle = dy

    def __init__(self, window, movement_speed=16, mouse_sensitivity=0.01):
        self.pivot_pos = np.array([0, 0, -32])
        self.camera_pos = np.array([12, 0, np.pi/2])

        self.forward_dir = np.array([0, 0, 1])
        self.right_dir = np.array([0, 0, 0])
        self.global_up_dir = np.array([0, -1, 0])
        self.local_up_dir = np.array([0, -1, 0])

        self.input_handler = Camera3D.InputHandler()

        window.push_handlers(self.input_handler)

        self.movement_speed = movement_speed
        self.mouse_sensitivity = mouse_sensitivity

    def move_horizontal(self, distance):
        self.pivot_pos = self.pivot_pos + self.right_dir * distance

    def move_vertical(self, distance):
        self.pivot_pos = self.pivot_pos + self.local_up_dir * distance

    def update(self, delta_time):
        self.move_horizontal(self.input_handler.dx_left * self.movement_speed)
        self.input_handler.dx_left = 0

        self.move_vertical(self.input_handler.dy_left * self.movement_speed)
        self.input_handler.dy_left = 0

        self.camera_pos[1] -= self.input_handler.dx_middle * self.mouse_sensitivity
        self.input_handler.dx_middle = 0

        self.camera_pos[2] += self.input_handler.dy_middle * self.mouse_sensitivity

        # Prevent the camera from flipping over
        epsilon = 0.1
        self.camera_pos[2] = max(epsilon, min(np.pi - epsilon, self.camera_pos[2]))
        self.input_handler.dy_middle = 0

        min_zoom = 0.05
        self.camera_pos[0] -= self.input_handler.scroll_y
        self.camera_pos[0] = max(min_zoom, self.camera_pos[0])
        self.input_handler.scroll_y = 0

    def draw(self):
        # Convert to Cartesian coordinates from spherical coordinates
        camera_trans = self.camera_pos[0] * np.array([
            np.sin(self.camera_pos[2]) * np.cos(self.camera_pos[1]),
            np.cos(self.camera_pos[2]),
            np.sin(self.camera_pos[2]) * np.sin(self.camera_pos[1])
        ])

        # Calculate unit vectors
        self.forward_dir = camera_trans / self.camera_pos[0]
        self.right_dir = np.cross(self.forward_dir, self.global_up_dir)
        self.local_up_dir = np.cross(self.forward_dir, self.right_dir)

        # Rotate and translate the camera
        gluLookAt(0, 0, 0, *self.forward_dir, *self.global_up_dir)
        glTranslatef(*camera_trans)
        glTranslatef(*self.pivot_pos)
