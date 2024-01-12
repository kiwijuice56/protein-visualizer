from input_handler import InputHandler

import numpy as np
import pyglet
from pyglet.gl import *


class Camera3D(object):
    def __init__(self, window, movement_speed=0.45, mouse_sensitivity=0.01, scroll_sensitivity=0.1):
        self.pivot_pos = np.array([0, 0, -32])
        self.camera_pos = np.array([12, 0, np.pi/2])

        self.forward_dir = np.array([0, 0, 1])
        self.right_dir = np.array([0, 0, 0])
        self.global_up_dir = np.array([0, -1, 0])
        self.local_up_dir = np.array([0, -1, 0])
        self.scroll_sensitivity = scroll_sensitivity

        self.input_handler = InputHandler()

        window.push_handlers(self.input_handler)

        self.movement_speed = movement_speed
        self.mouse_sensitivity = mouse_sensitivity

    def move_horizontal(self, distance):
        self.pivot_pos = self.pivot_pos + self.right_dir * distance

    def move_vertical(self, distance):
        self.pivot_pos = self.pivot_pos + self.local_up_dir * distance

    def update(self):
        self.move_horizontal(self.input_handler.dx_left * self.movement_speed * self.mouse_sensitivity * self.camera_pos[0])
        self.input_handler.dx_left = 0

        self.move_vertical(-self.input_handler.dy_left * self.movement_speed * self.mouse_sensitivity * self.camera_pos[0])
        self.input_handler.dy_left = 0

        self.camera_pos[1] -= self.input_handler.dx_middle * self.mouse_sensitivity
        self.input_handler.dx_middle = 0

        self.camera_pos[2] += self.input_handler.dy_middle * self.mouse_sensitivity

        # Prevent the camera from flipping over
        epsilon = 0.1
        self.camera_pos[2] = max(epsilon, min(np.pi - epsilon, self.camera_pos[2]))
        self.input_handler.dy_middle = 0

        min_zoom = 0.05
        self.camera_pos[0] -= self.input_handler.scroll_y * self.scroll_sensitivity * self.camera_pos[0]
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


class PDBRenderer:
    def __init__(self, protein, window, bounding_box=None, point_size=8, outline=True):
        """
        @param protein: Reference to Protein object to render
        @param window: Reference to the parent pyglet.window.Window
        @param bounding_box: Rendering region as an array of format [bottom_left_x, bottom_left_y, width, height]
        @param point_size: Initial area of plotted points
        @param outline: Whether plotted points have an outline
        """
        self.window = window
        self.protein = protein

        self.point_size = 0
        self.set_point_size(point_size)
        self.outline = outline

        if bounding_box is None:
            self.bounding_box = (0, -window.height, window.width, window.height)
        else:
            self.bounding_box = bounding_box

        self.residue_type = {}

        # Create vertex data for OpenGL
        point_coordinates = np.zeros(len(self.protein.atoms) * 3)
        for i, atom in enumerate(self.protein.atoms):
            point_coordinates[i * 3] = atom.bio_atom.get_coord()[0]
            point_coordinates[i * 3 + 1] = atom.bio_atom.get_coord()[1]
            point_coordinates[i * 3 + 2] = atom.bio_atom.get_coord()[2]

        self.atom_vertices = pyglet.graphics.vertex_list(
            len(self.protein.atoms),
            ('v3f', point_coordinates),
            ('c3B', np.zeros(len(self.protein.atoms) * 3, dtype=np.byte)))
        self.update_colors()

        self.outline_vertices = pyglet.graphics.vertex_list(
            len(self.protein.atoms),
            ('v3f', point_coordinates),
            ('c3B', np.full(len(self.protein.atoms) * 3, 32, dtype=np.byte)))

        self.camera = Camera3D(window)

    def color_atom(self, atom):
        """
        Returns the color of a particular atom depending on the renderer settings
        @param atom: The Atom to color
        @return: A length-3 RGB array, [0, 255]
        """
        return atom.residue.color

    def update_colors(self, start=0, end=-1):
        """
        Updates the rendered colors within the given range
        @param start: Beginning atom index, inclusive
        @param end: Final atom index, exclusive
        """
        if end == -1:
            end = len(self.protein.atoms)
        for i in range(start, end):
            self.atom_vertices.colors[i * 3: i * 3 + 3] = self.color_atom(self.protein.atoms[i])

    def set_point_size(self, new_size):
        if new_size <= 0:
            new_size = 1
        if new_size > 10:
            new_size = 10
        self.point_size = new_size

    def set_bounding_box(self, bounding_box):
        """
        Updates the bounding box and scales rendered vertices to fit
        @param bounding_box: Rendering region as an array of format [bottom_left_x, bottom_left_y, width, height]
        """
        self.bounding_box = bounding_box
        self.camera.input_handler.bounding_box = bounding_box

    def draw(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glEnable(GL_POINT_SMOOTH)
        glEnable(GL_SCISSOR_TEST)

        glLoadIdentity()
        glMatrixMode(GL_PROJECTION)
        gluPerspective(65, self.window.width / float(self.window.height), 0.8, 512)

        self.camera.draw()

        glClearColor(1.0, 1.0, 1.0, 1.0)
        glScissor(*self.bounding_box)

        if self.outline:
            glDisable(GL_DEPTH_TEST)
            glPointSize(self.point_size * 2)
            self.outline_vertices.draw(pyglet.gl.GL_POINTS)

        glPointSize(self.point_size)
        glEnable(GL_DEPTH_TEST)
        self.atom_vertices.draw(pyglet.gl.GL_POINTS)

        # Clean up
        glDisable(GL_SCISSOR_TEST)
        glDisable(GL_POINT_SMOOTH)
        glDisable(GL_DEPTH_TEST)

    def update(self):
        self.camera.update()