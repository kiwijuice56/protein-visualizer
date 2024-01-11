import numpy as np
import pyglet
from pyglet.gl import *

from input_handler import InputHandler


class Camera2D(object):
    def __init__(self, window, movement_speed=0.5, mouse_sensitivity=0.01):
        self.pos = np.zeros(2)

        self.input_handler = InputHandler()

        window.push_handlers(self.input_handler)

        self.movement_speed = movement_speed
        self.mouse_sensitivity = mouse_sensitivity

    def update(self):
        self.pos[0] += self.input_handler.dx_left * self.movement_speed
        self.input_handler.dx_left = 0

        self.pos[1] += self.input_handler.dy_left * self.movement_speed
        self.input_handler.dy_left = 0


class EmbeddingRenderer:
    def __init__(self, protein, window, bounding_box=None, point_size=6):
        """
        @param protein: Reference to Protein object to render
        @param window: Reference to the parent pyglet.window.Window
        @param bounding_box: Rendering region as an array of format [bottom_left_x, bottom_left_y, width, height]
        @param point_size: Initial area of plotted points
        """
        self.protein = protein
        self.window = window
        self.camera = Camera2D(window)

        # Find the region the points lie in, with extra padding to prevent points from clipping against the screen
        padding_amount = 4
        min_point = [min(self.protein.embedding_points[0::2]),
                     min(self.protein.embedding_points[1::2])]
        max_point = [max(self.protein.embedding_points[0::2]) + padding_amount,
                     max(self.protein.embedding_points[1::2]) + padding_amount]

        # Define a 2D space where the data is contained
        self.data_bounding_box = [min_point[0], min_point[1], max_point[0], max_point[1]]
        self.data_height = self.data_bounding_box[3] - self.data_bounding_box[1]
        self.data_width = self.data_bounding_box[2] - self.data_bounding_box[1]

        # Normalize the points to the domain [0, 1] for more convenient plotting
        self.norm_points = np.zeros(len(self.protein.embedding_points))
        for i in range(0, len(self.protein.embedding_points), 2):
            x = self.protein.embedding_points[i]
            y = self.protein.embedding_points[i + 1]
            dist_x = x - self.data_bounding_box[0]
            dist_y = y - self.data_bounding_box[1]

            self.norm_points[i] = dist_x / self.data_width
            self.norm_points[i + 1] = dist_y / self.data_height

        # Define a list of points scaled to the rendering bounding box
        self.scaled_points = np.zeros(len(self.protein.embedding_points))

        # Define the pyglet vertex list
        self.vertices = None

        # Default to full-screen render
        self.bounding_box = []
        if bounding_box is None:
            self.set_bounding_box([0, 0, window.width, window.height])
        else:
            self.set_bounding_box(bounding_box)

        self.point_size = point_size

    def set_bounding_box(self, bounding_box):
        """
        Updates the bounding box and scales rendered vertices to fit
        @param bounding_box: Rendering region as an array of format [bottom_left_x, bottom_left_y, width, height]
        """
        self.bounding_box = bounding_box
        self.camera.input_handler.bounding_box = bounding_box

        for i in range(0, len(self.norm_points), 2):
            x = self.norm_points[i]
            y = self.norm_points[i + 1]

            self.scaled_points[i] = self.camera.pos[0] + x * self.bounding_box[2] + self.bounding_box[0]
            self.scaled_points[i + 1] = self.camera.pos[1] + y * self.bounding_box[2] + self.bounding_box[1]
        self.vertices = pyglet.graphics.vertex_list(
            len(self.protein.embedding_points) // 2,
            ('v2f', self.scaled_points),
            ('c3B', [0] * (len(self.protein.embedding_points) // 2 * 3)))
        self.update_colors()

    def color_residue(self, residue):
        """
        Returns the color of a particular residue depending on the renderer settings
        @param residue: The Residue to color
        @return: A length-3 RGB array, [0, 255]
        """
        return residue.color

    def update_colors(self, start=0, end=-1):
        """
        Updates the rendered colors within the given range
        @param start: Beginning residue index, inclusive
        @param end: Final residue index, exclusive
        """
        if end == -1:
            end = len(self.protein.residues)
        for i in range(start, end):
            self.vertices.colors[i * 3: i * 3 + 3] = self.color_residue(self.protein.residues[i])

    def draw(self):
        glEnable(GL_SCISSOR_TEST)
        glEnable(GL_POINT_SMOOTH)

        glScissor(*self.bounding_box)
        glPointSize(self.point_size)

        glLoadIdentity()
        glDisable(GL_DEPTH_TEST)
        glEnable(GL_BLEND)
        glOrtho(0, self.window.width, 0, self.window.height, 0, 1000)

        background = pyglet.shapes.Rectangle(*self.bounding_box, color=(255, 255, 255))
        background.opacity = 210
        background.draw()
        self.vertices.draw(pyglet.gl.GL_POINTS)

        # Clean up
        glDisable(GL_SCISSOR_TEST)
        glDisable(GL_POINT_SMOOTH)
        glDisable(GL_DEPTH_TEST)
        glDisable(GL_BLEND)
