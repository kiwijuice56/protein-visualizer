import numpy as np
import pyglet
from pyglet.gl import *

from gui.input_handler import InputHandler


class Camera2D(object):
    ZOOM_RANGE = (0.01, 8)

    def __init__(self, window, movement_speed=0.003, mouse_sensitivity=0.5, scroll_sensitivity=0.07):
        self.pos = np.zeros(2)
        self.scale = 1

        self.input_handler = InputHandler()

        window.push_handlers(self.input_handler)

        self.movement_speed = movement_speed
        self.mouse_sensitivity = mouse_sensitivity
        self.scroll_sensitivity = scroll_sensitivity

        self.updated = False

    def update(self):
        self.updated = self.input_handler.dx_left != 0 or self.input_handler.dy_left != 0 or self.input_handler.scroll_y != 0

        self.pos[0] -= self.input_handler.dx_left * self.movement_speed * self.mouse_sensitivity * self.scale
        self.input_handler.dx_left = 0

        self.pos[1] -= self.input_handler.dy_left * self.movement_speed * self.mouse_sensitivity * self.scale
        self.input_handler.dy_left = 0

        self.scale -= self.scroll_sensitivity * self.scale * self.input_handler.scroll_y
        self.scale = max(self.ZOOM_RANGE[0], min(self.ZOOM_RANGE[1], self.scale))
        self.input_handler.scroll_y = 0


# Renders the 2D embeddings of a protein
class EmbeddingRenderer:
    POINT_SIZE_RANGE = (1, 16)
    POINT_OPACITY = 200
    BACKGROUND_COLOR = (0, 0, 0, 140)
    GRID_LINE_COUNT = 32
    GRID_LINE_COLOR = (255, 255, 255, 32)

    def __init__(self, protein, window, bounding_box=None, point_size=8):
        """
        @param protein: Reference to Protein object to render
        @param window: Reference to the parent pyglet.window.Window
        @param bounding_box: Rendering region as an array of format [bottom_left_x, bottom_left_y, width, height]
        @param point_size: Initial area of plotted points
        """
        self.protein = protein
        self.window = window
        self.camera = Camera2D(window)

        # Find the region the points lie in
        min_point = [min(self.protein.embedding_points[0::2]),
                     min(self.protein.embedding_points[1::2])]
        max_point = [max(self.protein.embedding_points[0::2]),
                     max(self.protein.embedding_points[1::2])]

        # Define a 2D space where the data is contained
        self.data_bounding_box = [min_point[0], min_point[1], max_point[0], max_point[1]]
        self.data_height = self.data_bounding_box[3] - self.data_bounding_box[1]
        self.data_width = self.data_bounding_box[2] - self.data_bounding_box[0]

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
        self.outline_vertices = None
        self.id_vertices = None

        self.hovered_residue = -1

        # Default to full-screen render
        self.bounding_box = []
        self.grid_list = None
        if bounding_box is None:
            self.set_bounding_box([0, 0, window.width, window.height])
        else:
            self.set_bounding_box(bounding_box)

        self.point_size = 0
        self.set_point_size(point_size)

    def set_point_size(self, new_size):
        self.point_size = max(self.POINT_SIZE_RANGE[0], min(self.POINT_SIZE_RANGE[1], new_size))

    def set_bounding_box(self, bounding_box):
        """
        Updates the bounding box and scales rendered vertices to fit
        @param bounding_box: Rendering region as an array of format [bottom_left_x, bottom_left_y, width, height]
        """
        self.bounding_box = bounding_box
        self.camera.input_handler.bounding_box = bounding_box

        world_box = [0, 0, 1, 1]

        world_box[0] += self.camera.pos[0]
        world_box[1] += self.camera.pos[1]

        world_box[0] += 0.5 - self.camera.scale / 2
        world_box[1] += 0.5 - self.camera.scale / 2

        world_box[2] *= self.camera.scale
        world_box[3] *= self.camera.scale

        for i in range(0, len(self.norm_points), 2):
            x = self.norm_points[i]
            y = self.norm_points[i + 1]

            u = (x - world_box[0]) / world_box[2]
            v = (y - world_box[1]) / world_box[3]

            self.scaled_points[i] = u * self.bounding_box[2] + self.bounding_box[0]
            self.scaled_points[i + 1] = v * self.bounding_box[3] + self.bounding_box[1]

        grid_points = []
        for i in range(0, self.GRID_LINE_COUNT + 1):
            size = 1.0 / self.GRID_LINE_COUNT
            u1 = (0.0 - world_box[0]) / world_box[2] * self.bounding_box[2] + self.bounding_box[0]
            u2 = (1.0 - world_box[0]) / world_box[2] * self.bounding_box[2] + self.bounding_box[0]
            v1 = (i * size - world_box[1]) / world_box[3] * self.bounding_box[3] + self.bounding_box[1]

            grid_points.extend([u1, v1, u2, v1])

            u3 = (i * size - world_box[0]) / world_box[2] * self.bounding_box[2] + self.bounding_box[0]
            v2 = (0.0 - world_box[1]) / world_box[3] * self.bounding_box[3] + self.bounding_box[1]
            v3 = (1.0 - world_box[1]) / world_box[3] * self.bounding_box[3] + self.bounding_box[1]

            grid_points.extend([u3, v2, u3, v3])

        self.grid_list = pyglet.graphics.vertex_list(
            len(grid_points) // 2,
            ('v2f', grid_points),
            ('c4B', self.GRID_LINE_COLOR * (len(grid_points) // 2))
        )

        rgb_id = [[(i >> 16) & 255, (i >> 8) & 255, i & 255] for i in range(1, len(self.scaled_points) // 2 + 1)]
        self.vertices = pyglet.graphics.vertex_list(
            len(self.protein.embedding_points) // 2,
            ('v2f', self.scaled_points),
            ('c4B', np.zeros(len(self.protein.embedding_points) // 2 * 4, dtype=np.byte)))
        self.id_vertices = pyglet.graphics.vertex_list(
            len(self.protein.embedding_points) // 2,
            ('v2f', self.scaled_points),
            ('c3B', [x for xs in rgb_id for x in xs]))
        self.update_colors()

    def update_colors(self, start=0, end=-1):
        """
        Updates the rendered colors within the given range
        @param start: Beginning residue index, inclusive
        @param end: Final residue index, exclusive
        """
        if end == -1:
            end = len(self.protein.residues)
        for i in range(start, end):
            self.vertices.colors[i * 4: i * 4 + 3] = self.protein.residues[i].color
            self.vertices.colors[i * 4 + 3] = self.POINT_OPACITY

    def detect_mouse(self):
        glClearColor(*[0] * 4)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glEnable(GL_SCISSOR_TEST)
        glDisable(GL_POINT_SMOOTH)

        glScissor(*self.bounding_box)

        glLoadIdentity()
        glDisable(GL_DEPTH_TEST)
        glEnable(GL_BLEND)
        glOrtho(0, self.window.width, 0, self.window.height, 0, 1000)

        if self.camera.updated:
            self.set_bounding_box(self.bounding_box)
            self.camera.updated = False

        glPointSize(int(self.point_size * 2.0))
        self.id_vertices.draw(pyglet.gl.GL_POINTS)

        read = (GLubyte * 3)(0)
        glReadPixels(self.camera.input_handler.x, self.camera.input_handler.y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, read)
        selected_id = sum([read[i] << (8 * (2 - i)) for i in range(3)])
        if not selected_id == 0:
            self.hovered_residue = selected_id - 1
        else:
            self.hovered_residue = -1

        # Clean up
        glDisable(GL_SCISSOR_TEST)
        glDisable(GL_DEPTH_TEST)
        glDisable(GL_BLEND)

    def draw(self):
        glEnable(GL_SCISSOR_TEST)
        glEnable(GL_POINT_SMOOTH)

        glScissor(*self.bounding_box)

        glLoadIdentity()
        glDisable(GL_DEPTH_TEST)
        glEnable(GL_BLEND)
        glOrtho(0, self.window.width, 0, self.window.height, 0, 1000)

        if self.camera.updated:
            self.set_bounding_box(self.bounding_box)
            self.camera.updated = False

        background = pyglet.shapes.Rectangle(*self.bounding_box, color=self.BACKGROUND_COLOR[0:3])
        background.opacity = self.BACKGROUND_COLOR[3]
        background.draw()

        self.grid_list.draw(pyglet.gl.GL_LINES)

        glPointSize(self.point_size)
        self.vertices.draw(pyglet.gl.GL_POINTS)

        # Clean up
        glDisable(GL_SCISSOR_TEST)
        glDisable(GL_POINT_SMOOTH)
        glDisable(GL_DEPTH_TEST)
        glDisable(GL_BLEND)
