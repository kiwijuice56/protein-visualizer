from gui.input_handler import InputHandler

import numpy as np
import pyglet
from pyglet.gl import *


class Camera3D(object):
    ZOOM_RANGE = (0.01, 512)
    VERTICAL_ROTATION_LIMIT = 0.1

    def __init__(self, window, movement_speed=0.45, mouse_sensitivity=0.003, scroll_sensitivity=0.1):
        """
        @param window: Reference to the parent pyglet.window.Window
        @param movement_speed: Translation speed
        @param mouse_sensitivity: Dampening amount to any mouse cursor movement
        @param scroll_sensitivity: Dampening amount to mouse scrolling
        """
        self.pivot_pos = np.array([0, 0, -32])
        self.camera_pos = np.array([12, 0, np.pi/2])

        self.forward_dir = np.array([0, 0, 1])
        self.right_dir = np.array([0, 0, 0])
        self.global_up_dir = np.array([0, 1, 0])
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
        self.move_horizontal(self.input_handler.dx_left * self.movement_speed * self.mouse_sensitivity * max(1.0, self.camera_pos[0]))
        self.input_handler.dx_left = 0

        self.move_vertical(-self.input_handler.dy_left * self.movement_speed * self.mouse_sensitivity * max(1.0, self.camera_pos[0]))
        self.input_handler.dy_left = 0

        self.camera_pos[1] += self.input_handler.dx_middle * self.mouse_sensitivity
        self.input_handler.dx_middle = 0

        self.camera_pos[2] -= self.input_handler.dy_middle * self.mouse_sensitivity

        # Prevent the camera from flipping over
        self.camera_pos[2] = max(self.VERTICAL_ROTATION_LIMIT, min(np.pi - self.VERTICAL_ROTATION_LIMIT, self.camera_pos[2]))
        self.input_handler.dy_middle = 0

        self.camera_pos[0] -= self.input_handler.scroll_y * self.scroll_sensitivity * self.camera_pos[0]
        self.camera_pos[0] = max(self.ZOOM_RANGE[0], min(self.ZOOM_RANGE[1], self.camera_pos[0]))
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


# Renders the 3D structure of a protein
class PDBRenderer:
    GRID_LINE_COUNT = 16
    BACKGROUND_COLOR = (98, 98, 98, 255)
    X_AXIS_COLOR = (255, 56, 89, 98)
    Z_AXIS_COLOR = (56, 109, 255, 98)
    GRID_LINE_COLOR = (255, 255, 255, 32)
    POINT_SIZE_RANGE = (1, 20)
    FOV = 65  # Degrees
    Z_NEAR = 2
    Z_FAR = 800

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
        self.hovered_residue = -1

        if bounding_box is None:
            self.bounding_box = (0, -window.height, window.width, window.height)
        else:
            self.bounding_box = bounding_box

        self.residue_type = {}

        # Create vertex data for OpenGL
        lowest_atom, offset = 0, self.protein.atoms[0].bio_atom.get_coord()
        for i, atom in enumerate(self.protein.atoms):
            if atom.bio_atom.get_coord()[1] < self.protein.atoms[lowest_atom].bio_atom.get_coord()[1]:
                lowest_atom = i
                offset = self.protein.atoms[lowest_atom].bio_atom.get_coord()
        point_coordinates = np.zeros(len(self.protein.atoms) * 3)
        for i, atom in enumerate(self.protein.atoms):
            point_coordinates[i * 3] = atom.bio_atom.get_coord()[0] - offset[0]
            point_coordinates[i * 3 + 1] = atom.bio_atom.get_coord()[1] - offset[1]
            point_coordinates[i * 3 + 2] = atom.bio_atom.get_coord()[2] - offset[2]

        rgb_id = [[(i >> 16) & 255, (i >> 8) & 255, i & 255] for i in range(1, len(self.protein.atoms) + 1)]
        self.atom_vertices = pyglet.graphics.vertex_list(
            len(self.protein.atoms),
            ('v3f', point_coordinates),
            ('c3B', np.zeros(len(self.protein.atoms) * 3, dtype=np.byte)))
        self.id_vertices = pyglet.graphics.vertex_list(
            len(self.protein.atoms),
            ('v3f', point_coordinates),
            ('c3B', [x for xs in rgb_id for x in xs]))
        self.update_colors()

        self.outline_vertices = pyglet.graphics.vertex_list(
            len(self.protein.atoms),
            ('v3f', point_coordinates),
            ('c3B', np.full(len(self.protein.atoms) * 3, 32, dtype=np.byte)))

        self.camera = Camera3D(window)

        grid_points, colors = [], []
        for i in range(0, self.GRID_LINE_COUNT + 1):
            size = 1.0 / self.GRID_LINE_COUNT

            grid_points.extend([-128, -8, i * size * 256 - 128, 128, -8, i * size * 256 - 128])
            grid_points.extend([i * size * 256 - 128, -8, -128, i * size * 256 - 128, -8, 128])
            if i == self.GRID_LINE_COUNT / 2:
                colors.extend(self.X_AXIS_COLOR * 2)
                colors.extend(self.Z_AXIS_COLOR * 2)
            else:
                colors.extend(self.GRID_LINE_COLOR * 4)

        self.grid_list = pyglet.graphics.vertex_list(
            len(grid_points) // 3,
            ('v3f', grid_points),
            ('c4B', colors)
        )

    def update_colors(self, start=0, end=-1):
        """
        Updates the rendered colors within the given range
        @param start: Beginning atom index, inclusive
        @param end: Final atom index, exclusive
        """
        if end == -1:
            end = len(self.protein.atoms)
        for i in range(start, end):
            self.atom_vertices.colors[i * 3: i * 3 + 3] = self.protein.atoms[i].color

    def set_point_size(self, new_size):
        self.point_size = max(self.POINT_SIZE_RANGE[0], min(self.POINT_SIZE_RANGE[1], new_size))

    def set_bounding_box(self, bounding_box):
        """
        Updates the bounding box and scales rendered vertices to fit
        @param bounding_box: Rendering region as an array of format [bottom_left_x, bottom_left_y, width, height]
        """
        self.bounding_box = bounding_box
        self.camera.input_handler.bounding_box = bounding_box

    def draw(self):
        glClearColor(*[0] * 4)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glEnable(GL_POINT_SMOOTH)
        glEnable(GL_SCISSOR_TEST)

        glLoadIdentity()
        glMatrixMode(GL_PROJECTION)
        gluPerspective(self.FOV, self.window.width / float(self.window.height), self.Z_NEAR, self.Z_FAR)

        self.camera.draw()

        glScissor(*self.bounding_box)

        glPointSize(int(1.5 * self.point_size))
        glDisable(GL_POINT_SMOOTH)
        glEnable(GL_DEPTH_TEST)
        self.id_vertices.draw(GL_POINTS)

        read = (GLubyte * 3)(0)
        glReadPixels(self.camera.input_handler.x, self.camera.input_handler.y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, read)
        selected_id = sum([read[i] << (8 * (2 - i)) for i in range(3)])
        if not selected_id == 0:
            self.hovered_residue = self.protein.atoms[selected_id - 1].residue.index
        else:
            self.hovered_residue = -1

        glClearColor(*[c / 255.0 for c in self.BACKGROUND_COLOR])
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glEnable(GL_LINE_SMOOTH)
        self.grid_list.draw(GL_LINES)

        glEnable(GL_POINT_SMOOTH)
        if self.outline:
            glDisable(GL_DEPTH_TEST)
            glPointSize(int(self.point_size * 1.5))
            self.outline_vertices.draw(GL_POINTS)

        glPointSize(self.point_size)
        glEnable(GL_DEPTH_TEST)
        self.atom_vertices.draw(GL_POINTS)

        # Clean up
        glDisable(GL_SCISSOR_TEST)
        glDisable(GL_LINE_SMOOTH)
        glDisable(GL_POINT_SMOOTH)
        glDisable(GL_DEPTH_TEST)

    def update(self):
        self.camera.update()
