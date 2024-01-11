import pyglet
from pyglet.gl import *


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
        point_coordinates = []
        for atom in self.protein.atoms:
            point_coordinates.extend(atom.bio_atom.get_coord())

        self.atom_vertices = pyglet.graphics.vertex_list(
            len(self.protein.atoms),
            ('v3f', point_coordinates),
            ('c3B', [0] * len(self.protein.atoms) * 3))
        self.update_colors()

        self.outline_vertices = pyglet.graphics.vertex_list(
            len(self.protein.atoms),
            ('v3f', point_coordinates),
            ('c3B', [32] * len(self.protein.atoms) * 3))

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

    def draw(self):
        glEnable(GL_POINT_SMOOTH)
        glEnable(GL_SCISSOR_TEST)

        glClearColor(1.0, 1.0, 1.0, 1.0)
        glScissor(self.bounding_box[0], self.bounding_box[1] + self.window.height, self.bounding_box[2], self.bounding_box[3])

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
