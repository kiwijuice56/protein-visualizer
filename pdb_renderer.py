import pyglet
from pyglet.gl import *

from enum import Enum


class PDBRenderer:
    class ColorMode(Enum):
        CPK = 0
        CHAINBOW = 1
        CONTRAST = 2

    def __init__(self, protein, window, bounding_box=None, color_mode=ColorMode.CPK, point_size=8, outline=True):
        """
        @param protein: Reference to Protein object to render
        @param window: Reference to the parent pyglet.window.Window
        @param color_mode: Initial color mode of plotted points
        @param bounding_box: Rendering region as an array of format [bottom_left_x, bottom_left_y, width, height]
        @param point_size: Initial area of plotted points
        @param outline: Whether plotted points have an outline
        """
        self.window = window
        self.protein = protein

        self.highlighted_index = 0
        self.color_mode = color_mode
        self.point_size = point_size
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

    def color_atom(self, atom, highlight=False):
        # Based on CPK coloring (https://en.wikipedia.org/wiki/CPK_coloring)
        cpk_colors = {"C": (64, 58, 64), "O": (219, 73, 70), "N": (70, 110, 219), "S": (235, 208, 56),
                      "P": (235, 145, 56),
                      "_": (255, 255, 255)}
        residue_colors = [(191, 53, 15), (204, 135, 24), (212, 209, 25), (110, 204, 33), (33, 204, 130), (33, 178, 204),
                          (36, 51, 212), (133, 36, 212), (212, 36, 145)]
        poisson_colors = [(187, 176, 148), (128, 118, 101), (89, 82, 70), (51, 51, 51), (25, 31, 34), (47, 68, 67),
                          (59, 94, 88), (90, 140, 108), (139, 180, 141), (192, 208, 165), (247, 239, 199),
                          (161, 205, 176), (112, 147, 149), (74, 120, 123), (56, 49, 64), (115, 77, 92),
                          (167, 103, 114), (204, 134, 125), (224, 186, 139), (195, 130, 82), (161, 86, 60),
                          (111, 52, 45), (68, 39, 31)]
        base_color = (128, 128, 128)

        match self.color_mode:
            case self.ColorMode.CPK:
                bio_id = atom.bio_atom.get_id()
                base_color = cpk_colors[bio_id[0]] if bio_id[0] in cpk_colors else cpk_colors['_']
            case self.ColorMode.CHAINBOW:
                base_color = residue_colors[atom.residue.index % len(residue_colors)]
            case self.ColorMode.CONTRAST:
                base_color = (255, 213, 25) if highlight else (46, 29, 115)

        return tuple([min(255, b + 128) for b in base_color]) if highlight else base_color

    def update_colors(self, start=0, end=-1):
        if end == -1:
            end = len(self.protein.atoms)
        for i in range(start, end):
            is_highlighted = self.protein.atoms[i] in self.protein.residues[self.highlighted_index].atoms
            self.atom_vertices.colors[i * 3: i * 3 + 3] = self.color_atom(self.protein.atoms[i], is_highlighted)

    def set_color_mode(self, new_color_mode):
        self.color_mode = new_color_mode
        self.update_colors()

    def set_highlighted_index(self, new_index):
        old_index = self.highlighted_index
        self.highlighted_index = new_index

        self.update_colors(self.protein.residues[old_index].atoms[0].index, self.protein.residues[old_index].atoms[-1].index + 1)
        self.update_colors(self.protein.residues[new_index].atoms[0].index, self.protein.residues[new_index].atoms[-1].index + 1)

    def set_point_size(self, new_size):
        if new_size <= 0:
            new_size = 1
        if new_size > 10:
            new_size = 10
        self.point_size = new_size

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
