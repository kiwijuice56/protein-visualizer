import Bio.PDB
import Bio.SeqRecord

import pyglet
from pyglet.gl import *

from enum import Enum


# Wrapper class for collections of atoms
class Residue:
    def __init__(self, atoms, bio_residue, index):
        self.atoms = atoms
        self.bio_residue = bio_residue
        self.index = index


# Wrapper class for single atoms
class Atom:
    def __init__(self, bio_atom, residue, index):
        self.bio_atom = bio_atom
        self.residue = residue
        self.index = index


class PDBRenderer:
    class ColorMode(Enum):
        CPK = 0
        RESIDUE = 1
        CONTRAST = 2

    def __init__(self, pdb_path, window, color_mode=ColorMode.CPK, point_size=8, outline=True):
        self.window = window

        pdbparser = Bio.PDB.PDBParser(QUIET=True)
        self.bio_structure = pdbparser.get_structure("id", pdb_path)
        self.residues = []
        self.atoms = []

        self.highlighted_index = 0
        self.color_mode = color_mode
        self.point_size = point_size
        self.outline = outline

        # Initialize data from PDB file
        atom_index = 0
        for i, bio_residue in enumerate(self.bio_structure.get_residues()):
            residue_atoms = []
            residue = Residue(residue_atoms, bio_residue, i)
            for bio_atom in bio_residue.get_atoms():
                residue_atoms.append(Atom(bio_atom, residue, atom_index))
                atom_index += 1
            self.residues.append(residue)
            self.atoms.extend(residue_atoms)

        # Create vertex data for OpenGL
        point_coordinates = []
        for atom in self.atoms:
            point_coordinates.extend(atom.bio_atom.get_coord())

        self.atom_vertices = pyglet.graphics.vertex_list(
            len(self.atoms),
            ('v3f', point_coordinates),
            ('c3B', [0] * len(self.atoms) * 3))
        self.update_colors()

        self.outline_vertices = pyglet.graphics.vertex_list(
            len(self.atoms),
            ('v3f', point_coordinates),
            ('c3B', [32] * len(self.atoms) * 3))

    def color_atom(self, atom, highlight=False):
        # Based on CPK coloring (https://en.wikipedia.org/wiki/CPK_coloring)
        cpk_colors = {"C": (64, 58, 64), "O": (219, 73, 70), "N": (70, 110, 219), "S": (235, 208, 56),
                      "P": (235, 145, 56),
                      "_": (255, 255, 255)}
        residue_colors = [(191, 53, 15), (204, 135, 24), (212, 209, 25), (110, 204, 33), (33, 204, 130), (33, 178, 204),
                          (36, 51, 212), (133, 36, 212), (212, 36, 145)]
        base_color = (128, 128, 128)

        match self.color_mode:
            case self.ColorMode.CPK:
                bio_id = atom.bio_atom.get_id()
                base_color = cpk_colors[bio_id[0]] if bio_id[0] in cpk_colors else cpk_colors['_']
            case self.ColorMode.RESIDUE:
                base_color = residue_colors[atom.residue.index % len(residue_colors)]
            case self.ColorMode.CONTRAST:
                base_color = (255, 213, 25) if highlight else (46, 29, 115)

        return tuple([min(255, b + 128) for b in base_color]) if highlight else base_color

    def update_colors(self, start=0, end=-1):
        if end == -1:
            end = len(self.atoms)
        for i in range(start, end):
            is_highlighted = self.atoms[i] in self.residues[self.highlighted_index].atoms
            self.atom_vertices.colors[i * 3: i * 3 + 3] = self.color_atom(self.atoms[i], is_highlighted)

    def set_color_mode(self, new_color_mode):
        self.color_mode = new_color_mode
        self.update_colors()

    def set_highlighted_index(self, new_index):
        old_index = self.highlighted_index
        self.highlighted_index = new_index

        self.update_colors(self.residues[old_index].atoms[0].index, self.residues[old_index].atoms[-1].index + 1)
        self.update_colors(self.residues[new_index].atoms[0].index, self.residues[new_index].atoms[-1].index + 1)

    def set_point_size(self, new_size):
        if new_size <= 0:
            new_size = 1
        if new_size > 10:
            new_size = 10
        self.point_size = new_size

    def draw(self):
        glClearColor(1.0, 1.0, 1.0, 1.0)

        glEnable(GL_VERTEX_PROGRAM_POINT_SIZE)
        glEnable(GL_POINT_SMOOTH)

        if self.outline:
            glDisable(GL_DEPTH_TEST)
            glPointSize(self.point_size * 2)
            self.outline_vertices.draw(pyglet.gl.GL_POINTS)

        glPointSize(self.point_size)
        glEnable(GL_DEPTH_TEST)
        self.atom_vertices.draw(pyglet.gl.GL_POINTS)
