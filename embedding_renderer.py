import h5py
import pyglet
from sklearn.manifold import TSNE
from pyglet.gl import *


class EmbeddingRenderer:
    def __init__(self, embedding_path, sequence_path, window, bounding_box=None, point_size=8):
        """
        @param embedding_path: Path to a `.h5` file containing a single database with the protein's embeddings
        @param sequence_path: Path to a `.fa` file (FASTA format) containing the protein's amino acid sequence
        @param window: Reference to the parent pyglet.window.Window
        @param bounding_box: Rendering region as an array of format [bottom_left_x, bottom_left_y, width, height]
        @param point_size: Initial area of plotted points
        """
        self.window = window

        embedding_file = h5py.File(embedding_path, 'r')
        sequence_file = open(sequence_path, 'r')

        title = sequence_file.readline()[1:].strip()
        sequence = sequence_file.readline().strip()

        embeddings = embedding_file[title][()]

        # Use the t-SNE algorithm to transform the embeddings into 2D vectors
        self.transform = TSNE(n_components=2, learning_rate='auto', init='random', perplexity=3).fit_transform(embeddings)
        self.points = self.transform.flatten().tolist()

        # Find the region the points lie in, with some added padding
        min_point, max_point = self.transform.min(axis=0) - [4.0, 4.0], self.transform.max(axis=0) + [4.0, 4.0]

        # Define a 2D space where the data is contained
        self.data_bounding_box = [min_point[0], min_point[1], max_point[0], max_point[1]]
        self.data_height = self.data_bounding_box[3] - self.data_bounding_box[1]
        self.data_width = self.data_bounding_box[2] - self.data_bounding_box[1]

        # Normalize the points to the domain [0, 1] for more convenient plotting
        self.norm_points = []
        for x, y in zip(self.points[::2], self.points[1::2]):
            dist_x = x - self.data_bounding_box[0]
            dist_y = y - self.data_bounding_box[1]

            norm_x = dist_x / self.data_width
            norm_y = dist_y / self.data_height

            self.norm_points.extend([norm_x, norm_y])

        # Define a list of points scaled to the rendering bounding box
        self.scaled_points = []

        # Define the pyglet vertex list
        self.vertices = None

        # Default to full-screen render
        self.bounding_box = []
        if bounding_box is None:
            self.set_bounding_box([0, -window.height, window.width, window.height])
        else:
            self.set_bounding_box(bounding_box)

        self.point_size = point_size

    def set_bounding_box(self, new_box):
        self.bounding_box = new_box

        self.scaled_points = []
        for x, y in zip(self.norm_points[::2], self.norm_points[1::2]):
            scaled_x = x * self.bounding_box[2] + self.bounding_box[0]
            scaled_y = y * self.bounding_box[2] + self.bounding_box[1]
            self.scaled_points.extend([scaled_x, scaled_y])
        self.vertices = pyglet.graphics.vertex_list(
            len(self.points) // 2,
            ('v2f', self.scaled_points),
            ('c3B', [64] * (len(self.points) // 2 * 3)))

    def draw(self):
        glEnable(GL_SCISSOR_TEST)
        glEnable(GL_POINT_SMOOTH)

        glScissor(self.bounding_box[0], self.bounding_box[1] + self.window.height, self.bounding_box[2], self.bounding_box[3])
        glClearColor(1.0, 1.0, 1.0, 1.0)
        glPointSize(self.point_size)

        glLoadIdentity()
        glDisable(GL_DEPTH_TEST)
        glEnable(GL_BLEND)
        glOrtho(0, self.window.width, -self.window.height, 0, 0, 1000)

        self.vertices.draw(pyglet.gl.GL_POINTS)

        # Clean up
        glDisable(GL_SCISSOR_TEST)
        glDisable(GL_POINT_SMOOTH)
        glDisable(GL_DEPTH_TEST)
        glDisable(GL_BLEND)
