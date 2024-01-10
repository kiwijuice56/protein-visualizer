import h5py
import pyglet
from sklearn.manifold import TSNE
from pyglet.gl import *


class EmbeddingRenderer:
    def __init__(self, embedding_path, sequence_path, window, point_size=8):
        self.window = window

        embedding_file = h5py.File(embedding_path, 'r')
        sequence_file = open(sequence_path, 'r')
        title = sequence_file.readline()[1:].strip()
        sequence = sequence_file.readline().strip()
        embeddings = embedding_file[title][()]

        self.transform = TSNE(n_components=2, learning_rate='auto', init='random', perplexity=3).fit_transform(embeddings)
        self.points = self.transform.flatten().tolist()

        # Padding
        min_point, max_point = self.transform.min(axis=0) - [4.0, 4.0], self.transform.max(axis=0) + [4.0, 4.0]

        self.data_bounding_box = [min_point[0], min_point[1], max_point[0], max_point[1]]
        self.data_height = self.data_bounding_box[3] - self.data_bounding_box[1]
        self.data_width = self.data_bounding_box[2] - self.data_bounding_box[1]

        self.norm_points = []
        for x, y in zip(self.points[::2], self.points[1::2]):
            dist_x = x - self.data_bounding_box[0]
            dist_y = y - self.data_bounding_box[1]

            norm_x = dist_x / self.data_width
            norm_y = dist_y / self.data_height

            self.norm_points.extend([norm_x, norm_y])

        self.scaled_points = []
        self.bounding_box = []
        self.vertices = None

        self.set_bounding_box([0, -window.height, window.width, window.height])

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
        glScissor(self.bounding_box[0], self.bounding_box[1] + self.window.height, self.bounding_box[2], self.bounding_box[3])
        glEnable(GL_SCISSOR_TEST)

        glClearColor(1.0, 1.0, 1.0, 1.0)
        glEnable(GL_POINT_SMOOTH)
        glPointSize(self.point_size)

        glLoadIdentity()
        glDisable(GL_DEPTH_TEST)
        glEnable(GL_BLEND)
        glOrtho(0, self.window.width, -self.window.height, 0, 0, 1000)

        self.vertices.draw(pyglet.gl.GL_POINTS)
        glDisable(GL_SCISSOR_TEST)
