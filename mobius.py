from collections import defaultdict, namedtuple, deque
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve, lsqr
import sys

Point = namedtuple('Point', ['x', 'y', 'z'])

class Mesh:
    def __init__(self, vertices=None, faces=None):
        if not vertices: vertices = []
        if not faces: faces = []
        self.vertices = vertices.copy()
        self.faces = faces.copy()

        self.u = None
        self.ul = None
        self.me_vertices = None
        self.me_faces = None
        self.me_vertex_to_edge = None
        self.me_neighbors = None
        self.edge_to_me_vertex = None
        self.cut_face = None
        self.embedding = None

        self._calculate_edge_mesh()

    def _calculate_edge_mesh(self):
        self.me_vertices = []
        self.me_faces = []

        self.edge_to_me_vertex = {}

        self.me_neighbors = []

        for face in self.faces:
            me_face = []
            for edge in self.get_face_edges(face):
                if edge not in self.edge_to_me_vertex:
                    self.edge_to_me_vertex[edge] = len(self.edge_to_me_vertex)
                    self.me_vertices.append(self.get_mid_edge_point(edge))
                    self.me_neighbors.append([])
                me_face.append(self.edge_to_me_vertex[edge])

            self.me_faces.append(me_face)
            for u, v in zip(me_face, me_face[1:] + [me_face[0]]):
                self.me_neighbors[u].append((v, 1))
                self.me_neighbors[v].append((u, -1))

        self.me_vertex_to_edge = {}
        for key, value in self.edge_to_me_vertex.items():
            self.me_vertex_to_edge[value] = key


    def get_mid_edge_point(self, edge):
        u, v = edge
        p = self.vertices[u]
        q = self.vertices[v]
        x = (p.x + q.x) / 2
        y = (p.y + q.y) / 2
        z = (p.z + q.z) / 2
        return Point(x, y, z)
        
    def calculate_planar_embedding(self):
        self.calculate_u()
        self.calculate_ul()

        self.embedding = [None for _ in self.me_vertices]
        for r in range(len(self.me_vertices)):
            i, j = self.me_vertex_to_edge[r]
            x = (self.u[i] + self.u[j]) / 2
            y = self.ul[r]
            self.embedding[r] = (x, y)

        self.display_embedding()

    def display_embedding(self):
        LENGTH_THRESHOLD = 30000
        fig, ax = plt.subplots()
        for me_face in self.me_faces[:]:
            if self.is_cut_me_face(me_face) or self.embedded_triangle_length(me_face) > LENGTH_THRESHOLD:
                continue
            x, y = zip(*(self.embedding[vr] for vr in me_face+[me_face[0]]))
            ax.plot(x, y, '-')
        plt.show()

    def embedded_triangle_length(self, me_face):
        triangle = list(self.embedding[vr] for vr in me_face)
        length = 0
        for (px, py), (qx, qy) in zip(triangle, triangle[1:] + triangle[:1]):
            x = px - qx
            y = py - qy
            length += np.sqrt(x**2 + y**2)
        return length


    def is_cut_me_face(self, me_face):
        return sorted(me_face) == sorted(self.edge_to_me_vertex[edge] for edge in self.get_face_edges(self.cut_face))

    def calculate_u(self):
        """ Calculate u (harmonic function on vertices) by solving
            sparse linear system. """
        n = len(self.vertices)
        b = np.zeros(n-2)
        A_dict = defaultdict(float)


        # TODO: implement better cut_face choice algorithm.
        self.cut_face = self.faces[3000]
        x, y, z = self.cut_face

        # We will set u[cm1] = -1 and u[c1] = 1.
        # Without loss of generality, we select cm1 as smaller than c1.
        if x < y:
            cm1, c1 = x, y
        elif y < z:
            cm1, c1 = y, z
        else:
            cm1, c1 = z, x


        for face in self.faces:
            for (u, v), angle in zip(self.get_face_edges(face), self.get_face_angles(face)):
                cot = 1 / np.tan(angle)
                A_dict[(u, u)] += cot
                A_dict[(v, u)] -= cot
                A_dict[(u, v)] -= cot
                A_dict[(v, v)] += cot

        rows = []
        cols = []
        data = []
        # We remove rows cm1 and c1, and we
        # move columns cm1 and c1 to b[].
        for (i, j), val in A_dict.items():
            if i in (cm1, c1):
                continue
            i -= (i >= cm1) + (i >= c1)
            if j == cm1:
                b[i] += val
            elif j == c1:
                b[i] -= val
            else:
                j -= (j >= cm1) + (j >= c1)
                rows.append(i)
                cols.append(j)
                data.append(val)

        A = csc_matrix((data, (rows, cols)), shape=(n-2, n-2))
        u = spsolve(A, b)

        # Insert values of cm1 and c1 into to the solution vector.
        self.u = np.zeros(n)
        self.u[:cm1] = u[:cm1]
        self.u[cm1] = -1
        self.u[cm1+1:c1] = u[cm1:c1-1]
        self.u[c1] = 1
        self.u[c1+1:] = u[c1-1:]

    def calculate_ul(self):
        """ Calculate ul (harmonic function on mid-edge vertices) 
            in a graph traversal. """
        source = 0
        self.ul = np.array([None for _ in self.me_vertices])
        self.ul[source] = 0

        visited = set()
        visited.add(source)

        frontier = deque()
        frontier.appendleft(source)

        while frontier:
            node = frontier.pop()
            for neighbor, orientation in self.me_neighbors[node]:
                if neighbor not in visited:
                    self.calculate_one_ul(node, neighbor, orientation)
                    visited.add(neighbor)
                    frontier.appendleft(neighbor)


    def calculate_one_ul(self, vs, vr, orientation):
        er = self.me_vertex_to_edge[vr]
        es = self.me_vertex_to_edge[vs]
        vj = next(iter((set(er) & set(es))))
        vi = next(iter((set(er) - set(es))))
        vk = next(iter((set(es) - set(er))))
        ang_k = self.get_angle(vk, vi, vj)
        ang_i = self.get_angle(vi, vj, vk)
        ui = self.u[vi]
        uj = self.u[vj]
        uk = self.u[vk]

        diff = (((ui - uj) / np.tan(ang_k))  +  ((uk - uj) / np.tan(ang_i))) / 2
        self.ul[vr] = self.ul[vs] + orientation * diff


    def get_face_edges(self, face): 
        """ Get edges in order. Note that the convention used for an edge here
            is that the first element is always smaller than the second one. """
        edges = []
        for u, v in zip(face, face[1:] + face[:1]):
            edge = (u, v) if u < v else (v, u)
            edges.append(edge)
        return edges

    def get_face_angles(self, face):
        """ Get opposite angles to edges, in order. """
        angles = []
        u, v, w = face
        angles.append(self.get_angle(w, u, v))
        angles.append(self.get_angle(u, w, v))
        angles.append(self.get_angle(v, w, u))
        return angles

    def get_angle(self, u, v, w):
        """ Returns angle between vectors uv and uw. """
        norm = lambda p : np.sqrt(sum(x**2 for x in p))
        p = self.vertices[u]
        q = self.vertices[v]
        r = self.vertices[w]
        vq = Point(q.x - p.x, q.y - p.y, q.z - p.z)
        vr = Point(r.x - p.x, r.y - p.y, r.z - p.z)
        return np.arccos(np.dot(vq, vr) / (norm(vq) * norm(vr)))


def read_mesh(filename):
    vertices = []
    faces = []
    with open(filename) as f:
        for line in f.readlines():
            if not line.split():
                continue
            tag, *nums = line.split() 
            if tag == 'v':
                vertices.append(Point(*map(float, nums)))
            elif tag == 'f':
                faces.append(tuple(map(lambda x : int(x)-1, nums)))
    return Mesh(vertices, faces)



def main():
    mesh = read_mesh('./datasets/non-rigid-world/cat0.obj')
    # mesh = read_mesh('./datasets/simple/tetra.obj')
    # mesh = read_mesh('./datasets/simple/reg_tetra.obj')
    mesh.calculate_planar_embedding()

if __name__ == '__main__':
    main()
