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
        self._calculate_neighbors()
        self._calculate_edge_mesh()

    def _calculate_neighbors(self):
        self.neighbors = [[] for _ in self.vertices]
        for face in self.faces:
            for u in face:
                for v in face:
                    if u != v:
                        self.neighbors[u].append(v)

    def _calculate_edge_mesh(self):
        self.me_vertices = []
        self.me_faces = []

        self.me_vertex_idx = {}

        # Initialize bigger list, prune excess later.
        self.me_neighbors = [[] for _ in range(3*len(self.faces))]

        for face in self.faces:
            me_face = []
            for edge in self.get_face_edges(face):
                if edge not in self.me_vertex_idx:
                    self.me_vertex_idx[edge] = len(self.me_vertex_idx)
                    self.me_vertices.append(self.get_mid_edge_point(edge))
                me_face.append(self.me_vertex_idx[edge])

            self.me_faces.append(me_face)
            for u in me_face:
                for v in me_face:
                    if u != v:
                        self.me_neighbors[u].append(v)

        self.inv_me_vertex_idx = {}
        for key, value in self.me_vertex_idx.items():
            self.inv_me_vertex_idx[value] = key

        self.me_neighbors = self.me_neighbors[:len(self.me_vertices)]

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
            i, j = self.inv_me_vertex_idx[r]
            x = (self.u[i] + self.u[j]) / 2
            y = self.ul[r]
            self.embedding[r] = (x, y)

        self.display_embedding()

    def display_embedding(self):
        fig, ax = plt.subplots()
        for me_face in self.me_faces[:1000]:
            if self.is_cut_me_face(me_face):
                continue
            x, y = zip(*(self.embedding[vr] for vr in me_face+[me_face[0]]))
            ax.plot(x, y, '-')
        plt.show()

    def is_cut_me_face(self, me_face):
        return sorted(me_face) == sorted(self.me_vertex_idx[edge] for edge in self.get_face_edges(self.cut_face))

    def calculate_u(self):
        """ Calculate u (harmonic function on vertices) by solving
            sparse linear system. """
        A = defaultdict(float)
        b = defaultdict(float)

        # The embedding should be quasi-invariant mod Mobius transformation,
        # so selecting the cut-face smartly is not very important.
        self.cut_face = self.faces[0]

        # We implicitly set u[u1] = 1 and u[um1] = -1
        um1, u1, _ = sorted(self.cut_face)
        print(um1, u1)

        for face in self.faces:
            for (u, v), angle in zip(self.get_face_edges(face), self.get_face_angles(face)):
                cot = 1 / np.tan(angle)
                for othr, vtx in ((u, v), (v, u)):
                    A[(vtx, vtx)] += cot
                    A[(othr, vtx)] -= cot


        n = len(self.vertices)
        rows = []
        cols = []
        data = []
        # um1 = u1 = 100000
        for (i, j), val in A.items():
            # print((i, j), val)
            if i in (um1, u1):
                continue
            i -= (i >= um1) + (i >= u1)
            if j == um1:
                b[i] += val
            elif j == u1:
                b[i] -= val
            else:
                j -= (j >= um1) + (j >= u1)
                rows.append(i)
                cols.append(j)
                data.append(val)

        """
        Ans = np.zeros((n, n))
        for (i, j), val in A.items():
            Ans[i, j] = val
        print(Ans)
        """
        A = csc_matrix((data, (rows, cols)), shape=(n-2, n-2))

        rows = []
        cols = []
        data = []
        bl = np.zeros(n-2)
        for j, val in b.items():
            bl[j] = val
            rows.append(j)
            cols.append(0)
            data.append(val)
        b = bl

        print('A')
        print(A.toarray())
        print('b')
        print(b)
        u = spsolve(A, b)
        print(u)
        # u[um1] = -1
        # u[u1] = 1
        # print(lsqr_data)
        self.u = np.zeros(n)
        self.u[:um1] = u[:um1]
        self.u[um1] = -1
        self.u[um1+1:u1] = u[um1:u1-1]
        self.u[u1] = 1
        self.u[u1+1:] = u[u1-1:]
        # The intuition is that those two values would
        # be somewhat equal.
        print(self.u)
        print(len(self.u[self.u > 0]))
        print(len(self.u[self.u < 0]))
        print('----------')
        print(np.sum(abs(A @ u - b)))

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
            for neighbor in self.me_neighbors[node]:
                if neighbor not in visited:
                    self.calculate_one_ul(node, neighbor)
                    visited.add(neighbor)
                    frontier.appendleft(neighbor)

        print(self.ul)
        print(len(self.ul[self.ul > 0]))
        print(len(self.ul[self.ul < 0]))



    def calculate_one_ul(self, vs, vr):
        er = self.inv_me_vertex_idx[vr]
        es = self.inv_me_vertex_idx[vs]
        vj = next(iter((set(er) & set(es))))
        vi = next(iter((set(er) - set(es))))
        vk = next(iter((set(es) - set(er))))
        ang_k = self.get_angle(vk, vi, vj)
        ang_i = self.get_angle(vi, vj, vk)
        ui = self.u[vi]
        uj = self.u[vj]
        uk = self.u[vk]
        diffs = abs(ui - uj), abs(ui - uk), abs(uj - uk)
        # if max(np.log(diffs)) > -2:
            # print(*('{:.2f}'.format(diff) for diff in diffs))

        diff = (((ui - uj) / np.tan(ang_k))  +  ((uk - uj) / np.tan(ang_i))) / 2
        self.ul[vr] = self.ul[vs] + diff


    def get_face_edges(self, face): 
        """ Get ordered opposite edges. """
        u, v, w = sorted(face)
        return [(v, w), (u, w), (u, v)]

    def get_face_angles(self, face):
        """ Get ordered angles. """
        angles = []
        dface = sorted(face) + sorted(face)
        for i in range(len(face)): 
            angles.append(self.get_angle(*dface[i:i+3]))
        return angles

    def get_angle(self, u, v, w):
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
    mesh = read_mesh('./datasets/non-rigid-world/cat1.obj')
    # mesh = read_mesh('./datasets/simple/tetra.obj')
    mesh.calculate_planar_embedding()

if __name__ == '__main__':
    main()
