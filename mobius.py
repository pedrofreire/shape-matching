from collections import defaultdict, namedtuple, deque
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve, lsqr
import heapq
import sys
import pyximport; pyximport.install()

Point = namedtuple('Point', ['x', 'y', 'z'])

class Mesh:
    def __init__(self, vertices=None, faces=None):
        if not vertices: vertices = []
        if not faces: faces = []
        self.vertices = vertices.copy()
        self.faces = faces.copy()

        self.neighbors = None
        self.u = None
        self.ul = None
        self.me_vertices = None
        self.me_faces = None
        self.me_vertex_to_edge = None
        self.me_neighbors = None
        self.edge_to_me_vertex = None
        self.cut_face = None
        self.embedding = None

        self._calculate_neighbors()
        self._calculate_edge_mesh()

    def _calculate_neighbors(self):
        self.neighbors = [[] for _ in self.vertices]
        for face in self.faces:
            for u, v in zip(face, face[1:] + face[:1]):
                self.neighbors[u].append(v)
                self.neighbors[v].append(u)

        self.vertex_faces = [[] for _ in self.vertices]
        for face_idx, face in enumerate(self.faces):
            for v in face:
                self.vertex_faces[v].append(face_idx)

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

        # self.display_embedding()

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

    def calculate_sample(self, N):
        self.calculate_geodesic_distances()
        self.sample = np.array(self.find_curvature_extrema())
        if len(self.sample) > N:
            self.sample = np.random.choice(self.sample, N, replace=False)
        self.sampleset = set(self.sample)

        while len(self.sample) < N:
            farthest_point = self.get_farthest_point()
            if not farthest_point:
                break
            self.sample.append(farthest_point)
            self.sampleset.add(farthest_point)

    def calculate_geodesic_distances(self):
        n = len(self.vertices)
        self.dists = [[float('inf') for _ in range(n)] for _ in range(n)]

        for source in range(n):
            # self.single_source_geodesic_distances(source)
            self.approx_single_source_geo_dists(source)

    def approx_single_source_geo_dists(self, source): 
        visited = [False for _ in self.vertices]
        dists = self.dists[source]
        frontier = deque()

        visited[source] = True
        dists[source] = 0
        frontier.appendleft(source)

        while frontier:
            node = frontier.pop()
            for neighbor in self.neighbors[node]:
                if not visited[neighbor]:
                    dists[neighbor] = dists[node] + 1
                    visited[neighbor] = True
                    frontier.appendleft(neighbor)

    def single_source_geodesic_distances(self, source):
        visited = set()
        dists = self.dists[source]
        frontier = []

        visited.add(source)
        dists[source] = 0
        frontier.append((dists[source], source))

        while frontier:
            dist, node = heapq.heappop(frontier)
            # Lazy removal of node.
            if dist != dists[node]:
                continue
            visited.add(node)

            for neighbor in self.neighbors[node]:
                if neighbor not in visited:
                    dist = self.get_edge_dist((node, neighbor))
                    if dists[neighbor] > dists[node] + dist:
                        dists[neighbor] = dists[node] + dist
                        heapq.heappush(frontier, (dists[neighbor], neighbor))

    def get_edge_dist(self, edge):
        p, q = (self.vertices[v] for v in edge)
        return np.sqrt(sum((px - qx)**2 for px, qx in zip(p, q)))

    def find_curvature_extrema(self):
        self.find_curvatures()
        extrema = []
        for vertex in range(len(self.vertices)):
            diff = self.curvature[vertex] - max(self.curvature[neighbor] for neighbor in self.neighbors[vertex])
            if diff > 0:
                extrema.append(vertex)
        return extrema

    def find_curvatures(self):
        self.curvature = [self.find_vertex_curvature(v) for v in range(len(self.vertices))]

    def find_vertex_curvature(self, u):
        area = 0
        sum_angles = 0
        for face_idx in self.vertex_faces[u]:
            face = self.faces[face_idx]
            v, w = set(face) - set([u])
            sum_angles += self.get_angle(u, v, w)
            area = self.face_area(face) / 3
        return (2*np.pi - sum_angles) / area

    def face_area(self, face):
        p, q, r = [self.vertices[v] for v in face]
        z1 = np.array([
         q.x - p.x,
         q.y - p.y,
         q.z - p.z,
        ])
        z2 = np.array([
         r.x - p.x,
         r.y - p.y,
         r.z - p.z,
        ])
        return np.linalg.norm(np.cross(z1, z2)) / 2

    def get_farthest_point(self):
        farthest_point = None
        max_dist = -1
        for point in range(len(self.vertices)):
            if point in self.sampleset:
                continue
            dist = min(self.dists[point][q] for q in self.sample)
            if dist > max_dist:
                max_dist = dist
                farthest_point = point
        return farthest_point 

    def project(self, points):
        return np.array([self.point_project(point) for point in points])

    def point_project(self, point):
        closest = min(self.neighbors[point], key=lambda q : self.get_edge_dist((point, q)))
        edge = tuple(sorted([point, closest]))
        me_point = self.edge_to_me_vertex[edge]
        x, y = self.embedding[me_point]
        return x + 1j*y
    

def mobius_voting(M1, M2, num_it=None):
    Z = M1.project(M1.sample)
    W = M2.project(M2.sample)
    N = len(Z)

    if num_it is None:
        num_it = N**3
    votes = np.zeros((N, N))
    for _ in range(num_it):
        Z_triplet = np.random.choice(Z, 3, replace=False)
        W_triplet = np.random.choice(W, 3, replace=False)

        Z_mobius = get_mobius_transform(Z_triplet)
        W_mobius = get_mobius_transform(W_triplet)

        Z_bar = Z_mobius(Z)
        W_bar = W_mobius(W)

        mutual_pairs = get_mutually_closest_pairs(Z_bar, W_bar)

        EPS = 1e-5
        THRESHOLD = 3 # 0.4 * len(Z) 

        if len(mutual_pairs) <= THRESHOLD:
            continue

        deformation_error = get_deformation_error(Z_bar, W_bar, mutual_pairs)
        for z, w in mutual_pairs:
            votes[z, w] += 1 / (EPS + deformation_error)

    correspondence = extract_correspondence(votes)
    return correspondence

def get_mobius_transform(Z):
    Y = np.array([np.exp(2j*k*np.pi/3) for k in range(3)])
    Ty = get_mobius_matrix(Y)
    Tz = get_mobius_matrix(Z)
    T = np.linalg.inv(Ty) @ Tz
    a, b, c, d = T.flatten()
    return np.vectorize(lambda x : (a*x + b) / (c*x + d))

def get_mobius_matrix(Z):
    z1, z2, z3 = Z
    return np.array([
        [z2 - z3, z1*z3 - z1*z2],
        [z2 - z1, z1*z3 - z3*z2],
    ])


def get_mutually_closest_pairs(Z, W):
    N = len(Z)
    # TODO: check if precomputing dists increase performance.
    dists = np.array([[abs(z - w) for w in W] for z in Z])
    
    Z_closest = np.argmin(dists, axis=1)
    W_closest = np.argmin(dists, axis=0)

    pairs = []
    for i, _ in enumerate(Z):
        j = Z_closest[i]
        if W_closest[j] == i:
            pairs.append((i, j))
    return pairs

def get_deformation_error(Z, W, correspondence):
    return sum(abs(Z[i] - W[j]) for i, j in correspondence) / len(correspondence)

def extract_correspondence(votes):
    correspondence = []
    if np.max(votes) > 0:
        votes /= np.max(votes)
    N = votes.shape[0]
    EPS = 1e-5
    THRESHOLD = 0
    for _ in range(N):
        row, col = np.unravel_index(np.argmax(votes), votes.shape)
        confidence = votes[row, col]
        if confidence < THRESHOLD:
            break
        correspondence.append((confidence, row, col))
        votes[row, :] = np.zeros(N)
        votes[:, col] = np.zeros(N)
    return correspondence



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

def print_correspondence(meshes, correspondence):
    print(*(confidence for confidence, *_ in correspondence), sep='\n', end='\n\n')
    for i, mesh in enumerate(meshes):
        for _, *v in correspondence:
            print(*mesh.vertices[v[i]])
        print('')

def main():
    file1 = 'cat0'
    file2 = 'cat1'
    meshes = []
    for filename in (file1, file2):
        if filename == file2: sys.exit(0)
        mesh = read_mesh(f'./datasets/non-rigid-world/{filename}.obj')
        mesh.calculate_planar_embedding()
        # mesh.display_embedding()
        mesh.calculate_sample(30)
        meshes.append(mesh)
        # for v in mesh.sample:
        #     print(*mesh.vertices[v])
        # print('')
    correspondence = mobius_voting(*meshes)
    print_correspondence(meshes, correspondence)

if __name__ == '__main__':
    main()
