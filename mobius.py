from collections import defaultdict, namedtuple, deque
import datetime
import heapq
import random
import re
import sys
import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve, lsqr

from shortest import voting, short_paths

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

        self.calculate_neighbors()
        self.calculate_edge_mesh()

    def calculate_neighbors(self):
        self.neighbors = [[] for _ in self.vertices]
        for face in self.faces:
            for u, v in zip(face, face[1:] + face[:1]):
                self.neighbors[u].append(v)
                self.neighbors[v].append(u)

        self.vertex_faces = [[] for _ in self.vertices]
        for face_idx, face in enumerate(self.faces):
            for v in face:
                self.vertex_faces[v].append(face_idx)

    def calculate_edge_mesh(self):
        self.me_vertices = []
        self.me_faces = []

        self.edge_to_me_vertex = {}

        self.me_neighbors = []
        self.face_to_me_face = {}

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

            self.face_to_me_face[face] = me_face

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


        # The choice of cut face should not impact the result,
        # so we just pick the first face.
        self.cut_face = self.faces[0]
        self.me_cut_face = self.face_to_me_face[self.cut_face]
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
                b[i] -= (-1) * val
            elif j == c1:
                b[i] -=  (1) * val
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
                if set([node, neighbor]).issubset(set(self.me_cut_face)):
                    continue
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
        angles.append(self.get_angle(u, v, w))
        angles.append(self.get_angle(v, w, u))
        return angles

    def get_angle(self, u, v, w):
        """ Returns angle between vectors uv and uw. """
        p = self.vertices[u]
        q = self.vertices[v]
        r = self.vertices[w]
        vq = Point(q.x - p.x, q.y - p.y, q.z - p.z)
        vr = Point(r.x - p.x, r.y - p.y, r.z - p.z)
        return np.arccos(np.dot(vq, vr) / (np.linalg.norm(vq) * np.linalg.norm(vr)))

    def calculate_sample(self, N):
        self.dists = self.fast_calculate_geodesic_distances()
        self.sample = self.find_curvature_extrema()
        if len(self.sample) > N:
            self.sample = np.random.choice(self.sample, N, replace=False)
        self.sampleset = set(self.sample)

        while len(self.sample) < N:
            farthest_point = self.get_farthest_point()
            self.sample.append(farthest_point)
            self.sampleset.add(farthest_point)
        self.sample = np.array(self.sample)

    def fast_calculate_geodesic_distances(self):
        return short_paths.w_all_pairs_shortest_paths(self.neighbors, self.vertices)

    def calculate_geodesic_distances(self):
        n = len(self.vertices)
        self.dists = [[float('inf') for _ in range(n)] for _ in range(n)]

        for source in range(n):
            self.single_source_geodesic_distances(source)

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
        good_extrema = set()
        for node in range(len(self.vertices)):
            is_boundary = any(len(set(self.neighbors[neighbor]) & set(self.neighbors[node])) < 2 for neighbor in self.neighbors[node])
            is_extrema = self.curvature[node] > max(self.curvature[neighbor] for neighbor in self.neighbors[node]) + 1

            # We remove an extremum if we find a better extremum in its region.
            is_near_better_extrema = False
            to_remove = []
            for extremum in good_extrema:
                is_close = self.dists[node][extremum] < 10
                if is_close: 
                    if self.curvature[node] > self.curvature[extremum]:
                        to_remove.append(extremum)
                    else:
                        is_near_better_extrema = True
                        break

            if not is_boundary and is_extrema and not is_near_better_extrema:
                good_extrema -= set(to_remove)
                good_extrema.add(node)
        return list(good_extrema)

    def find_curvatures(self):
        self.curvature = [self.find_vertex_curvature(v) for v in range(len(self.vertices))]

    def find_vertex_curvature(self, u):
        area = 0
        sum_angles = 0
        sqr_dist = lambda p,q : np.sqrt(sum((px - qx)**2 for px, qx in zip(p, q)))
        for face_idx in self.vertex_faces[u]:
            face = self.faces[face_idx]
            v, w = set(face) - set([u])
            sum_angles += self.get_angle(u, v, w)

            obtuse=any(angle > np.pi/2 for angle in self.get_face_angles(face))
            if not obtuse:
                a, b = self.get_angle(w, u, v), self.get_angle(v, w, u)
                area += 1/8 * (self.get_edge_dist((u, v)) / np.tan(a) + self.get_edge_dist((u, w)) / np.tan(b))
            else:                
                if(self.get_angle(u, v, w) > np.pi/2):             
                    area += self.face_area(face) / 2
                else:
                    area += self.face_area(face) / 4
                
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

    def total_area(self):
        return sum(self.face_area(face) for face in self.faces)
    
ITER_CONST = 10
MUTUAL_PAIRS_THRESHOLD = 0.4

def fast_mobius_voting(M1, M2, num_it=None):
    Z0 = M1.project(M1.sample)
    W0 = M2.project(M2.sample)

    votes = np.array(voting.mobius_voting(Z0, W0, ITER_CONST, MUTUAL_PAIRS_THRESHOLD))
    correspondence = extract_correspondence(votes)
    return correspondence

def mobius_voting(M1, M2, num_it=None):
    Z0 = M1.project(M1.sample)
    W0 = M2.project(M2.sample)
    N = len(Z0)

    EPS = 1e-8
    THRESHOLD = 0.40 * len(Z0) 

    if num_it is None:
        num_it = 10*N**3
    votes = np.zeros((N, N))
    for _ in range(num_it):
        Z_triplet = np.random.choice(Z0, 3, replace=False)
        W_triplet = np.random.choice(W0, 3, replace=False)

        Z_mobius = get_mobius_transform(Z_triplet)
        W_mobius = get_mobius_transform(W_triplet)

        Z = apply_mobius(Z0, Z_mobius)
        W = apply_mobius(W0, W_mobius)

        paired_z, paired_w = get_mutually_closest_pairs(Z, W)

        if len(paired_z) <= THRESHOLD:
            continue

        deformation_error = get_deformation_error(Z, W, paired_z, paired_w)
        votes[paired_z, paired_w] += 1 / (EPS + deformation_error) 

    correspondence = extract_correspondence(votes)
    return correspondence


def get_mobius_transform(Z):
    # Hardcode matrix to speed computation.
    # Y = np.array([np.exp(2j*k*np.pi/3) for k in range(3)])
    # Ty = get_mobius_matrix(Y)
    # yTy = np.inalginv(Ty)
    iTy = np.array([[ 0.16666667-2.88675135e-01j, -0.33333333-1.11022302e-16j],
                    [ 0.16666667+2.88675135e-01j, -0.33333333-1.12172823e-16j]])
    Tz = get_mobius_matrix(Z)
    a, b, c, d = (iTy @ Tz).flatten()
    return a, b, c, d

def get_mobius_matrix(Z):
    z1, z2, z3 = Z
    return np.array([
        [z2 - z3, z1*z3 - z1*z2],
        [z2 - z1, z1*z3 - z3*z2],
    ])

def apply_mobius(Z, Z_mobius):
    a, b, c, d = Z_mobius
    return (a*Z + b) / (c*Z + d)

def get_mutually_closest_pairs(Z, W):
    N = Z.shape[0]
    diffs = Z - W.reshape(-1, 1)
    dists_sqr = diffs.real**2 + diffs.imag**2
    
    Z_closest = np.argmin(dists_sqr, axis=0)
    W_closest = np.argmin(dists_sqr, axis=1)

    paired_z = np.nonzero(W_closest[Z_closest] == np.arange(N))[0]
    paired_w = Z_closest[paired_z]
    return paired_z, paired_w

def get_deformation_error(Z, W, paired_z, paired_w):
    return np.sum(np.abs(Z[paired_z] - W[paired_w])) / len(paired_z)

def extract_correspondence(votes):
    correspondence = []
    if np.max(votes) > 0:
        votes /= np.max(votes)
    N = votes.shape[0]
    THRESHOLD = 0 # 0.35
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
            print(*mesh.sample[v[i]])
        print('')

def output_pointset(points, filename):
    with open(filename, 'w') as f:
        for point in points:
            f.write('{:.2f} {:.2f} {:.2f}\n'.format(*point))

def evaluate(meshes, filenames, correspondence):
    Z_truth = read_groundtruth(filenames[0])
    W_truth = read_groundtruth(filenames[1])

    Z_pts = np.array([meshes[0].vertices[meshes[0].sample[i]] for _, i, _ in correspondence])
    W_pts = np.array([meshes[1].vertices[meshes[1].sample[j]] for _, _, j in correspondence])

    dists = np.array([[np.linalg.norm(pt - pt_truth) for pt in Z_pts] for pt_truth in Z_truth])
    
    Z_idxs = np.argmin(dists, axis=0)

    
    diffs = W_truth[Z_idxs] - W_pts
    dists = np.linalg.norm(diffs, axis=1)
    err = np.sum(dists) / (np.sqrt(meshes[0].total_area()) * len(correspondence))
    return err


OBJS_FOLDER = 'datasets/non-rigid-world'
GROUNDTRUTH_FOLDER = 'datasets/tosca'

def read_groundtruth(filename):
    points = []
    with open(f'{GROUNDTRUTH_FOLDER}/{filename}.obj') as f:
        for line in f.readlines():
            tag, *nums = line.split() 
            if tag == 'v':
                points.append(list(map(float, nums)))
    return np.array(points)

def symlink_force(target, link_name):
    os.symlink(target, 'tmp')
    os.rename('tmp', link_name)

def get_timestamp():
    return re.sub(r'[: \.]', '-', str(datetime.datetime.now()))

def output_error(output_folder, err, corresp_size):
    with open(f'{output_folder}/eval', 'w') as f:
        f.write(f'{err}\n')
        f.write(f'{corresp_size}\n')

def log_execution(meshes, filenames, correspondence, err):
    inter_folder = f'{filenames[0]}_{filenames[1]}_{get_timestamp()}'
    output_folder = f'outputs/{inter_folder}'
    os.makedirs(output_folder, exist_ok=True)

    for i, (mesh, filename) in enumerate(zip(meshes, filenames)):
        symlink_force(f'../../{OBJS_FOLDER}/{filename}.obj', f'{output_folder}/mesh{i}.obj')
        output_pointset((mesh.vertices[i] for i in mesh.sample), f'{output_folder}/sample{i}')
        output_pointset((mesh.vertices[mesh.sample[pair[i]]] for _, *pair in correspondence), f'{output_folder}/correspondence{i}')

    output_error(output_folder, err, len(correspondence))
    symlink_force(inter_folder, 'outputs/last')

    num_saved_runs = len(os.listdir('outputs/all_runs'))
    symlink_force(f'../{inter_folder}', f'outputs/all_runs/{num_saved_runs}')

def run_pair(filenames, evaluate=True):
    
    meshes = []
    for i, filename in enumerate(filenames):
        mesh = read_mesh(f'{OBJS_FOLDER}/{filename}.obj')
        meshes.append(mesh)

        mesh.calculate_planar_embedding()
        # mesh.display_embedding()
        mesh.calculate_sample(60)

    correspondence = fast_mobius_voting(*meshes)
    err = evaluate(meshes, filenames, correspondence)
    log_execution(meshes, filenames, correspondence, err)

def run_group(filenames):
    N = len(filenames)
    errors = [[None for _ in range(N)] for _ in range(N)]
    meshes = []
    for i, filename in enumerate(filenames):
        mesh = read_mesh(f'{OBJS_FOLDER}/{filename}.obj')
        meshes.append(mesh)
        mesh.calculate_planar_embedding()
        mesh.calculate_sample(40)

    for i, (first_mesh, first_name) in enumerate(zip(meshes, filenames)):
        for j, (second_mesh, second_name) in enumerate(zip(meshes, filenames)):
            if j < i:
                continue

            correspondence = fast_mobius_voting(first_mesh, second_mesh)
            # err = evaluate([first_mesh, second_mesh], [first_name, second_name], correspondence)
            err = 100
            errors[i][j] = err
            errors[j][i] = err

            log_execution([first_mesh, second_mesh], [first_name, second_name], correspondence, err)

    # print(*errors, sep='\n')


def run_experiments():
    filenames = []
    with open('datasets/good_models.txt') as f:
        for line in f.readlines():
            filenames.append(line.strip())

    filenames = filenames[:3]
    meshes = []
    for i, filename in enumerate(filenames):
        mesh = read_mesh(f'{OBJS_FOLDER}/{filename}.obj')
        meshes.append(mesh)
        mesh.calculate_planar_embedding()
        mesh.calculate_sample(40)

    full_samples = [mesh.sample for mesh in meshes]
    num_runs = 5
    for _ in range(num_runs):
        i = random.randrange(len(filenames))
        j = random.randrange(len(filenames))
        N = random.randrange(30, 40+1)
        first_mesh = meshes[i]
        first_name = filenames[i]
        second_mesh = meshes[j]
        second_name = filenames[j]

        first_mesh.sample = full_samples[i][:N]
        second_mesh.sample = full_samples[j][:N]

        correspondence = fast_mobius_voting(first_mesh, second_mesh)
        if first_name[:3] == second_name[:3]:
            err = evaluate([first_mesh, second_mesh], [first_name, second_name], correspondence)
        else:
            err = 100000.0
        log_execution([first_mesh, second_mesh], [first_name, second_name], correspondence, err)

def main():
    # files = ['cat0', 'cat1', 'cat10', 'cat2'] #, 'cat3', 'cat4', 'cat5', 'cat7', 'cat8']

    run_experiments()
    return

    filenames = [
        'victoria10',
        'victoria12',
        'victoria17',
        'victoria2',
    ]

    filenames = [
        'dog0',
        'dog1',
        'horse0',
        'wolf0',
    ]

    filenames = [
        'david0',
        'david1',
        'david10',
        'michael1',
        'michael10',
        'michael11',
            ]

    run_group(filenames)
    return
    filenames = ('horse0', 'horse14')
    run_pair(filenames)


if __name__ == '__main__':
    main()
