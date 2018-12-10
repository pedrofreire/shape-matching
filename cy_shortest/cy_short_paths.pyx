from collections import deque

def all_pairs_shortest_paths(mesh):
    n = len(mesh.vertices)
    dists = [[float('inf') for _ in range(n)] for _ in range(n)]

    for source in range(n):
        approx_single_source_geo_dists(dists, mesh, source)
    return dists


def approx_single_source_geo_dists(dists_mat, mesh, source): 
    visited = [False for _ in mesh.vertices]
    dists = dists_mat[source]
    frontier = deque()

    visited[source] = True
    dists[source] = 0
    frontier.appendleft(source)

    cdef int node
    cdef int neighbor
    while frontier:
        node = frontier.pop()
        for neighbor in mesh.neighbors[node]:
            if not visited[neighbor]:
                dists[neighbor] = dists[node] + 1
                visited[neighbor] = True
                frontier.appendleft(neighbor)
