#include <cmath>
#include <set>
#include <queue>
#include <unordered_set>
#include <utility>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
using namespace std;

typedef pair<float, int> fi;

namespace py = pybind11;

const float INF = 1e9;
const float EPS = 1e-9;

void unweighted_single_shortest_paths(vector<vector<float>>& dists_mat, const vector<vector<int>>& adj_lists, int source) {
    set<int> visited;
    queue<int> frontier;

    auto& dists = dists_mat[source];
    dists[source] = 0;
    visited.insert(source);
    frontier.push(source);

    while(!frontier.empty()) {
        int node = frontier.front();
        frontier.pop();
        for(int neighbor : adj_lists[node]) {
            if(visited.count(neighbor) == 0) {
                dists[neighbor] = dists[node] + 1;
                visited.insert(neighbor);
                frontier.push(neighbor);
            }
        }
    }
}

vector<vector<float>> all_pairs_shortest_paths(const vector<vector<int>>& adj_lists) {
    int n = adj_lists.size();
    vector<vector<float>> dists(n, vector<float>(n, INF));
    for(int i = 0; i != n; ++i)
        unweighted_single_shortest_paths(dists, adj_lists, i);
    return dists;
}

float calculate_edge_dist(const vector<vector<float>>& coords, int u, int v) {
    auto& p = coords[u];
    auto& q = coords[v];
    float x = p[0] - q[0];
    float y = p[1] - q[1];
    float z = p[2] - q[2];
    return sqrt(x*x + y*y + z*z);
}

void weighted_single_shortest_paths(vector<vector<float>>& dists_mat, vector<vector<float>>& edge_dists, 
                                    const vector<vector<float>>& coords, const vector<vector<int>>& adj_lists, int source) {
    unordered_set<int> visited;
    priority_queue<fi, vector<fi>, greater<fi>> frontier;

    auto& dists = dists_mat[source];
    dists[source] = 0;
    frontier.push({dists[source], source});

    while(!frontier.empty()) {
        auto top = frontier.top();
        frontier.pop();
        float dist = top.first;
        int node = top.second;

        // Lazy node removal.
        if(dist != dists[node])
            continue;
        if(visited.count(node) == 1)
            continue;

        visited.insert(source);

        for(int neighbor : adj_lists[node]) {
            if(visited.count(neighbor) == 0) {
                if(edge_dists[node][neighbor] == INF)
                    edge_dists[node][neighbor] = calculate_edge_dist(coords, node, neighbor);
                int weight = edge_dists[node][neighbor];
                if(dists[neighbor] > dists[node] + weight) {
                    dists[neighbor] = dists[node] + weight;
                    frontier.push({dists[neighbor], neighbor});
                }
            }
        }
    }
}


vector<vector<float>> w_all_pairs_shortest_paths(const vector<vector<int>>& adj_lists, const vector<vector<float>>& coords) {
    int n = adj_lists.size();
    vector<vector<float>> edge_dists(n, vector<float>(n, INF));
    vector<vector<float>> dists(n, vector<float>(n, INF));
    for(int i = 0; i != n; ++i)
        weighted_single_shortest_paths(dists, edge_dists, coords, adj_lists, i);
    return dists;
}

PYBIND11_MODULE(short_paths, m) {
    m.doc() = "Shortest paths functions";
    m.def("all_pairs_shortest_paths", &all_pairs_shortest_paths, "Compute all pairs shortest paths");
    m.def("w_all_pairs_shortest_paths", &w_all_pairs_shortest_paths, "Compute all pairs shortest paths");
}
