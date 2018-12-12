#include <algorithm>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <unordered_set>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
using namespace std;

typedef vector<complex<float>> vcf;

const float EPS = 1e-8;

float rand_float() {
    return rand() / (RAND_MAX + EPS);
}

int rand_int(int a, int b) {
    float r = rand_float();
    return a + r * (b - a);
}

vector<int> rand_triplet(int a, int b) {
    unordered_set<int> elems;
    vector<int> triplet;
    while(triplet.size() < 3) {
        int r = rand_int(a, b);
        if(elems.count(r) == 0) {
            elems.insert(r);
            triplet.push_back(r);
        }
    }
    return triplet;
}

complex<float> t0(0.1666, -0.28867);
complex<float> t1(-0.3333, 0);
complex<float> t2(0.1666, +0.28867);
complex<float> t3(-0.3333, 0);

void transform(const vector<int>& x_triplet, const vcf& X0, vcf& X) {
    const auto& x1 = X0[x_triplet[0]];
    const auto& x2 = X0[x_triplet[1]];
    const auto& x3 = X0[x_triplet[2]];


    auto r0 = x2 - x3;
    auto r1 = x1*x3 - x1*x2;
    auto r2 = x2 - x1;
    auto r3 = x1*x3 - x3*x2;

    auto a = t0*r0 + t1*r2;
    auto b = t0*r1 + t1*r3;

    auto c = t2*r0 + t3*r2;
    auto d = t2*r1 + t3*r3;

    for(int i = 0; i != (int)X0.size(); ++i) {
        X[i] = (a * X0[i] + b) / (c * X0[i] + d);
    }

}


vector<vector<float>> mobius_voting(const vcf& Z0, const vcf& W0, const float alpha_iter=10, float match_threshold=0.4) {
    
    int N = Z0.size();
    int num_its = alpha_iter*N*N*N;

    vector<int> z_triplet(3);
    vector<int> w_triplet(3);

    vector<int> matched_z;
    vector<int> matched_w;

    vector<int> z_closest(N);
    vector<int> w_closest(N);

    vector<vector<float>> dists(N, vector<float>(N));
    vector<vector<float>> votes(N, vector<float>(N));

    vcf Z(N);
    vcf W(N);
    
    for(int t = 0; t != num_its; ++t) {
        z_triplet = rand_triplet(0, N);
        w_triplet = rand_triplet(0, N);

        transform(z_triplet, Z0, Z);
        transform(w_triplet, W0, W);

        for(int i = 0; i != N; ++i) {
            for(int j = 0; j != N; ++j) {
                complex<float> diff = Z[i] - W[j];
                float re = diff.real();
                float im = diff.imag();
                dists[i][j] = re*re + im*im;
            }
        }

        for(int i = 0; i != N; ++i) {
            int closest = -1;
            float min_dist = 1e20;
            for(int j = 0; j != N; ++j) {
                if(min_dist > dists[i][j]) {
                    min_dist = dists[i][j];
                    closest = j;
                }
            }
            z_closest[i] = closest;
        }

        for(int j = 0; j != N; ++j) {
            int closest = -1;
            float min_dist = 1e20;
            for(int i = 0; i != N; ++i) {
                if(min_dist > dists[i][j]) {
                    min_dist = dists[i][j];
                    closest = i;
                }
            }
            w_closest[j] = closest;
        }

        matched_z.clear();
        matched_w.clear();

        for(int i = 0; i != N; ++i) {
            if(w_closest[z_closest[i]] == i) {
                matched_z.push_back(i);
                matched_w.push_back(z_closest[i]);
            }
        }
        int n = matched_z.size();

        if(n < match_threshold * N)
            continue;

        float energy = 0;
        for(int i = 0; i != n; ++i) {
            int z = matched_z[i];
            int w = matched_w[i];
            energy += abs(Z[z] - W[w]);
        }

        for(int i = 0; i != n; ++i) {
            int z = matched_z[i];
            int w = matched_w[i];
            votes[z][w] += 1.0 / (EPS + energy/n);
        }
    }

    return votes;
}

PYBIND11_MODULE(voting, m) {
    m.doc() = "Mobius voting";
    m.def("mobius_voting", &mobius_voting, "Find pointset correspondences using Mobius transformations");
}
