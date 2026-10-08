// Closed-form ground-state energy of N spinless fermions on an open chain of
// L sites with hopping t: the N lowest levels -2t cos(pi k / (L+1)).
#include <cmath>

double open_chain_energy(int L, double t, int N) {
    const double pi = std::acos(-1.0);
    double e = 0.0;
    for (int k = 1; k <= N; ++k) e += -2.0 * t * std::cos(pi * k / (L + 1));
    return e;
}
