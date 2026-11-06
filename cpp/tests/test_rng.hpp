// Portable random numbers for the tests. std::mt19937 output is fixed by the C++ standard, but
// std::uniform_real_distribution / std::normal_distribution are not (libc++, libstdc++ and MSVC differ), so the
// simulated images would change between platforms. These replacements give the same values everywhere.
#pragma once
#include <cmath>
#include <random>

namespace test_rng {
// uniform in [a, b) from 53 random bits
struct Uniform {
    double a, b;
    Uniform(double lo, double hi) : a(lo), b(hi) {}
    double operator()(std::mt19937 &g) const {
        unsigned long long hi = g() >> 5, lo = g() >> 6; // 27 + 26 bits
        double u = (double(hi) * 67108864.0 + double(lo)) / 9007199254740992.0;
        return a + (b - a) * u;
    }
};
// normal (mean, sigma) by Box–Muller; keeps the second deviate for the next call
struct Normal {
    double mu, sd;
    bool have = false;
    double spare = 0;
    Normal(double m, double s) : mu(m), sd(s) {}
    double operator()(std::mt19937 &g) {
        if (have) {
            have = false;
            return mu + sd * spare;
        }
        Uniform u01(0, 1);
        double u1 = 0;
        while (u1 <= 0)
            u1 = u01(g);
        double u2 = u01(g), r = std::sqrt(-2 * std::log(u1)), t = 6.283185307179586 * u2;
        spare = r * std::sin(t);
        have = true;
        return mu + sd * r * std::cos(t);
    }
};
} // namespace test_rng
