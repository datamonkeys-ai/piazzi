// Validation of offline identification (api_orbits.cpp) against JPL Horizons ephemerides.
// Orbits: cpp/tests/data/orbits_mpcorb.dat and orbits_comets.dat (80+ column MPCORB / COMET format), with heliocentric
// osculating elements from JPL Horizons at epoch 2020-12-17.0 TDB (see cpp/tests/data/README.md).
// Reference positions: https://ssd.jpl.nasa.gov/api/horizons.api (QUANTITIES=1,9: ICRF astrometric RA/Dec and
// APmag; CENTER I41@399 = Palomar/ZTF for "topo", 500 = geocenter). Offline propagation is two-body, so the
// tolerance grows with distance from the osculation epoch.
// Usage: test_orbits [data-dir]
#include "piazzi/astro_c.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

static int fails = 0;
static double field(const std::string &js, const std::string &packed, const char *key) {
    size_t p = js.find("\"packed\":\"" + packed + "\"");
    if (p == std::string::npos)
        return NAN;
    size_t e = js.find('}', p), k = js.find(std::string("\"") + key + "\":", p);
    if (k == std::string::npos || k > e)
        return NAN;
    return std::atof(js.c_str() + k + std::strlen(key) + 3);
}
static void check(bool ok, const char *what) {
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what);
    if (!ok)
        ++fails;
}
// Longitude no .cfg: positiva = Oeste. Palomar (I41): 243,1402° E -> 116,8598 Oeste
static void site(AstroSession *s, bool topo) {
    astro_set(s, "longitude", topo ? "116.8598" : "0");
    astro_set(s, "latitude", topo ? "33.35728" : "0");
    astro_set(s, "altitude_m", topo ? "1665.5" : "0");
}

int main(int argc, char **argv) {
    std::string dir = argc > 1 ? argv[1] : "cpp/tests/data/";
    AstroSession *s = astro_new();
    int n = astro_load_orbits(s, (dir + "orbits_mpcorb.dat").c_str(), (dir + "orbits_comets.dat").c_str());
    std::printf("orbits: %d\n", n);
    check(n == 6, "read MPCORB (4) + COMET (2)");

    // (a) (5157) Hindemith in the ZTF example (mid 1st exposure, Palomar)
    site(s, true);
    double mra = 107.934503553, mde = 22.404402578, jd0 = 2459251.7159606;
    std::string js = astro_identify_offline(s, jd0, mra, mde, 0.5);
    double sep = field(js, "05157", "sep_arcsec"), rate = field(js, "05157", "rate_arcsec_min"), pa = field(js, "05157", "pa_deg");
    char b[240];
    std::snprintf(b, sizeof b, "Hindemith: sep %.2f\" rate %.3f\"/min PA %.2f (Horizons 0.353\"/min, 276.7)", sep, rate, pa);
    check(sep < 1.0 && std::fabs(rate - 0.353) < 0.01 && std::fabs(pa - 276.7) < 1, b);

    // (b) JPL Horizons references
    struct R {
        const char *packed;
        double jd;
        bool topo;
        double ra, dec, V; // V < 0: do not compare (comets: Horizons uses a different magnitude law)
        double tol_as, tol_v;
    } ref[] = {{"05157", jd0, true, mra, mde, 17.654, 1.0, 0.1},
               {"00001", 2459251.7, false, 357.846711658, -9.861057672, 9.250, 2, 0.1},
               {"00433", 2459400.5, true, 258.298716652, -36.043082231, 11.474, 30, 0.1},
               {"00003", 2462000.5, false, 164.505715401, 6.642763109, 10.645, 2 * 3600.0, 0.2}, // sanity, 7.7 yr
               {"0002P", 2459030.5, false, 119.116196091, 18.550758700, -1, 200, 0},
               {"0009P", 2459640.5, true, 281.162375677, -22.640950666, -1, 300, 0}};
    for (auto &r : ref) {
        site(s, r.topo);
        js = astro_identify_offline(s, r.jd, r.ra, r.dec, std::max(0.1, 2 * r.tol_as / 3600));
        double sp = field(js, r.packed, "sep_arcsec"), V = field(js, r.packed, "vmag");
        std::snprintf(b, sizeof b, "%s JD %.5f vs Horizons: %.2f\" (tol %.0f\")  dV %.3f", r.packed, r.jd, sp, r.tol_as, V - r.V);
        check(sp < r.tol_as && (r.V < 0 || std::fabs(V - r.V) < r.tol_v), b);
    }

    // performance: time per orbit, extrapolated to ~1.4 million orbits
    site(s, true);
    auto t0 = std::chrono::steady_clock::now();
    int reps = 20000;
    for (int i = 0; i < reps; ++i)
        astro_identify_offline(s, jd0 + i * 0.01, mra, mde, 1.0);
    double us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count() / (reps * double(n));
    std::printf("time: %.3f us/orbit -> %.2f s for 1.4M\n", us, us * 1.4e6 / 1e6);
    check(us * 1.4e6 / 1e6 < 3, "performance");
    astro_free(s);
    std::printf(fails ? "%d failures\n" : "all ok\n", fails);
    return fails != 0;
}
