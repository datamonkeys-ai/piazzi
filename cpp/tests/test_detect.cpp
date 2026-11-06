// Tests of astro::detect_sources on the 3 ZTF frames of the bundled example (cpp/gui/assets/example):
// plausible detection count, bright stars of the example Gaia DR3 catalog (gaia_dr3.tsv)
// recovered within < 1 px after reduction, determinism and consistency of measure_at.
// Uso: test_detect [raiz-do-repo]
#include "piazzi/core.hpp"
#include "piazzi/reduce.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace astro;

static int fails = 0;
static void check(bool ok, const char *what) {
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what);
    if (!ok)
        ++fails;
}

// RA/Dec → pixel via the plate (linear part + Newton iterations with the full polynomial)
static void sky2pix(const Plate &p, double ra, double dec, double &x, double &y) {
    double xi, eta;
    radec_to_std(ra, dec, p.ra0, p.dec0, xi, eta);
    double a = p.cx[1], b = p.cx[2], c = p.cy[1], d = p.cy[2], det = a * d - b * c;
    x = (d * (xi - p.cx[0]) - b * (eta - p.cy[0])) / det;
    y = (-c * (xi - p.cx[0]) + a * (eta - p.cy[0])) / det;
    for (int it = 0; it < 10; ++it) {
        auto t = p.terms(x, y);
        double fx = 0, fy = 0;
        for (size_t i = 0; i < t.size(); ++i) {
            fx += t[i] * p.cx[i];
            fy += t[i] * p.cy[i];
        }
        x += (d * (xi - fx) - b * (eta - fy)) / det;
        y += (-c * (xi - fx) + a * (eta - fy)) / det;
    }
}

int main(int argc, char **argv) {
    std::string root = argc > 1 ? std::string(argv[1]) + "/" : "";
    Settings s = Settings::defaults();
    const std::string ex = root + "cpp/gui/assets/example/";
    s.load(ex + "example.cfg");
    std::vector<Image> ims;
    char b[200];
    for (int k = 0; k < 3; ++k) {
        Image im =
            load_image(ex + "ztf_hindemith_" + std::to_string(k) + ".fits", s.num("flip_horizontal") != 0, s.num("flip_vertical") != 0);
        auto t0 = std::chrono::steady_clock::now();
        std::vector<Detection> D = detect_sources(im, detect_params(s));
        double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        std::snprintf(b, sizeof b, "image %d: %zu detections in %.1f ms (expected 60–400)", k + 1, D.size(), ms);
        check(D.size() >= 60 && D.size() <= 400, b);
        // determinism: second run identical, bit for bit
        std::vector<Detection> E = detect_sources(im, detect_params(s));
        bool same = D.size() == E.size();
        for (size_t i = 0; same && i < D.size(); ++i)
            same = D[i].x == E[i].x && D[i].y == E[i].y && D[i].flux == E[i].flux && D[i].fwhm == E[i].fwhm && D[i].snr == E[i].snr;
        check(same, "determinism");
        // sane fields and flux ordering
        bool sane = true;
        for (size_t i = 0; i < D.size(); ++i) {
            const Detection &d = D[i];
            sane = sane && d.x >= 0 && d.y >= 0 && d.x < im.w && d.y < im.h && d.flux > 0 && d.fwhm > 0.5 && d.snr_tot > 0 &&
                   (i == 0 || D[i - 1].flux >= d.flux);
        }
        check(sane, "fields and ordering");
        // measure_at on the pixel of a bright star reproduces the centroid
        if (!D.empty()) {
            Detection m;
            bool ok = measure_at(im, D[0].x + 0.7, D[0].y - 0.6, detect_params(s), m);
            std::snprintf(b, sizeof b, "measure_at on the brightest star: d = %.3f px", ok ? std::hypot(m.x - D[0].x, m.y - D[0].y) : -1.0);
            check(ok && std::hypot(m.x - D[0].x, m.y - D[0].y) < 0.05, b);
        }
        ims.push_back(im);
    }
    // bright catalog stars projected through the plate → detection within < 1 px
    double t = ims[0].jd_mid;
    Catalog cat = load_vizier_tsv(ex + "gaia_dr3.tsv", 2000.0 + (t - 2451545.0) / 365.25, s.num("catalog_faint_mag"));
    check(cat.size() > 500, "example catalog");
    Reduction R = reduce(ims, cat, 107.931392, 22.404939, s, false);
    for (size_t k = 0; k < R.frames.size(); ++k) {
        const Frame &f = R.frames[k];
        if (!f.solved) {
            check(false, "plate solved");
            continue;
        }
        int n = 0, found = 0;
        double worst = 0;
        for (size_t i = 0; i < cat.size(); ++i) {
            if (cat.mag[i] > 16.0 || cat.mag[i] < s.num("catalog_bright_mag"))
                continue;
            double x, y;
            sky2pix(f.plate, cat.ra[i], cat.dec[i], x, y);
            if (x < 8 || y < 8 || x > f.im.w - 9 || y > f.im.h - 9)
                continue;
            if (f.im.at(int(std::lround(x)), int(std::lround(y))) >= s.num("saturation"))
                continue;
            ++n;
            double bd = 1e9;
            for (auto &d : f.det)
                bd = std::min(bd, std::hypot(d.x - x, d.y - y));
            if (bd < 1.0) {
                ++found;
                worst = std::max(worst, bd);
            }
        }
        std::snprintf(b, sizeof b, "image %zu: %d/%d stars RP ≤ 16 detected within < 1 px (worst %.2f px)", k + 1, found, n, worst);
        check(n >= 10 && found >= n - 1, b);
    }
    std::printf(fails ? "%d failures\n" : "OK\n", fails);
    return fails ? 1 : 0;
}
