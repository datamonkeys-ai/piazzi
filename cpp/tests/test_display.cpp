// Offline tests: alternative catalogs (VizieR-format samples), URLs, levels/histogram and FITS with WCS.
// Usage: test_display <repository root> [out.fits ...]  (writes FITS of the 3 ZTF example frames and a .txt with
// the plate pix_to_radec at corners/center, for external WCS verification)
#include "piazzi/astro_c.h"
#include "piazzi/reduce.hpp"
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>
using namespace astro;

static int fails = 0;
#define CHECK(c)                                                                                                                           \
    do {                                                                                                                                   \
        if (!(c)) {                                                                                                                        \
            std::printf("FAILED %s:%d: %s\n", __FILE__, __LINE__, #c);                                                                     \
            ++fails;                                                                                                                       \
        }                                                                                                                                  \
    } while (0)

int main(int argc, char **argv) {
    std::string root = argc > 1 ? argv[1] : ".", d = root + "/cpp/tests/data/";
    double ep = 2000.04;
    // UCAC4: rmag, pm epoch 2000, line without mag discarded, sorted by mag
    Catalog c = load_vizier_tsv(d + "cat_ucac4.tsv", ep, 18);
    CHECK(c.size() == 2 && c.mag[0] == 13.1 && c.mag[1] == 15.2);
    CHECK(std::fabs((c.dec[1] - 29.54) * 3.6e6 - (-20.0) * 0.04) < 1e-3);
    c = load_vizier_tsv(d + "cat_cmc15.tsv", ep, 18);
    CHECK(c.size() == 2 && c.mag[0] == 12.5 && c.ra[0] == 125.57);
    c = load_vizier_tsv(d + "cat_ppmxl.tsv", ep, 18);
    CHECK(c.size() == 2 && c.mag[0] == 11.0 && std::fabs((c.dec[1] - 29.54) * 3.6e6 - 8.0 * 0.04) < 1e-3);
    for (const char *f : {"cat_atlas.tsv", "cat_atlas.txt"}) { // ATLAS: RA_ICRS/DE_ICRS, pm since 2015.5
        c = load_vizier_tsv(d + f, ep, 18);
        CHECK(c.size() == 2 && c.mag[0] == 13.1);
        CHECK(std::fabs((c.dec[1] - 29.54) * 3.6e6 - (-20.0) * (ep - 2015.5)) < 1e-3);
    }
    c = load_vizier_tsv(d + "cat_gaia.tsv", ep, 18); // Gaia: RPmag, cut at maglim
    CHECK(c.size() == 2 && c.mag[0] == 13.0);
    const std::string ex = root + "/cpp/gui/assets/example/";
    c = load_vizier_tsv(ex + "gaia_dr3.tsv", ep, 18);
    CHECK(c.size() > 500);

    AstroSession *s = astro_new();
    astro_load_cfg(s, (ex + "example.cfg").c_str());
    std::string u = astro_catalog_url_for(s, 7, "08 22 17.0", "+29 32 31");
    CHECK(u.find("I/355/gaiadr3") != std::string::npos && u.find("c.r=30.0") != std::string::npos);
    u = astro_catalog_url_for(s, 4, "08 22 17.0", "+29 32 31");
    CHECK(u.find("I/322A/out") != std::string::npos && u.find("rmag=%3C") != std::string::npos);
    u = astro_catalog_url_for(s, 5, "08 22 17.0", "+29 32 31");
    CHECK(u.find("I/327/cmc15") != std::string::npos);
    u = astro_catalog_url_for(s, 6, "08 22 17.0", "+29 32 31");
    CHECK(u.find("I/317/sample") != std::string::npos);
    u = astro_catalog_url_for(s, 8, "08 22 17.0", "+29 32 31");
    CHECK(u.find("refcat2") != std::string::npos && u.find("RA_ICRS") != std::string::npos);
    for (int k = 0; k < 3; ++k)
        astro_add_image(s, (ex + "ztf_hindemith_" + std::to_string(k) + ".fits").c_str());
    u = astro_catalog_url_for(s, 7, "07 11 43.534", "+22 24 17.78");
    CHECK(u.find("c.r=30.0") == std::string::npos);
    std::printf("%s\n", u.c_str());
    if (argc > 2)
        astro_set(s, "plate_order", argv[2]);
    int nm = astro_reduce(s, (ex + "gaia_dr3.tsv").c_str(), "07 11 43.534", "+22 24 17.78", "5157");
    CHECK(nm >= 1);
    int w = 0, h = 0;
    CHECK(astro_frame_size(s, 0, &w, &h));
    std::vector<unsigned char> rgba(size_t(w) * h * 4);
    astro_render_levels(s, 0, 50, 99.5, 0, 0, rgba.data());
    int sat = 0, blk = 0;
    for (size_t i = 0; i < rgba.size(); i += 4) {
        sat += rgba[i] == 255;
        blk += rgba[i] == 0;
    }
    CHECK(sat > 0 && blk >= int(w * h * 0.4));
    astro_render_levels(s, 0, 50, 99.5, 2, 1, rgba.data());
    CHECK(rgba[3] == 255);
    std::vector<double> hist(64);
    double mm[2];
    CHECK(astro_histogram(s, 0, 64, hist.data(), mm) == 1 && mm[1] > mm[0]);
    double tot = 0;
    for (double v : hist)
        tot += v;
    CHECK(tot > 0.99 * w * h);
    const std::string tmp = (std::filesystem::temp_directory_path() / "astro_display_").string();
    FILE *t = std::fopen((tmp + "wcs_check.txt").c_str(), "w");
    CHECK(t != nullptr);
    if (!t)
        return 1;
    for (int k = 0; k < astro_frame_count(s); ++k) {
        std::string p = tmp + "frame" + std::to_string(k) + ".fits";
        CHECK(astro_save_fits(s, k, p.c_str()) == 1);
        double pts[5][2] = {{0, 0}, {w - 1.0, 0}, {0, h - 1.0}, {w - 1.0, h - 1.0}, {w / 2.0, h / 2.0}};
        for (auto &q : pts) {
            double ra, de;
            astro_pix_to_sky(s, k, q[0], q[1], &ra, &de);
            std::fprintf(t, "%s %.6f %.6f %.12f %.12f\n", p.c_str(), q[0], q[1], ra, de);
        }
    }
    std::fclose(t);
    astro_free(s);
    std::printf(fails ? "%d failures\n" : "OK\n", fails);
    return fails ? 1 : 0;
}
