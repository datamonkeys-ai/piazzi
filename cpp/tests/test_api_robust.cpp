// Robustness tests of the C API and the reports: calls with s == NULL, huge/odd designations (no exceptions),
// line breaks in the object name and JSON escaping.
#include "../src/session_internal.hpp"
#include "piazzi/astro_c.h"
#include "piazzi/reduce.hpp"
#include "piazzi/version.h"
#include <cstdio>
#include <cstring>
#include <set>
#include <string>
using namespace astro;

static int fails = 0;
#define CHECK(c)                                                                                                                           \
    do {                                                                                                                                   \
        if (!(c)) {                                                                                                                        \
            std::printf("FAILED %s:%d: %s\n", __FILE__, __LINE__, #c);                                                                     \
            ++fails;                                                                                                                       \
        }                                                                                                                                  \
    } while (0)

static int count_lines(const std::string &x) {
    int n = 0;
    for (char c : x)
        n += c == '\n';
    return n;
}

// minimal JSON validator (RFC 8259): rejects nan/inf and any invalid syntax
struct JsonCheck {
    const char *p;
    void ws() {
        while (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t')
            ++p;
    }
    bool str() {
        if (*p != '"')
            return false;
        for (++p; *p && *p != '"'; ++p) {
            if ((unsigned char)*p < 0x20)
                return false;
            if (*p == '\\' && !*++p)
                return false;
        }
        return *p++ == '"';
    }
    bool num() {
        const char *b = p;
        if (*p == '-')
            ++p;
        if (!(*p >= '0' && *p <= '9'))
            return false;
        while ((*p >= '0' && *p <= '9') || *p == '.' || *p == 'e' || *p == 'E' || *p == '+' || *p == '-')
            ++p;
        return p > b;
    }
    bool lit(const char *w) {
        size_t n = std::strlen(w);
        if (std::strncmp(p, w, n))
            return false;
        p += n;
        return true;
    }
    bool val() {
        ws();
        bool ok;
        if (*p == '{') {
            ++p;
            ws();
            if (*p == '}')
                return ++p, true;
            do {
                ws();
                if (!str())
                    return false;
                ws();
                if (*p++ != ':' || !val())
                    return false;
                ws();
            } while (*p == ',' && ++p);
            ok = *p++ == '}';
        } else if (*p == '[') {
            ++p;
            ws();
            if (*p == ']')
                return ++p, true;
            do {
                if (!val())
                    return false;
                ws();
            } while (*p == ',' && ++p);
            ok = *p++ == ']';
        } else if (*p == '"')
            ok = str();
        else
            ok = lit("true") || lit("false") || lit("null") || num();
        return ok;
    }
    static bool valid(const std::string &j) {
        JsonCheck c{j.c_str()};
        if (!c.val())
            return false;
        c.ws();
        return *c.p == 0;
    }
};

int main() {
    // --- s == NULL: every function returns the failure value, without crashing
    AstroSession *z = nullptr;
    int w = 0, h = 0;
    double a = 0, b = 0, out[12], counts[4], mm[2];
    unsigned char px[4];
    astro_free(z);
    astro_load_cfg(z, "x");
    astro_set(z, "a", "b");
    CHECK(std::string(astro_get(z, "a")).empty());
    CHECK(astro_add_image(z, "x") == 0);
    astro_clear_images(z);
    CHECK(std::string(astro_catalog_url(z, "0", "0")).empty());
    CHECK(astro_reduce(z, "x", "0", "0", "o") == -1);
    CHECK(std::string(astro_result_json(z)).empty());
    CHECK(std::string(astro_report(z, 0)).empty());
    CHECK(astro_frame_count(z) == 0);
    CHECK(astro_frame_size(z, 0, &w, &h) == 0);
    astro_render(z, 0, px);
    CHECK(astro_sky_to_pix(z, 0, 1, 1, &a, &b) == 0);
    CHECK(astro_frame_center(z, 0, &a, &b) == 0);
    CHECK(astro_pix_to_sky(z, 0, 1, 1, &a, &b) == 0);
    astro_set_mover_desig(z, 0, "x");
    CHECK(astro_set_dark(z, "x") == 0);
    CHECK(astro_set_flat(z, "x") == 0);
    CHECK(astro_image_count(z) == 0);
    CHECK(std::string(astro_image_header(z, 0)).empty());
    CHECK(astro_image_jd(z, 0) == 0);
    astro_set_image_jd(z, 0, 1);
    CHECK(astro_measure(z, 0, 1, 1, out) == 0);
    CHECK(astro_add_manual(z, -1, 0, 1, 1) == -1);
    astro_set_note(z, 0, 0, "I");
    astro_set_excluded(z, 0, 0, 1);
    CHECK(astro_stack(z, 1, 1, 0) == -1);
    CHECK(std::string(astro_catalog_url_for(z, 7, "0", "0")).empty());
    astro_render_levels(z, 0, 1, 99, 0, 0, px);
    CHECK(astro_histogram(z, 0, 4, counts, mm) == 0);
    CHECK(astro_save_fits(z, 0, "/nonexistent/x.fits") == 0);
    CHECK(astro_load_orbits(z, "x", "y") == 0);
    CHECK(std::string(astro_identify_offline(z, 2451545, 0, 0, 1)) == "[]");
    CHECK(std::string(astro_known_offline(z, 0, 20)) == "[]");
    CHECK(std::string(astro_image_catalog(z, 0)) == "[]");

    // --- session without a result: reports/frames fail without exceptions; an invalid setting does not crash
    AstroSession *s = astro_new();
    CHECK(std::string(astro_report(s, 0)).empty());
    CHECK(astro_frame_size(s, 0, &w, &h) == 0);
    astro_set(s, "catalog_faint_mag", "nao-numero");
    astro_add_image(s, "/nonexistent/a.fit");
    astro_catalog_url(s, "10 00 00", "+10 00 00");
    CHECK(astro_reduce(s, "/nonexistent/cat.tsv", "0", "0", "o\nbj") == -1);
    CHECK(std::string(astro_result_json(s)).find("error") != std::string::npos);
    astro_free(s);

    // --- designation packing: no std::stol on a huge number
    CHECK(packed_from("(433)") == "00433");
    CHECK(packed_from("100000") == "A0000");
    CHECK(packed_from("360017") == "a0017");
    CHECK(packed_from("619999") == "z9999");
    CHECK(packed_from("620000") == "620000");
    std::string big(400, '9');
    CHECK(packed_from(big) == big);
    CHECK(packed_from("2001 SP75") == "K01S75P");
    CHECK(packed_from("2000 ACx") == "2000 AC");
    CHECK(packed_from("2000 AC99999999999999999") == "2000 AC");

    // --- reports with a line break in the object name: one line per measurement
    Reduction R;
    R.frames.resize(1);
    Track T;
    T.desig = "Obj\nnovo\r|x";
    Measurement m{};
    m.frame = 0;
    m.jd = 2451545.0;
    m.snr = 10;
    m.mag = 15;
    T.meas.push_back(m);
    Settings st = Settings::defaults();
    std::string ades = ades_report(R, {T}, st), ph = phot_report(R, {T}, st), mpc = mpc_report(R, {T}, st);
    CHECK(ades.find("Obj nov") != std::string::npos && count_lines(ades) == count_lines(ades_report(R, {}, st)) + 1);
    CHECK(ph.find("Obj novo") != std::string::npos);
    CHECK(count_lines(mpc) == 7); // 6 header lines + 1 measurement
    Track G;
    G.desig = big;
    G.meas.push_back(m);
    CHECK(!mpc_report(R, {G}, st).empty() && !ades_report(R, {G}, st).empty());

    // --- JSON escaping: quotes, backslash and control characters
    CHECK(session_util::json_esc("a\"b\\c\nd\te\x01") == "a\\\"b\\\\c\\nd\\te\\u0001");

    // --- real reduction without an object (bundled ZTF example, cpp/gui/assets/example): empty designation;
    // ADES without provID "Piazzi", only trkSub; MPC with a distinct temporary designation per track;
    // valid JSON (no nan/inf) with the new fields
    CHECK(JsonCheck::valid("{\"a\":[1,-2.5e3,null,true,\"x\"]}") && !JsonCheck::valid("{\"a\":nan}") && !JsonCheck::valid("{\"a\":inf}"));
    CHECK(std::string(astro_version()) == PIAZZI_VERSION);
    {
        AstroSession *r = astro_new();
        astro_load_cfg(r, "cpp/gui/assets/example/example.cfg");
        astro_set(r, "observer", "Test");
        astro_set(r, "measurer", "Test");
        for (int k = 0; k < 3; ++k)
            astro_add_image(r, ("cpp/gui/assets/example/ztf_hindemith_" + std::to_string(k) + ".fits").c_str());
        int nmov = astro_reduce(r, "cpp/gui/assets/example/gaia_dr3.tsv", "07 11 43.534", "+22 24 17.78", nullptr);
        CHECK(nmov >= 1);
        // second track: manual measurement on a reference star in frames 0 and 1
        const auto &f0 = r->R.frames[0];
        CHECK(!f0.plate.ref_idx.empty());
        if (!f0.plate.ref_idx.empty()) {
            const auto &d = f0.det[f0.plate.ref_idx[0]];
            int t = astro_add_manual(r, -1, 0, d.x, d.y);
            CHECK(t >= 0);
            if (t >= 0)
                astro_add_manual(r, t, 1, d.x, d.y);
        }
        std::string js = astro_result_json(r);
        CHECK(JsonCheck::valid(js));
        CHECK(js.find("nan") == std::string::npos && js.find("inf") == std::string::npos);
        CHECK(js.find("\"trksub\":\"P") != std::string::npos && js.find("\"sig_ra\":") != std::string::npos &&
              js.find("\"sig_dec\":") != std::string::npos && js.find("\"trail_len\":") != std::string::npos &&
              js.find("\"remarks\":") != std::string::npos);
        CHECK(js.find("\"desig\":\"Piazzi\"") == std::string::npos);
        std::string ades = astro_report(r, 0);
        CHECK(!ades.empty() && ades.find("Piazzi|") == std::string::npos && ades.find("| Piazzi") == std::string::npos);
        std::string mpc = astro_report(r, 2);
        std::set<std::string> ids;
        size_t nobs = 0;
        for (size_t a = 0, e; a < mpc.size(); a = e + 1) {
            e = mpc.find('\n', a);
            if (e == std::string::npos)
                e = mpc.size();
            std::string ln = mpc.substr(a, e - a);
            if (ln.size() >= 80) { // observation lines (80 columns); the header is short
                ++nobs;
                ids.insert(ln.substr(0, 12));
            }
        }
        CHECK(nobs >= 5 && ids.size() >= 2); // 3 from the asteroid + 2 manual
        for (auto &x : ids)
            CHECK(x.find('P') != std::string::npos); // temporary designation 'P' + base 36
        astro_free(r);
    }

    if (fails)
        std::printf("%d failure(s)\n", fails);
    else
        std::printf("test_api_robust: OK\n");
    return fails ? 1 : 0;
}
