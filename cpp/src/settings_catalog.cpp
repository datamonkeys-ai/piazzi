#include "piazzi/reduce.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <initializer_list>
#include <numeric>
#include <sstream>

namespace astro {

Settings Settings::defaults() {
    Settings s;
    const char *d[][2] = {{"aperture_radius", "4"},
                          {"detect_sigma", "4.0"},
                          {"plate_max_rms", "0.20"},
                          {"catalog_radius", "0.75"},
                          {"background_mode", "1"},
                          {"catalog", "7"},
                          {"catalog_faint_mag", "18.0"},
                          {"catalog_bright_mag", "10.0"},
                          {"max_position_residual", "0.50"},
                          {"max_mag_residual", "0.50"},
                          {"plate_order", "0"},
                          {"match_star_count", "50"},
                          {"match_radius", "2.00"},
                          {"align_star_count", "30"},
                          {"pixel_width_mm", "0.0038"},
                          {"pixel_height_mm", "0.0038"},
                          {"saturation", "60000"},
                          {"report_band", "R"},
                          {"focal_length_mm", "1000.0"},
                          {"focal_length_tolerance", "1.0"},
                          {"position_angle", "0.0"},
                          {"position_angle_tolerance", "10.0"},
                          {"pointing_tolerance", "5.0"},
                          {"flip_horizontal", "0"},
                          {"flip_vertical", "0"},
                          {"time_reference", "0"},
                          {"min_fwhm", "0.70"},
                          {"report_mags", "1"},
                          {"extra_time_digit", "0"},
                          {"mpc_code", "XXX"},
                          {"ades_six_digit_positions", "0"},
                          {"ades_two_digit_mags", "0"},
                          {"telescope_aperture_m", "0.0"},
                          {"telescope_design", ""},
                          {"telescope_detector", "CCD"},
                          {"observer", ""},
                          {"measurer", ""},
                          {"submitter", ""},
                          {"telescope", ""}};
    for (auto &kv : d)
        s.kv[kv[0]] = kv[1];
    return s;
}

static std::string trim(const std::string &x) {
    size_t a = x.find_first_not_of(" \t\r\n"), b = x.find_last_not_of(" \t\r\n");
    return a == std::string::npos ? "" : x.substr(a, b - a + 1);
}

void Settings::load(const std::string &path) {
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '[' || line[0] == ';' || line[0] == '#')
            continue;
        auto p = line.find('=');
        if (p == std::string::npos)
            continue;
        kv[trim(line.substr(0, p))] = trim(line.substr(p + 1));
    }
}

double Settings::num(const std::string &k) const {
    auto it = kv.find(k);
    if (it == kv.end() || it->second.empty())
        return 0.0;
    try {
        return std::stod(it->second);
    } catch (...) {
        return 0.0;
    }
}
std::string Settings::str(const std::string &k) const {
    auto it = kv.find(k);
    return it == kv.end() ? "" : it->second;
}

// ---------------------------------------------------------------- catalogue (VizieR asu-tsv)
Catalog load_vizier_tsv(const std::string &path, double epoch_jyear, double maglim, double ruwe_max) {
    std::ifstream f(path);
    std::string line;
    std::vector<std::string> hdr;
    std::vector<std::vector<std::string>> rows;
    int nl = 0;
    std::vector<std::string> lines;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        lines.push_back(line);
    }
    // asu-txt (fixed columns): names in the description lines, widths in the dash line
    int dash = -1;
    for (size_t i = 0; i < lines.size(); ++i) {
        const auto &l = lines[i];
        if (l.empty() || l[0] != '-')
            continue;
        bool ok = true;
        int groups = 0;
        bool in = false;
        for (char ch : l) {
            if (ch != '-' && ch != ' ') {
                ok = false;
                break;
            }
            if (ch == '-' && !in) {
                ++groups;
                in = true;
            }
            if (ch == ' ')
                in = false;
        }
        if (ok && groups >= 2) {
            dash = int(i);
            break;
        }
    }
    if (dash >= 0) {
        for (int i = 0; i < dash; ++i) {
            const auto &l = lines[i];
            if (l.size() > 5 && l.compare(0, 4, "    ") == 0 && l[4] != ' ' && l.find(" (") != std::string::npos) {
                std::string n = l.substr(4, l.find(' ', 4) - 4);
                hdr.push_back(n);
            }
        }
        std::vector<std::pair<size_t, size_t>> span;
        const auto &d = lines[dash];
        for (size_t i = 0; i < d.size();) {
            if (d[i] == '-') {
                size_t j = i;
                while (j < d.size() && d[j] == '-')
                    ++j;
                span.push_back({i, j - i});
                i = j;
            } else
                ++i;
        }
        for (size_t i = dash + 1; i < lines.size(); ++i) {
            const auto &l = lines[i];
            if (l.empty() || l[0] == '#') {
                if (!rows.empty())
                    break;
                continue;
            }
            std::vector<std::string> cols;
            for (auto &sp : span)
                cols.push_back(sp.first < l.size() ? trim(l.substr(sp.first, sp.second)) : "");
            rows.push_back(cols);
        }
    } else
        for (auto &ln : lines) {
            line = ln;
            if (line.empty() || line[0] == '#')
                continue;
            std::vector<std::string> cols;
            std::stringstream ss(line);
            std::string c;
            while (std::getline(ss, c, '\t'))
                cols.push_back(trim(c));
            if (nl == 0)
                hdr = cols;
            else if (nl >= 3 && !(line.size() && line[0] == '-'))
                rows.push_back(cols);
            ++nl;
        }
    auto col = [&](const std::string &n) {
        for (size_t i = 0; i < hdr.size(); ++i)
            if (hdr[i] == n)
                return int(i);
        return -1;
    };
    // columns by name (Gaia DR3 unchanged); other supported catalogues:
    // UCAC4 rmag, CMC-15 r'mag, PPMXL r1mag (pm epoch 2000), ATLAS-RefCat2 rmag with RA_ICRS/DE_ICRS (pm epoch 2015.5)
    auto first = [&](std::initializer_list<const char *> ns) {
        for (auto n : ns) {
            int i = col(n);
            if (i >= 0)
                return i;
        }
        return -1;
    };
    int iR = first({"RAJ2000", "_RAJ2000", "RA_ICRS", "RAdeg"}), iD = first({"DEJ2000", "_DEJ2000", "DE_ICRS", "DEdeg"});
    int ipr = col("pmRA"), ipd = col("pmDE");
    int im = first({"RPmag", "rmag", "r'mag", "r1mag", "R1mag", "Rmag"});
    int ibp = col("BPmag"), irp = col("RPmag");
    double ep0 = (col("RAJ2000") < 0 && col("_RAJ2000") < 0 && col("RA_ICRS") >= 0) ? 2015.5 : 2000.0;
    auto val = [](const std::vector<std::string> &r, int i) -> double {
        if (i < 0 || i >= int(r.size()) || r[i].empty())
            return NAN;
        try {
            return std::stod(r[i]);
        } catch (...) {
            return NAN;
        }
    };
    int iera = first({"e_RA_ICRS", "e_RAJ2000", "ra_error"}), iede = first({"e_DE_ICRS", "e_DEJ2000", "dec_error"});
    int iepr = first({"e_pmRA", "pmra_error"}), iepd = first({"e_pmDE", "pmdec_error"}), iruwe = first({"RUWE", "ruwe"});
    struct Row {
        double ra, dec, mag, err, color;
    };
    std::vector<Row> out;
    double dt = epoch_jyear - ep0;
    for (auto &r : rows) {
        double ra = val(r, iR), de = val(r, iD), pr = val(r, ipr), pd = val(r, ipd), m = val(r, im);
        if (std::isnan(pr))
            pr = 0;
        if (std::isnan(pd))
            pd = 0;
        de = de + pd * dt / 3.6e6;
        ra = ra + pr * dt / 3.6e6 / std::cos(de * PI / 180);
        // RUWE > 1.4: unreliable 5-parameter astrometric solution (binaries, extended sources)
        double ru = val(r, iruwe);
        if (ruwe_max > 0 && std::isfinite(ru) && ru > ruwe_max)
            continue;
        // position σ (″, mean of axes) propagated from the Gaia DR3 reference epoch (2016.0) to the image
        double err = NAN;
        if (iera >= 0 && iede >= 0) {
            double ea = val(r, iera), ed = val(r, iede), ep = val(r, iepr), eq = val(r, iepd), t = epoch_jyear - 2016.0;
            if (std::isfinite(ea) && std::isfinite(ed)) {
                double va = ea * ea + (std::isfinite(ep) ? ep * ep * t * t : 0), vd = ed * ed + (std::isfinite(eq) ? eq * eq * t * t : 0);
                err = std::sqrt(0.5 * (va + vd)) / 1000.0;
            }
        }
        if (std::isfinite(ra) && std::isfinite(m) && m <= maglim)
            out.push_back({ra, de, m, err, val(r, ibp) - val(r, irp)});
    }
    std::vector<size_t> o(out.size());
    std::iota(o.begin(), o.end(), 0);
    std::stable_sort(o.begin(), o.end(), [&](size_t a, size_t b) { return out[a].mag < out[b].mag; });
    Catalog c;
    bool has_err = iera >= 0 && iede >= 0;
    for (size_t i : o) {
        c.ra.push_back(out[i].ra);
        c.dec.push_back(out[i].dec);
        c.mag.push_back(out[i].mag);
        if (has_err)
            c.err.push_back(out[i].err);
        if (ibp >= 0 && irp >= 0)
            c.color.push_back(out[i].color);
    }
    return c;
}

std::string gaia_box_url(double ra, double dec, double w_arcmin, double maglim, const std::string &host) {
    char b[512];
    std::snprintf(b, sizeof b,
                  "https://%s/viz-bin/asu-tsv?-source=I/355/gaiadr3&-c=%.6f%%20%+.6f&-out=RAJ2000,DEJ2000,pmRA,pmDE,Gmag,BPmag,RPmag,Plx,"
                  "e_RA_ICRS,e_DE_ICRS,e_pmRA,e_pmDE,RUWE&-out.max=20000&-oc.form=d&RPmag=%%3C%.1f&-c.geom=b&-c.r=%.1f&-c.u=arcmin",
                  host.c_str(), ra, dec, maglim + 0.1, w_arcmin);
    return b;
}

std::string catalog_box_url(int catalog, double ra, double dec, double w_arcmin, double maglim, const std::string &host) {
    if (catalog == 7)
        return gaia_box_url(ra, dec, w_arcmin, maglim, host);
    const char *src, *out, *mag; // fonte, colunas, banda (extras.py VIZ; "%27" = ')
    switch (catalog) {
    case 4:
        src = "I/322A/out";
        out = "RAJ2000,DEJ2000,pmRA,pmDE,rmag";
        mag = "rmag";
        break;
    case 5:
        src = "I/327/cmc15";
        out = "RAJ2000,DEJ2000,r%27mag";
        mag = "r%27mag";
        break;
    case 6:
        src = "I/317/sample";
        out = "RAJ2000,DEJ2000,pmRA,pmDE,r1mag";
        mag = "r1mag";
        break;
    case 8:
        src = "J/ApJ/867/105/refcat2";
        out = "RA_ICRS,DE_ICRS,pmRA,pmDE,rmag";
        mag = "rmag";
        break;
    default:
        return "";
    }
    char b[640];
    std::snprintf(b, sizeof b,
                  "https://%s/viz-bin/"
                  "asu-tsv?-source=%s&-c=%.6f%%20%+.6f&-out=%s&-out.max=20000&-oc.form=d&%s=%%3C%.1f&-c.geom=b&-c.r=%.1f&-c.u=arcmin",
                  host.c_str(), src, ra, dec, out, mag, maglim + 0.1, w_arcmin);
    return b;
}

} // namespace astro
