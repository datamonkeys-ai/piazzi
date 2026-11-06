// Online astronomy services: SkyBoT, VizieR (catalog, VSX, Gaia), SIMBAD and Sesame.
#include "app.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// trim spaces (and \r) from both ends of a field
static std::string trim_field(const std::string &x) {
    size_t a = x.find_first_not_of(' '), b = x.find_last_not_of(" \r");
    return a == std::string::npos ? "" : x.substr(a, b - a + 1);
}

bool skybot(double jd, double ra, double dec, double rd, std::string obs, std::vector<Known> &out, const std::atomic<bool> *cancel) {
    out.clear();
    if (obs.empty() || obs == "XXX")
        obs = "500";
    char u[512];
    std::snprintf(
        u, sizeof u,
        "https://ssp.imcce.fr/webservices/skybot/api/conesearch.php?-ep=%.6f&-ra=%.6f&-dec=%.6f&-rd=%.5f&-mime=text&-output=object&"
        "-observer=%s&-objFilter=110&-from=Piazzi",
        jd, ra, dec, rd, url_encode(obs).c_str());
    std::string t;
    if (!http_get(u, t, cancel, 25, 2))
        return false;
    std::stringstream ss(t);
    std::string l;
    while (std::getline(ss, l)) {
        if (l.empty() || l[0] == '#' || l.find('|') == std::string::npos)
            continue;
        std::vector<std::string> c;
        std::stringstream ls(l);
        std::string x;
        while (std::getline(ls, x, '|'))
            c.push_back(trim_field(x));
        if (c.size() < 8)
            continue;
        try {
            out.push_back({c[0] == "-" ? "" : c[0], c[1], c[4], sexa(c[2]) * 15, sexa(c[3]), std::stod(c[5]), std::stod(c[7])});
        } catch (...) {
        }
    }
    std::sort(out.begin(), out.end(), [](const Known &a, const Known &b) { return a.sep < b.sep; });
    return true;
}

bool download_catalog() {
    int cat = int(std::lround(std::stod(get("catalog").empty() ? "7" : get("catalog"))));
    std::string base = astro_catalog_url_for(A.s, cat, A.ra.c_str(), A.dec.c_str());
    for (const char *h : {"vizier.idia.ac.za", "vizier.cds.unistra.fr", "vizier.cfa.harvard.edu"}) {
        if (A.cancel)
            return false;
        status(std::string(tr("Downloading catalog from ")) + h + "…");
        std::string url = base;
        auto p = url.find("vizier.cds.unistra.fr");
        if (p != std::string::npos)
            url.replace(p, 21, h);
        std::string t;
        if (!http_get(url, t, &A.cancel, 40, 2, 200))
            continue;
        int rows = 0;
        std::stringstream ss(t);
        std::string l;
        while (std::getline(ss, l))
            if (!l.empty() && (std::isdigit((unsigned char)l[0]) || l[0] == ' '))
                ++rows;
        if (rows < 5)
            continue;
        std::string dst = prefs_dir() + "catalog.tsv";
        std::ofstream(dst, std::ios::binary) << t;
        std::lock_guard<std::recursive_mutex> lk(A.mu);
        A.catalog = dst;
        status(std::string(tr("Catalog downloaded from ")) + h + ".");
        return true;
    }
    status(tr("No VizieR mirror answered (no internet?). Load a catalog file (Gaia DR3 TSV)."));
    return false;
}

// TSV cone query to VizieR; ok tells whether the network responded.
static std::vector<std::vector<std::string>> vizier_rows(const char *src, const char *cols, double ra, double dec, double r, bool &ok,
                                                         const std::atomic<bool> *cancel) {
    std::vector<std::vector<std::string>> rows;
    std::string t;
    char u[400];
    std::snprintf(u, sizeof u,
                  "https://vizier.idia.ac.za/viz-bin/asu-tsv?-source=%s&-c=%.6f%%20%s%.6f&-c.rs=%.1f&-out=%s&-oc.form=d&-sort=_r", src, ra,
                  dec >= 0 ? "%2B" : "", dec, r, cols);
    ok = http_get(u, t, cancel, 15, 1);
    if (!ok)
        return rows;
    std::stringstream ss(t);
    std::string l;
    int n = 0;
    while (std::getline(ss, l)) {
        if (l.empty() || l[0] == '#')
            continue;
        if (++n <= 3)
            continue; // header, units and dashes
        std::vector<std::string> c;
        std::stringstream ls(l);
        std::string v;
        while (std::getline(ls, v, '\t'))
            c.push_back(trim_field(v));
        rows.push_back(c);
    }
    return rows;
}

static std::vector<std::vector<std::string>> simbad_rows(double ra, double dec, bool &ok, const std::atomic<bool> *cancel) {
    char q[600];
    std::snprintf(q, sizeof q,
                  "SELECT TOP 2 main_id, otype, DISTANCE(POINT('ICRS',ra,dec),POINT('ICRS',%.6f,%.6f))*3600 AS d FROM basic WHERE "
                  "CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',%.6f,%.6f,0.00167))=1 ORDER BY d",
                  ra, dec, ra, dec);
    std::vector<std::vector<std::string>> rows;
    std::string sb;
    ok = http_get("https://simbad.cfa.harvard.edu/simbad/sim-tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=" + url_encode(q), sb,
                  cancel, 15, 1);
    if (!ok)
        return rows;
    std::stringstream ss(sb);
    std::string l;
    std::getline(ss, l);
    while (std::getline(ss, l)) {
        std::vector<std::string> c;
        std::stringstream ls(l);
        std::string v;
        while (std::getline(ls, v, ',')) {
            v.erase(std::remove(v.begin(), v.end(), '"'), v.end());
            c.push_back(v);
        }
        if (c.size() >= 3)
            rows.push_back(c);
    }
    return rows;
}

// "What is this?": queries SkyBoT, VSX, SIMBAD and Gaia in the background for point (x, y) of the current plate.
// The thread builds a local Inspect and publishes it under A.mu only if it is still the latest query (insp_req).
void inspect_at(double x, double y) {
    double ra, dec;
    if (!astro_pix_to_sky(A.s, A.k, x, y, &ra, &dec))
        return;
    double jd = A.s->R.frames[A.k].im.jd_mid;
    std::string code = get("mpc_code");
    unsigned long req;
    auto cancel = std::make_shared<std::atomic<bool>>(false);
    {
        std::lock_guard<std::recursive_mutex> lk(A.mu);
        A.insp = {};
        A.insp.on = true;
        A.insp.loading = true;
        A.insp.x = x;
        A.insp.y = y;
        A.insp.ra = ra;
        A.insp.dec = dec;
        A.blink = false;
        req = ++A.insp_req;
        if (A.insp_cancel)
            *A.insp_cancel = true; // the previous query is stale: abort its curl
        A.insp_cancel = cancel;
    }
    A.live_threads++;
    std::thread([ra, dec, jd, code, x, y, req, cancel] {
        std::vector<Known> kn;
        bool ok_sb = skybot(jd, ra, dec, 10.0 / 3600, code, kn, cancel.get()), ok_vsx = false, ok_gaia = false, ok_simbad = false;
        auto vsx = vizier_rows("B/vsx/vsx", "Name,Type,max,min,Period", ra, dec, 8, ok_vsx, cancel.get());
        auto gaia = vizier_rows("I/355/gaiadr3", "_r,Gmag,BP-RP,Plx", ra, dec, 6, ok_gaia, cancel.get());
        auto simbad = simbad_rows(ra, dec, ok_simbad, cancel.get());
        bool offline = !ok_sb && !ok_vsx && !ok_gaia && !ok_simbad;
        App::Inspect I;
        I.on = true;
        I.loading = false;
        I.x = x;
        I.y = y;
        I.ra = ra;
        I.dec = dec;
        bool has_gaia = !gaia.empty() && gaia[0].size() > 1 && !gaia[0][1].empty();
        if (!kn.empty()) {
            I.title = (kn[0].number.empty() ? "" : "(" + kn[0].number + ") ") + kn[0].name;
            I.detail = tr("Asteroid · ") + kn[0].cls;
            I.col = T::amber;
        } else if (!vsx.empty() && !vsx[0].empty() && !vsx[0][0].empty()) {
            I.title = vsx[0][0];
            I.detail = tr("Variable star ") + (vsx[0].size() > 1 ? vsx[0][1] : "");
            I.col = T::red;
        } else if (!simbad.empty()) {
            I.title = simbad[0][0];
            I.detail = "SIMBAD · " + simbad[0][1];
            I.col = T::cyan;
        } else if (has_gaia) {
            I.title = tr("Gaia DR3 star");
            I.detail = "G " + gaia[0][1];
            I.col = T::green;
        } else if (offline) {
            I.title = tr("No network answer");
            I.detail = tr("SkyBoT, VizieR and SIMBAD did not answer (no internet?).");
            I.col = T::red;
        } else {
            I.title = tr("Nothing known here");
            I.detail = tr("No catalogued object within 6\".");
            I.col = T::dim;
        }
        if (!kn.empty())
            I.rows.push_back({"SkyBoT", kn[0].name});
        for (auto &v : vsx)
            if (!v.empty() && !v[0].empty())
                I.rows.push_back({"VSX", v[0] + " · " + (v.size() > 1 ? v[1] : "")});
        for (auto &s : simbad)
            I.rows.push_back({"SIMBAD", s[0] + " · " + s[1]});
        if (has_gaia)
            I.rows.push_back({"Gaia DR3", "G " + gaia[0][1] + (gaia[0].size() > 2 && !gaia[0][2].empty() ? " · BP-RP " + gaia[0][2] : "")});
        {
            std::lock_guard<std::recursive_mutex> lk(A.mu);
            if (!A.quit && A.insp.on && A.insp_req == req)
                A.insp = std::move(I);
        }
        A.live_threads--;
    }).detach();
}

bool sesame_resolve(const std::string &name, double &ra, double &dec, bool &net_ok, const std::atomic<bool> *cancel) {
    std::string q = url_encode(name), t;
    net_ok = false;
    for (const char *b : {"https://vizier.cfa.harvard.edu/viz-bin/nph-sesame/-oI/A?", "https://cds.unistra.fr/cgi-bin/nph-sesame/-oI/A?"}) {
        if (!http_get(std::string(b) + q, t, cancel, 20, 1))
            continue;
        net_ok = true;
        auto p = t.find("%J ");
        double a, d;
        if (p != std::string::npos && std::sscanf(t.c_str() + p + 3, "%lf %lf", &a, &d) == 2) {
            ra = a;
            dec = d;
            return true;
        }
    }
    return false;
}
