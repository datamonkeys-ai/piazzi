// Session and automatic pipeline: images -> catalog -> plates -> motion -> identification (on the worker thread).
#include "app.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

App A;

void status(const std::string &t) {
    std::lock_guard<std::recursive_mutex> l(A.mu);
    A.status = t;
}
static void step(int i, int st, const std::string &d = "") {
    std::lock_guard<std::recursive_mutex> l(A.mu);
    A.steps[i].st = st;
    A.steps[i].detail = d;
    int e = (int)std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - A.t0).count();
    char b[16];
    std::snprintf(b, sizeof b, "%02d:%02d", e / 60, e % 60);
    A.steps[i].at = b;
}
void fresh_steps() {
    std::lock_guard<std::recursive_mutex> l(A.mu);
    A.steps = {{"Images"}, {"Reference catalog"}, {"Astrometric plate"}, {"Moving objects"}, {"Identification"}};
    A.t0 = std::chrono::steady_clock::now();
}

void clear_session() {
    std::lock_guard<std::recursive_mutex> l(A.mu);
    astro_clear_images(A.s);
    A.images.clear();
    A.catalog.clear();
    A.ra.clear();
    A.dec.clear();
    A.object.clear();
    A.known.clear();
    A.identified.clear();
    A.insp = {};
    A.meas = {};
    A.manual_track = -1;
}
void add_file(const std::string &p) {
    std::string e = lower(fs::path(p).extension().string());
    if (e == ".cfg")
        astro_load_cfg(A.s, p.c_str());
    else if (e == ".txt" || e == ".tsv")
        A.catalog = p;
    else if (e == ".dat") {
        if (lower(fs::path(p).filename().string()).find("comet") != std::string::npos)
            A.comets = p;
        else
            A.mpcorb = p;
        A.norbits = astro_load_orbits(A.s, A.mpcorb.c_str(), A.comets.c_str()); // as in Settings: offline orbits
    } else {
        A.images.push_back(p);
        astro_add_image(A.s, p.c_str());
    }
}

// FITS header text ("KEY = value / comment") -> map key -> unquoted value
static std::map<std::string, std::string> parse_header(const std::string &txt) {
    auto tr = [](const std::string &x) {
        size_t a = x.find_first_not_of(" '"), b = x.find_last_not_of(" '\r");
        return a == std::string::npos ? std::string() : x.substr(a, b - a + 1);
    };
    std::map<std::string, std::string> h;
    std::stringstream ss(txt);
    std::string l;
    while (std::getline(ss, l)) {
        auto eq = l.find('=');
        if (eq == std::string::npos)
            continue;
        std::string k = l.substr(0, eq), v = l.substr(eq + 1);
        auto sl = v.find(" /");
        if (v.find('\'') == std::string::npos && sl != std::string::npos)
            v = v.substr(0, sl);
        if (auto q = v.find('\''); q != std::string::npos) {
            auto q2 = v.find('\'', q + 1);
            v = v.substr(q + 1, q2 == std::string::npos ? std::string::npos : q2 - q - 1);
        }
        h[tr(k)] = tr(v);
    }
    return h;
}
void apply_header() {
    if (A.images.empty())
        return;
    auto h = parse_header(astro_image_header(A.s, 0));
    try {
        if (A.ra.empty() || A.dec.empty()) {
            if (h.count("OBJCTRA") && h.count("OBJCTDEC")) {
                A.ra = h["OBJCTRA"];
                A.dec = h["OBJCTDEC"];
            } else if (h.count("CRVAL1") && h.count("CRVAL2")) {
                A.ra = fmt_hms(std::stod(h["CRVAL1"]) / 15, false);
                A.dec = fmt_hms(std::stod(h["CRVAL2"]), true);
            }
        }
        if (h.count("FOCALLEN") && std::stod(h["FOCALLEN"]) > 0)
            set("focal_length_mm", h["FOCALLEN"]);
        if (h.count("XPIXSZ") && std::stod(h["XPIXSZ"]) > 0) {
            setn("pixel_width_mm", std::stod(h["XPIXSZ"]) / 1000);
            setn("pixel_height_mm", std::stod(h.count("YPIXSZ") ? h["YPIXSZ"] : h["XPIXSZ"]) / 1000);
        }
    } catch (...) { // header with an invalid number: continue without those values
    }
}

static std::string known_label(const Known &a) {
    char b[200];
    std::snprintf(b, sizeof b, "%s%s · %s · V %.1f · %.0f\"", a.number.empty() ? "" : ("(" + a.number + ") ").c_str(), a.name.c_str(),
                  a.cls.c_str(), a.v, a.sep);
    return b;
}
// numeric value after `key` in the simple JSON returned by the C API, searching from p
static double json_num(const std::string &j, const char *key, size_t p) {
    auto q = j.find(key, p);
    return q == std::string::npos ? 0.0 : std::atof(j.c_str() + q + std::strlen(key));
}
// known asteroids in the field from local orbits (MPCORB/COMET.DAT)
static std::vector<Known> known_offline(int frame) {
    std::vector<Known> kn;
    std::string j = astro_known_offline(A.s, frame, 21.0);
    size_t p = 0;
    while ((p = j.find("\"name\":\"", p)) != std::string::npos) {
        Known k;
        p += 8;
        k.name = j.substr(p, j.find('"', p) - p);
        k.ra = json_num(j, "\"ra\":", p);
        k.dec = json_num(j, "\"dec\":", p);
        k.v = json_num(j, "\"vmag\":", p);
        k.cls = tr("local orbit");
        kn.push_back(k);
    }
    return kn;
}
static void identify() {
    auto &R = A.s->R;
    if (R.frames.empty())
        return;
    std::string code = get("mpc_code");
    // field from the first frame with a solved plate; if none, skip the SkyBoT query
    double cra = 0, cdec = 0;
    int fc = -1;
    for (int k = 0; k < (int)R.frames.size(); ++k)
        if (astro_frame_center(A.s, k, &cra, &cdec)) {
            fc = k;
            break;
        }
    std::vector<Known> kn;
    bool net = false;
    if (fc >= 0) {
        const astro::Frame &f0 = R.frames[fc];
        net = skybot(f0.im.jd_mid, cra, cdec, std::max(f0.im.w, f0.im.h) * 1.1 / 3600 * 0.75, code, kn, &A.cancel);
        if (!net)
            status(A.norbits > 0 ? tr("SkyBoT did not answer (no internet?) — using the local orbits.")
                                 : tr("SkyBoT did not answer (no internet?) — no identification."));
    } else
        status(tr("No plate solved — SkyBoT not queried."));
    if (A.cancel)
        return;
    bool online = !kn.empty() || A.norbits == 0;
    if (kn.empty() && A.norbits > 0)
        kn = known_offline(fc >= 0 ? fc : 0);
    {
        std::lock_guard<std::recursive_mutex> l(A.mu);
        A.known = kn;
    }
    for (size_t i = 0; i < R.movers.size(); ++i) {
        if (A.cancel)
            return;
        const auto &q = R.movers[i].meas[0];
        std::string lab;
        std::vector<Known> l;
        if (online && net)
            skybot(q.jd, q.ra, q.dec, 2.0 / 60, code, l, &A.cancel);
        if (!l.empty()) {
            lab = known_label(l[0]);
            std::string d = A.object.empty() ? (l[0].number.empty() ? l[0].name : l[0].number) : A.object;
            astro_set_mover_desig(A.s, (int)i, d.c_str());
        } else if (A.norbits > 0) {
            std::string j = astro_identify_offline(A.s, q.jd, q.ra, q.dec, 2.0 / 60);
            auto p = j.find("\"name\":\"");
            if (p != std::string::npos) {
                p += 8;
                std::string n = j.substr(p, j.find('"', p) - p);
                double sep = json_num(j, "\"sep_arcsec\":", 0);
                char b[160];
                std::snprintf(b, sizeof b, tr("%s · local orbit · %.0f\""), n.c_str(), sep);
                lab = b;
                astro_set_mover_desig(A.s, (int)i, (A.object.empty() ? n : A.object).c_str());
            }
        }
        std::lock_guard<std::recursive_mutex> lk(A.mu);
        A.identified[(int)i] = lab; // empty = unidentified
    }
}

void pipeline() {
    if (A.images.empty()) {
        step(0, 3, tr("No images — add at least two exposures of the same field"));
        return;
    }
    step(0, 2,
         A.images.size() >= 3 ? std::to_string(A.images.size()) + tr(" images")
                              : std::to_string(A.images.size()) + tr(" image(s) — motion needs 3+"));
    if (A.ra.empty() || A.dec.empty())
        apply_header();
    if (A.ra.empty() || A.dec.empty()) {
        step(1, 3, tr("Enter the field centre"));
        A.need_center = true;
        return;
    }
    astro_set_dark(A.s, A.dark.c_str());
    astro_set_flat(A.s, A.flat.c_str());
    if (A.catalog.empty()) {
        step(1, 1);
        if (!download_catalog()) {
            step(1, 3, A.status);
            return;
        }
    }
    if (A.cancel)
        return;
    step(1, 2, fs::path(A.catalog).filename().string());
    step(2, 1);
    status(tr("Reducing the plates…"));
    auto t0 = std::chrono::steady_clock::now();
    int n = astro_reduce(A.s, A.catalog.c_str(), A.ra.c_str(), A.dec.c_str(), A.object.c_str());
    if (n < 0) {
        std::string j = astro_result_json(A.s);
        step(2, 3, core_error_message(j));
        return;
    }
    auto &R = A.s->R;
    double rms = 0;
    int nref = 0;
    for (auto &f : R.frames) {
        rms += (f.plate.dra_mean + f.plate.dde_mean) / 2;
        nref += f.plate.nref;
    }
    char b[160];
    std::snprintf(b, sizeof b, tr("%d refs per image · %.2f\" · %.2f s"), nref / std::max<int>(1, (int)R.frames.size()),
                  rms / std::max<size_t>(1, R.frames.size()), std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    step(2, 2, b);
    step(3, 2, R.movers.empty() ? tr("none") : std::to_string(R.movers.size()) + tr(" object(s)"));
    if (A.cancel)
        return;
    step(4, 1);
    status(tr("Identifying…"));
    identify();
    if (A.cancel)
        return;
    std::string lab = A.identified.count(0) ? A.identified[0].substr(0, A.identified[0].find(" ·")) : "";
    step(4, 2, R.movers.empty() ? std::to_string(A.known.size()) + tr(" known asteroids in the field") : lab);
    std::lock_guard<std::recursive_mutex> l(A.mu);
    A.tex_dirty = true;
    A.k = 0;
    A.blink = true;
    A.phase = Phase::Results;
}

void stop_worker() {
    if (A.worker.joinable()) {
        A.cancel = true;
        A.worker.join();
    }
}
void run_async(std::function<void()> fn) {
    stop_worker();
    A.cancel = false;
    A.busy = true;
    A.worker = std::thread([fn = std::move(fn)] {
        auto fail = [](const std::string &msg) {
            status(msg);
            // mark the running step (or the first unfinished one) as failed so "Try again" appears
            std::lock_guard<std::recursive_mutex> l(A.mu);
            int idx = -1;
            for (size_t i = 0; i < A.steps.size() && idx < 0; ++i)
                if (A.steps[i].st == 1)
                    idx = (int)i;
            for (size_t i = 0; i < A.steps.size() && idx < 0; ++i)
                if (A.steps[i].st != 2)
                    idx = (int)i;
            if (idx >= 0)
                step(idx, 3, msg);
        };
        try {
            fn();
        } catch (std::exception &e) {
            fail(std::string(tr("Error: ")) + e.what());
        } catch (...) {
            fail(tr("Unexpected error in the worker thread."));
        }
        A.busy = false;
    });
}
void start_session(std::function<void()> prep) {
    // first wait for the old thread (without the lock, which it also uses); only then swap steps and phase, under A.mu,
    // so it does not write to A.steps/A.phase after they are recreated
    stop_worker();
    {
        std::lock_guard<std::recursive_mutex> l(A.mu);
        fresh_steps();
        A.phase = Phase::Session;
    }
    run_async([prep = std::move(prep)] {
        prep();
        pipeline();
    });
}
// restore the settings example.cfg changed (call with the worker thread stopped or from inside it)
static void restore_user_settings() {
    for (auto &kv : A.user_settings)
        set(kv.first, kv.second);
    A.user_settings.clear();
}
void start_with_files(const std::vector<std::string> &files) {
    start_session([files] {
        if (A.example.exchange(false))
            restore_user_settings();
        clear_session();
        for (auto &p : files)
            add_file(p);
        apply_header();
    });
}
// built-in example (assets/example): images, catalog and center come from example.cfg — nothing is downloaded
void start_example() {
    std::string d = assets_dir() + "example/", cfg = d + "example.cfg";
    if (!fs::exists(cfg)) {
        status(tr("Example files not found: ") + cfg);
        return;
    }
    start_session([d, cfg] {
        if (!A.example.exchange(true)) { // save the user settings the example will overwrite
            std::ifstream in(cfg);
            std::string ln;
            while (std::getline(in, ln)) {
                auto eq = ln.find('=');
                if (eq == std::string::npos || ln.empty() || ln[0] == ';' || ln[0] == '[')
                    continue;
                std::string k = ln.substr(0, eq);
                k.erase(k.find_last_not_of(" \t") + 1);
                A.user_settings[k] = get(k);
            }
        }
        clear_session();
        astro_load_cfg(A.s, cfg.c_str());
        std::string imgs = get("example_images");
        std::stringstream ss(imgs);
        std::string f;
        while (std::getline(ss, f, ';'))
            if (!f.empty())
                add_file(d + f);
        std::lock_guard<std::recursive_mutex> l(A.mu);
        A.catalog = d + get("example_catalog");
        A.ra = get("example_ra");
        A.dec = get("example_dec");
    });
}
void rerun_pipeline() {
    stop_worker(); // see start_session
    {
        std::lock_guard<std::recursive_mutex> l(A.mu);
        fresh_steps();
        A.phase = Phase::Session;
    }
    run_async(pipeline);
}
void go_home() {
    stop_worker();
    if (A.example.exchange(false))
        restore_user_settings();
    clear_session();
    A.phase = Phase::Home;
}
