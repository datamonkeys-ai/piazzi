// Internal session state of the C API — shared by the files implementing astro_c.h.
#pragma once
#include "piazzi/astro_c.h"
#include "piazzi/reduce.hpp"
#include <cstdio>
#include <map>
#include <string>
#include <vector>

struct AstroSession {
    astro::Settings s = astro::Settings::defaults();
    std::vector<std::string> images; // paths, in the order they were added
    astro::Reduction R;
    bool ok = false;                     // result of the last astro_reduce
    std::string buf, json, report, getv; // string return buffers
    // calibration (empty = none)
    std::string dark_path, flat_path;
    // forced epoch per image path (UTC JD), overrides the header one
    std::map<std::string, double> jd_override;
    // offline orbits (MPCORB / COMET.DAT), filled by astro_load_orbits
    struct Orbit {
        std::string packed, name;
        double H = 0, G = 0.15, epoch = 0, M = 0, w = 0, node = 0, inc = 0, e = 0, a = 0, q = 0, Tp = 0, n = 0;
        bool comet = false;
        double M1 = 0, K1 = 0;
    };
    std::vector<Orbit> orbits;
    // manual measurements/verification: notes and exclusions per (object, image)
    std::map<std::pair<int, int>, std::string> notes;
    std::map<std::pair<int, int>, bool> excluded;
    std::string buf2, buf3;
    // medida manual / track & stack (api_measure.cpp)
    astro::Catalog cat;               // catalogue used in the last reduction (filtered by catalog_bright_mag, as inside reduce)
    std::vector<astro::Track> manual; // hand-measured objects; global index = R.movers.size() + i
    std::string object;               // designation given to astro_reduce (empty: track without provID, trkSub only)
    std::string hdr;                  // astro_image_header buffer
};

// Loads image `idx` of s->images with flip, time adjustments (ExpUnit, time_reference, TimeOffsetH, TimeDelay, time_offset_s,
// forced epoch) and, if calibrate, dark/flat — like extras.load_calibrated. Throws on read error.
astro::Image session_load_image(AstroSession *s, size_t idx, bool calibrate);
// Objects for the reports: movers + manual ones, with notes applied and excluded measurements removed.
std::vector<astro::Track> session_report_tracks(AstroSession *s);
// Rebuilds s->json from s->R / s->manual (astro_c.cpp).
void session_rebuild_json(AstroSession *s);

// Plate/statistics utilities shared by api_measure.cpp and api_stack.cpp.
namespace session_util {
// RA/Dec -> pixel pela placa (parte linear + Newton), como astro_sky_to_pix / extras.radec_to_pix.
void sky2pix(const astro::Plate &p, double ra, double dec, double &px, double &py);
// mean plate scale (cx/cy units per pixel)
double plate_scale(const astro::Plate &p);
// mediana ignorando NaN (como numpy.nanmedian); reordena v
double nanmedian(std::vector<float> &v);
// escapes a string for inclusion in JSON quotes: " and \\, and controls U+0000–U+001F (\\n, \\t, ... or \\u00XX)
inline std::string json_esc(const std::string &x) {
    std::string o;
    for (char c : x) {
        unsigned char u = (unsigned char)c;
        switch (c) {
        case '"':
            o += "\\\"";
            break;
        case '\\':
            o += "\\\\";
            break;
        case '\n':
            o += "\\n";
            break;
        case '\r':
            o += "\\r";
            break;
        case '\t':
            o += "\\t";
            break;
        case '\b':
            o += "\\b";
            break;
        case '\f':
            o += "\\f";
            break;
        default:
            if (u < 0x20) {
                char b[8];
                std::snprintf(b, sizeof b, "\\u%04x", u);
                o += b;
            } else
                o += c;
        }
    }
    return o;
}
} // namespace session_util
