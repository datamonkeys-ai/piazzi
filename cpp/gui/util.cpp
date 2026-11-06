// App utilities: sexagesimal formatting, paths and files.
#include "app.hpp"
#include <SDL.h>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

// Sexagesimal "HH MM SS.s" (sign: "+DD MM SS.s"). Rounds to tenths of a second BEFORE splitting, so it never
// yields "59 60.0"; unsigned (RA) is reduced modulo 24; signed, the "-" also applies between -1° and 0° and is
// dropped if the rounded value is zero.
std::string fmt_hms(double v, bool sign) {
    if (!std::isfinite(v))
        v = 0;
    bool neg = v < 0;
    long long t = std::llround(std::fabs(v) * 36000.0); // tenths of a second
    if (!sign) {
        const long long day = 24LL * 36000;
        t = neg ? (day - t % day) % day : t % day;
        neg = false;
    }
    if (t == 0)
        neg = false;
    long long h = t / 36000, m = (t / 600) % 60, ds = t % 600;
    char b[48];
    std::snprintf(b, sizeof b, "%s%02lld %02lld %02lld.%lld", sign ? (neg ? "-" : "+") : "", h, m, ds / 10, ds % 10);
    return b;
}
double sexa(const std::string &s) {
    double a = 0, b = 0, c = 0;
    int n = std::sscanf(s.c_str(), "%lf %lf %lf", &a, &b, &c);
    if (n <= 1)
        return a;
    double sg = s.find('-') != std::string::npos ? -1 : 1;
    return sg * (std::fabs(a) + b / 60 + c / 3600);
}
std::string lower(std::string x) {
    for (auto &c : x)
        c = (char)std::tolower((unsigned char)c);
    return x;
}
std::string url_encode(const std::string &s) {
    std::string enc;
    for (unsigned char ch : s) {
        if (std::isalnum(ch))
            enc += (char)ch;
        else {
            char h[4];
            std::snprintf(h, sizeof h, "%%%02X", ch);
            enc += h;
        }
    }
    return enc;
}
// On first run (new "Piazzi" folder empty), copies preferences from the old "Astrometria" folder (profiles.txt,
// downloaded catalog, etc.), if it exists. The old folder is left intact.
static void migrate_old_prefs(const std::string &dst) {
    std::error_code ec;
    if (!fs::is_directory(dst, ec) || !fs::is_empty(dst, ec) || ec)
        return;
    char *op = SDL_GetPrefPath("datamonkeys", "Astrometria"); // creates the folder if missing: removed below if empty
    if (!op)
        return;
    std::string old = op;
    SDL_free(op);
    if (fs::is_empty(old, ec)) {
        fs::remove(old, ec);
        return;
    }
    for (const auto &e : fs::directory_iterator(old, ec)) {
        std::string n = e.path().filename().string();
        if (n.rfind("http_", 0) == 0) // skip download temp files
            continue;
        fs::copy(e.path(), fs::path(dst) / n, fs::copy_options::recursive | fs::copy_options::skip_existing, ec);
    }
}
std::string prefs_dir() {
    char *p = SDL_GetPrefPath("datamonkeys", "Piazzi");
    bool ok = p != nullptr;
    std::string r = ok ? p : "./";
    SDL_free(p);
    static std::once_flag once;
    if (ok)
        std::call_once(once, [&r] { migrate_old_prefs(r); });
    return r;
}
std::string assets_dir() {
    char *p = SDL_GetBasePath();
    std::string r = p ? p : "./";
    SDL_free(p);
    return r + "assets/";
}
std::string output_dir() {
    std::string d = A.outdir.empty() ? prefs_dir() : A.outdir;
    if (d.back() != '/' && d.back() != '\\')
        d += '/';
    return d;
}
std::string read_file(const std::string &p) {
    std::ifstream f(p, std::ios::binary);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}
