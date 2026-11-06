// Reading of SBIG (ST-x, compressed or not) and FITS images — mirrors piazzi/sbig.py and pipeline.load.
#include "piazzi/core.hpp"

#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <stdexcept>

namespace astro {
namespace {

std::vector<unsigned char> read_file(const std::string &path) {
    std::ifstream f(path, std::ios::binary);
    if (!f)
        throw std::runtime_error("cannot open " + path);
    return std::vector<unsigned char>((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

std::string trim(const std::string &s) {
    size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
}

// Gregorian -> Julian day number (noon), Fliegel & Van Flandern algorithm.
long jdn(long y, long m, long d) {
    long a = (14 - m) / 12, yy = y + 4800 - a, mm = m + 12 * a - 3;
    return d + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32045;
}

// UTC JD from civil date + seconds of day.
double jd_utc(int y, int m, int d, double sec_of_day) {
    return double(jdn(y, m, d)) - 0.5 + sec_of_day / 86400.0;
}

double parse_hms(const std::string &t) { // "HH:MM:SS(.sss)"
    int h = 0, mi = 0;
    double s = 0;
    const char *p = t.c_str();
    char *e;
    h = int(std::strtol(p, &e, 10));
    if (*e == ':')
        p = e + 1;
    else
        return h * 3600.0;
    mi = int(std::strtol(p, &e, 10));
    if (*e == ':')
        p = e + 1;
    else
        return h * 3600.0 + mi * 60.0;
    s = std::strtod(p, &e);
    return h * 3600.0 + mi * 60.0 + s;
}

inline unsigned rd16le(const unsigned char *p) {
    return unsigned(p[0]) | (unsigned(p[1]) << 8);
}

// Validates dimensions (0 < w,h <= 65535) and returns w*h*bpp in 64 bits, rejecting what does not fit in size_t
// (wasm32 has a 32-bit size_t) — avoids multiplication overflow before allocating.
size_t checked_size(long w, long h, int bpp, const char *fmt) {
    if (w <= 0 || h <= 0 || w > 65535 || h > 65535)
        throw std::runtime_error(std::string(fmt) + ": bad size");
    uint64_t n = uint64_t(w) * uint64_t(h) * uint64_t(bpp);
    if (n > uint64_t(SIZE_MAX) / 2)
        throw std::runtime_error(std::string(fmt) + ": image too large");
    return size_t(n);
}

void flip(Image &im, bool fh, bool fv) {
    if (fh)
        for (int y = 0; y < im.h; ++y)
            for (int x = 0; x < im.w / 2; ++x)
                std::swap(im.pix[size_t(y) * im.w + x], im.pix[size_t(y) * im.w + im.w - 1 - x]);
    if (fv)
        for (int y = 0; y < im.h / 2; ++y)
            for (int x = 0; x < im.w; ++x)
                std::swap(im.pix[size_t(y) * im.w + x], im.pix[size_t(im.h - 1 - y) * im.w + x]);
}

Image load_sbig(const std::vector<unsigned char> &raw, bool time_is_start) {
    if (raw.size() < 2048)
        throw std::runtime_error("SBIG: short file");
    size_t hl = 0;
    while (hl < 2048 && raw[hl])
        ++hl;
    std::string txt(raw.begin(), raw.begin() + hl);
    std::vector<std::string> lines;
    std::string cur;
    for (char c : txt) {
        if (c == '\r' || c == '\n') {
            if (!trim(cur).empty())
                lines.push_back(trim(cur));
            cur.clear();
        } else
            cur += c;
    }
    if (!trim(cur).empty())
        lines.push_back(trim(cur));
    if (lines.empty())
        throw std::runtime_error("SBIG: empty header");
    std::map<std::string, std::string> hdr;
    for (size_t i = 1; i < lines.size(); ++i) {
        size_t e = lines[i].find('=');
        if (e != std::string::npos)
            hdr[trim(lines[i].substr(0, e))] = trim(lines[i].substr(e + 1));
    }
    Image im;
    long sw = std::strtol(hdr["Width"].c_str(), nullptr, 10), sh = std::strtol(hdr["Height"].c_str(), nullptr, 10);
    size_t nbytes = checked_size(sw, sh, 2, "SBIG");
    im.w = int(sw);
    im.h = int(sh);
    const unsigned char *data = raw.data() + 2048;
    size_t n_data = raw.size() - 2048;
    bool compressed = lines[0].find("Compressed") != std::string::npos;
    if (!compressed && n_data < nbytes)
        throw std::runtime_error("SBIG: short data");
    im.pix.assign(nbytes / 2, 0);
    if (!compressed) {
        for (size_t i = 0; i < im.pix.size(); ++i)
            im.pix[i] = uint16_t(rd16le(data + 2 * i));
    } else {
        size_t p = 0;
        for (int y = 0; y < im.h; ++y) {
            if (p + 2 > n_data)
                throw std::runtime_error("SBIG: truncated");
            size_t n = rd16le(data + p);
            p += 2;
            if (p + n > n_data)
                throw std::runtime_error("SBIG: truncated");
            const unsigned char *row = data + p;
            p += n;
            uint16_t *out = &im.pix[size_t(y) * im.w];
            if (n == size_t(im.w) * 2) {
                for (int x = 0; x < im.w; ++x)
                    out[x] = uint16_t(rd16le(row + 2 * x));
                continue;
            }
            if (n < 2)
                throw std::runtime_error("SBIG: linha truncada");
            size_t q = 0;
            unsigned v = rd16le(row);
            q += 2;
            out[0] = uint16_t(v);
            for (int x = 1; x < im.w; ++x) {
                if (q + 1 > n)
                    throw std::runtime_error("SBIG: linha truncada");
                unsigned b = row[q++];
                if (b == 0x80) {
                    if (q + 2 > n)
                        throw std::runtime_error("SBIG: linha truncada");
                    v = rd16le(row + q);
                    q += 2;
                } else
                    v = (v + unsigned(int(b > 127 ? int(b) - 256 : int(b)))) & 0xFFFFu;
                out[x] = uint16_t(v);
            }
        }
    }
    im.exptime = std::atoi(hdr["Exposure"].c_str()) / 100.0;
    int mm = 0, dd = 0, yy = 0;
    std::sscanf(hdr["Date"].c_str(), "%d/%d/%d", &mm, &dd, &yy);
    yy += yy < 50 ? 2000 : 1900;
    // SBIG: "Time" is the exposure start; jd_mid = mid-exposure (as in FITS with time_reference=start)
    im.jd_mid = jd_utc(yy, mm, dd, parse_hms(hdr["Time"]) + (time_is_start ? im.exptime / 2 : 0.0));
    return im;
}

Image load_fits(const std::vector<unsigned char> &raw, bool time_is_start) {
    std::map<std::string, std::string> hdr;
    size_t pos = 0;
    bool end = false;
    while (!end) {
        if (pos + 80 > raw.size())
            throw std::runtime_error("FITS: no END");
        std::string card(raw.begin() + pos, raw.begin() + pos + 80);
        pos += 80;
        std::string key = trim(card.substr(0, 8));
        if (key == "END") {
            end = true;
            break;
        }
        if (card.size() >= 10 && card[8] == '=') {
            std::string v = card.substr(10);
            std::string val;
            size_t a = v.find_first_not_of(' ');
            if (a != std::string::npos && v[a] == '\'') {
                size_t b = a + 1;
                while (b < v.size()) {
                    if (v[b] == '\'') {
                        if (b + 1 < v.size() && v[b + 1] == '\'') {
                            val += '\'';
                            b += 2;
                            continue;
                        }
                        break;
                    }
                    val += v[b++];
                }
                val = trim(val);
            } else {
                size_t sl = v.find('/');
                val = trim(v.substr(0, sl));
            }
            if (!hdr.count(key))
                hdr[key] = val;
        }
    }
    pos = (pos + 2879) / 2880 * 2880;
    auto num = [&](const char *k, double def) {
        auto it = hdr.find(k);
        return it == hdr.end() ? def : std::atof(it->second.c_str());
    };
    int bitpix = int(num("BITPIX", 0)), naxis = int(num("NAXIS", 0));
    if (naxis != 2 && !(naxis == 3 && int(num("NAXIS3", 1)) == 1))
        throw std::runtime_error("FITS: unsupported NAXIS");
    Image im;
    double d1 = num("NAXIS1", 0), d2 = num("NAXIS2", 0);
    if (!(d1 > 0 && d1 <= 65535 && d2 > 0 && d2 <= 65535)) // also rejects NaN before the int conversion
        throw std::runtime_error("FITS: bad size");
    im.w = int(d1);
    im.h = int(d2);
    double bzero = num("BZERO", 0), bscale = num("BSCALE", 1);
    if (!(bitpix == 8 || bitpix == 16 || bitpix == 32 || bitpix == -32 || bitpix == -64))
        throw std::runtime_error("FITS: BITPIX");
    int bpp = std::abs(bitpix) / 8;
    size_t nbytes = checked_size(im.w, im.h, bpp, "FITS");
    size_t n = size_t(im.w) * size_t(im.h);
    if (pos > raw.size() || nbytes > raw.size() - pos)
        throw std::runtime_error("FITS: short data");
    im.pix.resize(n);
    const unsigned char *p = raw.data() + pos;
    for (size_t i = 0; i < n; ++i, p += bpp) {
        double v;
        if (bitpix == 8)
            v = p[0];
        else if (bitpix == 16)
            v = int16_t(uint16_t((p[0] << 8) | p[1]));
        else if (bitpix == 32)
            v = int32_t((uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) | (uint32_t(p[2]) << 8) | p[3]);
        else if (bitpix == -32) {
            uint32_t u = (uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) | (uint32_t(p[2]) << 8) | p[3];
            float f;
            std::memcpy(&f, &u, 4);
            v = f;
        } else {
            uint64_t u = 0;
            for (int k = 0; k < 8; ++k)
                u = (u << 8) | p[k];
            std::memcpy(&v, &u, 8);
        }
        // pipeline: astype(float32); detect: clip(rint(d), 0, 65535)
        float f = float(v * bscale + bzero);
        double r = std::nearbyint(double(f));
        if (!(r > 0))
            r = 0;
        if (r > 65535)
            r = 65535;
        im.pix[i] = uint16_t(r);
    }
    im.exptime = num("EXPTIME", num("EXPOSURE", 0));
    im.gain = num("GAIN", num("EGAIN", 0));
    if (!(im.gain > 0 && im.gain < 1000))
        im.gain = 0;
    std::string dob = hdr["DATE-OBS"], tob;
    size_t T = dob.find('T');
    if (T != std::string::npos) {
        tob = dob.substr(T + 1);
        dob = dob.substr(0, T);
    } else
        tob = hdr.count("TIME-OBS") ? hdr["TIME-OBS"] : hdr["UT"];
    int y = 0, m = 0, d = 0;
    std::sscanf(dob.c_str(), "%d-%d-%d", &y, &m, &d);
    im.jd_mid = jd_utc(y, m, d, parse_hms(tob) + (time_is_start ? im.exptime / 2 : 0.0));
    if (y == 0) { // no DATE-OBS (e.g. ZTF/IRSA cutouts): OBSJD or MJD-OBS = exposure start
        double jd = num("OBSJD", 0);
        if (jd == 0 && num("MJD-OBS", 0) != 0)
            jd = num("MJD-OBS", 0) + 2400000.5;
        if (jd != 0)
            im.jd_mid = jd + (time_is_start ? im.exptime / 2 / 86400.0 : 0.0);
    }
    return im;
}

bool ends_with_ci(const std::string &s, const char *suf) {
    size_t n = std::strlen(suf);
    if (s.size() < n)
        return false;
    for (size_t i = 0; i < n; ++i)
        if (std::tolower((unsigned char)s[s.size() - n + i]) != suf[i])
            return false;
    return true;
}

} // namespace

Image load_image(const std::string &path, bool flip_h, bool flip_v, bool time_is_start) {
    std::vector<unsigned char> raw = read_file(path);
    Image im = (ends_with_ci(path, ".fit") || ends_with_ci(path, ".fits") || ends_with_ci(path, ".fts")) ? load_fits(raw, time_is_start)
                                                                                                         : load_sbig(raw, time_is_start);
    flip(im, flip_h, flip_v);
    im.path = path;
    return im;
}

} // namespace astro
