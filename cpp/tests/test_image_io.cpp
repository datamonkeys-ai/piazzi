// Robustness of the SBIG/FITS readers against malformed files: each case must throw, never crash or read out of bounds.
// Also checks reading a valid uncompressed SBIG (pixels, exposure and header instant).
// Usage: test_image_io (writes synthetic files to the system temporary directory)
#include "piazzi/core.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

static int fails = 0;

static std::string write_tmp(const std::string &name, const std::vector<unsigned char> &bytes) {
    std::string p = (fs::temp_directory_path() / ("astro_test_image_io_" + name)).string();
    std::ofstream(p, std::ios::binary).write(reinterpret_cast<const char *>(bytes.data()), std::streamsize(bytes.size()));
    return p;
}

// 2048-byte SBIG header (zero-terminated text) followed by the data
static std::vector<unsigned char> sbig(const std::string &first, const std::string &w, const std::string &h,
                                       const std::vector<unsigned char> &data) {
    std::string t = first + "\r\nWidth = " + w + "\r\nHeight = " + h + "\r\nEnd\r\n";
    std::vector<unsigned char> v(2048, 0);
    std::copy(t.begin(), t.end(), v.begin());
    v.insert(v.end(), data.begin(), data.end());
    return v;
}

static std::vector<unsigned char> fits(const std::vector<std::string> &cards, size_t ndata) {
    std::string h;
    for (auto c : cards) {
        c.resize(80, ' ');
        h += c;
    }
    std::string end = "END";
    end.resize(80, ' ');
    h += end;
    h.resize((h.size() + 2879) / 2880 * 2880, ' ');
    std::vector<unsigned char> v(h.begin(), h.end());
    v.resize(v.size() + ndata, 0);
    return v;
}

static void expect_throw(const char *what, const std::string &name, const std::vector<unsigned char> &bytes) {
    std::string p = write_tmp(name, bytes);
    try {
        astro::load_image(p, false, false);
        std::printf("FAIL %s: no exception\n", what);
        ++fails;
    } catch (const std::exception &e) {
        std::printf("ok    %s: %s\n", what, e.what());
    }
    std::error_code ec;
    fs::remove(p, ec);
}

static void expect_ok(const char *what, const std::string &name, const std::vector<unsigned char> &bytes, int w, int h) {
    std::string p = write_tmp(name, bytes);
    try {
        astro::Image im = astro::load_image(p, false, false);
        if (im.w != w || im.h != h) {
            std::printf("FAIL %s: size %dx%d\n", what, im.w, im.h);
            ++fails;
        } else
            std::printf("ok    %s\n", what);
    } catch (const std::exception &e) {
        std::printf("FAIL %s: %s\n", what, e.what());
        ++fails;
    }
    std::error_code ec;
    fs::remove(p, ec);
}

int main() {
    const std::string C = "ST-7 Compressed Image", U = "ST-7 Image";
    // compressed SBIG: the line declares 3 bytes (initial value + 1 delta) but the width requires 8 pixels
    expect_throw("SBIG truncated line", "t1.st7", sbig(C, "8", "1", {3, 0, 10, 0, 1}));
    // escape delta 0x80 without the 2 following bytes
    expect_throw("SBIG truncated escape", "t2.st7", sbig(C, "4", "1", {3, 0, 10, 0, 0x80}));
    // line shorter than 2 bytes
    expect_throw("SBIG line < 2 bytes", "t3.st7", sbig(C, "4", "1", {1, 0, 10}));
    // line length past end of file
    expect_throw("SBIG length past end", "t4.st7", sbig(C, "4", "1", {200, 0, 10, 0}));
    expect_throw("SBIG negative width", "t5.st7", sbig(U, "-5", "10", std::vector<unsigned char>(100)));
    expect_throw("SBIG zero height", "t6.st7", sbig(U, "10", "0", std::vector<unsigned char>(100)));
    expect_throw("SBIG huge", "t7.st7", sbig(U, "2000000000", "2000000000", std::vector<unsigned char>(100)));
    expect_throw("SBIG 65536", "t8.st7", sbig(U, "65536", "1", std::vector<unsigned char>(100)));
    expect_throw("SBIG short data", "t9.st7", sbig(U, "100", "100", std::vector<unsigned char>(100)));
    expect_throw("SBIG short file", "t10.st7", std::vector<unsigned char>(100));
    // valid compressed: 4 pixels = 10, 11, 9, 300
    expect_ok("SBIG valid compressed", "t11.st7", sbig(C, "4", "1", {7, 0, 10, 0, 1, 0xFE, 0x80, 0x2C, 0x01}), 4, 1);
    // valid uncompressed SBIG with Exposure (hundredths of s), Date (mm/dd/yy) and Time (start): 3×2 pixels LE
    {
        std::vector<unsigned char> px;
        for (int i = 0; i < 6; ++i) {
            px.push_back((unsigned char)((1000 + 7 * i) & 0xFF));
            px.push_back((unsigned char)((1000 + 7 * i) >> 8));
        }
        std::string p = write_tmp("t12.st6", sbig("ST-6 Image\r\nExposure = 6000\r\nDate = 01/14/00\r\nTime = 19:06:11", "3", "2", px));
        try {
            astro::Image im = astro::load_image(p, false, false, true);
            // 2000-01-14 19:06:11 UTC + 30 s = JD 2451558.29630787
            bool ok = im.w == 3 && im.h == 2 && im.pix[0] == 1000 && im.pix[5] == 1035 && im.exptime == 60.0 &&
                      std::fabs(im.jd_mid - 2451558.2963079) < 2e-7;
            std::printf("%s SBIG uncompressed: %dx%d pix0 %d exp %.1f s jd_mid %.7f\n", ok ? "ok   " : "FAIL ", im.w, im.h,
                        im.w > 0 ? int(im.pix[0]) : -1, im.exptime, im.jd_mid);
            if (!ok)
                ++fails;
        } catch (const std::exception &e) {
            std::printf("FAIL SBIG uncompressed: %s\n", e.what());
            ++fails;
        }
        std::error_code ec;
        fs::remove(p, ec);
    }

    const std::vector<std::string> base = {"SIMPLE  =                    T", "BITPIX  =                   16",
                                           "NAXIS   =                    2"};
    auto with = [&](const std::string &n1, const std::string &n2) {
        auto c = base;
        c.push_back("NAXIS1  = " + n1);
        c.push_back("NAXIS2  = " + n2);
        return c;
    };
    expect_throw("FITS negative NAXIS1", "f1.fits", fits(with("-10", "10"), 2880));
    expect_throw("FITS zero NAXIS2", "f2.fits", fits(with("10", "0"), 2880));
    expect_throw("FITS huge", "f3.fits", fits(with("4000000000", "4000000000"), 2880));
    expect_throw("FITS 65536", "f4.fits", fits(with("65536", "1"), 2880));
    expect_throw("FITS short data", "f5.fits", fits(with("100", "100"), 2880));
    expect_throw("FITS without END", "f6.fits", std::vector<unsigned char>(100, ' '));
    {
        auto c = with("10", "10");
        c[1] = "BITPIX  =                   12";
        expect_throw("FITS invalid BITPIX", "f7.fits", fits(c, 2880));
    }
    expect_ok("FITS valid", "f8.fits", fits(with("10", "10"), 2880), 10, 10);

    std::printf(fails ? "%d failure(s)\n" : "OK\n", fails);
    return fails ? 1 : 0;
}
