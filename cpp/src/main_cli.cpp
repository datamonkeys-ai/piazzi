// piazzi: native command line (macOS, Windows, Linux).
//   piazzi --cfg config.cfg --catalog gaia.tsv --center 107.9314 22.4049 --object "5157" --out dir img1.fits img2.fits ...
// Without --catalog, prints the VizieR box query URL to download.
// --inject N: injection-recovery test. Injects N synthetic moving objects (each frame's real PSF convolved
// with the exposure trail, flux from the frame zero point, Poisson noise from the gain) into the real images,
// reruns the reduction and reports which were recovered and the O−C of each measurement against truth. --inject-seed S,
// --inject-mag LO HI (default 18 21.5), --inject-rate LO HI ("/min, default 0.15 0.8), --inject-snr LO HI (SNR/frame
// instead of magnitude), --inject-color DMAG (odd frames fainter), --inject-cross (over a star),
// --inject-edge (leaves the crop). Tool: tools/benchmark/run.py inject.
#include "piazzi/reduce.hpp"
#include "piazzi/version.h"
#include "reduce_internal.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <random>
#include <string>
#include <vector>

using namespace astro;

static double parse_angle(const std::string &s, bool hours) {
    // accepts decimal degrees ("125.57") or sexagesimal ("08 22 17.0" / "+29 32 31")
    double a = 0, b = 0, c = 0;
    int n = std::sscanf(s.c_str(), "%lf %lf %lf", &a, &b, &c);
    if (n <= 1)
        return a;
    double sg = (s.find('-') != std::string::npos) ? -1 : 1;
    double v = std::fabs(a) + b / 60 + c / 3600;
    return sg * v * (hours ? 15 : 1);
}

// Injection options. Scenarios (combinable): rate = rate range (long trails when v·t_exp ≫ FWHM);
// snr > 0: the flux in each frame gives this optimal (PSF) SNR in the frame sky noise, instead of magnitude; cross: the
// object passes over a star in a random frame; color: odd frames get +color mag (color differs between filters);
// edge: the object enters and leaves the crop mid-sequence.
struct InjectOpts {
    int n = 0;
    unsigned seed = 1;
    double mlo = 18, mhi = 21.5, vlo = 0.15, vhi = 0.8, snrlo = 0, snrhi = 0, color = 0;
    bool cross = false, edge = false;
    int probe = -1; // --inject-probe I: link_probe/stack_probe on object I (diagnostics)
};

// Synthetic object injection and recovery (see header). Each object is drawn with the frame's real PSF (σx, σy
// from isolated references, frame_psf) convolved with the segment traversed during the exposure (the trail model of
// Vereš et al. 2012, psf_pixel), total flux = zero-point aperture flux / aperture correction, Poisson noise
// from the gain. Truth is exact (mid-exposure position from the 1st-reduction plate); the output has one line per
// object ("inj|...") and per recovered measurement ("injm|..." with O−C vs. truth and ADES σ), for run.py inject.
static int inject_run(std::vector<Image> ims, const Catalog &cat, double ra0, double dec0, const Settings &s, const InjectOpts &io) {
    Reduction R0 = reduce(ims, cat, ra0, dec0, s, false);
    int ref = -1;
    double tsum = 0;
    int ns = 0;
    double t0 = INFINITY, t1 = -INFINITY;
    for (size_t k = 0; k < R0.frames.size(); ++k)
        if (R0.frames[k].solved && R0.frames[k].seeing > 0) {
            if (ref < 0)
                ref = int(k);
            tsum += R0.frames[k].im.jd_mid;
            t0 = std::min(t0, R0.frames[k].im.jd_mid);
            t1 = std::max(t1, R0.frames[k].im.jd_mid);
            ++ns;
        }
    if (ref < 0 || ns < 2) {
        std::fprintf(stderr, "injection: plates not solved\n");
        return 1;
    }
    const double tref = tsum / ns, span = (t1 - t0) * 1440; // min
    const Frame &F0 = R0.frames[ref];
    const double sc = 0.5 * (F0.plate.scale_x() + F0.plate.scale_y());
    double rc, dc;
    F0.plate.pix_to_radec((F0.im.w - 1) / 2.0, (F0.im.h - 1) / 2.0, rc, dc);
    const double dsat = s.str("saturation").empty() ? 60000 : s.num("saturation");
    const int ap = std::max(1, int(s.num("aperture_radius")));
    // PSF, aperture correction and sky noise of each frame
    std::vector<internal::PsfShape> psf(R0.frames.size());
    std::vector<double> apc(R0.frames.size(), 1), sky(R0.frames.size(), 0);
    std::vector<char> okf(R0.frames.size(), 0);
    for (size_t k = 0; k < R0.frames.size(); ++k) {
        const Frame &f = R0.frames[k];
        if (!f.solved || !(f.seeing > 0))
            continue;
        const double sg = f.seeing / (0.5 * (f.plate.scale_x() + f.plate.scale_y())) / 2.3548;
        if (!internal::frame_psf(f, dsat, ap, psf[k], apc[k]) || !(apc[k] > 0)) {
            psf[k] = internal::PsfShape{sg, sg, 0, 1, 0};
            apc[k] = 1;
        }
        std::vector<double> sv;
        for (auto &d : f.det)
            if (d.sigma > 0)
                sv.push_back(d.sigma);
        if (sv.empty())
            continue;
        std::nth_element(sv.begin(), sv.begin() + sv.size() / 2, sv.end());
        sky[k] = sv[sv.size() / 2];
        okf[k] = 1;
    }
    struct Obj {
        double xi, eta, vx, vy, mag, rate, pa, snr = 0;
    };
    std::mt19937 rng(io.seed);
    std::uniform_real_distribution<double> U(0, 1);
    std::normal_distribution<double> nd(0, 1);
    auto pos_t = [&](const Obj &o, const Frame &f, double jd, double &x, double &y) {
        double dt = (jd - tref) * 1440, ra, dec;
        std_to_radec(o.xi + o.vx * dt, o.eta + o.vy * dt, rc, dc, ra, dec);
        x = f.im.w / 2.0 - f.offx;
        y = f.im.h / 2.0 - f.offy;
        return plate_radec_to_pix(f.plate, ra, dec, x, y);
    };
    auto pos = [&](const Obj &o, const Frame &f, double &x, double &y) { return pos_t(o, f, f.im.jd_mid, x, y); };
    auto inside = [&](const Obj &o, size_t k, double marg) {
        const Frame &f = R0.frames[k];
        double x, y;
        return pos(o, f, x, y) && x >= marg && y >= marg && x < f.im.w - marg && y < f.im.h - marg;
    };
    // maximum rate that fits in the frame over the arc (for scenarios that do not leave the crop)
    const double fsz = std::min(F0.im.w, F0.im.h) * sc;
    double vmax = io.vhi;
    if (!io.edge && span > 0)
        vmax = std::min(vmax, 0.75 * fsz / span);
    if (vmax < io.vlo) {
        std::printf("injsum|0|0|0|infeasible: arc of %.1f min, max rate %.2f\"/min\n", span, vmax);
        return 0;
    }
    // stars for the cross scenario (unsaturated, SNR 30–2000)
    std::vector<std::pair<size_t, size_t>> stars;
    if (io.cross)
        for (size_t k = 0; k < R0.frames.size(); ++k)
            if (okf[k])
                for (size_t i = 0; i < R0.frames[k].det.size(); ++i) {
                    const Detection &d = R0.frames[k].det[i];
                    if (d.peak < 0.8 * dsat && d.snr_tot > 30 && d.snr_tot < 2000 && d.trail_len == 0)
                        stars.push_back({k, i});
                }
    std::vector<Obj> O;
    const double hw = 0.42 * F0.im.w * sc, hh = 0.42 * F0.im.h * sc;
    for (int tries = 0; int(O.size()) < io.n && tries < 20000; ++tries) {
        Obj o;
        // log-uniform rate when the range exceeds a factor of 5 (equal weight per decade)
        o.rate = vmax > 5 * io.vlo && io.vlo > 0 ? io.vlo * std::pow(vmax / io.vlo, U(rng)) : io.vlo + (vmax - io.vlo) * U(rng);
        o.pa = 360 * U(rng);
        o.vx = o.rate * std::sin(o.pa * PI / 180);
        o.vy = o.rate * std::cos(o.pa * PI / 180);
        o.mag = io.mlo + (io.mhi - io.mlo) * U(rng);
        if (io.snrhi > 0)
            o.snr = io.snrlo + (io.snrhi - io.snrlo) * U(rng);
        o.xi = (2 * U(rng) - 1) * hw;
        o.eta = (2 * U(rng) - 1) * hh;
        if (io.cross && !stars.empty()) {
            // trail center over the star in the random frame, within 0.5 FWHM of it perpendicularly
            auto st = stars[size_t(U(rng) * stars.size()) % stars.size()];
            const Frame &f = R0.frames[st.first];
            const Detection &d = f.det[st.second];
            double ra, dec, xi, eta;
            f.plate.pix_to_radec(d.x, d.y, ra, dec);
            radec_to_std(ra, dec, rc, dc, xi, eta);
            double off = (2 * U(rng) - 1) * 0.5 * f.seeing, dt = (f.im.jd_mid - tref) * 1440;
            o.xi = xi - o.vx * dt + off * std::cos(o.pa * PI / 180);
            o.eta = eta - o.vy * dt - off * std::sin(o.pa * PI / 180);
        }
        int nin = 0;
        for (size_t k = 0; k < R0.frames.size(); ++k)
            nin += okf[k] && inside(o, k, 20 + 0.5 * o.rate * R0.frames[k].im.exptime / 60 / sc);
        if (io.edge) {
            // enters or leaves the crop during the sequence: inside in ≥ 2 frames and outside in ≥ 1
            if (nin < 2 || nin == ns)
                continue;
        } else if (nin < ns)
            continue;
        bool far = true;
        for (auto &q : O)
            far &= std::hypot(q.xi - o.xi, q.eta - o.eta) > 60;
        if (far)
            O.push_back(o);
    }
    for (size_t k = 0; k < ims.size() && k < R0.frames.size(); ++k) {
        const Frame &f = R0.frames[k];
        if (!okf[k])
            continue;
        Image &im = ims[k];
        const double g = R0.gain > 0 ? R0.gain : 1;
        const double sgx = psf[k].sx, sgy = psf[k].sy, sgm = std::sqrt(sgx * sgy);
        for (auto &o : O) {
            double x0, y0, xa, ya, xb, yb;
            if (!pos(o, f, x0, y0))
                continue;
            const double te = f.im.exptime / 2 / 86400;
            if (!pos_t(o, f, f.im.jd_mid - te, xa, ya) || !pos_t(o, f, f.im.jd_mid + te, xb, yb))
                continue;
            const double L = std::hypot(xb - xa, yb - ya);
            internal::PsfShape P{sgx, sgy, L, L > 1e-6 ? (xb - xa) / L : 1, L > 1e-6 ? (yb - ya) / L : 0};
            double mag = o.mag + ((io.color != 0 && (k % 2)) ? io.color : 0);
            double Fl = std::pow(10.0, 0.4 * (f.plate.zp - mag)) / apc[k];
            if (o.snr > 0) {
                // optimal (PSF) SNR in the frame noise: SNR² (F/g + n_eff σ²) = F², n_eff = 4πσ² (1 + L/(2√π σ))
                const double neff = 4 * PI * sgm * sgm * (1 + L / (2 * std::sqrt(PI) * sgm)), S2 = o.snr * o.snr;
                Fl = 0.5 * (S2 / g + std::sqrt(S2 * S2 / (g * g) + 4 * S2 * neff * sky[k] * sky[k]));
                if (int(k) == ref)
                    o.mag = f.plate.zp - 2.5 * std::log10(Fl * apc[k]);
            }
            const int Rr = int(std::ceil(4 * std::max(sgx, sgy) + L / 2)) + 1;
            for (int y = int(y0) - Rr; y <= int(y0) + Rr; ++y)
                for (int x = int(x0) - Rr; x <= int(x0) + Rr; ++x) {
                    if (x < 0 || y < 0 || x >= im.w || y >= im.h)
                        continue;
                    double e = Fl * internal::psf_pixel(P, x - x0, y - y0); // expected ADU in the pixel
                    if (!(e > 1e-4))
                        continue;
                    double v = im.pix[size_t(y) * im.w + x] + e + std::sqrt(e / g) * nd(rng);
                    im.pix[size_t(y) * im.w + x] = uint16_t(std::clamp(std::lround(v), 0L, 65535L));
                }
        }
    }
    Settings s2 = s;
    if (io.probe >= 0 && io.probe < int(O.size())) {
        double ra, dec;
        char pv[200];
        std_to_radec(O[size_t(io.probe)].xi, O[size_t(io.probe)].eta, rc, dc, ra, dec);
        std::snprintf(pv, sizeof pv, "%.7f %.7f %.7f %.4f %.4f", ra, dec, tref, O[size_t(io.probe)].vx * 60, O[size_t(io.probe)].vy * 60);
        s2.set("link_probe", pv);
        s2.set("stack_probe", pv);
    }
    Reduction R = reduce(ims, cat, ra0, dec0, s2, true);
    if (io.probe >= 0)
        for (auto &l : R.log)
            if (l.find("probe") != std::string::npos || l.rfind("trail", 0) == 0)
                std::printf("  · %s\n", l.c_str());
    std::vector<char> hit(R.movers.size(), 0);
    int nrec = 0;
    for (size_t i = 0; i < O.size(); ++i) {
        int best = -1;
        size_t nm = 0;
        for (size_t t = 0; t < R.movers.size(); ++t) {
            const Track &T = R.movers[t];
            int ok = 0;
            for (auto &m : T.meas) {
                double x, y;
                if (pos(O[i], R.frames[m.frame], x, y) && std::hypot(x - m.x, y - m.y) * sc < 2.0)
                    ++ok;
            }
            if (ok >= 2 && 2 * ok > int(T.meas.size()) && T.meas.size() > nm) {
                best = int(t);
                nm = T.meas.size();
            }
        }
        int nin = 0;
        for (size_t k = 0; k < R0.frames.size(); ++k)
            nin += okf[k] && inside(O[i], k, 0);
        const double lpx = O[i].rate * F0.im.exptime / 60 / sc;
        std::string rem = best >= 0 && !R.movers[size_t(best)].meas.empty() ? R.movers[size_t(best)].meas[0].remark : "";
        std::printf("inj|%zu|%.3f|%.3f|%.1f|%.2f|%.2f|%d|%d|%zu|%s\n", i, O[i].mag, O[i].rate, O[i].pa, O[i].snr, lpx, nin, best >= 0, nm,
                    rem.c_str());
        if (best < 0)
            continue;
        hit[size_t(best)] = 1;
        ++nrec;
        const Track &T = R.movers[size_t(best)];
        for (auto &m : T.meas) {
            const Frame &f = R.frames[m.frame];
            double dt = (m.jd - tref) * 1440, ra, dec;
            std_to_radec(O[i].xi + O[i].vx * dt, O[i].eta + O[i].vy * dt, rc, dc, ra, dec);
            double dra = (m.ra - ra) * 3600 * std::cos(dec * PI / 180), dde = (m.dec - dec) * 3600;
            if (dra > 648000)
                dra -= 1296000;
            if (dra < -648000)
                dra += 1296000;
            double era, ede, cra, cde;
            measurement_sigma(R, T, m, s, era, ede, &cra, &cde);
            int di = m.frame < int(T.idx.size()) ? T.idx[m.frame] : -1;
            double tl = di >= 0 && di < int(f.det.size()) ? f.det[size_t(di)].trail_len : 0;
            std::printf("injm|%zu|%d|%+.4f|%+.4f|%.4f|%.4f|%.4f|%.4f|%.2f|%.1f|%s\n", i, m.frame, dra, dde, era, ede, cra, cde, tl,
                        m.snr_tot, m.remark.c_str());
        }
    }
    int other = 0;
    for (char h : hit)
        other += !h;
    std::printf("injsum|%d|%zu|%d|\n", nrec, O.size(), other);
    std::printf("injection: %d of %zu recovered; %d other track(s) (real or false objects)\n", nrec, O.size(), other);
    return 0;
}

static void write(const std::string &path, const std::string &text) {
    std::ofstream(path, std::ios::binary) << text;
}

int main(int argc, char **argv) {
    Settings s = Settings::defaults();
    std::string catalog, out = ".", object = "";
    double ra0 = NAN, dec0 = NAN;
    bool moving = true;
    InjectOpts io;
    std::vector<std::string> images;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "missing value for %s\n", a.c_str());
                std::exit(2);
            }
            return std::string(argv[++i]);
        };
        if (a == "--cfg")
            s.load(next());
        else if (a == "--catalog")
            catalog = next();
        else if (a == "--center") {
            ra0 = parse_angle(next(), true);
            dec0 = parse_angle(next(), false);
        } else if (a == "--center-deg") {
            ra0 = std::stod(next());
            dec0 = std::stod(next());
        } else if (a == "--object")
            object = next();
        else if (a == "--out")
            out = next();
        else if (a == "--set") {
            auto kv = next();
            auto p = kv.find('=');
            s.set(kv.substr(0, p), kv.substr(p + 1));
        } else if (a == "--inject")
            io.n = std::stoi(next());
        else if (a == "--inject-seed")
            io.seed = unsigned(std::stoul(next()));
        else if (a == "--inject-mag") {
            io.mlo = std::stod(next());
            io.mhi = std::stod(next());
        } else if (a == "--inject-rate") {
            io.vlo = std::stod(next());
            io.vhi = std::stod(next());
        } else if (a == "--inject-snr") {
            io.snrlo = std::stod(next());
            io.snrhi = std::stod(next());
        } else if (a == "--inject-color")
            io.color = std::stod(next());
        else if (a == "--inject-cross")
            io.cross = true;
        else if (a == "--inject-edge")
            io.edge = true;
        else if (a == "--inject-probe")
            io.probe = std::stoi(next());
        else if (a == "--no-moving")
            moving = false;
        else if (a == "--version") {
            std::puts("piazzi " PIAZZI_VERSION);
            return 0;
        } else if (a == "-h" || a == "--help") {
            std::puts("usage: piazzi --cfg F --catalog F.tsv --center RA DEC [--object DESIG] [--set K=V] [--out DIR] [--inject N "
                      "[--inject-seed S] "
                      "[--inject-mag LO HI] [--inject-rate LO HI] [--inject-snr LO HI] [--inject-color DMAG] [--inject-cross] "
                      "[--inject-edge]] images...");
            return 0;
        } else
            images.push_back(a);
    }
    if (images.empty() || std::isnan(ra0)) {
        std::fprintf(stderr, "usage: piazzi --cfg F --catalog F.tsv --center RA DEC images...\n");
        return 2;
    }

    std::vector<Image> ims;
    for (auto &p : images)
        ims.push_back(load_image(p, s.num("flip_horizontal") != 0, s.num("flip_vertical") != 0, s.num("time_reference") == 0));
    double sx = s.num("pixel_width_mm") / s.num("focal_length_mm") * 206264.8,
           sy = s.num("pixel_height_mm") / s.num("focal_length_mm") * 206264.8;
    double boxw = ims[0].w * sx / 60 + 2 * s.num("pointing_tolerance");
    if (catalog.empty()) {
        std::printf("download the catalog and pass it with --catalog:\n  curl -o gaia.tsv '%s'\n",
                    gaia_box_url(ra0, dec0, boxw, s.num("catalog_faint_mag")).c_str());
        return 1;
    }
    double t = ims[0].jd_mid;
    for (auto &im : ims)
        t = std::min(t, im.jd_mid);
    double epoch = 2000.0 + (t - 2451545.0) / 365.25;
    Catalog cat = load_vizier_tsv(catalog, epoch, s.num("catalog_faint_mag"), s.str("plate_ruwe").empty() ? 1.4 : s.num("plate_ruwe"));
    std::printf("catalog: %zu stars\n", cat.size());

    if (io.n > 0)
        return inject_run(ims, cat, ra0, dec0, s, io);
    Reduction R = reduce(ims, cat, ra0, dec0, s, moving);
    // the "image"/"moving object" lines below are parsed by tools/benchmark/run.py
    for (size_t k = 0; k < R.frames.size(); ++k) {
        const Frame &f = R.frames[k];
        const Plate &p = f.plate;
        std::printf("image %zu: %zu detections | %d of %d refs, dRA %.2f\" dDe %.2f\" | photometry %d, ZP %.3f | seeing %.1f\"\n", k + 1,
                    f.det.size(), p.nref, p.nmatch, p.dra_mean, p.dde_mean, p.nphot, p.zp, f.seeing);
    }
    // filters match log prefixes emitted by the reduction code (cpp/src); keep them in sync
    for (auto &l : R.log)
        if (l.find("stack") != std::string::npos || l.find("2-detection") != std::string::npos || s.num("plate_debug") != 0 ||
            (s.num("trail_debug") != 0 && l.rfind("trail", 0) == 0) || l.rfind("link probe", 0) == 0)
            std::printf("  · %s\n", l.c_str());
    for (auto &T : R.movers)
        std::printf("moving object: %.2f\"/min PA %.1f°\n", T.rate, T.pa);
    std::vector<Track> tr = R.movers;
    for (auto &T : tr)
        T.desig = object; // empty: no provID (ADES carries only the trkSub)
    if (!tr.empty()) {
        write(out + "/ADESReport.psv", ades_report(R, tr, s));
        write(out + "/PhotReport.txt", phot_report(R, tr, s));
        write(out + "/MPCReport.txt", mpc_report(R, tr, s));
        std::printf("reports in %s\n", out.c_str());
    }
    return 0;
}
