// Track & Stack (extras.track_and_stack / session.stack) and the plate's Image Catalog.
#include "piazzi/astro_c.h"
#include "piazzi/reduce.hpp"
#include "session_internal.hpp"
#include <exception>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>

using namespace astro;

using namespace astro;
using namespace session_util;

extern "C" {

int astro_stack(AstroSession *s, double rate, double pa, int mode) {
    try {
        if (!s || !s->ok)
            return -1;
        std::vector<int> use;
        for (int k = 0; k < int(s->R.frames.size()); ++k)
            if (!s->R.frames[k].stacked && s->R.frames[k].solved)
                use.push_back(k);
        if (use.empty())
            return -1;
        const Frame &f0 = s->R.frames[use[0]];
        const int W = f0.im.w, H = f0.im.h;
        const size_t N = size_t(W) * H;
        double cra, cde;
        f0.plate.pix_to_radec(W / 2.0, H / 2.0, cra, cde);
        double x0, y0;
        sky2pix(f0.plate, cra, cde, x0, y0);
        const float NaN = std::numeric_limits<float>::quiet_NaN();
        std::vector<std::vector<float>> layers;
        double exps = 0;
        for (int k : use) {
            const Frame &f = s->R.frames[k];
            exps += f.im.exptime;
            double dt = (f.im.jd_mid - f0.im.jd_mid) * 1440.0;
            double dxi = rate * dt * std::sin(pa * PI / 180), deta = rate * dt * std::cos(pa * PI / 180), rt, dtt, xt, yt;
            std_to_radec(dxi, deta, cra, cde, rt, dtt);
            sky2pix(f.plate, rt, dtt, xt, yt);
            double sx = x0 - xt, sy = y0 - yt; // ndimage.shift(order=1, mode='constant', cval=nan)
            std::vector<float> L(N, NaN);
            for (int y = 0; y < H; ++y)
                for (int x = 0; x < W; ++x) {
                    double u = x - sx, v = y - sy;
                    if (u < 0 || v < 0 || u > f.im.w - 1 || v > f.im.h - 1)
                        continue;
                    int iu = std::min(int(u), f.im.w - 2), iv = std::min(int(v), f.im.h - 2);
                    double fu = u - iu, fv = v - iv;
                    if (f.im.w < 2 || f.im.h < 2) {
                        iu = iv = 0;
                        fu = fv = 0;
                    }
                    auto P = [&](int a, int b) { return double(f.im.at(a, b)); };
                    L[size_t(y) * W + x] = float((1 - fu) * (1 - fv) * P(iu, iv) + fu * (1 - fv) * P(iu + 1, iv) +
                                                 (1 - fu) * fv * P(iu, iv + 1) + fu * fv * P(iu + 1, iv + 1));
                }
            layers.push_back(std::move(L));
        }
        std::vector<double> ped(layers.size(), 0);
        if (mode == 0)
            for (size_t l = 0; l < layers.size(); ++l) {
                std::vector<float> t = layers[l];
                ped[l] = nanmedian(t);
            }
        double pmed = 0;
        if (mode == 0) {
            std::vector<float> t(ped.begin(), ped.end());
            pmed = nanmedian(t);
        }
        std::vector<float> out(N, NaN), col;
        col.reserve(layers.size());
        for (size_t i = 0; i < N; ++i) {
            col.clear();
            double sum = 0;
            int n = 0;
            for (size_t l = 0; l < layers.size(); ++l) {
                float q = layers[l][i];
                if (std::isnan(q))
                    continue;
                col.push_back(q);
                sum += q - (mode == 0 ? ped[l] : 0);
                ++n;
            }
            if (mode == 2)
                out[i] = float(nanmedian(col));
            else if (mode == 0)
                out[i] = float(sum + pmed); // nansum(cube - ped) + median(ped): zero valid samples -> pmed only
            else if (n)
                out[i] = float(sum / n);
        }
        {
            std::vector<float> t = out;
            double fill = nanmedian(t);
            if (std::isnan(fill))
                fill = 0;
            for (auto &q : out)
                if (std::isnan(q))
                    q = float(fill);
        }
        Frame st;
        st.im.w = W;
        st.im.h = H;
        st.im.pix.resize(N);
        for (size_t i = 0; i < N; ++i) {
            double r = std::nearbyint(double(out[i]));
            st.im.pix[i] = uint16_t(r > 0 ? std::min(r, 65535.0) : 0.0);
        }
        st.im.jd_mid = f0.im.jd_mid;
        st.im.exptime = exps;
        st.im.path = "stack(" + std::to_string(use.size()) + ")";
        st.det = detect_sources(st.im, detect_params(s->s));
        std::stable_sort(st.det.begin(), st.det.end(), [](const Detection &a, const Detection &b) { return a.flux > b.flux; });
        st.plate = f0.plate;
        st.plate.ref_idx.clear();
        st.plate.ref_ok.clear();
        st.plate.cat_idx.clear(); // indices referred to the detections of f0
        if (mode == 0)
            st.plate.zp = f0.plate.zp + 2.5 * std::log10(double(use.size()));
        st.solved = true;
        st.stacked = true;
        st.offx = f0.offx;
        st.offy = f0.offy;
        std::vector<double> fw;
        for (auto &d : st.det)
            fw.push_back(d.fwhm);
        if (!fw.empty()) {
            std::nth_element(fw.begin(), fw.begin() + fw.size() / 2, fw.end());
            st.seeing = fw[fw.size() / 2] * (st.plate.scale_x() + st.plate.scale_y()) / 2;
        } else
            st.seeing = f0.seeing;
        s->R.frames.push_back(std::move(st));
        session_rebuild_json(s);
        return int(s->R.frames.size()) - 1;
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return -1;
    } catch (...) {
        return -1;
    }
}

const char *astro_image_catalog(AstroSession *s, int k) {
    try {
        if (!s)
            return "[]";
        s->buf3 = "[]";
        if (!s->ok || k < 0 || k >= int(s->R.frames.size()))
            return s->buf3.c_str();
        const Frame &f = s->R.frames[k];
        const Plate &p = f.plate;
        std::vector<char> ty(f.det.size(), '?');
        for (size_t i = 0; i < f.typ.size() && i < ty.size(); ++i)
            if (f.typ[i] == 1 || f.typ[i] == 2)
                ty[i] = 'S';
            else if (f.typ[i] == 3) // point fragment of a trail (deferred search)
                ty[i] = 'T';
        for (size_t i = 0; i < ty.size(); ++i)
            if (f.det[i].trail_len > 0 && ty[i] != 'r' && ty[i] != 'R') // attached trail (elongated source)
                ty[i] = 'L';
        for (size_t i = 0; i < p.ref_idx.size(); ++i)
            if (p.ref_idx[i] < int(ty.size()))
                ty[p.ref_idx[i]] = p.ref_ok[i] ? 'R' : 'r';
        if (!f.stacked)
            for (auto &T : s->R.movers)
                if (k < int(T.idx.size()) && T.idx[k] >= 0 && T.idx[k] < int(ty.size()))
                    ty[T.idx[k]] = 'M';
        double sc = f.solved ? plate_scale(p) : 0;
        std::ostringstream j;
        j.precision(10);
        j << "[";
        for (size_t i = 0; i < f.det.size(); ++i) {
            const Detection &d = f.det[i];
            double ra = 0, dec = 0;
            if (f.solved)
                p.pix_to_radec(d.x, d.y, ra, dec);
            double mag = (f.solved && d.flux > 0) ? p.zp - 2.5 * std::log10(d.flux) : 0;
            j << (i ? "," : "") << "{\"x\":" << d.x << ",\"y\":" << d.y << ",\"ra\":" << ra << ",\"dec\":" << dec << ",\"mag\":" << mag
              << ",\"snr\":" << d.snr << ",\"fwhm\":" << d.fwhm * sc << ",\"fit\":" << d.fit << ",\"flux\":" << d.flux << ",\"type\":\""
              << ty[i] << "\"}";
        }
        j << "]";
        s->buf3 = j.str();
        return s->buf3.c_str();
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return "[]";
    } catch (...) {
        return "[]";
    }
}

} // extern "C"
