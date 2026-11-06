// Results panel tools: Track & Stack, contrast, detection list, header/timing, save FITS.
#include "app.hpp"
#include "imgui.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <filesystem>
#include <set>
#include <string>
#include <vector>

static void tool_stack(const astro::Reduction &R) {
    static float rate = -1, pa = 0;
    static int mode = 1;
    if (rate < 0) {
        rate = R.movers.empty() ? 1.0f : (float)R.movers[0].rate;
        pa = R.movers.empty() ? 0.0f : (float)R.movers[0].pa;
    }
    ImGui::InputFloat("″/min", &rate, 0.1f, 1.0f, "%.2f");
    ImGui::InputFloat("AP °", &pa, 1, 10, "%.1f");
    const char *modes[] = {tr("sum"), tr("mean"), tr("median")};
    ImGui::Combo(tr("mode"), &mode, modes, 3);
    if (ImGui::Button(tr("Stack"), ImVec2(-1, 0))) {
        int f = astro_stack(A.s, rate, pa, mode);
        if (f >= 0) {
            A.tex_dirty = true;
            A.k = f;
            A.blink = false;
            status(tr("Stacked: plate ") + std::to_string(f + 1) + ".");
        }
    }
    ImGui::TextColored(T::faint, "%s", tr("Measurements on the stacked plate get note K."));
}

static void tool_contrast() {
    static std::vector<float> hist;
    static int hk = -1;
    if (hk != A.k) {
        std::vector<double> c(128);
        double mm[2];
        if (astro_histogram(A.s, A.k, 128, c.data(), mm)) {
            hist.assign(c.begin(), c.end());
            for (auto &v : hist)
                v = std::log10(1 + v);
        }
        hk = A.k;
    }
    if (!hist.empty())
        ImGui::PlotHistogram("##h", hist.data(), (int)hist.size(), 0, nullptr, 0, FLT_MAX, ImVec2(-1, 70));
    bool ch = ImGui::SliderFloat(tr("black %"), &A.lo_pct, 0, 99, "%.1f");
    ch |= ImGui::SliderFloat(tr("white %"), &A.hi_pct, 50, 100, "%.2f");
    const char *str[] = {"linear", tr("square root"), "log"};
    ch |= ImGui::Combo(tr("stretch"), &A.stretch, str, 3);
    ch |= ImGui::Checkbox(tr("invert"), &A.invert);
    if (ch)
        A.tex_dirty = true;
}

static void tool_detections(const astro::Reduction &R) {
    const astro::Frame &f = R.frames[A.k];
    std::set<int> refs(f.plate.ref_idx.begin(), f.plate.ref_idx.end());
    ImGui::TextDisabled("%s", tr("R reference · S star · T trail fragment · L trail · ? loose"));
    if (!ImGui::BeginTable("##cat", 6, ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY | ImGuiTableFlags_BordersInnerH, ImVec2(0, 260)))
        return;
    for (const char *h : {"T", "x", "y", "mag", "SNR", "FWHM\""})
        ImGui::TableSetupColumn(h);
    ImGui::TableHeadersRow();
    double sc = std::sqrt(std::fabs(f.plate.cx[1] * f.plate.cy[2] - f.plate.cx[2] * f.plate.cy[1])); // "/px
    for (size_t i = 0; i < f.det.size(); ++i) {
        auto &d = f.det[i];
        bool ref = refs.count((int)i) > 0;
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        int ty = i < f.typ.size() ? f.typ[i] : 0;
        const char *code = ref ? "R" : d.trail_len > 0 ? "L" : ty == 3 ? "T" : (ty == 1 || ty == 2) ? "S" : "?";
        ImGui::TextColored(ref ? T::green : T::dim, "%s", code);
        ImGui::TableNextColumn();
        ImGui::Text("%.2f", d.x);
        ImGui::TableNextColumn();
        ImGui::Text("%.2f", d.y);
        ImGui::TableNextColumn();
        ImGui::Text("%.2f", f.plate.zp - 2.5 * std::log10(std::max(d.flux, 1e-9)));
        ImGui::TableNextColumn();
        ImGui::Text("%.1f", d.snr);
        ImGui::TableNextColumn();
        ImGui::Text("%.1f", d.fwhm * sc);
    }
    ImGui::EndTable();
}

static void tool_header() {
    static double off[64];
    int n = astro_image_count(A.s);
    for (int i = 0; i < n; ++i) {
        ImGui::PushID(i);
        double jd = astro_image_jd(A.s, i);
        ImGui::TextColored(T::dim, "%s  ·  %s", fs::path(A.images[i]).filename().string().c_str(), astro::jd_to_isot(jd).c_str());
        if (ImGui::TreeNode(tr("header"))) {
            ImGui::TextUnformatted(astro_image_header(A.s, i));
            ImGui::TreePop();
        }
        ImGui::SetNextItemWidth(160);
        ImGui::InputDouble(tr("correct (s)"), &off[i % 64], 1, 10, "%.1f");
        ImGui::SameLine();
        if (ImGui::SmallButton(tr("apply"))) {
            astro_set_image_jd(A.s, i, jd + off[i % 64] / 86400);
            off[i % 64] = 0;
            status(tr("Time changed — reduce again to apply."));
        }
        ImGui::PopID();
    }
    if (ImGui::Button(tr("Reduce again"), ImVec2(-1, 0)))
        rerun_pipeline();
}

static void tool_save_fits() {
    if (!ImGui::Button(tr("Save FITS with WCS (current plate)"), ImVec2(-1, 0)))
        return;
    std::string src = A.images.empty() ? tr("image") : A.images[std::min<size_t>(A.k, A.images.size() - 1)];
    std::string p = output_dir() + fs::path(src).stem().string() + "_wcs.fits";
    status(astro_save_fits(A.s, A.k, p.c_str()) ? tr("FITS saved to ") + p : tr("Could not save the FITS."));
}

void sheet_tools() {
    auto &R = A.s->R;
    ImGui::PushFont(T::mono, 13.0f);
    if (ImGui::CollapsingHeader("Track & Stack"))
        tool_stack(R);
    if (ImGui::CollapsingHeader(tr("Contrast")))
        tool_contrast();
    if (ImGui::CollapsingHeader(tr("Detection list (Image Catalog)")))
        tool_detections(R);
    if (ImGui::CollapsingHeader(tr("Header and time")))
        tool_header();
    if (ImGui::CollapsingHeader(tr("Save image")))
        tool_save_fits();
    ImGui::PopFont();
}
