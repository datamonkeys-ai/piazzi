// Results screen: plate viewer (blink, zoom, layers, click to identify/measure) + side panel.
#include "app.hpp"
#include "imgui.h"
#include <mutex>
#include <string>
#include <vector>
#ifdef __EMSCRIPTEN__
#include <GLES3/gl3.h>
#else
#include <SDL_opengl.h>
#endif
#include <algorithm>
#include <cmath>

void make_textures() {
    for (GLuint t : A.tex)
        glDeleteTextures(1, &t);
    A.tex.clear();
    auto &R = A.s->R;
    for (size_t k = 0; k < R.frames.size(); ++k) {
        int w = R.frames[k].im.w, h = R.frames[k].im.h;
        std::vector<unsigned char> rgba(size_t(w) * h * 4);
        astro_render_levels(A.s, (int)k, A.lo_pct, A.hi_pct, A.stretch, A.invert, rgba.data());
        GLuint id;
        glGenTextures(1, &id);
        glBindTexture(GL_TEXTURE_2D, id);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgba.data());
        A.tex.push_back(id);
    }
    A.tex_dirty = false;
}

// top bar: back, plates, blink, click mode, layers. False if the user went back to the home screen.
static bool viewer_toolbar(const astro::Reduction &R) {
    if (ImGui::Button("<", ImVec2(36, 36))) {
        go_home();
        return false;
    }
    ImGui::SameLine();
    ImGui::PushFont(T::mono, 13.0f);
    for (size_t i = 0; i < R.frames.size(); ++i) {
        ImGui::SameLine();
        bool on = (int)i == A.k;
        if (on) {
            ImGui::PushStyleColor(ImGuiCol_Button, T::amber);
            ImGui::PushStyleColor(ImGuiCol_Text, T::ink);
        }
        if (ImGui::Button((std::to_string(i + 1) + (R.frames[i].stacked ? "Σ" : "")).c_str(), ImVec2(34, 36))) {
            A.k = (int)i;
            A.blink = false;
        }
        if (on)
            ImGui::PopStyleColor(2);
    }
    ImGui::SameLine();
    if (ImGui::Button(A.blink ? tr("pause") : tr("blink"), ImVec2(80, 36)))
        A.blink = !A.blink;
    ImGui::SameLine();
    ImGui::SetNextItemWidth(200);
    {
        const char *modes[] = {tr("click: identify"), tr("click: measure")};
        ImGui::Combo("##mode", &A.click_mode, modes, 2);
    }
    ImGui::SameLine();
    if (ImGui::Button(tr("layers")))
        ImGui::OpenPopup("##lay");
    if (ImGui::BeginPopup("##lay")) {
        ImGui::Checkbox(tr("Object track"), &A.lay_track);
        ImGui::Checkbox(tr("Known asteroids"), &A.lay_known);
        ImGui::Checkbox(tr("Reference stars"), &A.lay_refs);
        ImGui::Checkbox(tr("All detections"), &A.lay_dets);
        ImGui::EndPopup();
    }
    ImGui::PopFont();
    return true;
}

// overlays on the image; P converts a plate pixel to a screen coordinate
template <class F> static void draw_overlays(ImDrawList *dl, const astro::Reduction &R, const astro::Frame &f, float s, F P) {
    if (A.lay_dets)
        for (auto &d : f.det)
            dl->AddCircle(P(d.x, d.y), 3 * s, T::u(T::cyan, 0.35f));
    if (A.lay_refs)
        for (size_t i = 0; i < f.plate.ref_idx.size(); ++i) {
            auto &d = f.det[f.plate.ref_idx[i]];
            reticle(dl, P(d.x, d.y), 5 * s, f.plate.ref_ok[i] ? T::u(T::green, 0.7f) : T::u(T::red, 0.7f), 1, 0.5f);
        }
    if (A.lay_known)
        for (auto &a : A.known) {
            bool same = false;
            for (auto &[i, lab] : A.identified)
                if (lab.find(a.name) != std::string::npos)
                    same = true;
            if (same)
                continue;
            double x, y;
            if (!astro_sky_to_pix(A.s, A.k, a.ra, a.dec, &x, &y))
                continue;
            ImVec2 c = P(x, y);
            reticle(dl, c, 7 * s, T::u(T::cyan, 0.7f), 1, 0.35f);
            ImGui::PushFont(T::mono, 10.0f);
            dl->AddText(ImVec2(c.x - ImGui::CalcTextSize(a.name.c_str()).x / 2, c.y - 12 * s - 10), T::u(T::cyan, 0.85f), a.name.c_str());
            ImGui::PopFont();
        }
    for (size_t i = 0; i < R.movers.size(); ++i) {
        auto &m = R.movers[i];
        if (A.lay_track && m.meas.size() > 1)
            for (size_t j = 1; j < m.meas.size(); ++j)
                dl->AddLine(P(m.meas[j - 1].x, m.meas[j - 1].y), P(m.meas[j].x, m.meas[j].y), T::u(T::amber, 0.5f), 1.2f);
        for (auto &q : m.meas)
            if (q.frame == A.k) {
                ImVec2 c = P(q.x, q.y);
                reticle(dl, c, 10 * s, T::u(T::amber), 1.8f);
                std::string lab = A.identified.count((int)i) && !A.identified[(int)i].empty()
                                      ? A.identified[(int)i].substr(0, A.identified[(int)i].find(" ·"))
                                      : tr("new object?");
                ImGui::PushFont(T::mono, 12.0f);
                dl->AddText(ImVec2(c.x - ImGui::CalcTextSize(lab.c_str()).x / 2, c.y + 14 * s), T::u(T::amber), lab.c_str());
                ImGui::PopFont();
            }
    }
    {
        std::lock_guard<std::recursive_mutex> lk(A.mu); // A.insp is written by the "What is this?" thread
        if (A.insp.on)
            reticle(dl, P(A.insp.x, A.insp.y), 8 * s, IM_COL32_WHITE, 1.4f, 0.6f);
    }
    if (A.meas.on && A.meas.frame == A.k)
        reticle(dl, P(A.meas.v[5], A.meas.v[6]), 8 * s, T::u(T::green), 1.4f, 0.6f);
}

// zoomable image, coordinate tooltip and click (identify the nearest detection or measure)
static void viewer_image(const astro::Reduction &R, const astro::Frame &f) {
    ImGui::BeginChild("##img", ImVec2(0, 0), 0, ImGuiWindowFlags_HorizontalScrollbar);
    ImVec2 area = ImGui::GetContentRegionAvail();
    float base = std::min(area.x / f.im.w, area.y / f.im.h), s = base * A.zoom;
    const ImGuiIO &io = ImGui::GetIO();
    if (ImGui::IsWindowHovered() && io.MouseWheel != 0 && (io.KeyCtrl || io.KeySuper))
        A.zoom = std::clamp(A.zoom * (io.MouseWheel > 0 ? 1.15f : 0.87f), 1.0f, 8.0f);
    ImVec2 o = ImGui::GetCursorScreenPos();
    if (A.k < (int)A.tex.size())
        ImGui::Image((ImTextureID)(intptr_t)A.tex[A.k], ImVec2(f.im.w * s, f.im.h * s));
    bool hovered = ImGui::IsItemHovered();
    auto P = [&](double x, double y) { return ImVec2(o.x + float((x + 0.5) * s), o.y + float((y + 0.5) * s)); };
    draw_overlays(ImGui::GetWindowDrawList(), R, f, s, P);
    if (hovered) {
        ImVec2 mp = ImGui::GetMousePos();
        double x = (mp.x - o.x) / s - 0.5, y = (mp.y - o.y) / s - 0.5, ra, dec;
        if (astro_pix_to_sky(A.s, A.k, x, y, &ra, &dec)) {
            ImGui::PushFont(T::mono, 12.0f);
            ImGui::SetTooltip("%s  %s\nx %.1f  y %.1f  ·  ADU %d", fmt_hms(ra / 15, false).c_str(), fmt_hms(dec, true).c_str(), x, y,
                              f.im.at(int(x + .5), int(y + .5)));
            ImGui::PopFont();
        }
        // the 2nd click of a double click only toggles zoom (does not redo the query/measurement)
        if (ImGui::IsMouseClicked(0) && ImGui::GetIO().MouseClickedCount[0] == 1) {
            if (A.click_mode == 0) {
                double bx = x, by = y, bd = 8;
                for (auto &d : f.det) {
                    double dd = std::hypot(d.x - x, d.y - y);
                    if (dd < bd) {
                        bd = dd;
                        bx = d.x;
                        by = d.y;
                    }
                }
                inspect_at(bx, by);
            } else {
                A.meas.on = astro_measure(A.s, A.k, x, y, A.meas.v) != 0;
                A.meas.frame = A.k;
                A.blink = false;
            }
        }
        if (ImGui::IsMouseDoubleClicked(0))
            A.zoom = A.zoom > 1 ? 1 : 2.5f;
    }
    ImGui::EndChild();
}

void ui_results() {
    auto &R = A.s->R;
    if (R.frames.empty()) {
        go_home(); // called without A.mu: go_home waits for the worker thread, which also uses the lock
        return;
    }
    if (A.tex_dirty)
        make_textures();
    A.k = std::clamp(A.k, 0, (int)R.frames.size() - 1);
    if (A.blink && ImGui::GetTime() - A.last_blink > 0.6) {
        int nb = 0;
        for (auto &fr : R.frames)
            nb += !fr.stacked;
        A.k = (A.k + 1) % std::max(1, nb);
        A.last_blink = ImGui::GetTime();
    }
    float panelW = 430;
    ImVec2 avail = ImGui::GetContentRegionAvail();
    ImGui::BeginChild("##viewer", ImVec2(avail.x - panelW - 12, 0), 0, ImGuiWindowFlags_NoScrollbar);
    if (!viewer_toolbar(R)) {
        ImGui::EndChild();
        return;
    }
    const astro::Frame &f = R.frames[A.k];
    tag(tr("Plate ") + std::to_string(A.k + 1) + " · " + astro::jd_to_isot(f.im.jd_mid).substr(11), T::amber);
    viewer_image(R, f);
    ImGui::EndChild();
    ImGui::SameLine();
    ImGui::BeginChild("##sheet", ImVec2(panelW, 0));
    sheet();
    ImGui::EndChild();
}
