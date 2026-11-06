// Home and session-progress screens (observing log).
#include "app.hpp"
#include "imgui.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <mutex>
#include <string>

static void home_card(const char *n, const char *tg, const char *ttl, const char *desc, bool primary, const std::function<void()> &act) {
    ImGui::PushID(n);
    ImVec2 p = ImGui::GetCursorScreenPos();
    float w = ImGui::GetContentRegionAvail().x, h = 104;
    bool clicked = ImGui::InvisibleButton("##c", ImVec2(w, h));
    bool hov = ImGui::IsItemHovered();
    ImDrawList *dl = ImGui::GetWindowDrawList();
    dl->AddRectFilled(p, ImVec2(p.x + w, p.y + h), T::u(T::panel, 0.92f), 8);
    dl->AddRect(p, ImVec2(p.x + w, p.y + h), primary || hov ? T::u(T::amber, 0.45f) : T::u(T::line), 8);
    dl->AddRectFilled(ImVec2(p.x + 16, p.y + 16), ImVec2(p.x + 50, p.y + 50), primary ? T::u(T::amber) : T::u(T::amberDim), 4);
    ImGui::PushFont(T::mono, 14.0f);
    dl->AddText(ImVec2(p.x + 23, p.y + 24), primary ? T::u(T::ink) : T::u(T::amber), n);
    ImGui::PopFont();
    std::string tu = upper_spaced(tg, false) + ' ';
    ImGui::PushFont(T::mono, 11.0f);
    dl->AddText(ImVec2(p.x + 66, p.y + 16), T::u(T::dim), tu.c_str());
    ImGui::PopFont();
    ImGui::PushFont(T::serif, 22.0f);
    dl->AddText(ImVec2(p.x + 66, p.y + 34), T::u(T::text), ttl);
    ImGui::PopFont();
    ImGui::PushFont(T::sans, 14.5f);
    dl->AddText(nullptr, 0, ImVec2(p.x + 66, p.y + 66), T::u(T::dim), desc, nullptr, w - 110);
    ImGui::PopFont();
    ImGui::PushFont(T::sans, 16.0f);
    dl->AddText(ImVec2(p.x + w - 30, p.y + 16), primary ? T::u(T::amber) : T::u(T::faint), "↗");
    ImGui::PopFont();
    if (clicked)
        act();
    ImGui::Dummy(ImVec2(0, 4));
    ImGui::PopID();
}

void ui_home() {
    ImVec2 avail = ImGui::GetContentRegionAvail();
    float W = std::min(620.0f, avail.x);
    ImGui::SetCursorPosX((ImGui::GetWindowWidth() - W) / 2);
    ImGui::BeginChild("##home", ImVec2(W, 0), 0, ImGuiWindowFlags_NoScrollbar);
    tag(tr("Pocket observatory"));
    ImGui::SameLine(W - 260);
    ImGui::PushStyleColor(ImGuiCol_Text, T::amber);
    ImGui::PushFont(T::mono, 12.0f);
    if (ImGui::Button((std::string(tr(A.profile.c_str())) + tr("  ·  settings")).c_str(), ImVec2(260, 0)))
        A.show_settings = true;
    ImGui::PopFont();
    ImGui::PopStyleColor();
    ImGui::Dummy(ImVec2(0, 14));
    title("Piazzi", 50);
    ImGui::PushFont(T::serif, 20.0f);
    ImGui::TextColored(T::dim, "%s", tr("Measure the sky. Find what moves."));
    ImGui::PopFont();
    ImGui::Dummy(ImVec2(0, 10));
    // star chart with the "object" blinking between two positions
    ImVec2 p = ImGui::GetCursorScreenPos();
    ImDrawList *dl = ImGui::GetWindowDrawList();
    dl->AddRectFilled(p, ImVec2(p.x + W, p.y + 150), IM_COL32(5, 6, 8, 255), 4);
    star_chart(dl, p, ImVec2(W, 150), 11);
    bool ph = std::fmod(ImGui::GetTime(), 1.4) < 0.7;
    ImVec2 c(p.x + W * (ph ? .58f : .55f), p.y + 150 * (ph ? .43f : .47f));
    dl->AddLine(ImVec2(p.x + W * .49f, p.y + 82), ImVec2(p.x + W * .64f, p.y + 52), T::u(T::amber, 0.4f));
    reticle(dl, c, 15, T::u(T::amber));
    dl->AddCircleFilled(c, 1.6f, IM_COL32_WHITE);
    dl->AddRect(p, ImVec2(p.x + W, p.y + 150), T::u(T::line), 4);
    ImGui::PushFont(T::mono, 11.0f);
    dl->AddText(ImVec2(p.x + 10, p.y + 130), T::u(T::amber, 0.8f), "0 7  1 1  4 4 . 2 8  + 2 2  2 4  1 6");
    dl->AddText(ImVec2(p.x + W - 110, p.y + 130), T::u(T::amber, 0.8f), "0 . 3 5 \" / M I N");
    ImGui::PopFont();
    ImGui::Dummy(ImVec2(W, 164));
    home_card("01", tr("Your plates"), tr("Open my images"), tr("FITS or SBIG from your telescope. Or drag the files onto the window."),
              true, [] { BR.show(tr("Choose the images (and, optionally, the .cfg and the catalog)"), true, start_with_files); });
    home_card("02", tr("Example"), tr("Open example"), tr("Asteroid (5157) Hindemith in three ZTF images. Runs offline."), false,
              [] { start_example(); });
    ImGui::Dummy(ImVec2(0, 10));
    tag("open source   ·   Gaia DR3   ·   SkyBoT   ·   ADES / MPC", T::faint);
    ImGui::Dummy(ImVec2(0, 6));
    ImGui::PushFont(T::mono, 11.0f);
    ImGui::TextColored(T::faint, "Piazzi %s · open source · MIT", astro_version());
    ImGui::PopFont();
    ImGui::EndChild();
}

void ui_session() {
    float W = std::min(620.0f, ImGui::GetContentRegionAvail().x);
    ImGui::SetCursorPosX((ImGui::GetWindowWidth() - W) / 2);
    ImGui::BeginChild("##sess", ImVec2(W, 0));
    // the lock is released before go_home/rerun_pipeline: they wait for the worker thread, which also uses A.mu
    std::unique_lock<std::recursive_mutex> l(A.mu);
    int done = 0;
    bool fail = false;
    for (auto &s : A.steps) {
        done += s.st == 2;
        fail |= s.st == 3;
    }
    if (ImGui::Button("x", ImVec2(36, 36))) {
        l.unlock();
        go_home();
        ImGui::EndChild();
        return;
    }
    ImGui::SameLine(W - 110);
    tag(tr("Step ") + std::to_string(std::min(done + 1, 5)) + "/5", T::amber);
    ImVec2 p = ImGui::GetCursorScreenPos();
    ImDrawList *dl = ImGui::GetWindowDrawList();
    dl->AddRectFilled(p, ImVec2(p.x + W, p.y + 2), T::u(T::line));
    dl->AddRectFilled(p, ImVec2(p.x + W * done / 5.0f, p.y + 2), fail ? T::u(T::red) : T::u(T::amber));
    ImGui::Dummy(ImVec2(0, 20));
    title(fail ? tr("Needs attention.") : done == 5 ? tr("Done.") : tr("Reducing the plates…"), 34);
    ImGui::PushFont(T::sans, 15.0f);
    ImGui::TextColored(T::dim, "%s", A.status.c_str());
    ImGui::PopFont();
    ImGui::Dummy(ImVec2(0, 10));
    ImGui::PushStyleColor(ImGuiCol_ChildBg, T::panel);
    ImGui::BeginChild("##log", ImVec2(0, 300), ImGuiChildFlags_Borders);
    ImGui::Dummy(ImVec2(0, 4));
    tag(tr("Observing log"));
    for (auto &s : A.steps) {
        ImGui::PushFont(T::mono, 12.0f);
        ImGui::TextColored(T::faint, "%s", s.at.empty() ? "--:--" : s.at.c_str());
        ImGui::PopFont();
        ImGui::SameLine(70);
        const char *g = s.st == 0 ? "○" : s.st == 1 ? "◐" : s.st == 2 ? "●" : "✕";
        ImGui::PushFont(T::mono, 14.0f);
        ImGui::TextColored(s.st == 3 ? T::red : s.st == 0 ? T::faint : T::amber, "%s", g);
        ImGui::SameLine(96);
        ImGui::TextColored(s.st == 0 ? T::faint : T::text, "%s", tr(s.title.c_str()));
        ImGui::PopFont();
        if (!s.detail.empty()) {
            ImGui::SetCursorPosX(96);
            ImGui::PushFont(T::mono, 12.5f);
            ImGui::TextColored(s.st == 3 ? T::red : T::dim, "%s", s.detail.c_str());
            ImGui::PopFont();
        }
    }
    ImGui::EndChild();
    ImGui::PopStyleColor();
    ImGui::Dummy(ImVec2(0, 14));
    bool retry = false, home = false;
    if (fail) {
        // the worker thread may still be finishing (busy): accept "Try again" only once it has stopped
        bool busy = A.busy;
        ImGui::BeginDisabled(busy);
        retry = primary_button(tr("TRY AGAIN")) && !busy;
        ImGui::EndDisabled();
        home = ghost_button(tr("BACK TO START"));
    } else
        home = ghost_button(tr("CANCEL"));
    l.unlock();
    if (retry)
        rerun_pipeline();
    else if (home)
        go_home();
    ImGui::EndChild();
}
