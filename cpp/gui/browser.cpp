// Portable file browser (ImGui modal): picks files or a folder.
#include "app.hpp"
#include "imgui.h"
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <string>
#include <utility>
#include <vector>

Browser BR;

void Browser::show(const std::string &t, bool m, std::function<void(std::vector<std::string>)> cb, bool d) {
    title = t;
    multi = m;
    done = std::move(cb);
    dirs_only = d;
    sel.clear();
    open = true;
    if (cwd.empty()) {
        const char *h = std::getenv("HOME");
        if (!h)
            h = std::getenv("USERPROFILE");
        cwd = h ? h : ".";
    }
    pending_open = true; // popup ID depends on the current window: open it in draw(), same window as BeginPopupModal
}

void Browser::draw() {
    if (pending_open) {
        ImGui::OpenPopup("##browser");
        pending_open = false;
    }
    ImGui::SetNextWindowSize(ImVec2(720, 520), ImGuiCond_Appearing);
    if (!ImGui::BeginPopupModal("##browser", nullptr, ImGuiWindowFlags_NoTitleBar))
        return;
    tag(title, T::amber);
    if (ImGui::Button(tr("↑ up"))) {
        cwd = fs::path(cwd).parent_path().string();
        sel.clear();
    }
    ImGui::SameLine();
    ImGui::PushFont(T::mono, 13.0f);
    ImGui::TextWrapped("%s", cwd.c_str());
    ImGui::PopFont();
    // two-line footer (hint, then buttons): reserve the height of both
    ImGui::BeginChild("##list", ImVec2(0, -(54 + ImGui::GetTextLineHeightWithSpacing() + 6)), ImGuiChildFlags_Borders);
    std::error_code ec;
    std::vector<fs::directory_entry> es;
    for (auto &e : fs::directory_iterator(cwd, fs::directory_options::skip_permission_denied, ec)) {
        std::string n = e.path().filename().string();
        if (!n.empty() && n[0] != '.')
            es.push_back(e);
    }
    std::sort(es.begin(), es.end(), [](auto &a, auto &b) {
        if (a.is_directory() != b.is_directory())
            return a.is_directory();
        return a.path().filename() < b.path().filename();
    });
    for (auto &e : es) {
        std::string n = e.path().filename().string(), p = e.path().string();
        bool dir = e.is_directory(ec);
        if (dirs_only && !dir)
            continue;
        ImGui::PushFont(T::mono, 13.0f);
        if (ImGui::Selectable(((dir ? "▸ " : "  ") + n).c_str(), sel.count(p) > 0, ImGuiSelectableFlags_AllowDoubleClick)) {
            if (dir && ImGui::IsMouseDoubleClicked(0)) {
                cwd = p;
                sel.clear();
                ImGui::PopFont();
                break;
            }
            if (!dir || dirs_only) {
                const ImGuiIO &io = ImGui::GetIO();
                if (!multi || (!io.KeyCtrl && !io.KeySuper && !io.KeyShift))
                    sel.clear();
                sel.insert(p);
            }
        }
        ImGui::PopFont();
    }
    ImGui::EndChild();
    // hint on its own line above the buttons (it used to be drawn under Cancel/Open in narrow windows)
    tag(multi ? tr("Ctrl/⌘ + click selects several · double-click opens a folder") : tr("double-click opens a folder"), T::faint);
    ImGui::SetCursorPosX(std::max(ImGui::GetCursorPosX(), ImGui::GetWindowWidth() - 250));
    if (ImGui::Button(tr("Cancel"), ImVec2(110, 0))) {
        open = false;
        ImGui::CloseCurrentPopup();
    }
    ImGui::SameLine();
    if (ImGui::Button(dirs_only ? tr("Use folder") : tr("Open"), ImVec2(110, 0)) && (dirs_only || !sel.empty())) {
        std::vector<std::string> v(sel.begin(), sel.end());
        if (dirs_only && v.empty())
            v.push_back(cwd);
        open = false;
        ImGui::CloseCurrentPopup();
        if (done)
            done(v);
    }
    ImGui::EndPopup();
}
