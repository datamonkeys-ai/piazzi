// "Pocket observatory" theme and shared drawing widgets.
#include "app.hpp"
#include "imgui.h"
#include <cctype>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>

namespace T {
ImFont *serif = nullptr, *mono = nullptr, *sans = nullptr;
} // namespace T

void apply_theme() {
    ImGuiStyle &s = ImGui::GetStyle();
    s.WindowRounding = 0;
    s.ChildRounding = 8;
    s.FrameRounding = 6;
    s.PopupRounding = 8;
    s.GrabRounding = 4;
    s.TabRounding = 6;
    s.WindowBorderSize = 0;
    s.ChildBorderSize = 1;
    s.FrameBorderSize = 1;
    s.ItemSpacing = ImVec2(10, 9);
    s.FramePadding = ImVec2(10, 7);
    s.WindowPadding = ImVec2(22, 18);
    s.ScrollbarSize = 10;
    s.ScrollbarRounding = 6;
    ImVec4 *c = s.Colors;
    c[ImGuiCol_WindowBg] = T::ink;
    c[ImGuiCol_ChildBg] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_PopupBg] = T::panel;
    c[ImGuiCol_Border] = T::line;
    c[ImGuiCol_Text] = T::text;
    c[ImGuiCol_TextDisabled] = T::faint;
    c[ImGuiCol_FrameBg] = T::panel;
    c[ImGuiCol_FrameBgHovered] = ImVec4(0.09f, 0.10f, 0.12f, 1);
    c[ImGuiCol_FrameBgActive] = ImVec4(0.11f, 0.12f, 0.14f, 1);
    c[ImGuiCol_Button] = T::panel;
    c[ImGuiCol_ButtonHovered] = ImVec4(0.10f, 0.11f, 0.13f, 1);
    c[ImGuiCol_ButtonActive] = T::amberDim;
    c[ImGuiCol_Header] = T::amberDim;
    c[ImGuiCol_HeaderHovered] = ImVec4(1, 0.71f, 0.28f, 0.22f);
    c[ImGuiCol_HeaderActive] = ImVec4(1, 0.71f, 0.28f, 0.30f);
    c[ImGuiCol_CheckMark] = T::amber;
    c[ImGuiCol_SliderGrab] = T::amber;
    c[ImGuiCol_SliderGrabActive] = T::amber;
    c[ImGuiCol_Separator] = T::line;
    c[ImGuiCol_TitleBg] = T::ink;
    c[ImGuiCol_TitleBgActive] = T::panel;
    c[ImGuiCol_ScrollbarBg] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_ScrollbarGrab] = T::line;
    c[ImGuiCol_PlotHistogram] = T::amber;
    c[ImGuiCol_TableHeaderBg] = T::panel;
    c[ImGuiCol_TableBorderLight] = T::line;
    c[ImGuiCol_TableBorderStrong] = T::line;
    c[ImGuiCol_TableRowBgAlt] = ImVec4(1, 1, 1, 0.025f);
    c[ImGuiCol_ModalWindowDimBg] = ImVec4(0, 0, 0, 0.6f);
    c[ImGuiCol_TextSelectedBg] = T::amberDim;
}

// uppercase and spacing per UTF-8 character (not per byte); with accents, common Portuguese accents are uppercased
std::string upper_spaced(const std::string &t, bool accents) {
    static const std::map<std::string, std::string> up = {{"á", "Á"}, {"à", "À"}, {"â", "Â"}, {"ã", "Ã"}, {"é", "É"}, {"ê", "Ê"},
                                                          {"í", "Í"}, {"ó", "Ó"}, {"ô", "Ô"}, {"õ", "Õ"}, {"ú", "Ú"}, {"ç", "Ç"}};
    std::string u;
    for (size_t i = 0; i < t.size();) {
        unsigned char c = t[i];
        size_t n = c < 0x80 ? 1 : (c >> 5) == 6 ? 2 : (c >> 4) == 14 ? 3 : 4;
        std::string ch = t.substr(i, n);
        i += n;
        if (n == 1)
            ch = std::string(1, (char)std::toupper(c));
        else if (accents && up.count(ch))
            ch = up.at(ch);
        u += ch;
        u += ' ';
    }
    if (!u.empty())
        u.pop_back();
    return u;
}

// instrument-style text: uppercase, mono, letter-spaced
void tag(const std::string &t, ImVec4 col) {
    std::string u = upper_spaced(t, true);
    ImGui::PushFont(T::mono, 11.0f);
    ImGui::TextColored(col, "%s", u.c_str());
    ImGui::PopFont();
}
void title(const char *t, float size) {
    ImGui::PushFont(T::serif, size);
    ImGui::TextUnformatted(t);
    ImGui::PopFont();
}
bool primary_button(const char *t, float w) {
    ImGui::PushStyleColor(ImGuiCol_Button, T::amber);
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(1, 0.77f, 0.40f, 1));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.9f, 0.62f, 0.22f, 1));
    ImGui::PushStyleColor(ImGuiCol_Text, T::ink);
    ImGui::PushFont(T::mono, 14.0f);
    bool r = ImGui::Button(t, ImVec2(w < 0 ? ImGui::GetContentRegionAvail().x : w, 44));
    ImGui::PopFont();
    ImGui::PopStyleColor(4);
    return r;
}
bool ghost_button(const char *t, float w) {
    ImGui::PushFont(T::mono, 13.0f);
    bool r = ImGui::Button(t, ImVec2(w < 0 ? ImGui::GetContentRegionAvail().x : w, 38));
    ImGui::PopFont();
    return r;
}
void reticle(ImDrawList *dl, ImVec2 c, float r, ImU32 col, float th, float gap) {
    float l = r * gap;
    for (int sx : {-1, 1})
        for (int sy : {-1, 1}) {
            ImVec2 k(c.x + r * sx, c.y + r * sy);
            dl->AddLine(ImVec2(k.x, k.y - l * sy), k, col, th);
            dl->AddLine(k, ImVec2(k.x - l * sx, k.y), col, th);
        }
}
void star_chart(ImDrawList *dl, ImVec2 p0, ImVec2 sz, unsigned seed) {
    for (int i = 0; i <= 8; ++i) {
        float x = p0.x + sz.x * i / 8;
        dl->AddBezierQuadratic(ImVec2(x, p0.y), ImVec2(x + (x - p0.x - sz.x / 2) * .03f, p0.y + sz.y / 2),
                               ImVec2(x + (x - p0.x - sz.x / 2) * .12f, p0.y + sz.y), IM_COL32(255, 255, 255, 9), 0.6f);
    }
    for (int j = 0; j <= 12; ++j) {
        float y = p0.y + sz.y * j / 12;
        dl->AddBezierQuadratic(ImVec2(p0.x, y), ImVec2(p0.x + sz.x / 2, y - 18), ImVec2(p0.x + sz.x, y), IM_COL32(255, 255, 255, 9), 0.6f);
    }
    uint64_t s = seed;
    auto rnd = [&]() {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return float(s >> 33) / 2147483648.0f;
    };
    for (int i = 0; i < 260; ++i) {
        float x = p0.x + rnd() * sz.x, y = p0.y + rnd() * sz.y, m = rnd();
        float r = 0.35f + std::pow(m, 6) * 2.4f;
        dl->AddCircleFilled(ImVec2(x, y), r * 0.5f, IM_COL32(255, 255, 255, (int)(255 * (0.15f + m * 0.5f))));
    }
}
// "label ........ value" row with separator
void row(const char *k, const std::string &v, ImVec4 col) {
    tag(k);
    ImGui::SameLine(ImGui::GetContentRegionAvail().x - ImGui::CalcTextSize(v.c_str()).x * 1.0f);
    ImGui::PushFont(T::mono, 13.5f);
    ImGui::TextColored(col, "%s", v.c_str());
    ImGui::PopFont();
    ImGui::Separator();
}
