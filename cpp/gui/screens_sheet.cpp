// Results side panel: tr("What is this?"), manual measurement, object cards, plate quality.
#include "app.hpp"
#include "imgui.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

void save_reports() {
    std::string d = output_dir();
    std::ofstream(d + "ADESReport.psv", std::ios::binary) << astro_report(A.s, 0);
    std::ofstream(d + "PhotReport.txt", std::ios::binary) << astro_report(A.s, 1);
    std::ofstream(d + "MPCReport.txt", std::ios::binary) << astro_report(A.s, 2);
    status(tr("Reports saved to ") + d);
    open_folder(d);
}

static void inspect_card() {
    App::Inspect I;
    {
        std::lock_guard<std::recursive_mutex> lk(A.mu); // copy under lock: the query thread may replace A.insp
        I = A.insp;
    }
    ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(1, 1, 1, 0.35f));
    ImGui::PushStyleColor(ImGuiCol_ChildBg, T::panel);
    ImGui::BeginChild("##insp", ImVec2(0, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_AutoResizeY);
    tag(tr("What is this?"), T::text);
    ImGui::SameLine(ImGui::GetContentRegionAvail().x - 10);
    if (ImGui::SmallButton("✕")) {
        std::lock_guard<std::recursive_mutex> lk(A.mu);
        A.insp.on = false;
    }
    ImGui::PushFont(T::mono, 12.5f);
    ImGui::TextColored(T::dim, "%s   %s", fmt_hms(I.ra / 15, false).c_str(), fmt_hms(I.dec, true).c_str());
    ImGui::PopFont();
    if (I.loading) {
        ImGui::PushFont(T::mono, 12.5f);
        ImGui::TextColored(T::dim, "%s", tr("Querying SIMBAD, VSX, Gaia and SkyBoT…"));
        ImGui::PopFont();
    } else {
        title(I.title.c_str(), 24);
        ImGui::TextColored(I.col, "%s", I.detail.c_str());
        for (auto &[k, v] : I.rows)
            row(k.c_str(), v);
    }
    ImGui::EndChild();
    ImGui::PopStyleColor(2);
    ImGui::Dummy(ImVec2(0, 6));
}

static void measure_card() {
    ImGui::PushStyleColor(ImGuiCol_ChildBg, T::panel);
    ImGui::BeginChild("##meas", ImVec2(0, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_AutoResizeY);
    tag(tr("Manual measurement · plate ") + std::to_string(A.meas.frame + 1), T::amber);
    ImGui::SameLine(ImGui::GetContentRegionAvail().x - 10);
    if (ImGui::SmallButton("✕##m"))
        A.meas.on = false;
    auto &v = A.meas.v;
    char b[64];
    row(tr("RA / Dec"), fmt_hms(v[0] / 15, false) + "  " + fmt_hms(v[1], true));
    std::snprintf(b, sizeof b, "%.2f", v[2]);
    row("Magnitude", b);
    std::snprintf(b, sizeof b, "%.1f / %.1f", v[3], v[9]);
    row(tr("SNR peak / total"), b);
    std::snprintf(b, sizeof b, "%.1f\"", v[4]);
    row("FWHM", b);
    std::snprintf(b, sizeof b, "%.2f, %.2f", v[5], v[6]);
    row("x, y", b);
    std::snprintf(b, sizeof b, "%.3f", v[7]);
    row("Fit RMS", b);
    if (ghost_button(A.manual_track < 0 ? tr("CREATE OBJECT FROM THIS MEASUREMENT") : tr("ADD TO THE MANUAL OBJECT"))) {
        A.manual_track = astro_add_manual(A.s, A.manual_track, A.meas.frame, v[5], v[6]);
        status(tr("Measurement added to the manual object."));
    }
    ImGui::EndChild();
    ImGui::PopStyleColor();
    ImGui::Dummy(ImVec2(0, 6));
}

// split "a · b · c" into its parts
static std::vector<std::string> split_dot(const std::string &id) {
    std::vector<std::string> parts;
    size_t a = 0, b;
    while ((b = id.find(" · ", a)) != std::string::npos) {
        parts.push_back(id.substr(a, b - a));
        a = b + 4; // " · " is 4 bytes in UTF-8
    }
    parts.push_back(id.substr(a));
    return parts;
}

// card for a moving object: identification, rate, measurements (MPC note, exclusion) and reports
static void object_card(size_t i, const astro::Track &m) {
    std::string id = A.identified.count((int)i) ? A.identified[(int)i] : "";
    bool known = !id.empty();
    ImGui::PushID((int)i);
    ImGui::PushStyleColor(ImGuiCol_Border, known ? ImVec4(1, 0.71f, 0.28f, 0.45f) : ImVec4(1, 0.36f, 0.31f, 0.45f));
    ImGui::PushStyleColor(ImGuiCol_ChildBg, T::panel);
    ImGui::BeginChild("##obj", ImVec2(0, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_AutoResizeY);
    tag(known ? tr("Object sheet · identified") : tr("Object sheet · not identified"), known ? T::amber : T::red);
    std::vector<std::string> parts = split_dot(id);
    title(known ? parts[0].c_str() : (m.desig.empty() ? tr("Possible new object") : m.desig.c_str()), 28);
    char b[64];
    if (m.desig.empty()) { // no designation: the report uses the track's temporary designation; here a display label
        std::snprintf(b, sizeof b, tr("object %zu"), i + 1);
        row(tr("Label"), b);
    }
    if (known && parts.size() > 1)
        row(tr("Class"), parts[1]);
    if (known && parts.size() > 2)
        row(tr("Detail"), parts[2]);
    std::snprintf(b, sizeof b, "%.2f\"/min", m.rate);
    row(tr("Rate"), b);
    std::snprintf(b, sizeof b, "%.1f°", m.pa);
    row(tr("Position angle"), b);
    tag(tr("Measurements · MPC note · exclude"));
    for (size_t j = 0; j < m.meas.size(); ++j) {
        auto &q = m.meas[j];
        ImGui::PushID((int)j);
        ImGui::PushFont(T::mono, 12.0f);
        ImGui::TextColored(T::amber, "%zu", j + 1);
        ImGui::SameLine();
        ImGui::TextColored(T::dim, "%s", astro::jd_to_isot(q.jd).substr(11).c_str());
        ImGui::SameLine();
        ImGui::Text("%s %s", astro::fmt_ra(q.ra).c_str(), astro::fmt_dec(q.dec).c_str());
        ImGui::SameLine();
        ImGui::TextColored(T::dim, "%.1f", q.mag);
        ImGui::SameLine();
        char nb[4] = {0};
        auto key = std::make_pair((int)i, q.frame);
        if (A.s->notes.count(key))
            std::snprintf(nb, sizeof nb, "%s", A.s->notes[key].c_str());
        ImGui::SetNextItemWidth(30);
        if (ImGui::InputText("##n", nb, 2))
            astro_set_note(A.s, (int)i, q.frame, nb);
        ImGui::SameLine();
        bool ex = A.s->excluded.count(key) && A.s->excluded[key];
        if (ImGui::Checkbox("##x", &ex))
            astro_set_excluded(A.s, (int)i, q.frame, ex);
        ImGui::PopFont();
        ImGui::PopID();
    }
    ImGui::Dummy(ImVec2(0, 4));
    if (primary_button(tr("SUBMIT TO MPC")))
        A.show_mpc = true;
    if (ghost_button(tr("SAVE ADES · PHOTOMETRIC · MPC 80")))
        save_reports();
    ImGui::EndChild();
    ImGui::PopStyleColor(2);
    ImGui::PopID();
    ImGui::Dummy(ImVec2(0, 6));
}

static void plate_quality(const astro::Reduction &R) {
    ImGui::Dummy(ImVec2(0, 6));
    tag(tr("Plate quality"));
    for (size_t i = 0; i < R.frames.size(); ++i) {
        auto &f = R.frames[i];
        double r = (f.plate.dra_mean + f.plate.dde_mean) / 2;
        ImVec4 c = f.plate.nref < 6 ? T::red : r < 0.25 ? T::green : r < 0.6 ? T::amber : T::red;
        const char *t = f.plate.nref < 6 ? tr("few stars") : r < 0.25 ? tr("excellent") : r < 0.6 ? tr("good") : tr("poor");
        ImGui::PushFont(T::mono, 12.5f);
        ImGui::Text(tr("Plate %zu%s"), i + 1, f.stacked ? tr(" (stacked)") : "");
        ImGui::SameLine();
        ImGui::TextColored(c, "%s", t);
        ImGui::SameLine(ImGui::GetContentRegionAvail().x - 170);
        ImGui::TextColored(T::dim, "%d refs · %.2f\" · %.1f\"", f.plate.nref, r, f.seeing);
        ImGui::PopFont();
        ImVec2 p = ImGui::GetCursorScreenPos();
        float w = ImGui::GetContentRegionAvail().x;
        ImDrawList *dl = ImGui::GetWindowDrawList();
        dl->AddRectFilled(p, ImVec2(p.x + w, p.y + 3), T::u(T::line));
        dl->AddRectFilled(p, ImVec2(p.x + w * (float)std::clamp(1 - r, 0.05, 1.0), p.y + 3), T::u(c));
        ImGui::Dummy(ImVec2(0, 6));
    }
}

void sheet() {
    auto &R = A.s->R;
    bool insp_on;
    {
        std::lock_guard<std::recursive_mutex> lk(A.mu);
        insp_on = A.insp.on;
    }
    if (insp_on)
        inspect_card();
    if (A.meas.on)
        measure_card();
    if (R.movers.empty()) {
        tag(tr("Result"), T::amber);
        title(tr("Nothing moved."), 26);
        ImGui::TextWrapped(tr("No object moved between the images above the threshold. %d known asteroid(s) in the field."),
                           (int)A.known.size());
    }
    for (size_t i = 0; i < R.movers.size(); ++i)
        object_card(i, R.movers[i]);
    sheet_tools();
    plate_quality(R);
    if (!A.known.empty()) {
        ImGui::Dummy(ImVec2(0, 6));
        tag(tr("Known asteroids in the field"));
        for (size_t i = 0; i < A.known.size() && i < 8; ++i) {
            const Known &k = A.known[i];
            char b[32];
            std::snprintf(b, sizeof b, "V %.1f", k.v);
            row(((k.number.empty() ? "" : "(" + k.number + ") ") + k.name).c_str(), b, k.v < 18.5 ? T::cyan : T::dim);
        }
    }
    ImGui::Dummy(ImVec2(0, 8));
    if (ghost_button(tr("NEW SESSION")))
        go_home();
}
