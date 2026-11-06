// Dialogs: field center (with name lookup) and settings (profiles, parameters, calibration, output).
#include "app.hpp"
#include "imgui.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

struct SesameJob {
    std::atomic<bool> cancel{false}; // own flag: cancelling this lookup does not touch the session (A.cancel)
    std::mutex mu;
    bool done = false, found = false, net_ok = false;
    double ra = 0, dec = 0;
};

void ui_center() {
    static std::string msg;
    // Sesame lookup in the background: state is shared via shared_ptr and read under the job's own lock
    static std::shared_ptr<SesameJob> job;
    if (A.need_center.exchange(false)) {
        // new popup: nothing from the previous lookup (message or lookup still running, which is cancelled)
        if (job)
            job->cancel = true;
        job.reset();
        msg.clear();
        ImGui::OpenPopup("##center");
    }
    ImGui::SetNextWindowSize(ImVec2(520, 0), ImGuiCond_Appearing);
    if (!ImGui::BeginPopupModal("##center", nullptr, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_AlwaysAutoResize))
        return;
    tag(tr("Field centre"), T::amber);
    title(tr("Where was the telescope pointing?"), 24);
    ImGui::TextWrapped("%s", tr("No coordinates in the header. Search for the target by name or type RA and Dec."));
    static char name[128] = "", ra[64] = "", dec[64] = "";
    ImGui::SetNextItemWidth(380);
    ImGui::InputTextWithHint("##nm", "M42, NGC 2264, Betelgeuse…", name, sizeof name);
    ImGui::SameLine();
    bool searching = false;
    if (job) {
        std::lock_guard<std::mutex> lk(job->mu);
        if (job->done) {
            if (job->found) {
                std::snprintf(ra, sizeof ra, "%s", fmt_hms(job->ra / 15, false).c_str());
                std::snprintf(dec, sizeof dec, "%s", fmt_hms(job->dec, true).c_str());
                msg = tr("found");
            } else
                msg = job->net_ok ? tr("not found") : tr("no answer from Sesame (no internet?)");
        } else
            searching = true;
    }
    if (job && !searching)
        job.reset();
    ImGui::BeginDisabled(searching);
    if (ImGui::Button(searching ? tr("Searching…") : tr("Search")) && name[0]) {
        job = std::make_shared<SesameJob>();
        msg.clear();
        A.live_threads++;
        std::thread([j = job, q = std::string(name)] {
            double a = 0, d = 0;
            bool net_ok = false;
            bool found = sesame_resolve(q, a, d, net_ok, &j->cancel);
            {
                std::lock_guard<std::mutex> lk(j->mu);
                j->found = found;
                j->net_ok = net_ok;
                j->ra = a;
                j->dec = d;
                j->done = true;
            }
            A.live_threads--;
        }).detach();
    }
    ImGui::EndDisabled();
    if (!msg.empty())
        tag(msg, msg == tr("found") ? T::green : T::red);
    ImGui::SetNextItemWidth(240);
    ImGui::InputText(tr("RA h m s"), ra, sizeof ra);
    ImGui::SetNextItemWidth(240);
    ImGui::InputText("Dec ° ' \"", dec, sizeof dec);
    if (ImGui::Button(tr("Cancel"), ImVec2(120, 0))) {
        if (job)
            job->cancel = true;
        ImGui::CloseCurrentPopup();
    }
    ImGui::SameLine();
    if (primary_button(tr("CONTINUE"), 200) && ra[0] && dec[0]) {
        A.ra = ra;
        A.dec = dec;
        ImGui::CloseCurrentPopup();
        rerun_pipeline();
    }
    ImGui::EndPopup();
}

static void settings_profiles() {
    tag(tr("Equipment"), T::amber);
    for (auto &p : all_profiles()) {
        char b[200];
        std::snprintf(b, sizeof b, "%s   ·   %.0f mm · %.1f×%.1f µm · %.2f\"/px · %s", p.name.c_str(), p.focal, p.pw * 1000, p.ph * 1000,
                      p.pw / p.focal * 206264.8, p.mpc.c_str());
        ImGui::PushFont(T::mono, 12.5f);
        if (ImGui::Selectable(b, p.name == A.profile))
            apply_profile(p);
        ImGui::PopFont();
    }
    static char pn[80] = "";
    ImGui::SetNextItemWidth(380);
    ImGui::InputTextWithHint("##pn", tr("Profile name (e.g. C11 + ASI533)"), pn, sizeof pn);
    ImGui::SameLine();
    if (ImGui::Button(tr("Save current")) && pn[0]) {
        try {
            Profile p{pn,
                      std::stod(get("focal_length_mm")),
                      std::stod(get("pixel_width_mm")),
                      std::stod(get("pixel_height_mm")),
                      std::stoi(get("flip_horizontal")),
                      get("mpc_code")};
            A.profiles.erase(std::remove_if(A.profiles.begin(), A.profiles.end(), [&](auto &x) { return x.name == p.name; }),
                             A.profiles.end());
            A.profiles.push_back(p);
            apply_profile(p);
            pn[0] = 0;
        } catch (...) {
            status(tr("Invalid equipment parameters — profile not saved."));
        }
    }
}

// ---------------------------------------------------------------- parameter screen generated from the schema
// minimal JSON (objects, arrays, strings, numbers, true/false/null) to read astro_settings_schema()
namespace {
struct J {
    enum Ty { Null, Bool, Num, Str, Arr, Obj } t = Null;
    bool b = false;
    double n = 0;
    std::string s;
    std::vector<J> a;
    std::vector<std::pair<std::string, J>> o;
    const J &operator[](const char *k) const {
        static const J none;
        for (auto &kv : o)
            if (kv.first == k)
                return kv.second;
        return none;
    }
};
struct Parser {
    const char *p;
    void ws() {
        while (*p == ' ' || *p == '\n' || *p == '\t' || *p == '\r')
            ++p;
    }
    std::string str() {
        std::string out;
        ++p; // quote
        while (*p && *p != '"') {
            if (*p == '\\' && p[1]) {
                ++p;
                char c = *p == 'n' ? '\n' : *p == 't' ? '\t' : *p;
                if (*p == 'u') { // \uXXXX (only ASCII appears in the schema)
                    out += char(std::strtol(std::string(p + 1, 4).c_str(), nullptr, 16));
                    p += 5;
                    continue;
                }
                out += c;
                ++p;
            } else
                out += *p++;
        }
        if (*p)
            ++p;
        return out;
    }
    J val() {
        ws();
        J v;
        if (*p == '{') {
            v.t = J::Obj;
            ++p;
            ws();
            while (*p && *p != '}') {
                ws();
                std::string k = str();
                ws();
                if (*p == ':')
                    ++p;
                v.o.push_back({k, val()});
                ws();
                if (*p == ',')
                    ++p;
                ws();
            }
            if (*p)
                ++p;
        } else if (*p == '[') {
            v.t = J::Arr;
            ++p;
            ws();
            while (*p && *p != ']') {
                v.a.push_back(val());
                ws();
                if (*p == ',')
                    ++p;
                ws();
            }
            if (*p)
                ++p;
        } else if (*p == '"') {
            v.t = J::Str;
            v.s = str();
        } else if (!std::strncmp(p, "true", 4)) {
            v.t = J::Bool;
            v.b = true;
            p += 4;
        } else if (!std::strncmp(p, "false", 5)) {
            v.t = J::Bool;
            p += 5;
        } else if (!std::strncmp(p, "null", 4)) {
            p += 4;
        } else {
            char *e = nullptr;
            v.t = J::Num;
            v.n = std::strtod(p, &e);
            p = e && e != p ? e : p + 1;
        }
        return v;
    }
};
const J &schema() {
    static const J S = [] {
        Parser ps{astro_settings_schema()};
        return ps.val();
    }();
    return S;
}
// text in the current language from a node with "en"/"pt" keys (or hint_en/hint_pt)
const std::string &loc(const J &j, const char *en = "en", const char *pt = "pt") {
    const J &v = j[lang_pt() ? pt : en];
    return v.s.empty() ? j[en].s : v.s;
}
bool parse_num(const std::string &x, double &v) {
    char *e = nullptr;
    v = std::strtod(x.c_str(), &e);
    return !x.empty() && e && *e == 0;
}
} // namespace

static void setting_field(const J &f) {
    const std::string &key = f["key"].s, &type = f["type"].s, &def = f["default"].s;
    std::string cur = get(key);
    const char *label = loc(f).c_str();
    ImGui::PushID(key.c_str());
    // in-progress numeric edit per field: written to the core only when the value is valid
    static std::map<std::string, std::string> draft;
    bool invalid = false;
    const char *why = "";
    ImGui::SetNextItemWidth(220);
    if (type == "bool") {
        bool b = cur.empty() ? def == "1" : std::atof(cur.c_str()) != 0;
        if (ImGui::Checkbox(label, &b))
            set(key, b ? "1" : "0");
    } else if (type == "enum") {
        const auto &opts = f["options"].a;
        int sel = -1;
        std::string v = cur.empty() ? def : cur;
        for (size_t i = 0; i < opts.size(); ++i) {
            double a, b;
            if (opts[i]["value"].s == v || (parse_num(opts[i]["value"].s, a) && parse_num(v, b) && a == b))
                sel = int(i);
        }
        const char *preview = sel >= 0 ? loc(opts[sel]).c_str() : v.c_str();
        if (ImGui::BeginCombo(label, preview)) {
            for (size_t i = 0; i < opts.size(); ++i)
                if (ImGui::Selectable(loc(opts[i]).c_str(), int(i) == sel))
                    set(key, opts[i]["value"].s);
            ImGui::EndCombo();
        }
    } else {
        auto it = draft.find(key);
        std::string shown = it != draft.end() ? it->second : cur;
        char b[256];
        std::snprintf(b, sizeof b, "%s", shown.c_str());
        bool edited = ImGui::InputText(label, b, sizeof b);
        if (type == "number" || type == "int") {
            double v = 0;
            std::string t = b;
            if (!t.empty() && !parse_num(t, v)) {
                invalid = true;
                why = tr("not a number");
            } else if (!t.empty() && (v < f["min"].n || v > f["max"].n || (type == "int" && v != std::floor(v)))) {
                invalid = true;
                why = tr("out of range");
            }
            if (edited) {
                if (invalid)
                    draft[key] = t;
                else {
                    draft.erase(key);
                    set(key, t);
                }
            }
        } else if (edited)
            set(key, b);
    }
    // hint: tooltip with the explanation, the default and the range
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort | ImGuiHoveredFlags_ForTooltip)) {
        std::string tip = loc(f, "hint_en", "hint_pt");
        tip += std::string("\n") + tr("default") + ": " + (def.empty() ? "—" : def);
        if (type == "number" || type == "int") {
            char r[96];
            std::snprintf(r, sizeof r, "  [%g … %g]", f["min"].n, f["max"].n);
            tip += r;
        }
        ImGui::SetTooltip("%s", tip.c_str());
    }
    ImGui::SameLine();
    bool differs = (type == "number" || type == "int") ? [&] {
        double a, b;
        return !(parse_num(cur, a) && parse_num(def, b) && a == b) && cur != def;
    }()
                                                       : cur != def;
    ImGui::BeginDisabled(!differs && draft.find(key) == draft.end());
    if (ImGui::SmallButton("↺")) {
        draft.erase(key);
        set(key, def);
    }
    ImGui::SetItemTooltip("%s", tr("Restore default"));
    ImGui::EndDisabled();
    if (invalid || draft.count(key)) {
        ImGui::SameLine();
        ImGui::TextColored(T::red, "%s", *why ? why : tr("out of range"));
    }
    ImGui::PopID();
}

static void settings_params() {
    tag(tr("Parameters"), T::amber);
    const J &S = schema();
    for (const J &g : S["groups"].a) {
        const std::string &gid = g["id"].s;
        // fixed id (###) so the open/closed state does not depend on the language; "advanced" starts collapsed
        std::string head = loc(g) + "###grp_" + gid;
        if (!ImGui::CollapsingHeader(head.c_str(), g["collapsed"].b ? 0 : ImGuiTreeNodeFlags_DefaultOpen))
            continue;
        ImGui::Indent(8);
        for (const J &f : S["settings"].a)
            if (f["group"].s == gid)
                setting_field(f);
        ImGui::Unindent(8);
    }
    if (ImGui::Button(tr("Restore defaults"))) {
        for (const J &f : S["settings"].a)
            set(f["key"].s, f["default"].s);
        status(tr("Defaults restored."));
    }
}

static void pick_file(const char *label, std::string &path, bool orbits) {
    ImGui::PushFont(T::mono, 12.0f);
    ImGui::TextColored(T::dim, "%s: %s", label, path.empty() ? tr("(none)") : fs::path(path).filename().string().c_str());
    ImGui::PopFont();
    ImGui::SameLine();
    if (ImGui::SmallButton((std::string(tr("choose")) + "##" + label).c_str()))
        BR.show(label, false, [&path, orbits](std::vector<std::string> v) {
            if (A.busy) { // the session may have started with the browser open: do not touch it while processing
                status(tr("Wait for the session to finish before changing this file."));
                return;
            }
            path = v[0];
            if (orbits) {
                // the full MPCORB takes seconds: read it on the worker thread (A.busy locks the session meanwhile),
                // without freezing the UI
                std::string mp = A.mpcorb, co = A.comets;
                status(tr("Loading orbits…"));
                run_async([mp, co] {
                    A.norbits = astro_load_orbits(A.s, mp.c_str(), co.c_str());
                    status(std::to_string(A.norbits.load()) + tr(" orbits loaded."));
                });
            }
        });
    if (!path.empty()) {
        ImGui::SameLine();
        if (ImGui::SmallButton((std::string(tr("remove")) + "##" + label).c_str()))
            path.clear();
    }
}

void ui_settings() {
    if (!A.show_settings)
        return;
    ImGui::SetNextWindowSize(ImVec2(680, 720), ImGuiCond_Appearing);
    std::string wt = std::string(tr("Settings")) + "###settings"; // fixed id: language switch keeps the window
    if (!ImGui::Begin(wt.c_str(), &A.show_settings, ImGuiWindowFlags_NoCollapse)) {
        ImGui::End();
        return;
    }
    // language: stays outside the lock (does not touch the session)
    int li = lang_pt() ? 1 : 0;
    const char *langs[] = {"English", "Português"};
    ImGui::SetNextItemWidth(220);
    if (ImGui::Combo(tr("Language"), &li, langs, 2))
        set_lang(li == 1 ? "pt" : "en");
    ImGui::Separator();
    // session running: the worker uses the session (astro_reduce etc.) without a lock — nothing here may touch it
    bool busy = A.busy;
    if (busy)
        ImGui::TextColored(T::amber, "%s", tr("Session running — settings are locked until it finishes."));
    ImGui::BeginDisabled(busy);
    settings_profiles();
    ImGui::Separator();
    settings_params();
    ImGui::Separator();
    tag(tr("Calibration and orbits"), T::amber);
    pick_file("Dark", A.dark, false);
    pick_file("Flat", A.flat, false);
    pick_file("MPCORB.DAT", A.mpcorb, true);
    pick_file("COMET.DAT", A.comets, true);
    ImGui::TextColored(T::faint, tr("%d orbits loaded (offline identification)."), A.norbits.load());
    ImGui::Separator();
    tag(tr("Output"), T::amber);
    ImGui::PushFont(T::mono, 12.0f);
    ImGui::TextColored(T::dim, "%s", (A.outdir.empty() ? prefs_dir() : A.outdir).c_str());
    ImGui::PopFont();
    ImGui::SameLine();
    if (ImGui::SmallButton(tr("choose folder")))
        BR.show(
            tr("Folder for reports and FITS"), false,
            [](std::vector<std::string> v) {
                if (!A.busy)
                    A.outdir = v[0];
            },
            true);
    ImGui::EndDisabled();
    ImGui::End();
}
