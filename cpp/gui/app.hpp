// Piazzi — declarations shared by the desktop app modules (Dear ImGui + SDL2 + OpenGL3).
// Global state `A` (App), theme `T`, and the functions each module exports to the others.
#pragma once
#include "../src/session_internal.hpp"
#include "imgui.h"
#include "piazzi/astro_c.h"
#include "piazzi/reduce.hpp"
#include "strings.hpp"
#include <SDL.h>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>
#ifdef __EMSCRIPTEN__
#include <GLES3/gl3.h>
#else
#include <SDL_opengl.h>
#endif

namespace fs = std::filesystem;

// ================================================================ "pocket observatory" theme (theme.cpp)
namespace T {
inline const ImVec4 ink(0.027f, 0.031f, 0.039f, 1), panel(0.055f, 0.063f, 0.075f, 1), line(1, 1, 1, 0.09f), amber(1.0f, 0.71f, 0.28f, 1),
    amberDim(1.0f, 0.71f, 0.28f, 0.16f), red(1.0f, 0.36f, 0.31f, 1), cyan(0.50f, 0.83f, 1.0f, 1), green(0.61f, 0.90f, 0.39f, 1),
    text(0.93f, 0.91f, 0.86f, 1), dim(0.93f, 0.91f, 0.86f, 0.55f), faint(0.93f, 0.91f, 0.86f, 0.32f);
inline ImU32 u(const ImVec4 &c, float a = -1) {
    return ImGui::ColorConvertFloat4ToU32(ImVec4(c.x, c.y, c.z, a < 0 ? c.w : a));
}
extern ImFont *serif, *mono, *sans;
} // namespace T

void apply_theme();
std::string upper_spaced(const std::string &t, bool accents); // uppercase, spaced per UTF-8 character
void tag(const std::string &t, ImVec4 col = T::dim);
void title(const char *t, float size);
bool primary_button(const char *t, float w = -1);
bool ghost_button(const char *t, float w = -1);
void reticle(ImDrawList *dl, ImVec2 c, float r, ImU32 col, float th = 1.5f, float gap = 0.45f);
void star_chart(ImDrawList *dl, ImVec2 p0, ImVec2 sz, unsigned seed = 7);
void row(const char *k, const std::string &v, ImVec4 col = T::text);

// ================================================================ utilities (util.cpp)
std::string fmt_hms(double v, bool sign);
double sexa(const std::string &s);
std::string lower(std::string x);
std::string url_encode(const std::string &s); // everything non-alphanumeric becomes %XX
std::string prefs_dir();
std::string assets_dir();
std::string output_dir(); // chosen output folder (or preferences), with trailing slash
std::string read_file(const std::string &p);

// ================================================================ processes and network (net.cpp)
// Runs argv[0] (looked up on PATH) with the given arguments, WITHOUT a shell, and waits for it to finish.
// Returns the exit code, or -1 if the process could not be started or was cancelled: with `cancel` (or
// A.quit) set, the child gets SIGTERM (TerminateProcess on Windows) within ~50 ms.
int run_process(const std::vector<std::string> &argv, const std::atomic<bool> *cancel = nullptr);
// Opens a folder in the system file manager (without waiting).
void open_folder(const std::string &dir);
// GET via the system curl (bundled with Windows 10+, macOS and Linux), HTTPS only. False if no attempt returned
// minbytes or if `cancel` (the job's own flag; null = none) or A.quit was set (aborts the running curl
// and the waits between attempts). The worker thread passes &A.cancel; Sesame and "What is this?" pass their own flags.
bool http_get(const std::string &url, std::string &out, const std::atomic<bool> *cancel, int timeout = 60, int tries = 2,
              size_t minbytes = 1);

struct Known {
    std::string number, name, cls;
    double ra = 0, dec = 0, v = 0, sep = 0;
};
// SkyBoT (IMCCE) cone search. False = network failure (out stays empty); true with empty out = nothing in the cone.
bool skybot(double jd, double ra, double dec, double rd, std::string obs, std::vector<Known> &out, const std::atomic<bool> *cancel);
bool download_catalog(); // worker thread: cancellable via A.cancel
void inspect_at(double x, double y);
// Resolves a name (Sesame/CDS) to RA/Dec degrees. False if not found or offline (net_ok tells which).
bool sesame_resolve(const std::string &name, double &ra, double &dec, bool &net_ok, const std::atomic<bool> *cancel);

// ================================================================ app state
struct Step {
    std::string title;
    int st = 0; // 0 pending 1 running 2 ok 3 failed
    std::string detail, at;
};
struct Profile {
    std::string name;
    double focal, pw, ph;
    int flip;
    std::string mpc;
};
enum class Phase { Home, Session, Results };

struct App {
    AstroSession *s = astro_new();
    std::recursive_mutex mu;
    std::atomic<Phase> phase{Phase::Home}; // read by the UI and also written by the worker thread
    std::vector<std::string> images;
    std::string catalog, ra, dec, object, status;
    std::vector<Step> steps;
    std::chrono::steady_clock::time_point t0;
    // busy: a session is running on the worker thread (set in run_async); cancel: asks the thread to stop;
    // quit: the app is shutting down (also stops the detached threads)
    std::atomic<bool> busy{false}, cancel{false}, quit{false};
    // detached threads still alive ("What is this?", open folder, Sesame): the end of main waits for them
    std::atomic<int> live_threads{0};
    std::thread worker;
    std::vector<Known> known;
    std::map<int, std::string> identified; // object index -> label "name · class · V · Δ"
    bool tex_dirty = false;
    std::vector<GLuint> tex;
    int k = 0;
    bool blink = true;
    double last_blink = 0;
    float zoom = 1;
    bool lay_track = true, lay_known = true, lay_refs = false, lay_dets = false;
    float lo_pct = 50, hi_pct = 99.5f;
    int stretch = 1;
    bool invert = false;
    int click_mode = 0; // 0 identify, 1 measure
    struct Inspect {
        bool on = false, loading = false;
        double x = 0, y = 0, ra = 0, dec = 0;
        std::string title, detail;
        std::vector<std::pair<std::string, std::string>> rows;
        ImVec4 col = T::dim;
    } insp;                     // guarded by `mu` (the query thread swaps the result under the lock)
    unsigned long insp_req = 0; // current query id; responses to older queries are discarded
    // cancel flag for the running "What is this?" query (a new query cancels the previous one); under `mu`
    std::shared_ptr<std::atomic<bool>> insp_cancel;
    struct Meas {
        bool on = false;
        double v[12]{};
        int frame = 0;
    } meas;
    int manual_track = -1;
    std::atomic<bool> need_center{false}; // set by the worker thread, consumed by the UI
    bool show_settings = false, show_mpc = false;
    std::string outdir;
    std::vector<Profile> profiles;
    std::string profile = "Default";
    std::atomic<bool> example{false};                 // built-in example session: MPC submission blocked
    std::map<std::string, std::string> user_settings; // user settings changed by example.cfg (restored afterwards)
    std::string dark, flat, mpcorb, comets;
    std::atomic<int> norbits{0}; // written by the thread that loads the orbits
};
extern App A;

// ======================================================== session and automatic pipeline (workflow.cpp, profiles.cpp)
void set(const std::string &k, const std::string &v);
void setn(const std::string &k, double v);
std::string get(const std::string &k);
void status(const std::string &t);
void fresh_steps();
void clear_session();
void add_file(const std::string &p);
void apply_header();
void pipeline();
void run_async(std::function<void()> fn);
void stop_worker(); // cancels and waits for the worker thread (call WITHOUT holding A.mu: it also uses the lock)
void start_session(std::function<void()> prep);
void start_with_files(const std::vector<std::string> &files); // new session from picked/dropped files
void rerun_pipeline();
void start_example(); // opens the built-in example (assets/example) and runs the full pipeline offline
void go_home();
std::vector<Profile> all_profiles();
void save_profiles();
void apply_profile(const Profile &p);
void load_profiles();

// ================================================================ file browser (browser.cpp)
struct Browser {
    bool open = false, multi = true, dirs_only = false;
    bool pending_open = false; // OpenPopup must be called in the same window that draws the modal (draw)
    std::string cwd, title, filter_hint;
    std::set<std::string> sel;
    std::function<void(std::vector<std::string>)> done;
    void show(const std::string &t, bool m, std::function<void(std::vector<std::string>)> cb, bool d = false);
    void draw();
};
extern Browser BR;

// ================================================================ screens
void ui_home();                         // screens_home.cpp
void ui_session();                      // screens_home.cpp
void ui_results();                      // screens_results.cpp
void make_textures();                   // screens_results.cpp
void sheet();                           // screens_sheet.cpp
void sheet_tools();                     // screens_tools.cpp
void save_reports();                    // screens_sheet.cpp
void ui_center();                       // screens_dialogs.cpp
void ui_settings();                     // screens_dialogs.cpp
void ui_mpc();                          // screens_mpc.cpp
bool valid_email(const std::string &e); // screens_mpc.cpp
