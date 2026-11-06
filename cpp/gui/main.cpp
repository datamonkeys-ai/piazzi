// Piazzi — native C++ desktop app (Dear ImGui + SDL2 + OpenGL3).
// A single codebase for Windows, Linux and macOS; same core as the iPhone and web apps.
// Main loop: window, fonts, events (file drag-and-drop) and screen dispatch.
#include "app.hpp"
#include "imgui.h"
#include <SDL.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <mutex>
#ifdef __EMSCRIPTEN__
#include <GLES3/gl3.h>
#else
#include <SDL_opengl.h>
#endif
#include "imgui_impl_opengl3.h"
#include "imgui_impl_sdl2.h"
#include <cstdio>
#include <string>
#include <vector>

// window icon drawn in code (no PNG decoder): the same artwork as packaging/icon/piazzi.svg —
// night-blue rounded square, amber (#f0a24b) four-corner reticle, white dot and short asteroid trail
static void set_window_icon(SDL_Window *win) {
    const int N = 64;
    SDL_Surface *sf = SDL_CreateRGBSurfaceWithFormat(0, N, N, 32, SDL_PIXELFORMAT_RGBA32);
    if (!sf)
        return;
    auto *px = static_cast<Uint8 *>(sf->pixels);
    auto seg = [](float x, float y, float ax, float ay, float bx, float by) { // distance to segment
        float dx = bx - ax, dy = by - ay, t = ((x - ax) * dx + (y - ay) * dy) / (dx * dx + dy * dy);
        t = t < 0 ? 0 : t > 1 ? 1 : t;
        float ex = x - ax - t * dx, ey = y - ay - t * dy;
        return std::sqrt(ex * ex + ey * ey);
    };
    // coordinates in the SVG's 1024 space
    const float L[][4] = {{232, 400, 232, 232}, {232, 232, 400, 232}, {624, 232, 792, 232}, {792, 232, 792, 400},
                          {792, 624, 792, 792}, {792, 792, 624, 792}, {400, 792, 232, 792}, {232, 792, 232, 624}};
    for (int y = 0; y < N; ++y)
        for (int x = 0; x < N; ++x) {
            float X = (x + 0.5f) * 1024 / N, Y = (y + 0.5f) * 1024 / N, r = 224;
            float cx = X < r ? r - X : X > 1024 - r ? X - (1024 - r) : 0, cy = Y < r ? r - Y : Y > 1024 - r ? Y - (1024 - r) : 0;
            float inside = std::sqrt(cx * cx + cy * cy) <= r ? 1.0f : 0.0f;
            float g = std::sqrt((X - 512) * (X - 512) + (Y - 390) * (Y - 390)) / 760; // radial sky gradient
            g = g > 1 ? 1 : g;
            float R = 0x13 + (0x06 - 0x13) * g, G = 0x23 + (0x0b - 0x23) * g, B = 0x4a + (0x1a - 0x4a) * g;
            float dmin = 1e9f;
            for (auto &l : L)
                dmin = std::min(dmin, seg(X, Y, l[0], l[1], l[2], l[3]));
            if (dmin <= 32)
                R = 0xf0, G = 0xa2, B = 0x4b;
            else if (seg(X, Y, 388, 604, 474, 546) <= 20)
                R = R * 0.4f + 0xf0 * 0.6f, G = G * 0.4f + 0xa2 * 0.6f, B = B * 0.4f + 0x4b * 0.6f;
            if ((X - 530) * (X - 530) + (Y - 508) * (Y - 508) <= 62 * 62)
                R = G = B = 255;
            Uint8 *p = px + y * sf->pitch + x * 4;
            p[0] = Uint8(R), p[1] = Uint8(G), p[2] = Uint8(B), p[3] = Uint8(255 * inside);
        }
    SDL_SetWindowIcon(win, sf);
    SDL_FreeSurface(sf);
}

// Without a working OpenGL 3 driver the app cannot draw: say so instead of crashing.
static int gl_fail(SDL_Window *win, const char *what) {
    std::string m = std::string(what) + ": " + SDL_GetError() +
                    "\n\nPiazzi needs OpenGL 3.0 or newer. Update the graphics driver, or use the web version.";
    std::fprintf(stderr, "%s\n", m.c_str());
    SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "Piazzi", m.c_str(), win);
    if (win)
        SDL_DestroyWindow(win);
    SDL_Quit();
    return 1;
}

int main(int, char **) {
    SDL_SetHint(SDL_HINT_VIDEO_HIGHDPI_DISABLED, "0");
    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) != 0) {
        std::fprintf(stderr, "SDL: %s\n", SDL_GetError());
        return 1;
    }
#if defined(__APPLE__)
    const char *glsl = "#version 150";
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);
#else
    const char *glsl = "#version 130";
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
#endif
    SDL_Window *win = SDL_CreateWindow("Piazzi", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1440, 900,
                                       SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
    if (!win)
        return gl_fail(nullptr, "Could not create the window");
    set_window_icon(win);
    SDL_GLContext gl = SDL_GL_CreateContext(win);
    if (!gl)
        return gl_fail(win, "Could not create an OpenGL context");
    SDL_GL_MakeCurrent(win, gl);
    SDL_GL_SetSwapInterval(1);
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO &io = ImGui::GetIO();
    io.IniFilename = nullptr;
    std::string fd = assets_dir() + "fonts/";
    T::sans = io.Fonts->AddFontFromFileTTF((fd + "mono.ttf").c_str(), 15.0f); // general UI in mono (instrument look)
    T::mono = T::sans;
    T::serif = io.Fonts->AddFontFromFileTTF((fd + "serif.ttf").c_str(), 22.0f);
    if (!T::sans)
        T::sans = T::mono = io.Fonts->AddFontDefault();
    if (!T::serif)
        T::serif = T::sans;
    io.FontDefault = T::sans;
    apply_theme();
    ImGui_ImplSDL2_InitForOpenGL(win, gl);
    if (!ImGui_ImplOpenGL3_Init(glsl))
        return gl_fail(win, "Could not initialise OpenGL");
    SDL_EventState(SDL_DROPFILE, SDL_ENABLE);
    load_lang();
    load_profiles();
    std::vector<std::string> dropped;
    Uint32 last_drop = 0;
    bool running = true;
    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            ImGui_ImplSDL2_ProcessEvent(&e);
            if (e.type == SDL_QUIT)
                running = false;
            if (e.type == SDL_DROPFILE) {
                dropped.push_back(e.drop.file);
                SDL_free(e.drop.file);
                last_drop = SDL_GetTicks();
            }
        }
        // multiple files arrive as separate events: wait for 250 ms of silence before opening the session
        if (!dropped.empty() && SDL_GetTicks() - last_drop > 250 && A.phase == Phase::Home && !BR.open && !A.busy) {
            start_with_files(dropped);
            dropped.clear();
        }
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplSDL2_NewFrame();
        ImGui::NewFrame();
        ImGuiViewport *vp = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(vp->WorkPos);
        ImGui::SetNextWindowSize(vp->WorkSize);
        ImGui::Begin("##root", nullptr,
                     ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoBringToFrontOnFocus |
                         ImGuiWindowFlags_NoSavedSettings);
        ImDrawList *bg = ImGui::GetWindowDrawList();
        if (A.phase != Phase::Results)
            star_chart(bg, vp->WorkPos, vp->WorkSize, 7);
        Phase ph;
        {
            std::lock_guard<std::recursive_mutex> l(A.mu);
            ph = A.phase;
        }
        if (ph == Phase::Home)
            ui_home();
        else if (ph == Phase::Session)
            ui_session();
        else
            ui_results();
        ui_center();
        ui_mpc();
        BR.draw();
        ImGui::End();
        ui_settings();
        ImGui::Render();
        int w, h;
        SDL_GL_GetDrawableSize(win, &w, &h);
        glViewport(0, 0, w, h);
        glClearColor(T::ink.x, T::ink.y, T::ink.z, 1);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        SDL_GL_SwapWindow(win);
    }
    // shutdown: stop the worker thread and the detached ones (curl is killed within ~50 ms) and wait for them
    A.quit = true;
    A.cancel = true;
    if (A.worker.joinable())
        A.worker.join();
    for (int i = 0; i < 300 && A.live_threads > 0; ++i)
        SDL_Delay(10);
    // textures before destroying the GL context; the session is freed only if no thread can still use it
    for (GLuint t : A.tex)
        glDeleteTextures(1, &t);
    A.tex.clear();
    if (A.live_threads == 0) {
        astro_free(A.s);
        A.s = nullptr;
    }
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplSDL2_Shutdown();
    ImGui::DestroyContext();
    SDL_GL_DeleteContext(gl);
    SDL_DestroyWindow(win);
    SDL_Quit();
    // if a detached thread is still running after 3 s, exit without static destructors (it would use a destroyed A)
    if (A.live_threads > 0)
        std::_Exit(0);
    return 0;
}
