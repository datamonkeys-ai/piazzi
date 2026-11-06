// Process execution WITHOUT a shell, and HTTP via the system curl.
// Nothing here goes through /bin/sh or cmd.exe: each argument is passed separately to the child process, so
// e-mails, paths and URLs typed by the user cannot inject commands.
#include "app.hpp"
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <thread>
#include <vector>
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#elif !defined(__EMSCRIPTEN__)
#include <signal.h>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
extern char **environ;
#endif

#ifdef _WIN32
// Quotes an argument per the CommandLineToArgvW / CRT rules: backslashes are special only before
// quotes; 2n backslashes + quote -> n backslashes + delimiting quote, 2n+1 -> n backslashes + literal quote.
static std::wstring quote_arg_windows(const std::wstring &a) {
    if (!a.empty() && a.find_first_of(L" \t\n\v\"") == std::wstring::npos)
        return a;
    std::wstring r = L"\"";
    for (size_t i = 0;; ++i) {
        size_t nbs = 0;
        while (i < a.size() && a[i] == L'\\') {
            ++i;
            ++nbs;
        }
        if (i == a.size()) {
            r.append(nbs * 2, L'\\'); // double the backslashes before the closing quote
            break;
        }
        if (a[i] == L'"')
            r.append(nbs * 2 + 1, L'\\'); // escape the backslashes and the quote
        else
            r.append(nbs, L'\\');
        r.push_back(a[i]);
    }
    r.push_back(L'"');
    return r;
}
static std::wstring widen(const std::string &s) {
    if (s.empty())
        return {};
    int n = MultiByteToWideChar(CP_UTF8, 0, s.data(), (int)s.size(), nullptr, 0);
    std::wstring w(n, L'\0');
    MultiByteToWideChar(CP_UTF8, 0, s.data(), (int)s.size(), w.data(), n);
    return w;
}
// Absolute path of a system executable: curl -> %SystemRoot%\\System32\\curl.exe (GetSystemDirectoryW),
// explorer -> %SystemRoot%\\explorer.exe (GetWindowsDirectoryW). Other names: empty (not executed).
static std::wstring system_exe(const std::string &name) {
    wchar_t dir[MAX_PATH];
    UINT n = 0;
    std::wstring file;
    if (name == "curl") {
        n = GetSystemDirectoryW(dir, MAX_PATH);
        file = L"curl.exe";
    } else if (name == "explorer") {
        n = GetWindowsDirectoryW(dir, MAX_PATH);
        file = L"explorer.exe";
    }
    if (n == 0 || n >= MAX_PATH)
        return {};
    std::wstring r(dir, n);
    if (r.back() != L'\\')
        r += L'\\';
    return r + file;
}
#endif

// Cancellation requested (session cancelled or app shutting down)?
static bool cancel_requested(const std::atomic<bool> *cancel) {
    return (cancel && cancel->load()) || A.quit.load();
}

int run_process(const std::vector<std::string> &argv, const std::atomic<bool> *cancel) {
    if (argv.empty())
        return -1;
#if defined(_WIN32)
    std::wstring cmd;
    for (const auto &a : argv) {
        if (!cmd.empty())
            cmd += L' ';
        cmd += quote_arg_windows(widen(a));
    }
    STARTUPINFOW si{};
    PROCESS_INFORMATION pi{};
    si.cb = sizeof si;
    // absolute executable path (against "binary planting": never look up curl.exe/explorer.exe in the current folder
    // or on PATH); curl and explorer ship with Windows. No console window.
    std::wstring app = system_exe(argv[0]);
    if (app.empty())
        return -1;
    if (!CreateProcessW(app.c_str(), cmd.data(), nullptr, nullptr, FALSE, CREATE_NO_WINDOW, nullptr, nullptr, &si, &pi))
        return -1;
    // wait in short slices to honor cancellation: terminate the child and return -1
    while (WaitForSingleObject(pi.hProcess, 50) == WAIT_TIMEOUT) {
        if (cancel_requested(cancel)) {
            TerminateProcess(pi.hProcess, 1);
            WaitForSingleObject(pi.hProcess, 2000);
            CloseHandle(pi.hProcess);
            CloseHandle(pi.hThread);
            return -1;
        }
    }
    DWORD code = 1;
    GetExitCodeProcess(pi.hProcess, &code);
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);
    return (int)code;
#elif defined(__EMSCRIPTEN__)
    return -1; // no processes in the browser
#else
    std::vector<char *> args;
    args.reserve(argv.size() + 1);
    for (const auto &a : argv)
        args.push_back(const_cast<char *>(a.c_str()));
    args.push_back(nullptr);
    pid_t pid = 0;
    if (posix_spawnp(&pid, args[0], nullptr, nullptr, args.data(), environ) != 0)
        return -1;
    // poll with WNOHANG to honor cancellation: SIGTERM the child, reap the zombie and return -1
    int st = 0;
    for (;;) {
        pid_t r = waitpid(pid, &st, WNOHANG);
        if (r == pid)
            break;
        if (r < 0) {
            if (errno == EINTR)
                continue;
            return -1;
        }
        if (cancel_requested(cancel)) {
            kill(pid, SIGTERM);
            // up to ~1 s to exit on its own; then SIGKILL (cannot be ignored) and reap the zombie
            for (int i = 0; i < 50 && waitpid(pid, &st, WNOHANG) == 0; ++i)
                std::this_thread::sleep_for(std::chrono::milliseconds(20));
            if (waitpid(pid, &st, WNOHANG) == 0) {
                kill(pid, SIGKILL);
                while (waitpid(pid, &st, 0) < 0 && errno == EINTR) {
                }
            }
            return -1;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
    }
    return WIFEXITED(st) ? WEXITSTATUS(st) : -1;
#endif
}

void open_folder(const std::string &dir) {
#if defined(_WIN32)
    std::vector<std::string> argv{"explorer", dir};
#elif defined(__APPLE__)
    std::vector<std::string> argv{"open", dir};
#else
    std::vector<std::string> argv{"xdg-open", dir};
#endif
    // in the background: some file managers only return when the window is closed
    // (counts as a live thread: the end of main waits for it before destroying the global state)
    A.live_threads++;
    std::thread([argv] {
        run_process(argv);
        A.live_threads--;
    }).detach();
}

// Unique temp file name per call (pid + atomic counter + thread id): concurrent requests to the same URL
// do not race on the same file.
static std::string unique_tmp() {
    static std::atomic<unsigned long> counter{0};
#ifdef _WIN32
    unsigned long pid = (unsigned long)GetCurrentProcessId();
#elif defined(__EMSCRIPTEN__)
    unsigned long pid = 0;
#else
    unsigned long pid = (unsigned long)getpid();
#endif
    size_t tid = std::hash<std::thread::id>{}(std::this_thread::get_id());
    return prefs_dir() + "http_" + std::to_string(pid) + "_" + std::to_string(counter++) + "_" + std::to_string(tid) + ".tmp";
}

bool http_get(const std::string &url, std::string &out, const std::atomic<bool> *cancel, int timeout, int tries, size_t minbytes) {
    std::string tmp = unique_tmp();
    for (int t = 0; t < tries; ++t) {
        if (cancel_requested(cancel))
            return false;
        // HTTPS only, including redirects
        // -q must be the FIRST argument: ignores the user's .curlrc (which could change proxy, protocols, etc.)
        int rc = run_process({"curl", "-q", "-s", "-L", "-f", "--proto", "=https", "--proto-redir", "=https", "--max-time",
                              std::to_string(timeout), "-A", "Piazzi-desktop", "-o", tmp, url},
                             cancel);
        if (rc == 0) {
            out = read_file(tmp);
            std::error_code ec;
            fs::remove(tmp, ec);
            if (out.size() >= minbytes)
                return true;
        }
        std::error_code ec;
        fs::remove(tmp, ec);
        if (t + 1 == tries)
            break; // no wait after the last attempt
        // wait between attempts in 50 ms slices, honoring cancellation
        for (int w = 0; w < 30 * (t + 1); ++w) {
            if (cancel_requested(cancel))
                return false;
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }
    }
    return false;
}
