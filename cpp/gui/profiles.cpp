// Session settings and equipment profiles (plain file in the preferences folder).
#include "app.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

void set(const std::string &k, const std::string &v) {
    astro_set(A.s, k.c_str(), v.c_str());
}
void setn(const std::string &k, double v) {
    char b[64];
    std::snprintf(b, sizeof b, "%.10g", v);
    set(k, b);
}
std::string get(const std::string &k) {
    return astro_get(A.s, k.c_str());
}

// built-in profiles (none for now) + user profiles
std::vector<Profile> all_profiles() {
    return A.profiles;
}

void save_profiles() {
    std::ofstream f(prefs_dir() + "profiles.txt");
    f << "current\t" << A.profile << "\n";
    for (auto &p : A.profiles)
        f << p.name << '\t' << p.focal << '\t' << p.pw << '\t' << p.ph << '\t' << p.flip << '\t' << p.mpc << '\n';
}
void apply_profile(const Profile &p) {
    setn("focal_length_mm", p.focal);
    setn("pixel_width_mm", p.pw);
    setn("pixel_height_mm", p.ph);
    setn("flip_horizontal", p.flip);
    set("mpc_code", p.mpc);
    A.profile = p.name;
    save_profiles();
}
void load_profiles() {
    A.profiles.clear();
    std::ifstream f(prefs_dir() + "profiles.txt");
    std::string line, cur;
    while (std::getline(f, line)) {
        std::stringstream ss(line);
        std::vector<std::string> c;
        std::string x;
        while (std::getline(ss, x, '\t'))
            c.push_back(x);
        if (c.size() == 2 && c[0] == "current")
            cur = c[1];
        else if (c.size() >= 6) {
            try {
                A.profiles.push_back({c[0], std::stod(c[1]), std::stod(c[2]), std::stod(c[3]), std::stoi(c[4]), c[5]});
            } catch (...) { // corrupt line: skip it instead of crashing the app at startup
            }
        }
    }
    for (auto &p : all_profiles())
        if (p.name == cur)
            apply_profile(p);
}
