// Submission of the ADES report to the Minor Planet Center (submit_psv) via curl, without a shell.
#include "app.hpp"
#include "imgui.h"
#include <cctype>
#include <fstream>
#include <string>
#include <vector>

// simple name@domain.tld format: no spaces, quotes or control characters; a single @; domain contains a dot
bool valid_email(const std::string &e) {
    auto at = e.find('@');
    if (at == std::string::npos || at == 0 || e.find('@', at + 1) != std::string::npos || e.size() > 254)
        return false;
    std::string local = e.substr(0, at), dom = e.substr(at + 1);
    for (unsigned char c : local)
        if (!std::isalnum(c) && std::string("._%+-").find((char)c) == std::string::npos)
            return false;
    for (unsigned char c : dom)
        if (!std::isalnum(c) && c != '.' && c != '-')
            return false;
    auto dot = dom.rfind('.');
    if (dot == std::string::npos || dot == 0 || dom.size() - dot - 1 < 2 || dom.find("..") != std::string::npos || dom.front() == '-')
        return false;
    for (size_t i = dot + 1; i < dom.size(); ++i)
        if (!std::isalpha((unsigned char)dom[i]))
            return false;
    return local.front() != '.' && local.back() != '.';
}

void ui_mpc() {
    if (A.show_mpc) {
        ImGui::OpenPopup("##mpc");
        A.show_mpc = false;
    }
    ImGui::SetNextWindowSize(ImVec2(560, 0), ImGuiCond_Appearing);
    if (!ImGui::BeginPopupModal("##mpc", nullptr, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_AlwaysAutoResize))
        return;
    static char email[128] = "";
    static bool test = true;
    static std::string res;
    tag(tr("Submission to the Minor Planet Center"), T::amber);
    title(tr("Submit the ADES report?"), 24);
    // built-in example: the data is not the user's — submission is blocked
    bool example = A.example;
    if (example) {
        tag(tr("Example data — not your observation"), T::red);
        ImGui::TextWrapped("%s", tr("These measurements come from the bundled ZTF example and must not be sent to the MPC. Open your "
                                    "own images to submit."));
    } else
        ImGui::TextWrapped("%s", tr("The file ADESReport.psv will be sent to the MPC through the submit_psv interface. Review the "
                                    "measurements first. In test mode the MPC only validates and answers, without recording."));
    ImGui::SetNextItemWidth(360);
    ImGui::BeginDisabled(example);
    ImGui::InputTextWithHint(tr("e-mail for confirmation"), tr("you@example.com"), email, sizeof email);
    bool ok = valid_email(email);
    if (email[0] && !ok)
        tag(tr("invalid e-mail"), T::red);
    ImGui::Checkbox(tr("Test mode (recommended the first time)"), &test);
    ImGui::EndDisabled();
    if (!res.empty()) {
        ImGui::PushFont(T::mono, 12.0f);
        ImGui::TextWrapped("%s", res.c_str());
        ImGui::PopFont();
    }
    if (ImGui::Button(tr("Close"), ImVec2(120, 0))) {
        res.clear();
        ImGui::CloseCurrentPopup();
    }
    ImGui::SameLine();
    ImGui::BeginDisabled(!ok || example);
    bool send = primary_button(test ? tr("SUBMIT (TEST)") : tr("SUBMIT FOR REAL"), 280);
    ImGui::EndDisabled();
    if (send && ok && !example) {
        std::string d = prefs_dir(), f = d + "ADESReport.psv", out = d + "mpc_resp.txt";
        std::ofstream(f, std::ios::binary) << astro_report(A.s, 0);
        std::string url = std::string("https://minorplanetcenter.net/submit_psv") + (test ? "_test" : "");
        // --form-string: the value is sent literally (curl does not interpret @file, <file or ;type=)
        int rc = run_process({"curl", "-q", "-s", "--max-time", "60", "-o", out, "--form-string", "ack=Piazzi " + get("mpc_code"),
                              "--form-string", std::string("ac2=") + email, "--form-string", "obj_type=Unclassified", "-F",
                              "source=@\"" + f + "\"", url});
        res = rc == 0 ? tr("MPC answer: ") + read_file(out).substr(0, 600) : std::string(tr("Submission failed (no internet?)."));
    }
    ImGui::EndPopup();
}
