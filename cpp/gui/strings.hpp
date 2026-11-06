// Desktop app strings: English by default, Portuguese optional. The key is the English text itself.
#pragma once
#include <string>

// text in the current language (returns the English itself if there is no translation)
const char *tr(const char *en);
// current language: "en" or "pt"
const std::string &lang();
void set_lang(const std::string &l); // saved to prefs_dir()/lang.txt
void load_lang();
bool lang_pt();
// core error message (JSON {"error","code","detail"} from astro_result_json) in the current language
std::string core_error_message(const std::string &json);
