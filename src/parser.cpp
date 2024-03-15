#include "parser.h"
#include <iostream>
#include <sstream>
#include <algorithm>

Settings::Settings(const char* filename) {
    input_file.open(filename);
    if (!input_file) {
        throw SettingsParserException("Error opening settings file.");
    }
}

Settings::~Settings() {
    input_file.close();
}

Settings::option::option(std::string n, std::string v) {
    name = n;
    if (v == "true" || v == "True") {
        type = 'b';
        value_d = 1;
    } else if (v == "false" || v == "False") {
        type = 'b';
        value_d = 0;
    } else { 
        try {
            value_d = std::stod(v);
            type = 'd';
        } catch (...) {
            v.erase(remove(v.begin(), v.end(), '\"'), v.end());
            value_s = v;
            type = 's';
        }
    }
}

bool Settings::ParseLine(std::string line) {
    std::istringstream lineStream(line);
    std::ostringstream optionName;
    std::ostringstream optionValue;
    char c;
    bool split = false;

    while (lineStream >> c) {
        if (c == '#')
            break;
        else if (c == '=') {
            split = true;
        } else {
            if (split)
                optionValue << c;
            else
                optionName << c;
        }
    }
    if (optionName.tellp() != std::streampos(0) || optionValue.tellp() != std::streampos(0)) {
        settings.push_back(option(optionName.str(), optionValue.str()));
        return true;
    }
    return false; 
}

int Settings::LoadSettings() {
    std::string lineBuffer;
    int counter = 0;
    while (input_file) {
        getline(input_file, lineBuffer);
        if (ParseLine(lineBuffer)) {
            counter++;
        }
    }
    return counter;
}

double Settings::GetOptionDouble(std::string name) {
    std::string error = name + " not found.";
    if (settings.empty()) {
        throw SettingsParserException (error.c_str());
    }
    for (option setting : settings) {
        if (setting.name == name) {
            if (setting.type != 'd') {
                throw SettingsParserException ("Invalid Option Type.");
            }
            return setting.value_d;
        }
    }
    throw SettingsParserException (error.c_str());
}


std::string Settings::GetOptionString(std::string name) {
    if (settings.empty()) {
        throw SettingsParserException ("Option Not Found.");
    }
    for (option setting : settings) {
        if (setting.name == name) {
            if (setting.type != 's') {
                throw SettingsParserException ("Invalid Option Type.");
            }
            return setting.value_s;
        }
    }
    throw SettingsParserException ("Option Not Found.");
}

bool Settings::GetOptionBool(std::string name) {
    if (settings.empty()) {
        throw SettingsParserException ("Option Not Found.");
    }
    for (option setting : settings) {
        if (setting.name == name) {
            if (setting.type != 'b') {
                throw SettingsParserException ("Invalid Option Type.");
            }
            return setting.value_d;
        }
    }
    throw SettingsParserException ("Option Not Found.");
}

