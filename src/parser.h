#ifndef PARSER_H
#define PARSER
#include <iostream>
#include <vector>
#include <fstream>


class SettingsParserException : public std::exception {
/**
 * @brief Custom exception. DO NOT CHANGE.
 * 
 * This will be thrown whenever the parser encounters an error.
 */
private:
    std::string message;
public:
    SettingsParserException(const char* msg) : message(msg) {}
    const char* what() const throw() {
        return message.c_str();
    }
};

class Settings {
private:
    struct option {
        char type;
        std::string name;
        double value_d;
        std::string value_s;
        option(std::string n, double v) : name(n), value_d(0), value_s("") {}
        option(std::string n, std::string v);
        option() : name(""), value_d(0), value_s("") {}
    };
    std::vector<option> settings;
    std::ifstream input_file;
    bool ParseLine(std::string line);
public:
    /**
     * @brief Construct a new Settings object
     * 
     * @param filename Name of the configuration file.
     */
    Settings(const char* filename);
    ~Settings();

    /**
     * @brief Run this after the constructor to load settings from the file.
     * 
     * @return int Number of settings loaded.
     */
    int LoadSettings();

    /**
     * @brief Get the option as a double. Throws exception if value is not double.
     * 
     * @param name Name of the option to get.
     * @return double Value of the option.
     * 
     */
    double GetOptionDouble(std::string name);
    /**
     * @brief Get the option as a string. Throws exception if value is not string.
     * 
     * @param name Name of the option to get.
     * @return string Value of the option.
     * 
     */
    std::string GetOptionString(std::string name);
    /**
     * @brief Get the option as a boolean. Throws exception if value is not boolean.
     * 
     * @param name Name of the option to get.
     * @return bool Value of the option.
     * 
     */
    bool GetOptionBool(std::string name);
};

#endif