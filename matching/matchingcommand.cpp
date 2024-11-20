//
// Created by Shixuan Sun on 2018/6/29.
//

#include "matchingcommand.h"

MatchingCommand::MatchingCommand(const int argc, char **argv) : CommandParser(argc, argv) {
    // Initialize options value
    options_key[OptionKeyword::QueryGraphFile] = "-q";
    options_key[OptionKeyword::DataGraphFile] = "-d";
    options_key[OptionKeyword::Dist] = "-D";
    options_key[OptionKeyword::QueryNums] = "-Q";
    options_key[OptionKeyword::Methe] = "-M";
    processOptions();
};

void MatchingCommand::processOptions() {
    // Query graph file path
    options_value[OptionKeyword::QueryGraphFile] = getCommandOption(options_key[OptionKeyword::QueryGraphFile]);

    // Data graph file path
    options_value[OptionKeyword::DataGraphFile] = getCommandOption(options_key[OptionKeyword::DataGraphFile]);

    // D
    options_value[OptionKeyword::Dist] = getCommandOption(options_key[OptionKeyword::Dist]);

    // QueryNums
    options_value[OptionKeyword::QueryNums] = getCommandOption(options_key[OptionKeyword::QueryNums]);

    options_value[OptionKeyword::Methe] = getCommandOption(options_key[OptionKeyword::Methe]);
}