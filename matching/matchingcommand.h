//
// Created by Shixuan Sun on 2018/6/29.
//

#ifndef SUBGRAPHMATCHING_MATCHINGCOMMAND_H
#define SUBGRAPHMATCHING_MATCHINGCOMMAND_H

#include "utility/commandparser.h"
#include <map>
#include <iostream>
enum OptionKeyword {
    QueryGraphFile = 0,     // -q, The query graph file path, compulsive parameter
    DataGraphFile = 1,      // -d, The data graph file path, compulsive parameter
    Dist = 2,
    QueryNums = 3,
    Methe = 4,
    Thread = 5
};

class MatchingCommand : public CommandParser{
private:
    std::map<OptionKeyword, std::string> options_key;
    std::map<OptionKeyword, std::string> options_value;

private:
    void processOptions();

public:
    MatchingCommand(int argc, char **argv);

    std::string getDataGraphFilePath() {
        return options_value[OptionKeyword::DataGraphFile];
    }

    std::string getQueryGraphFilePath() {
        return options_value[OptionKeyword::QueryGraphFile];
    }

    std::string getDist() {
        return options_value[OptionKeyword::Dist];
    }

    std::string getQueryNums() {
        return options_value[OptionKeyword::QueryNums];
    }

    std::string getMeth() {
        return options_value[OptionKeyword::Methe];
    }

    std::string getThread() {
        return options_value[OptionKeyword::Thread];
    }
};


#endif //SUBGRAPHMATCHING_MATCHINGCOMMAND_H
