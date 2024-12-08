#ifndef INPUTPARSER_H
#define INPUTPARSER_H

#include <string>
#include "TriangulationInstance.h"

class InputParser {
public:
    explicit InputParser(const std::string& filePath);
    TriangulationInstance parse();

private:
    std::string filePath;
};

#endif // INPUTPARSER_H
