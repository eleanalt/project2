#include "OutputGenerator.h"
#include <fstream>
#include <iostream>

OutputGenerator::OutputGenerator(const std::string& outputFilePath)
    : outputFilePath(outputFilePath) {}

void OutputGenerator::write(const TriangulationResult& result) {
    nlohmann::json outputJson;

    outputJson["content_type"] = result.content_type;
    outputJson["instance_uid"] = result.instance_uid;
    outputJson["steiner_points_x"] = result.steiner_points_x;
    outputJson["steiner_points_y"] = result.steiner_points_y;
    outputJson["edges"] = result.edges;
    outputJson["obtuse_count"] = result.obtuse_count;
    outputJson["method"] = result.method;
    outputJson["parameters"] = result.parameters;

    std::ofstream outputFile(outputFilePath);
    if (!outputFile.is_open()) {
        throw std::runtime_error("Could not open output file: " + outputFilePath);
    }

    outputFile << outputJson.dump(4);
    outputFile.close();

    std::cout << "Output successfully written to: " << outputFilePath << "\n";
}
