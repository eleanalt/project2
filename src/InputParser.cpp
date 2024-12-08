#include "InputParser.h"
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

InputParser::InputParser(const std::string& filePath) : filePath(filePath) {}

TriangulationInstance InputParser::parse() {
    std::ifstream inputFile(filePath);
    if (!inputFile.is_open()) {
        throw std::runtime_error("Could not open input file: " + filePath);
    }

    json inputJson;
    inputFile >> inputJson;

    TriangulationInstance instance;

    if (inputJson.contains("instance_uid") && inputJson["instance_uid"].is_string()) {
        instance.instance_uid = inputJson["instance_uid"];
    } else {
        throw std::runtime_error("Invalid or missing 'instance_uid'");
    }

    if (inputJson.contains("num_points") && inputJson["num_points"].is_number_integer()) {
        instance.num_points = inputJson["num_points"];
    } else {
        throw std::runtime_error("Invalid or missing 'num_points'");
    }

    if (inputJson.contains("points_x") && inputJson["points_x"].is_array()) {
        instance.points_x = inputJson["points_x"].get<std::vector<int>>();
    } else {
        throw std::runtime_error("Invalid or missing 'points_x'");
    }

    if (inputJson.contains("points_y") && inputJson["points_y"].is_array()) {
        instance.points_y = inputJson["points_y"].get<std::vector<int>>();
    } else {
        throw std::runtime_error("Invalid or missing 'points_y'");
    }

    if (inputJson.contains("region_boundary") && inputJson["region_boundary"].is_array()) {
        instance.region_boundary = inputJson["region_boundary"].get<std::vector<int>>();
    } else {
        throw std::runtime_error("Invalid or missing 'region_boundary'");
    }

    if (inputJson.contains("num_constraints") && inputJson["num_constraints"].is_number_integer()) {
        instance.num_constraints = inputJson["num_constraints"];
    } else {
        throw std::runtime_error("Invalid or missing 'num_constraints'");
    }

    if (inputJson.contains("additional_constraints") && inputJson["additional_constraints"].is_array()) {
        for (const auto& constraint : inputJson["additional_constraints"]) {
            instance.additional_constraints.emplace_back(constraint[0], constraint[1]);
        }
    }

    if (inputJson.contains("method") && inputJson["method"].is_string()) {
        instance.method = inputJson["method"];
    } else {
        throw std::runtime_error("Invalid or missing 'method'");
    }

    if (inputJson.contains("parameters") && inputJson["parameters"].is_object()) {
        instance.parameters = inputJson["parameters"];
    } else {
        instance.parameters = {}; 
    }

    if (inputJson.contains("delaunay") && inputJson["delaunay"].is_boolean()) {
        instance.delaunay = inputJson["delaunay"];
    } else {
        instance.delaunay = true; 
    }

    return instance;
}

