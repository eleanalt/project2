#pragma once

#include <string>
#include <vector>
#include <utility>
#include <nlohmann/json.hpp>

struct TriangulationResult {
    std::string content_type = "CG_SHOP_2025_Solution";
    std::string instance_uid;
    std::vector<std::string> steiner_points_x;
    std::vector<std::string> steiner_points_y;
    std::vector<std::pair<int, int>> edges;
    int obtuse_count = 0;
    std::string method;
    nlohmann::json parameters;
};

