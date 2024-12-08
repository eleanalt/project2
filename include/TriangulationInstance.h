#ifndef TRIANGULATIONINSTANCE_H
#define TRIANGULATIONINSTANCE_H

#include <vector>
#include <string>
#include <nlohmann/json.hpp>  // Για την χρήση του nlohmann::json

struct TriangulationInstance {
    std::string instance_uid;
    int num_points;
    std::vector<int> points_x;
    std::vector<int> points_y;
    std::vector<int> region_boundary;
    int num_constraints;
    std::vector<std::pair<int, int>> additional_constraints;
    std::string method;  // Παράδειγμα για τη μέθοδο
    nlohmann::json parameters;  // Χρησιμοποιούμε το nlohmann::json για τις παραμέτρους
    bool delaunay;  // Αν η τριγωνοποίηση αρχίζει από την Delaunay τριγωνοποίηση
};

#endif // TRIANGULATIONINSTANCE_H

