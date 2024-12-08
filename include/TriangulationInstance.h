#ifndef TRIANGULATIONINSTANCE_H
#define TRIANGULATIONINSTANCE_H

#include <vector>
#include <string>
#include <nlohmann/json.hpp>  // ��� ��� ����� ��� nlohmann::json

struct TriangulationInstance {
    std::string instance_uid;
    int num_points;
    std::vector<int> points_x;
    std::vector<int> points_y;
    std::vector<int> region_boundary;
    int num_constraints;
    std::vector<std::pair<int, int>> additional_constraints;
    std::string method;  // ���������� ��� �� ������
    nlohmann::json parameters;  // �������������� �� nlohmann::json ��� ��� �����������
    bool delaunay;  // �� � ������������� ������� ��� ��� Delaunay �������������
};

#endif // TRIANGULATIONINSTANCE_H

