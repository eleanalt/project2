#include "TriangulationOptimizer.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>
#include <random>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <nlohmann/json.hpp>

using namespace std;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef CDT::Point Point;

// ��������� ��� �� ���������� ���� ������ ������ ������ ��� �����
double randomDouble(double min, double max) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}
// ��������� ��� ������� ����� ������ ���� ������������� ����� a�������
int TriangulationOptimizer::countObtuseAngles(const CDT& cdt) {
    int obtuse_count = 0;
    
    // ��������� ���� ��� ����� (faces) ��� ��������������
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        if (cdt.is_infinite(face)) { // �� � ��� ����� ������, ������������
            continue;
        }

        // ������ ��� ����� ������� ��� ��������
        Point p1 = face->vertex(0)->point();
        Point p2 = face->vertex(1)->point();
        Point p3 = face->vertex(2)->point();

        // ����������� ����������� �����
        auto v1 = p2 - p1;
        auto v2 = p3 - p2;
        auto v3 = p1 - p3;

        // ����������� ������ �� ��������� ��������
        double cos_angle1 = CGAL::to_double(v1 * v3) / (CGAL::sqrt(CGAL::squared_distance(p1, p2)) * CGAL::sqrt(CGAL::squared_distance(p1, p3)));
        double cos_angle2 = CGAL::to_double(v2 * (-v1)) / (CGAL::sqrt(CGAL::squared_distance(p2, p3)) * CGAL::sqrt(CGAL::squared_distance(p2, p1)));
        double cos_angle3 = CGAL::to_double((-v2) * (-v3)) / (CGAL::sqrt(CGAL::squared_distance(p3, p1)) * CGAL::sqrt(CGAL::squared_distance(p3, p2)));

        // ������� �� ������ ����� ����� ���������� ��� 90 ������
        if (cos_angle1 < 0 || cos_angle2 < 0 || cos_angle3 < 0) {
            ++obtuse_count;
        }
    }

    return obtuse_count;
}
// ����������� ��� ��������������� ������ (circumradius) ��� ��� ��� ��� ��������������
double TriangulationOptimizer::calculateCircumradius(CDT::Face_handle face) {
    auto p0 = face->vertex(0)->point();
    auto p1 = face->vertex(1)->point();
    auto p2 = face->vertex(2)->point();

    // ����������� ��� ������� ��� ��������
    double a = CGAL::squared_distance(p0, p1);
    double b = CGAL::squared_distance(p1, p2);
    double c = CGAL::squared_distance(p2, p0);

    // ����������� ��� ��������������� ������
    double s = (a + b + c) / 2.0; // ����������� ��������� ��� �� �������
    double area = sqrt(s * (s - a) * (s - b) * (s - c)); // ������� ��������
    double circumradius = (a * b * c) / (4.0 * area); // ����������� ��� ��������������� ������

    return circumradius;
}

// ������ ��� ������������ ������� ��� CDT ��� ��� �������� ������
Point TriangulationOptimizer::findNearestVertex(CDT& cdt, const Point& point) {
    double min_distance = std::numeric_limits<double>::infinity();// ������������ �� ��� ���������� ������ ������
    Point nearest_point;

    // ��������� ���� ��� ������� ��� CDT ��� �� ���� ��� �����������
    for (auto vertex = cdt.finite_vertices_begin(); vertex != cdt.finite_vertices_end(); ++vertex) {
        double distance = CGAL::squared_distance(vertex->point(), point); // ����������� ���������
        if (distance < min_distance) {
            min_distance = distance; // �� � �������� ����� ���������, ������������ �� ����������� ������
            nearest_point = vertex->point(); 
        }
    }
    return nearest_point;
}


//��������� ��� �������������� ��������������
TriangulationResult TriangulationOptimizer::optimize(const TriangulationInstance& instance) {
    TriangulationResult result;
    result.instance_uid = instance.instance_uid;// ���������� ��� ��������� �������������� ��� ��������������
    result.method = instance.method;//���������� ������� ��� ���������������
    result.parameters = instance.parameters;  //���������� ����������

    // ���������� ��� �������������� �� ��� ����� ��� PSLG
    CDT cdt;
    std::vector<Point> points;

    // �������� ��� ������� ��� �� PSLG
    for (size_t i = 0; i < instance.num_points; ++i) {
        points.emplace_back(instance.points_x[i], instance.points_y[i]);
        cdt.insert(points.back());  // �������� ������� ���� �������������
    }

    // �������� ��� �����������
    for (size_t i = 0; i < instance.region_boundary.size(); ++i) {
        size_t next = (i + 1) % instance.region_boundary.size();
        cdt.insert_constraint(points[instance.region_boundary[i]], points[instance.region_boundary[next]]);
    }

    // �������� ��� �������� �����������
    for (const auto& constraint : instance.additional_constraints) {
        cdt.insert_constraint(points[constraint.first], points[constraint.second]);
    }

    // ������������ ��� ������� Steiner ��� ��� �����
    result.steiner_points_x.clear();
    result.steiner_points_y.clear();
    result.edges.clear();

    // ������� ��� ����� ��� ��������������
    for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge) {
        auto segment = cdt.segment(*edge);
        int source_index = std::distance(points.begin(), std::find(points.begin(), points.end(), segment.source()));
        int target_index = std::distance(points.begin(), std::find(points.begin(), points.end(), segment.target()));
        result.edges.emplace_back(source_index, target_index);
    }

    // ����������� ��� ������� ��� ����������� ��������
    result.obtuse_count = 0;
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        if (isObtuse(face)) {
            result.obtuse_count++;
        }
    }

    // �������� ��� ����������� ������� ���������������
    if (instance.method == "local") {
        localSearch(cdt, instance,result);
    } else if (instance.method == "sa") {
        simulatedAnnealing(cdt, instance,result);
    } else if (instance.method == "ant") {
        antColonyOptimization(cdt, instance,result);
    }

    return result;
}

// ����������� ������� ��������������
double TriangulationOptimizer::evaluateTriangulation(CDT& cdt) {
    int obtuse_count = 0;
    //������� ��� ���� �� ����� ��� �� �������� ���� ���������� ��������
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        if (isObtuse(face)) {
            ++obtuse_count;
        }
    }

    return static_cast<double>(obtuse_count);
}
//������� ��������� ������� ��� ���������� Steiner �� ����� ��������� ��� ��������� ����
Point TriangulationOptimizer::selectBestPoint(CDT& cdt, const std::map<Point, double>& pheromones, double alpha, double beta) {
    double total_probability = 0.0;
    std::map<Point, double> probabilities;
    Point best_point;
    double best_score = std::numeric_limits<double>::infinity();

    // ����������� ����������� ��� ���� ������
    for (auto vertex = cdt.finite_vertices_begin(); vertex != cdt.finite_vertices_end(); ++vertex) {
        Point point = vertex->point();
        double pheromone = pow(pheromones.at(point), alpha);  // ����������� ���������
        double heuristic = pow(calculateHeuristic(cdt, point), beta);// ����������� ��������� �����
        probabilities[point] = pheromone * heuristic;  // ����������� ��������� �����������
        total_probability += probabilities[point]; // �������� ���� ��� �����������
    }

    // ������ ������� �� ���� ��� �����������
    double random_value = (rand() / static_cast<double>(RAND_MAX)) * total_probability;
    double cumulative_probability = 0.0;
    
    // ��������� ��� ����������� ��� ���������� �� ������ �� ��� ���������� ����������
    for (const auto& entry : probabilities) {
        cumulative_probability += entry.second;
        if (cumulative_probability >= random_value) {
            return entry.first;
        }
    }

    return cdt.finite_vertices_begin()->point();  // Default �� �������� 
}

// ����������� ��� ��������� �����
double TriangulationOptimizer::calculateHeuristic(CDT& cdt, const Point& point) {
    double heuristic = 0.0;
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        // ������ ��� Vertex_handle ��� �� ������ ���� ��� �������
        bool found = false;
        for (int i = 0; i < 3; ++i) {
            if (face->vertex(i)->point() == point) {
                found = true;
                break;
            }
        }

        if (found) {
            // ����������� ��� ��������� ����� ��� �� ������
            heuristic += calculateCircumradius(face); // �������������� �� circumradius ��� ��� �������� ����
        }
    }
    return heuristic;
}

// ��������� ��� ��������� ��� �� ������
void TriangulationOptimizer::updatePheromone(std::map<Point, double>& pheromones, const Point& point, double weight) {
    pheromones[point] = std::max(pheromones[point] + weight, 1.0); // Avoid decay to zero
}



// ���������� ������� ������� Steiner
std::set<Point> steiner_points_set; // Maintain unique Steiner points

Point TriangulationOptimizer::generateRandomSteinerPoint(CDT& cdt) {
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();

    // ������ ��� ����� ��� ������� ���� �������
    for (auto vertex = cdt.finite_vertices_begin(); vertex != cdt.finite_vertices_end(); ++vertex) {
        Point p = vertex->point();
        minX = std::min(minX, p.x());
        maxX = std::max(maxX, p.x());
        minY = std::min(minY, p.y());
        maxY = std::max(maxY, p.y());
    }
    // ���������� ������� ������� ����� ��� �����
    for (int attempts = 0; attempts < 100; ++attempts) {
        double x = randomDouble(minX, maxX);
        double y = randomDouble(minY, maxY);

        if (isValidPoint(x, y, cdt)) {
            std::cout << "Generated valid Steiner point: (" << x << ", " << y << ")\n";
            return Point(x, y);
        } else {
            std::cout << "Invalid Steiner point: (" << x << ", " << y << ")\n";
        }
    }

    throw std::runtime_error("Failed to generate valid Steiner point after 100 attempts.");
}


// ������� �� ��� ������ ����� ������ ���� �������������
bool TriangulationOptimizer::isValidPoint(double x, double y, CDT& cdt) {
    Point p(x, y);

    //������� �� �� ������ ����� ���� ���� ������������� � ����� ������
    CDT::Face_handle face = cdt.locate(p);
    bool is_infinite = cdt.is_infinite(face);

    std::cout << "Point (" << x << ", " << y << ") -> "
              << (is_infinite ? "Infinite face" : "Valid face") << "\n";

    return !is_infinite; // ��������� true �� �� ������ ����� ������
}

// ������� �� ��� ��� ����� ����������
bool TriangulationOptimizer::isObtuse(CDT::Face_handle face) {
    auto p0 = face->vertex(0)->point();
    auto p1 = face->vertex(1)->point();
    auto p2 = face->vertex(2)->point();

    double a = CGAL::squared_distance(p0, p1);
    double b = CGAL::squared_distance(p1, p2);
    double c = CGAL::squared_distance(p2, p0);
    
    // ���������� true �� � ��� ����� ����������
    return (a + b <= c) || (b + c <= a) || (c + a <= b);
}
// ����������� ��� ��������� ��� ��� ������������� �� ���� ��� ������ ��� ����������� �������� ��� ��� ������� Steiner
double TriangulationOptimizer::calculateEnergy(CDT& cdt, int steiner_count, double alpha, double beta) {
    int obtuse_count = 0;

    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        if (isObtuse(face)) {
            ++obtuse_count;
        }
    }

    // ��������� ���������� ����������� �������� ��� ������� Steiner
    return alpha * obtuse_count + beta * steiner_count;
}

// �������� ���� ������� ��� ��� �������������
void TriangulationOptimizer::removePoint(CDT& cdt, const Point& point) {
    // ��������� ��� ����������� ����� �� ��� ������������ ��������
    auto nearest_vertex = cdt.finite_vertices_begin();
    double min_distance = std::numeric_limits<double>::infinity();
    
    // ������ ��� ������������ �������
    for (auto vertex = cdt.finite_vertices_begin(); vertex != cdt.finite_vertices_end(); ++vertex) {
        double distance = CGAL::squared_distance(vertex->point(), point); // ����������� ���������
        if (distance < min_distance) {
            min_distance = distance;
            nearest_vertex = vertex;
        }
    }

    // ��������� ��� ����������� �����
    if (nearest_vertex != cdt.finite_vertices_end()) {
        cdt.remove(nearest_vertex);
    }
}

// ������ ��������� ��� �������������� ��������������
void TriangulationOptimizer::localSearch(CDT& cdt, const TriangulationInstance& instance, TriangulationResult& result) {
    int L = instance.parameters.value("L", 50);  // �������� ������� �����������
    bool improvement = true; // ������� ��� ������� ��������
    
    // �������� ��� L ����������� � ����� �� ��� ������� ��������
    for (int iteration = 0; iteration < L && improvement; ++iteration) {
        improvement = false;
        double best_cost = std::numeric_limits<double>::infinity();//�������. �������
        Point best_point;

        // ������ ��� 5 ������������ ������
        for (int trial = 0; trial < 5; ++trial) {
            Point candidate_point = generateRandomSteinerPoint(cdt);  // ���������� ���� ��������� �������
            cdt.insert(candidate_point);// �������� ��� ������� ���� �������������

            // ����������� ��� ������� ��� ���� ��������������
            double current_cost = evaluateTriangulation(cdt);
            if (current_cost < best_cost) {// �� ����� �������� ��� �� �������� ����� ����
                best_cost = current_cost;
                best_point = candidate_point;
                improvement = true; // ������� ��������
            }

            removePoint(cdt, candidate_point);   // ��������� ��� ��������������
        }


        // �� ���������� � �������������, ������� �� �������� ������
        if (improvement) {
           cdt.insert(best_point); // �������� ��� ��������� �������
           steiner_points_set.insert(best_point); // �������� ��� ������ ��� Steiner points
           result.steiner_points_x.push_back(std::to_string(best_point.x()));
           result.steiner_points_y.push_back(std::to_string(best_point.y()));
        }

    }
}


// ������������� �������� ��� �������������� ��������������
void TriangulationOptimizer::simulatedAnnealing(CDT& cdt, const TriangulationInstance& instance, TriangulationResult& result) {
    double alpha = instance.parameters.value("alpha", 1.0);  // ���������� ����������� ��������
    double beta = instance.parameters.value("beta", 1.0); // ���������� ������� Steiner
    double temperature = instance.parameters.value("temperature", 1000.0); // ������ �����������
    double cooling_rate = instance.parameters.value("cooling_rate", 0.95); // ����������� �����
    int L = instance.parameters.value("L", 100);

    double current_energy = calculateEnergy(cdt, result.steiner_points_x.size(), alpha, beta);
    
        // ��������� ��� L ����������� � ����� �� ����� ������ ����
    for (int iteration = 0; iteration < L; ++iteration) {
       
        bool improvement = false;
        Point best_point;
        double best_energy = std::numeric_limits<double>::infinity();
        
        // ������� ��� 5 ������������ ������
        for (int trial = 0; trial < 5; ++trial) {
            Point candidate_point = generateRandomSteinerPoint(cdt);

            cdt.insert(candidate_point);
            double new_energy = calculateEnergy(cdt, result.steiner_points_x.size() + 1, alpha, beta);

            if (new_energy < best_energy) { //�� � ��� �������� ����� ��������, ������ ��������
                best_energy = new_energy;
                best_point = candidate_point;
                improvement = true;
            }

            removePoint(cdt, candidate_point);
        }

        if (improvement) {
            cdt.insert(best_point);
            result.steiner_points_x.push_back(std::to_string(best_point.x()));
            result.steiner_points_y.push_back(std::to_string(best_point.y()));
            current_energy = best_energy;
        } else {
            std::cout << "No improvement in iteration " << iteration << ".\n";
        }
        // ������ ������������
        temperature *= cooling_rate;

        if (temperature < 1e-5) {
            std::cout << "Temperature below threshold. Stopping.\n";
            break;
        }
    }
}





// ������� ���������� ��� �������������� ��������������
void TriangulationOptimizer::antColonyOptimization(CDT& cdt, const TriangulationInstance& instance, TriangulationResult& result) {
    // ������������ ����������
    int L = instance.parameters.value("L", 50);  
    int K = instance.parameters.value("K", 10);  
    double alpha = instance.parameters.value("alpha", 2.0);
    double beta = instance.parameters.value("beta", 1.0);  
    double evaporation_rate = instance.parameters.value("evaporation_rate", 0.5); 
    // ������������ ���������
    std::map<CDT::Face_handle, double> pheromone;
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        pheromone[face] = 1.0;  // ������������ ���������
    }
    //�������� ���������� �������� 
    for (int cycle = 0; cycle < L; ++cycle) {
        for (int ant = 0; ant < K; ++ant) {
            // ������� ������� ����
            auto face = cdt.finite_faces_begin();
            std::advance(face, rand() % cdt.number_of_faces());

            // ����������� circumradius ��� ��������� ����������
            double circumradius = CGAL::sqrt(CGAL::squared_distance(face->vertex(0)->point(), face->vertex(1)->point())) +
                                  CGAL::sqrt(CGAL::squared_distance(face->vertex(1)->point(), face->vertex(2)->point())) +
                                  CGAL::sqrt(CGAL::squared_distance(face->vertex(2)->point(), face->vertex(0)->point()));
            double heuristic = circumradius / 3.0; // ��������� �����������

            // ����������� ��� �����������
            double probability = pow(pheromone[face], alpha) * pow(heuristic, beta);

            // �������� Steiner point
            Point steiner_point = CGAL::centroid(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point());
            cdt.insert(steiner_point);
        }

        // �������� ��������� ��� ����������
        for (auto& [face, pheromone_value] : pheromone) {
            int obtuse_count_before = countObtuseAngles(cdt);
            pheromone_value *= (1 - evaporation_rate);  // ��������
            if (!cdt.is_infinite(face)) {
               // �� � �������� ���������, ���������� �� ��������
              int obtuse_count_after = countObtuseAngles(cdt);
            if (obtuse_count_after < obtuse_count_before) {
               pheromone_value += 1.0 / (1 + obtuse_count_after);  // �������� ��������� �� �������� �� obtuse angles
        }
    }
}

    }

    // ���������� ��� Steiner points
    for (auto vertex = cdt.finite_vertices_begin(); vertex != cdt.finite_vertices_end(); ++vertex) {
        if (std::find(instance.points_x.begin(), instance.points_x.end(), vertex->point().x()) == instance.points_x.end()) {
            result.steiner_points_x.push_back(std::to_string(vertex->point().x()));
            result.steiner_points_y.push_back(std::to_string(vertex->point().y()));
        }
    }

    // ����������� obtuse angles
    result.obtuse_count = countObtuseAngles(cdt);
}

