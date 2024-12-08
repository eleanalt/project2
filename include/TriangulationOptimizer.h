#ifndef TRIANGULATIONOPTIMIZER_H
#define TRIANGULATIONOPTIMIZER_H

#include "TriangulationInstance.h"  // Εισαγωγή του TriangulationInstance.h
#include "TriangulationResult.h"
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <nlohmann/json.hpp>  // Για την χρήση του nlohmann::json
#include <map>

// Ορισμός του τύπου K για τον ακριβή υπολογισμό των σημείων
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
// Ορισμός του τύπου CDT (Constrained Delaunay Triangulation)
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
// Ορισμός του τύπου Point για τα σημεία της τριγωνοποίησης
typedef CDT::Point Point;

class TriangulationOptimizer {
public:
    TriangulationResult optimize(const TriangulationInstance& instance);
    double calculateCircumradius(CDT::Face_handle face); 
    Point findNearestVertex(CDT& cdt, const Point& point);
    int countObtuseAngles(const CDT& cdt);

private:
    double evaluateTriangulation(CDT& cdt);
    double calculateHeuristic(CDT& cdt, const Point& point);
    double calculateEnergy(CDT& cdt, int steiner_count, double alpha, double beta);
    Point generateRandomSteinerPoint(CDT& cdt);
    bool isValidPoint(double x, double y, CDT& cdt);
    void localSearch(CDT& cdt, const TriangulationInstance& instance, TriangulationResult& result);
    void simulatedAnnealing(CDT& cdt, const TriangulationInstance& instance, TriangulationResult& result);
    void removePoint(CDT& cdt, const Point& point);
	Point selectBestPoint(CDT& cdt, const std::map<Point, double>& pheromones, double alpha, double beta);

    void antColonyOptimization(CDT& cdt, const TriangulationInstance& instance, TriangulationResult& result);

    void updatePheromone(std::map<Point, double>& pheromones, const Point& point, double pheromone_weight);
    bool isObtuse(CDT::Face_handle face);
};

#endif // TRIANGULATIONOPTIMIZER_H

