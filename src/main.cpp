#include <iostream>
#include "InputParser.h"
#include "TriangulationOptimizer.h"
#include "OutputGenerator.h"

int main(int argc, char* argv[]) {
    if (argc != 5 || std::string(argv[1]) != "-i" || std::string(argv[3]) != "-o") {
        std::cerr << "Usage: ./opt_triangulation -i input.json -o output.json\n";
        return 1;
    }

    std::string inputFilePath = argv[2];
    std::string outputFilePath = argv[4];

    try {
        // ���������� �������
        InputParser parser(inputFilePath);
        TriangulationInstance instance = parser.parse();

        // �������� ���������������
        TriangulationOptimizer optimizer;
        TriangulationResult result = optimizer.optimize(instance);

        // ���������� ������� ������
        OutputGenerator outputGen(outputFilePath);
        outputGen.write(result);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

