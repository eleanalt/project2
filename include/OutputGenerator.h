#ifndef OUTPUTGENERATOR_H
#define OUTPUTGENERATOR_H

#include <string>
#include "TriangulationResult.h"

class OutputGenerator {
public:
    // ������������� ��� ������� �� �������� ��� ������� ������
    OutputGenerator(const std::string& outputFilePath);  // ���������� ��� ������� �� �������� ��� �������

    // ������� ��� ��� ������� ��� ������
    void write(const TriangulationResult& result);

private:
    std::string outputFilePath;  // � �������� ��� ������� ������
};

#endif // OUTPUTGENERATOR_H

