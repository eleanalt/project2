#ifndef OUTPUTGENERATOR_H
#define OUTPUTGENERATOR_H

#include <string>
#include "TriangulationResult.h"

class OutputGenerator {
public:
    // Κατασκευαστής που παίρνει τη διαδρομή του αρχείου εξόδου
    OutputGenerator(const std::string& outputFilePath);  // Υποθέτουμε ότι παίρνει τη διαδρομή του αρχείου

    // Μέθοδος για την εγγραφή στο αρχείο
    void write(const TriangulationResult& result);

private:
    std::string outputFilePath;  // Η διαδρομή του αρχείου εξόδου
};

#endif // OUTPUTGENERATOR_H

