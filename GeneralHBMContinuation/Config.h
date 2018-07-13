
struct Config;

#pragma once
#include <string>

struct Config
{
public:
    // number of dofs (in physical domain, not HBM)
    int DOFCount;
    // number of harmonic waves (first harmonic is the DC component)
    // that means that the number of harmonic dofs will be DOFCount * (2 * HarmonicCount - 1)
    int HarmonicCount;
    // min frequency
    double FrequencyMin;
    // max frequency
    double FrequencyMax;
    
    std::string MassMatrixFilePath;
    std::string StiffnessMatrixFile;
    std::string DampingMatrixFile;
};
