
struct Config;

#pragma once
#include <string>

struct Config
{
public:
    // number of harmonic waves (first harmonic is the DC component)
    // that means that the number of harmonic dofs will be DOFCount * (2 * HarmonicCount - 1)
    int HarmonicWaveCount;
    // min frequency
    double FrequencyMin;
    // max frequency
    double FrequencyMax;
    // integration point count
    double IntPointCount;
    
    std::string ConfigFilePath;
    std::string ConfigFileName;
    // relative paths with respect to the config path
    std::string MassMatrixFilePath;
    std::string StiffnessMatrixFile;
    std::string DampingMatrixFile;
    std::string ExcitationForceFile;
};
