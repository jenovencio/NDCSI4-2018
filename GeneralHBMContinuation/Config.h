
struct Config;

#pragma once
#include <string>

struct Config
{
public:
    
    std::string ConfigFilePath;
    std::string ConfigFileName;
    // relative paths with respect to the config path
    std::string MassMatrixFile;
    std::string MassMatrixType; // full or sparse
    std::string DampingMatrixFile;
    std::string DampingMatrixType; // full or sparse
    std::string StiffnessMatrixFile;
    std::string StiffnessMatrixType; // full or sparse
    std::string ExcitationForceFile;
    std::string ContinuationSettingsFile;
    
    // number of harmonic waves (first harmonic is the DC component)
    // that means that the number of harmonic dofs will be DOFCount * (2 * HarmonicCount - 1)
    int HarmonicWaveCount;
    // min frequency
    double FrequencyMin;
    // max frequency
    double FrequencyMax;
    // integration point count
    double IntPointCount;
    // should the interface store the whole solutions? (all harmonic coefficients)
    bool SaveWholeSolutions;
};
