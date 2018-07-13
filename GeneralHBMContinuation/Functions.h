
#pragma once
#include "NOX_LAPACK_Matrix.H"
#include "NOX_LAPACK_Vector.H"
#include "Config.h"

Config LoadConfig(const std::string& aConfigFilePath);
void PrintConfig(const Config& aConfig);
// nox lapack matrices are stored as col major      
// Access the (i,j) entry of A
// T& operator()(int i, int j) { return entries[i + (p*j)]; }
NOX::LAPACK::Matrix<double> LoadSquareMatrix(const std::string& aFilePath);

// aHarmonicCount here means total number of waves, so cosines and sines together
// aHarmonicIndex goes between 0 and aHarmonicCount - 1
inline int GetHBMDofIndex(const int& aStandardDofIndex, const int& aHarmonicIndex, const int& aHarmonicCount)
{
    return aStandardDofIndex * aHarmonicCount + aHarmonicIndex;
}
// index of two harmonic functions and a time point into a serial index
// aHarmonicCount here means total number of waves, so cosines and sines together
inline int GetBProductIndex(const int& aBIndex1, const int& aBIndex2, const int& aIntPointIndex, const int& aHarmonicCount, const int& aIntPointCount)
{
    return (aBIndex1 * aHarmonicCount + aBIndex2) * aIntPointCount + aIntPointIndex;
}
