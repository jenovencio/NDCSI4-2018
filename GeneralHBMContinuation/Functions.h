
#pragma once
#include <fstream>
#include <map>

#include "NOX_LAPACK_Matrix.H"
#include "NOX_LAPACK_Vector.H"
#include "Config.h"
#include "MatrixDefinition.h"
#include "Nonlinearities/NonlinearityDefinition.h"

std::map<std::string, std::vector<std::string>> ParseGeneralConfigFile(const std::string& aFilePath);
// nox lapack matrices are stored as col major      
// Access the (i,j) entry of A
// T& operator()(int i, int j) { return entries[i + (p*j)]; }
NOX::LAPACK::Matrix<double> LoadSquareMatrixFull(const std::string& aFilePath, int& aDim);
// loads sparse matrix from a file
NOX::LAPACK::Matrix<double> LoadSquareMatrixSparse(const std::string& aFilePath, int& aDim);
NOX::LAPACK::Matrix<double> LoadSquareMatrix(const std::string& aFilePath, const std::string& aMatrixType, int& aDim);
// aCheckHarmCount - total number of harmonic coefficients per physical dof
std::vector<double> LoadExcitationForce(const std::string& aFilePath, int& aDim, int& aHarmCoeffCount);
// skips comment lines and returns the next valid (non comment) line in the file
std::string GetNextValidLine(std::ifstream& aFile, bool& aIsValid);

// aHarmonicCount here means total number of waves, so cosines and sines together
// aHarmonicIndex goes between 0 and aHarmonicCount - 1
inline int GetHBMDofIndex(const int& aStandardDofIndex, const int& aHarmonicIndex, const int& aHarmonicCount)
{
    return aStandardDofIndex * aHarmonicCount + aHarmonicIndex;
}
void CheckString(const std::string& aString, const std::vector<std::string>& aPossibilities, const std::string& aGroupName = "");
// the same check string, but uses the keys of the provided mam as the vector of possibilities
template <class T>
void CheckString(const std::string& aString, const std::map<std::string, T>& aPossibilities, const std::string& aGroupName = "");
std::vector<double> GetRelativeTimePoints(const int& aTimePointCount);
NOX::LAPACK::Vector GetAverage(const std::vector<NOX::LAPACK::Vector>& aVectors);
std::string SkipWhiteSpaces(const std::string& aString);
