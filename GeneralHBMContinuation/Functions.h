
#pragma once
#include <NOX_LAPACK_Matrix.H>
#include "Config.h"

Config LoadConfig(const std::string& aConfigFilePath);
void PrintConfig(const Config& aConfig);
// nox lapack matrices are stored as col major      
//! Access the (i,j) entry of A
// T& operator()(int i, int j) { return entries[i + (p*j)]; }
NOX::LAPACK::Matrix<double> LoadSquareMatrix(const std::string& aFilePath);
