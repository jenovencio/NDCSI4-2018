
#pragma once

#include "NOX_LAPACK_Matrix.H"

NOX::LAPACK::Matrix<double> Multiply(const NOX::LAPACK::Matrix<double>& aA, const NOX::LAPACK::Matrix<double>& aB);
NOX::LAPACK::Matrix<double> Add(const NOX::LAPACK::Matrix<double>& aA, const NOX::LAPACK::Matrix<double>& aB);
void Add(NOX::LAPACK::Matrix<double>& aTarget, const NOX::LAPACK::Matrix<double>& aAddedMatrix);
