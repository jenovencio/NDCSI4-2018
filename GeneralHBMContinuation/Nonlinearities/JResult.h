
struct JResult;

#pragma once
#include "NOX_LAPACK_Matrix.H"

struct JResult
{
public:
    // derivatives of forces by current time displacements
    NOX::LAPACK::Matrix<double> DFbyDX;
    // derivatives of forces by previous time displacements (corrected)
    NOX::LAPACK::Matrix<double> DFbyDXp;
    // derivatives of current corrected previous time displacements by current time displacements
    NOX::LAPACK::Matrix<double> DGbyDX;
    // derivatives of current corrected previous time displacements by previous time displacements
    NOX::LAPACK::Matrix<double> DGbyDXp;
};
