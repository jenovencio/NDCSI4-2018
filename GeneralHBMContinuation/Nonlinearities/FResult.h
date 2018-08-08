
struct FResult;

#pragma once
#include "NOX_LAPACK_Vector.H"

struct FResult
{
public:
    // forces
    NOX::LAPACK::Vector F;
    // "correction", parameter used for the next time step
    NOX::LAPACK::Vector C;
};
