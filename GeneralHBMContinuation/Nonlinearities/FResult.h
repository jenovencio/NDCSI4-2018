
struct FResult;

#pragma once
#include "NOX_LAPACK_Vector.H"

struct FResult
{
public:
    NOX::LAPACK::Vector FValues;
    // correction of the provided X
    NOX::LAPACK::Vector XCorr;
    bool XCorrSet = false;
};
