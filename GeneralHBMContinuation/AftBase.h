
class AftBase;

#pragma once
#include <vector>

#include "NOX_LAPACK_Matrix.H"
#include "NOX_LAPACK_Vector.H"

class AftBase
{
private:
    
public:
    virtual const std::vector<NOX::LAPACK::Vector>& FrequencyToTime(const NOX::LAPACK::Vector& aXFreq, const double& aFrequency) = 0;
    virtual const NOX::LAPACK::Vector& TimeToFrequency(const std::vector<NOX::LAPACK::Vector>& aXTime, const double& aFrequency) = 0;
};
