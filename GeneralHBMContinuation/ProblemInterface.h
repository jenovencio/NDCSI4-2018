
class ProblemInterface;

#pragma once
#include "LOCA_LAPACK_Interface.H"
#include "LOCA_Parameter_Vector.H"
#include "Config.h"

class ProblemInterface : public LOCA::LAPACK::Interface
{
private:
    // number of standard degrees of freedom (in physical domain, not HBM)
    const int cDOFCount;
    double mFrequency = 0;
    
public:
    ProblemInterface(const Config& aConfig);
    virtual void setParams(const LOCA::ParameterVector& aParams) override;
};
