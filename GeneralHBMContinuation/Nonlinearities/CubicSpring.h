
class CubicSpring;
struct CubicSpringDef;

#pragma once
#include "NonlinearBase.h"

struct CubicSpringDef
{
public:
    int DofIndex;
    double StiffnessCoeff;
};

class CubicSpring : public NonlinearBase
{
private:
    std::vector<CubicSpringDef> mSprings;
public:
    // frequency domain to frequency domain
    virtual NOX::LAPACK::Vector ComputeRHS(const NOX::LAPACK::Vector& aX, const double& aFrequency) const override;
    // frequency domain to frequency domain
    virtual NOX::LAPACK::Matrix<double> ComputeJacobian(const NOX::LAPACK::Vector& aX, const double& aFrequency) const override;
    
    void AddCubicSpring(const CubicSpringDef& aDef);
    void AddCubicSpring(const int& aDofIndex, const double& aStiffnessCoeff);
    void ClearSprings();
};
