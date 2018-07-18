
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
    
protected:
    // frequency domain to frequency domain
    virtual NOX::LAPACK::Vector ComputeFTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const override;
    // frequency domain to frequency domain
    virtual NOX::LAPACK::Matrix<double> ComputeJacobianTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const override;
    virtual int NumberOfPrepLoops() const override;
    // returns indices of dofs (in time domain) that are nonlinear, i.e.
    // that can potentially return a nonzero value in their F evaluation
    // this is just to speed up the F and jacobian computations (so we don't uselessly iterate over all dofs)
    virtual std::vector<int> NonzeroFPositions() const override;
    
public:
    void AddCubicSpring(const CubicSpringDef& aDef);
    void AddCubicSpring(const int& aDofIndex, const double& aStiffnessCoeff);
    void ClearSprings();
};
