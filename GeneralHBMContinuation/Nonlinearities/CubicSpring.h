
class CubicSpring;
struct CubicSpringDef;

#pragma once
#include "NonlinearBase.h"

struct CubicSpringDef
{
public:
    // a cubic spring between two dofs
    // if one of the dofs is set to -1, the spring is considered to be grounded on that side
    int Dof1Index;
    int Dof2Index;
    double StiffnessCoeff;
};

class CubicSpring : public NonlinearBase
{
    
private:
    std::vector<CubicSpringDef> mSprings;
    
protected:
    // frequency domain to frequency domain
    virtual FResult ComputeFTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const override;
    // frequency domain to frequency domain
    virtual NOX::LAPACK::Matrix<double> ComputeJacobianTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const override;
    virtual int NumberOfPrepLoops() const override;
    // returns indices of dofs (in time domain) that are nonlinear, i.e.
    // that can potentially return a nonzero value in their F evaluation
    // this is just to speed up the F and jacobian computations (so we don't uselessly iterate over all dofs)
    virtual std::vector<int> NonzeroFPositions() const override;
    virtual bool IsHistoryDependent() const override;
    
public:
    void AddCubicSpring(const CubicSpringDef& aDef);
    void AddCubicSpring(const int& aDofIndex, const int& aDof2Index, const double& aStiffnessCoeff);
    void ClearSprings();
    virtual void LoadFromFile(const std::string& aFilePath) override;
    virtual std::string ClassName() const override;
};
