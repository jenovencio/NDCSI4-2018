

class NonlinearBaseFD;

#pragma once
#include "NOX_LAPACK_Matrix.H"
#include "NOX_LAPACK_Vector.H"

#include "NonlinearBase.h"

class NonlinearBaseFD : public NonlinearBase
{
private:
    NOX::LAPACK::Matrix<double> mFiniteDifferenceSteps;
    
public:
    virtual void Init(AftBase* const aAft, const ProblemParams& aProblemParams) override;
    
protected:
    // time domain to time domain
    virtual NOX::LAPACK::Matrix<double> ComputeJacobianTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const override;
    
    void SetFiniteDifferenceSteps(const double& aStep);
    void SetFiniteDifferenceSteps(const NOX::LAPACK::Matrix<double>& aSteps);
    
};
