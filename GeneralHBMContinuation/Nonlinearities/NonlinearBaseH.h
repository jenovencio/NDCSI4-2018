
class NonlinearBaseH;

#pragma once
#include "NonlinearBase.h"
#include "FResult.h"

// "hybrid" nonlinearity
// evaluates forces in time domain
// evaluates derivatives of forces in time domain, but by harmonic coefficients
class NonlinearBaseH : public NonlinearBase
{
protected:
    // frequency domain to frequency domain
    virtual NOX::LAPACK::Vector ComputeFInner(const NOX::LAPACK::Vector& aX, const double& aFrequency) const override;
    // frequency domain to frequency domain
    virtual NOX::LAPACK::Matrix<double> ComputeJacobianInner(const NOX::LAPACK::Vector& aX, const double& aFrequency) const override;
    
    // compute forces in time domain
    virtual FResult ComputeFTD(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev, const int& aTimePointIndex) const = 0;
    // compute derivative of forces by harmonic coefficients in time domain
    // size of the return matrix: (dof) x (dof*harm)
    virtual NOX::LAPACK::Matrix<double> ComputeDFDH(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev, const NOX::LAPACK::Matrix<double>& aJPrev, const int& aTimePointIndex) const = 0;
    
    virtual NOX::LAPACK::Matrix<double> GetFirstJ() const 
    {
        return NOX::LAPACK::Matrix<double>(mProblemParams->DofCountPhysical, mProblemParams->DofCountHBM);
    }
    
    virtual int NumberOfPrepLoops() const = 0;
    virtual bool IsHistoryDependent() const { return true; }
public:
    
};
