
class NonlinearBaseH;

#pragma once
#include "NonlinearBase.h"

// "hybrid" nonlinearity
// evaluates forces in time domain
// evaluates derivatives of forces in time domain, but by harmonic coefficients
class NonlinearBaseH : public NonlinearBase
{
protected:
    // frequency domain to frequency domain
    virtual NOX::LAPACK::Vector ComputeFInner(const NOX::LAPACK::Vector& aX, const double& aFrequency) final;
    // frequency domain to frequency domain
    virtual NOX::LAPACK::Matrix<double> ComputeJacobianInner(const NOX::LAPACK::Vector& aX, const double& aFrequency) final;
    
    // compute forces in time domain
    virtual NOX::LAPACK::Vector ComputeFTD(const NOX::LAPACK::Vector& aX, const int& aTimePointIndex) = 0;
    // compute derivative of forces by harmonic coefficients in time domain
    // size of the return matrix: (dof) x (dof*harm)
    virtual NOX::LAPACK::Matrix<double> ComputeDFDH(const NOX::LAPACK::Vector& aX, const int& aTimePointIndex) = 0;
        
    virtual int NumberOfPrepLoops() const = 0;
    
    // signals that there will be a series of calls to ComputeFTD following (time step after time step, possibly multiple cycles)
    // with the given "value" in frequency domain
    virtual void InitFComputation(const NOX::LAPACK::Vector& aX) { }
    // signals that there will be a series of calls to ComputeDFDH following (time step after time step, possibly multiple cycles)
    // with the given "value" in frequency domain
    virtual void InitJComputation(const NOX::LAPACK::Vector& aX) { }
};
