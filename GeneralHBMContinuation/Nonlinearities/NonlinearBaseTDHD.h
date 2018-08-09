
class NonlinearBaseTDHD;

#pragma once
#include "NonlinearBaseH.h"
#include "FResult.h"
#include "JResult.h"

// Nonlinearity purely in time domain
// history dependent (HD)

class NonlinearBaseTDHD : public NonlinearBaseH
{
private:
    // "correction", history dependence parameter, which cycles through the forces evaluations
    NOX::LAPACK::Vector mC;
    NOX::LAPACK::Matrix<double> mDGbyDH;
    bool mIsFirst = false;
    
protected:
    
    // compute forces in time domain
    virtual NOX::LAPACK::Vector ComputeFTD(const NOX::LAPACK::Vector& aX, const int& aTimePointIndex) final;
    // compute derivative of forces by harmonic coefficients in time domain
    // size of the return matrix: (dof) x (dof*harm)
    virtual NOX::LAPACK::Matrix<double> ComputeDFDH(const NOX::LAPACK::Vector& aX, const int& aTimePointIndex) final;
    
    virtual bool IsHistoryDependent() const final { return true; }
    
    virtual FResult ComputeFAndC(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aC) const = 0;
    // 4 derivatives necessary to compute the DFDH for the history dependent nonlinearity
    // the aC parameter comes from the previous F evaluation (for previous time step)
    virtual JResult ComputeDerivatives(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aC) const = 0;
    virtual NOX::LAPACK::Vector InitC(const NOX::LAPACK::Vector& aX) const = 0;
    virtual NOX::LAPACK::Matrix<double> InitDGbyDH(const NOX::LAPACK::Vector& aX) const { return NOX::LAPACK::Matrix<double>(mProblemParams.DofCountPhysical, mProblemParams.DofCountHBM); }
    
    // signals that there will be a series of calls to ComputeFTD following (time step after time step, possibly multiple cycles)
    // with the given "value" in frequency domain
    virtual void InitFComputation(const NOX::LAPACK::Vector& aX) final;
    // signals that there will be a series of calls to ComputeDFDH following (time step after time step, possibly multiple cycles)
    // with the given "value" in frequency domain
    virtual void InitJComputation(const NOX::LAPACK::Vector& aX) final;
    
    
};
