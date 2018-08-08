
class NonlinearBaseTD;

#pragma once
#include <functional>

#include "fftw3.h"

#include "NOX_LAPACK_Matrix.H"
#include "NOX_LAPACK_Vector.H"

#include "NonlinearBaseH.h"
#include "../Aft/AftBase.h"
#include "../ProblemParams.h"

#define DEFAULT_FD_STEP 1e-8

// nonlinearity defined purely in time domain
// history independent
class NonlinearBaseTD : public NonlinearBaseH
{
    
private:
    // indices of nonlinear dofs in time domain, not to be touched outside the Init and Finalise methods
//     std::vector<int> mNonzeroFPositions;
    
public:
    // "locks" the nonlinearity to any changes (logically, not actually of course)
    // the point is, after this method is called, the nonlinearity should not change anymore (in terms of it's mathematical definition at least),
    // so for instance for a cubic spring nonlinearity, no springs should be added or removed, etc.
    // changes to the nonlinearity after this function call might not be taken into consideration when the nonlinearity is not eveluated, that's 
    // the point of this mechanism
    virtual void Finalise() override;
    
protected:
    // compute forces in time domain
    virtual NOX::LAPACK::Vector ComputeFTD(const NOX::LAPACK::Vector& aX, const int& aTimePointIndex) final 
    {
        return ComputeFTDInner(aX);
    }
    virtual NOX::LAPACK::Matrix<double> ComputeDFDH(const NOX::LAPACK::Vector& aX, const int& aTimePointIndex) final;
    // time domain to time domain
    virtual NOX::LAPACK::Vector ComputeFTDInner(const NOX::LAPACK::Vector& aX) = 0;
    // time domain to time domain
    virtual NOX::LAPACK::Matrix<double> ComputeJacobianTimeDomain(const NOX::LAPACK::Vector& aX) = 0;
    
    virtual int NumberOfPrepLoops() const final { return 0; }
    
    // signals that there will be a series of calls to ComputeFTD following (time step after time step, possibly multiple cycles)
    // with the given "value" in frequency domain
    virtual void InitFComputation(const NOX::LAPACK::Vector& aX) final { }
    // signals that there will be a series of calls to ComputeDFDH following (time step after time step, possibly multiple cycles)
    // with the given "value" in frequency domain
    virtual void InitJComputation(const NOX::LAPACK::Vector& aX) final { }
    // returns indices of elements in the F vector (in time domain) that are nonzero, i.e.
    // where some nonlinearity occurs
    // this is just to speed up the F and jacobian computations in the frequency domain (so we don't uselessly iterate over all elements)
//     virtual std::vector<int> NonzeroFPositions() const = 0;
    
//     // general jacobian evaluation functions, using arbitrary F evaluation functions
//     NOX::LAPACK::Matrix<double> ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)>& aFEval, double aStep = DEFAULT_FD_STEP) const;
//     NOX::LAPACK::Matrix<double> ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)>& aFEval, const NOX::LAPACK::Matrix<double>& aSteps) const;
//     
//     // Same functions as above, but using the ComputeF as the aFEval function (in frequency domain)
//     NOX::LAPACK::Matrix<double> ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const double& aFrequency, double aStep = DEFAULT_FD_STEP) const;
//     NOX::LAPACK::Matrix<double> ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const double& aFrequency, const NOX::LAPACK::Matrix<double>& aSteps) const;
//     
//     // Same functions as above, but using the ComputeF as the aFEval function (in time domain)
//     NOX::LAPACK::Matrix<double> ComputeJacobianFiniteDifferenceTD(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev, double aStep = DEFAULT_FD_STEP) const;
//     NOX::LAPACK::Matrix<double> ComputeJacobianFiniteDifferenceTD(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev, const NOX::LAPACK::Matrix<double>& aSteps) const;
    
};

