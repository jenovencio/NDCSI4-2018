
class NonlinearBaseTDHD;

#pragma once
#include "NonlinearBaseH.h"
#include "JResult.h"

// Nonlinearity purely in time domain
// history dependent (HD)

class NonlinearBaseTDHD : public NonlinearBaseH
{
protected:
    
    // compute derivative of forces by harmonic coefficients in time domain
    // size of the return matrix: (dof) x (dof*harm)
    virtual NOX::LAPACK::Matrix<double> ComputeDFDH(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev, const NOX::LAPACK::Matrix<double>& aJPrev, const int& aTimePointIndex) const final;
    
    virtual bool IsHistoryDependent() const final { return true; }
    
    // 4 derivatives necessary to compute the DFDH for the history dependent nonlinearity
    virtual JResult ComputeDerivatives(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const = 0;
    
};
