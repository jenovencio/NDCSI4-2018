
#include "NonlinearBaseTDHD.h"

// compute forces in time domain
NOX::LAPACK::Vector NonlinearBaseTDHD::ComputeFTD(const NOX::LAPACK::Vector& aX, const int& aTimePointIndex)
{
    FResult lResult = ComputeFAndC(aX, mC);
    mC = lResult.C;
    
    return lResult.F;
}

// compute derivative of forces by harmonic coefficients in time domain
// size of the return matrix: (dof) x (dof*harm)
NOX::LAPACK::Matrix<double> NonlinearBaseTDHD::ComputeDFDH(const NOX::LAPACK::Vector& aX, const int& aTimePointIndex)
{
    // check
//     if (aFrequency <= 0) throw "Frequency must be a positive value!";
    if (mProblemParams.HarmonicCount <= 0) throw "Number of harmonic coefficients must be a positive integer!";
    
    int lDofCount = mProblemParams.DofCountPhysical;
    
    NOX::LAPACK::Matrix<double> lReturnMatrix(lDofCount, mProblemParams.DofCountHBM);
    
    //TODO do the proper jacobian with history logic
    
    return lReturnMatrix;
}

// signals that there will be a series of calls to ComputeFTD following (time step after time step, possibly multiple cycles)
// with the given "value" in frequency domain
void NonlinearBaseTDHD::InitFComputation(const NOX::LAPACK::Vector& aX)
{
    mC = InitC(aX);
}
// signals that there will be a series of calls to ComputeDFDH following (time step after time step, possibly multiple cycles)
// with the given "value" in frequency domain
void NonlinearBaseTDHD::InitJComputation(const NOX::LAPACK::Vector& aX)
{
    mC = InitC(aX);
    mDGbyDH = InitDGbyDH(aX);
}
