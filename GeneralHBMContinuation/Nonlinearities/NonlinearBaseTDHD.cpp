
#include "NonlinearBaseTDHD.h"

// compute forces in time domain
NOX::LAPACK::Vector NonlinearBaseTDHD::ComputeFTD(const NOX::LAPACK::Vector& aX, const int& aTimePointIndex)
{
    FResult lResult = ComputeFAndC(aX, mC);
    
    if (lResult.F.length() != mProblemParams.DofCountPhysical)
        throw "Invalid size of the F vector! Expected size: " + std::to_string(mProblemParams.DofCountPhysical);
    
    if (lResult.C.length() != mProblemParams.DofCountPhysical)
        throw "Invalid size of the C vector! Expected size: " + std::to_string(mProblemParams.DofCountPhysical);
    
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
    
    if (mC.length() != mProblemParams.DofCountPhysical) throw "Invalid size of the C vector! Expected size: " + std::to_string(mProblemParams.DofCountPhysical);
}
// signals that there will be a series of calls to ComputeDFDH following (time step after time step, possibly multiple cycles)
// with the given "value" in frequency domain
void NonlinearBaseTDHD::InitJComputation(const NOX::LAPACK::Vector& aX)
{
    InitFComputation(aX);
    
    
    mDGbyDH = InitDGbyDH(aX);
    
    if (mDGbyDH.numRows() != mProblemParams.DofCountPhysical) throw "Invalid number of rows of the DGbyDH matrix! Expected size: " + std::to_string(mProblemParams.DofCountPhysical);
    
    if (mDGbyDH.numCols() != mProblemParams.DofCountHBM) throw "Invalid number of columns of the DGbyDH matrix! Expected size: " + std::to_string(mProblemParams.DofCountHBM);
}
