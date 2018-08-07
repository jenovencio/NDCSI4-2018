
#include "NonlinearBaseTDHD.h"


// compute derivative of forces by harmonic coefficients in time domain
// size of the return matrix: (dof) x (dof*harm)
NOX::LAPACK::Matrix<double> NonlinearBaseTDHD::ComputeDFDH(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev, const NOX::LAPACK::Matrix<double>& aJPrev, const int& aTimePointIndex) const
{
    // check
//     if (aFrequency <= 0) throw "Frequency must be a positive value!";
    if (mProblemParams.HarmonicCount <= 0) throw "Number of harmonic coefficients must be a positive integer!";
    
    int lDofCount = mProblemParams.DofCountPhysical;
    
    NOX::LAPACK::Matrix<double> lReturnMatrix(lDofCount, mProblemParams.DofCountHBM);
    
    
    return lReturnMatrix;
}
