
#include <algorithm>

#include "FricElem1D.h"
#include "../Misc.h"
#include "../Functions.h"

FResult FricElem1D::ComputeFAndC(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aC) const
{
    if (mDofID < 0 || mDofID >= aX.length()) throw "Dof ID is out of range!";
    
    FResult lReturnResult;
    
    lReturnResult.F = NOX::LAPACK::Vector(aX.length());
    lReturnResult.C = NOX::LAPACK::Vector(aX.length());
    
    double lUx = 0;
    
    double lUxr;
    
    double lCoul;
    
    double lTx;
    double lN;
    
    lUx = aX(mDofID);
    
    lUxr = aC(mDofID);
    
    lN = mN0;
    
    lCoul = mMu * lN;
    lTx = mKtx * (lUx - lUxr);
    
    if (std::abs(lTx) > lCoul)
    {
        lTx = Sgn(lTx) * lCoul;
        lUxr = lUx - lTx / mKtx;
    }
    
    lReturnResult.F(mDofID) = lTx;
    lReturnResult.C(mDofID) = lUxr;

    return lReturnResult;
}
JResult FricElem1D::ComputeDerivatives(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aC) const
{
    
}
NOX::LAPACK::Vector FricElem1D::InitC(const NOX::LAPACK::Vector& aX) const
{
    NOX::LAPACK::Vector lReturnVector(mProblemParams.DofCountPhysical);
    
    for (int iDof = 0; iDof < mProblemParams.DofCountPhysical; iDof++)
    {
        int lIndexHBM = GetHBMDofIndex(iDof, 0, mProblemParams.HarmonicCount);
        
        lReturnVector(iDof) = aX(lIndexHBM);
    }
    
    return lReturnVector;
}
NOX::LAPACK::Matrix<double> FricElem1D::InitDGbyDH(const NOX::LAPACK::Vector& aX) const
{
    return NOX::LAPACK::Matrix<double>(mProblemParams.DofCountPhysical, mProblemParams.DofCountHBM);
}

// std::vector<int> FricElem1D::NonzeroFPositions() const
// {
//     std::vector<int> lReturnVector;
//     lReturnVector.push_back(mDofID);
//     return lReturnVector;
// }
