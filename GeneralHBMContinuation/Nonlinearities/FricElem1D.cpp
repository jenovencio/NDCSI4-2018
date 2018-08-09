
#include <algorithm>
#include <sstream>

#include "FricElem1D.h"
#include "../Misc.h"
#include "../Functions.h"

FResult FricElem1D::ComputeFAndC(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aC) const
{
    FResult lReturnResult;
    
    lReturnResult.F = NOX::LAPACK::Vector(aX.length());
    lReturnResult.C = NOX::LAPACK::Vector(aX.length());
    
    for (int i = 0; i < mFrictionDefinitions.size(); i++)
    {
        const FricElem1DDef& lDef = mFrictionDefinitions.at(i);
        int lDofID = lDef.DofID;
        
        if (lDofID < 0 || lDofID >= aX.length()) throw "Dof ID is out of range!";
        
        double lUx = aX(lDofID);
        double lUxr = aC(lDofID);
        double lN = lDef.N;
        double lCoul = lDef.Mu * lN;
        
        double lTx = lDef.K * (lUx - lUxr);
        
        if (std::abs(lTx) > lCoul)
        {
            lTx = Sgn(lTx) * lCoul;
            lUxr = lUx - lTx / lDef.K;
        }
        
        lReturnResult.F(lDofID) = lTx;
        lReturnResult.C(lDofID) = lUxr;

    }
    return lReturnResult;
}
JResult FricElem1D::ComputeDerivatives(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aC) const
{
    JResult lReturnResult;
    
    lReturnResult.DFbyDX = NOX::LAPACK::Matrix<double>(aX.length(), aX.length());
    lReturnResult.DFbyDC = NOX::LAPACK::Matrix<double>(aX.length(), aX.length());
    lReturnResult.DGbyDX = NOX::LAPACK::Matrix<double>(aX.length(), aX.length());
    lReturnResult.DGbyDC = NOX::LAPACK::Matrix<double>(aX.length(), aX.length());
    
    for (int i = 0; i < mFrictionDefinitions.size(); i++)
    {
        const FricElem1DDef& lDef = mFrictionDefinitions.at(i);
        int lDofID = lDef.DofID;
        
        if (lDofID < 0 || lDofID >= aX.length()) throw "Friction dof ID is out of range!";
        
        double lUx = aX(lDofID);
        double lUxr = aC(lDofID);
        double lN = lDef.N;
        double lCoul = lDef.Mu * lN;
        
        double lTx = lDef.K * (lUx - lUxr);
        
        lReturnResult.DFbyDX(lDofID, lDofID) = lDef.K;
        lReturnResult.DFbyDC(lDofID, lDofID) = -lDef.K;
        lReturnResult.DGbyDX(lDofID, lDofID) = 0;
        lReturnResult.DGbyDC(lDofID, lDofID) = 1;
        
        if (std::abs(lTx) > lCoul)
        {
            // no need to set these, just to see the expressions so we can make the derivatives
    //         lTx = Sgn(lTx) * lCoul;
    //         lUxr = lUx - lTx / mKtx;
            
            lReturnResult.DFbyDX(lDofID, lDofID) = 0;
            lReturnResult.DFbyDC(lDofID, lDofID) = 0;
            lReturnResult.DGbyDX(lDofID, lDofID) = 1;
            lReturnResult.DGbyDC(lDofID, lDofID) = 0;
        }
    }
    return lReturnResult;
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
void FricElem1D::LoadFromFile(const std::string& aFilePath)
{
    std::vector<std::string> lAllLines = ReadAllValidLines(aFilePath);
    
    for (int iLine = 0; iLine < lAllLines.size(); iLine++)
    {
        std::string lLine = lAllLines[iLine];
        
        std::stringstream lStream(lLine);
        
        int lDofID;
        double lK;
        double lN;
        double lMu;
        
        lStream >> lDofID >> lK >> lN >> lMu;
        
        FricElem1DDef lDef;
        lDef.DofID = lDofID;
        lDef.K = lK;
        lDef.N = lN;
        lDef.Mu = lMu;
        
        AddFrictionDef(lDef);
    }
}
void FricElem1D::AddFrictionDef(const FricElem1DDef& aDef)
{
    if (IsFinalised()) throw "Can not add another friction element, the object is already finalised!";
    
    int lDofIndex = aDef.DofID;
    
    if (lDofIndex < 0) throw "Friction dof index must be non negative!";
    if (mFrictionDofs.find(lDofIndex) != mFrictionDofs.end()) throw "Dof " + std::to_string(lDofIndex) + " is already added to the system!";
    
    mFrictionDefinitions.push_back(aDef);
    mFrictionDofs.insert(lDofIndex);
}

// std::vector<int> FricElem1D::NonzeroFPositions() const
// {
//     std::vector<int> lReturnVector;
//     lReturnVector.push_back(mDofID);
//     return lReturnVector;
// }
