
#include <algorithm>

#include "FricElem3D.h"
#include "../Misc.h"

FResult FricElem3D::ComputeFTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const
{
    if (aX.length() % 3 != 0) throw "Number of dofs must be divisible by 3!";
    
    FResult lReturnResult;
    lReturnResult.FValues = NOX::LAPACK::Vector(aX.length());
    lReturnResult.XCorr = NOX::LAPACK::Vector(aX.length());
    
    lReturnResult.XCorrSet = true;
    
    // calculations start
    int lDpN = 3;
    int lNodeCount = aX.length() / lDpN;
    
    double lKn  = 1;
    double lKtx = 1;
    double lKty = 1;
    
    double lN0 = 0.04;
    double lMu = 0.2;
    
    NOX::LAPACK::Vector lUx(lNodeCount);
    NOX::LAPACK::Vector lUy(lNodeCount);
    NOX::LAPACK::Vector lV (lNodeCount);
    
    NOX::LAPACK::Vector lUxr(lNodeCount);
    NOX::LAPACK::Vector lUyr(lNodeCount);
    
    NOX::LAPACK::Vector lCoul(lNodeCount);
    
    NOX::LAPACK::Vector lTx(lNodeCount);
    NOX::LAPACK::Vector lTy(lNodeCount);
    NOX::LAPACK::Vector lN (lNodeCount);
    
    for (int iNode = 0; iNode < lNodeCount; iNode++)
    {
        lUx(iNode) = aX(lDpN * iNode);
        lUy(iNode) = aX(lDpN * iNode + 1);
        lV (iNode) = aX(lDpN * iNode + 2);
        
        lUxr(iNode) = aXPrev(lDpN * iNode);
        lUyr(iNode) = aXPrev(lDpN * iNode + 1);
        
        lN(iNode) = std::max(lN0 + lKn * lV(iNode), 0.0);
        
        lCoul(iNode) = lMu * lN(iNode);
        lTx(iNode) = lKtx * (lUx(iNode) - lUxr(iNode));
        lTy(iNode) = lKty * (lUy(iNode) - lUyr(iNode));
        
        if (std::abs(lTx(iNode)) > lCoul(iNode))
        {
            lTx(iNode) = Sgn(lTx(iNode)) * lCoul(iNode);
            lUxr(iNode) = lUx(iNode) - lTx(iNode) / lKtx;
        }
        if (std::abs(lTy(iNode)) > lCoul(iNode))
        {
            lTy(iNode) = Sgn(lTy(iNode)) * lCoul(iNode);
            lUyr(iNode) = lUy(iNode) - lTy(iNode) / lKty;
        }
        
        lReturnResult.FValues(iNode * lDpN) = lTx(iNode);
        lReturnResult.FValues(iNode * lDpN + 1) = lTy(iNode);
        lReturnResult.FValues(iNode * lDpN + 2) = lN(iNode) - lN0;
        
        lReturnResult.XCorr(iNode * lDpN)     = lUxr(iNode);
        lReturnResult.XCorr(iNode * lDpN + 1) = lUyr(iNode);
        lReturnResult.XCorr(iNode * lDpN + 2) = aXPrev(iNode * lDpN + 2); // not sure about this, but it probably doesn't matter
    }
    
//     std::cout << "Displacements      : " << aX << std::endl;
//     std::cout << "Displacements corr : " << lReturnResult.XCorr << std::endl;
//     std::cout << "Displacements prev : " << aXPrev << std::endl;
//     std::cout << "Forces             : " << lReturnResult.FValues << std::endl;
//     
//     STOP
    
    return lReturnResult;
}
// time domain to time domain
NOX::LAPACK::Matrix<double> FricElem3D::ComputeJacobianTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const
{
    if (aX.length() % 3 != 0) throw "Number of dofs must be divisible by 3!";
    
    NOX::LAPACK::Matrix<double> lReturnMatrix(aX.length(), aX.length());
    
    // calculations start
    int lDpN = 3;
    int lNodeCount = aX.length() / lDpN;
    
    double lKn  = 1;
    double lKtx = 1;
    double lKty = 1;
    
    double lN0 = 0.04;
    double lMu = 0.2;
    
    NOX::LAPACK::Vector lUx(lNodeCount);
    NOX::LAPACK::Vector lUy(lNodeCount);
    NOX::LAPACK::Vector lV (lNodeCount);
    
    NOX::LAPACK::Vector lUxr(lNodeCount);
    NOX::LAPACK::Vector lUyr(lNodeCount);
    
    NOX::LAPACK::Vector lCoul(lNodeCount);
    
    NOX::LAPACK::Vector lTx(lNodeCount);
    NOX::LAPACK::Vector lTy(lNodeCount);
    NOX::LAPACK::Vector lN (lNodeCount);
    
    for (int iNode = 0; iNode < lNodeCount; iNode++)
    {
        lUx(iNode) = aX(lDpN * iNode);
        lUy(iNode) = aX(lDpN * iNode + 1);
        lV (iNode) = aX(lDpN * iNode + 2);
        
        lUxr(iNode) = aXPrev(lDpN * iNode);
        lUyr(iNode) = aXPrev(lDpN * iNode + 1);
        
//         lN(iNode) = std::max(lN0 + lKn * lV(iNode), 0.0);
        
        lCoul(iNode) = lMu * lN(iNode);
        lTx(iNode) = lKtx * (lUx(iNode) - lUxr(iNode));
        lTy(iNode) = lKty * (lUy(iNode) - lUyr(iNode));
        
        if (std::abs(lTx(iNode)) > lCoul(iNode))
        {
            lReturnMatrix(iNode * lDpN, iNode * lDpN + 2) += lKn * lMu * Sgn(lTx(iNode));
        }
        else
        {
            lReturnMatrix(iNode * lDpN, iNode * lDpN) += lKtx;
        }
        
        if (std::abs(lTy(iNode)) > lCoul(iNode))
        {
            lReturnMatrix(iNode * lDpN + 1, iNode * lDpN + 2) += lKn * lMu * Sgn(lTy(iNode));
        }
        else
        {
            lReturnMatrix(iNode * lDpN + 1, iNode * lDpN + 1) += lKty;
        }
        
        if (lN0 + lKn * lV(iNode) > 0)
        {
            lReturnMatrix(iNode * lDpN + 2, iNode * lDpN + 2) += lKn;
        }
    }
    
//     std::cout << "Displacements      : " << aX << std::endl;
//     std::cout << "Displacements corr : " << lReturnResult.XCorr << std::endl;
//     std::cout << "Displacements prev : " << aXPrev << std::endl;
//     std::cout << "Forces             : " << lReturnResult.FValues << std::endl;
//     
//     STOP
    
//     NOX::LAPACK::Matrix<double> lJacTest = ComputeJacobianFiniteDifferenceTD(aX, aXPrev);
//     
//     std::cout << "Analytic jacobian: " << std::endl;
//     std::cout << lReturnMatrix << std::endl;
//     
//     std::cout << "Finite diff jacobian: " << std::endl;
//     std::cout << lJacTest << std::endl;
//     
//     STOP
    
    return lReturnMatrix;
}
int FricElem3D::NumberOfPrepLoops() const
{
    return 1;
}
std::vector<int> FricElem3D::NonzeroFPositions() const
{
    std::vector<int> lReturnVector;
    
    // all dofs can be nonlinear
    for (int i = 0; i < DofCountTimeDomain(); i++)
    {
        lReturnVector.push_back(i);
    }
    
    return lReturnVector;
}
bool FricElem3D::IsCorrectingX() const
{
    return true;
}
void FricElem3D::LoadFromFile(const std::string& aFilePath)
{
}
std::string FricElem3D::ClassName() const
{
    return "Friction Element 3D";
}
