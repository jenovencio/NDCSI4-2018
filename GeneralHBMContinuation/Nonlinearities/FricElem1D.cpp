
// #include <algorithm>
// 
// #include "FricElem1D.h"
// #include "../Misc.h"
// 
// FResult FricElem1D::ComputeFTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const
// {
//     if (mDofID < 0 || mDofID >= aX.length()) throw "Dof ID is out of range!";
//     
//     FResult lReturnResult;
//     lReturnResult.FValues = NOX::LAPACK::Vector(aX.length());
//     lReturnResult.XCorr = NOX::LAPACK::Vector(aX.length());
//     
//     lReturnResult.XCorrSet = true;
//                 
//     double lUx = 0;
//     
//     double lUxr;
//     
//     double lCoul;
//     
//     double lTx;
//     double lN;
//     
//     lUx = aX(mDofID);
//     
//     lUxr = aXPrev(mDofID);
//     
//     lN = mN0;
//     
//     lCoul = mMu * lN;
//     lTx = mKtx * (lUx - lUxr);
//     
//     if (std::abs(lTx) > lCoul)
//     {
//         lTx = Sgn(lTx) * lCoul;
//         lUxr = lUx - lTx / mKtx;
//     }
//     
//     lReturnResult.FValues(mDofID) = lTx;
//     lReturnResult.XCorr(mDofID) = lUxr;
// 
//     return lReturnResult;
// }
// // time domain to time domain
// NOX::LAPACK::Matrix<double> FricElem1D::ComputeJacobianTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const
// {
//     if (mDofID < 0 || mDofID >= aX.length()) throw "Dof ID is out of range!";
//     
//     NOX::LAPACK::Matrix<double> lReturnMatrix(aX.length(), aX.length());
//     
//     double lUx = 0;
//     
//     double lUxr;
//     
//     double lCoul;
//     
//     double lTx;
//     double lN;
//     
//     lUx = aX(mDofID);
//     
//     lUxr = aXPrev(mDofID);
//     
//     lN = mN0;
//     lCoul = mMu * lN;
//     
//     lTx = mKtx * (lUx - lUxr);
//     
//     lReturnMatrix(mDofID, mDofID) = mKtx;
//     
//     if (std::abs(lTx) > lCoul)
//     {
//         lReturnMatrix(mDofID, mDofID) = 0.0;
//     }
//     
//     return lReturnMatrix;
// }
// int FricElem1D::NumberOfPrepLoops() const
// {
//     return 1;
// }
// std::vector<int> FricElem1D::NonzeroFPositions() const
// {
//     std::vector<int> lReturnVector;
//     lReturnVector.push_back(mDofID);
//     return lReturnVector;
// }
// bool FricElem1D::IsHistoryDependent() const
// {
//     return true;
// }
// void FricElem1D::LoadFromFile(const std::string& aFilePath)
// {
//     
// }
// std::string FricElem1D::ClassName() const
// {
//     return "Friction Element 1D";
// }
