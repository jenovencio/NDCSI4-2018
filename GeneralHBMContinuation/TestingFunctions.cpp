
// #include "TestingFunctions.h"
// #include "Functions.h"
// #include "Misc.h"
// 
// NOX::LAPACK::Vector ComputeFTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev, const double& aNormalStaticLoad)
// {
//     if (aX.length() % 3 != 0) throw "Number of dofs must be divisible by 3!";
//     
//     FResult lReturnResult;
//     lReturnResult.FValues = NOX::LAPACK::Vector(aX.length());
//     lReturnResult.XCorr = NOX::LAPACK::Vector(aX.length());
//     
//     lReturnResult.XCorrSet = true;
//     
//     // calculations start
//     int lDpN = 3;
//     int lNodeCount = aX.length() / lDpN;
//     
//     double lKn  = 1;
//     double lKtx = 1;
//     double lKty = 1;
//     
//     double lN0 = aNormalStaticLoad;
//     double lMu = 0.2;
//     
//     NOX::LAPACK::Vector lUx(lNodeCount);
//     NOX::LAPACK::Vector lUy(lNodeCount);
//     NOX::LAPACK::Vector lV (lNodeCount);
//     
//     NOX::LAPACK::Vector lUxr(lNodeCount);
//     NOX::LAPACK::Vector lUyr(lNodeCount);
//     
//     NOX::LAPACK::Vector lCoul(lNodeCount);
//     
//     NOX::LAPACK::Vector lTx(lNodeCount);
//     NOX::LAPACK::Vector lTy(lNodeCount);
//     NOX::LAPACK::Vector lN (lNodeCount);
//     
//     for (int iNode = 0; iNode < lNodeCount; iNode++)
//     {
//         lUx(iNode) = aX(lDpN * iNode);
//         lUy(iNode) = aX(lDpN * iNode + 1);
//         lV (iNode) = aX(lDpN * iNode + 2);
//         
//         lUxr(iNode) = aXPrev(lDpN * iNode);
//         lUyr(iNode) = aXPrev(lDpN * iNode + 1);
//         
//         lN(iNode) = std::max(lN0 + lKn * lV(iNode), 0.0);
//         
//         lCoul(iNode) = lMu * lN(iNode);
//         lTx(iNode) = lKtx * (lUx(iNode) - lUxr(iNode));
//         lTy(iNode) = lKty * (lUy(iNode) - lUyr(iNode));
//         
//         if (std::abs(lTx(iNode)) > lCoul(iNode))
//         {
//             lTx(iNode) = Sgn(lTx(iNode)) * lCoul(iNode);
//             lUxr(iNode) = lUx(iNode) - lTx(iNode) / lKtx;
//         }
//         if (std::abs(lTy(iNode)) > lCoul(iNode))
//         {
//             lTy(iNode) = Sgn(lTy(iNode)) * lCoul(iNode);
//             lUyr(iNode) = lUy(iNode) - lTy(iNode) / lKty;
//         }
//         
//         lReturnResult.FValues(iNode * lDpN) = lTx(iNode);
//         lReturnResult.FValues(iNode * lDpN + 1) = lTy(iNode);
//         lReturnResult.FValues(iNode * lDpN + 2) = lN(iNode) - lN0;
//         
//         lReturnResult.XCorr(iNode * lDpN)     = lUxr(iNode);
//         lReturnResult.XCorr(iNode * lDpN + 1) = lUyr(iNode);
//         lReturnResult.XCorr(iNode * lDpN + 2) = aXPrev(iNode * lDpN + 2); // not sure about this, but it probably doesn't matter
//     }
//     
// //     std::cout << "Displacements      : " << aX << std::endl;
// //     std::cout << "Displacements corr : " << lReturnResult.XCorr << std::endl;
// //     std::cout << "Displacements prev : " << aXPrev << std::endl;
// //     std::cout << "Forces             : " << lReturnResult.FValues << std::endl;
// //     
// //     STOP
//     
//     return lReturnResult;
// }
// 
// std::vector<NOX::LAPACK::Vector> GetForcesForMotion(const std::vector<NOX::LAPACK::Vector> aX, const double& aNormalStaticLoad)
// {
//     std::vector<NOX::LAPACK::Vector> lReturnVector;
//     
//     NOX::LAPACK::Vector lXPrev = GetAverage(aX);
// //     NOX::LAPACK::Vector lXPrev = aX[0];
//     
//     std::cout << "Starting X: " << lXPrev << std::endl;
//     
//     for (int i = 0; i < aX.size(); i++)
//     {
//         FResult lResult = ComputeFTimeDomain(aX[i], lXPrev, aNormalStaticLoad);
//         
//         lReturnVector.push_back(lResult.FValues);
//         
//         lXPrev = lResult.XCorr;
//     }
//     
//     return lReturnVector;
// }
