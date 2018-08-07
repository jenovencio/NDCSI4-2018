
#include "DSBuilderSimple.h"

void DSBuilderSimple::Init(const Config& aConfig, const ProblemParams& aParams)
{
    DSBuilderBase::Init(aConfig, aParams);
    
    if (aConfig.MassMatrices.size() != 1) throw "Exactly 1 mass matrix muse be provided for this dynamic stiffness builder!";
    if (aConfig.DampingMatrices.size() != 1) throw "Exactly 1 damping matrix muse be provided for this dynamic stiffness builder!";
    if (aConfig.StiffnessMatrices.size() != 1) throw "Exactly 1 stiffness matrix muse be provided for this dynamic stiffness builder!";
    
    mMassMatrix = LoadSquareMatrix(aConfig.ConfigFilePath + "/" + aConfig.MassMatrices[0].File, aConfig.MassMatrices[0].Type);
    mDampingMatrix = LoadSquareMatrix(aConfig.ConfigFilePath + "/" + aConfig.DampingMatrices[0].File, aConfig.DampingMatrices[0].Type);
    mStiffnessMatrix = LoadSquareMatrix(aConfig.ConfigFilePath + "/" + aConfig.StiffnessMatrices[0].File, aConfig.StiffnessMatrices[0].Type);
    
    CheckMatrixSize(mDampingMatrix, aParams.DofCountPhysical, "Damping");
    CheckMatrixSize(mMassMatrix, aParams.DofCountPhysical, "Mass");
    CheckMatrixSize(mStiffnessMatrix, aParams.DofCountPhysical, "Stiffness");
}

NOX::LAPACK::Matrix<double>* DSBuilderSimple::CreateDynamicStiffnessMatrixInner(const double& aFrequency)
{
    if (aFrequency <= 0) throw "Frequency must be a positive value! (attempted to set " + std::to_string(aFrequency) + ")";
    
    NOX::LAPACK::Matrix<double>* lReturnMatrix = new NOX::LAPACK::Matrix<double>(mProblemParams.DofCountHBM, mProblemParams.DofCountHBM);
    NOX::LAPACK::Matrix<double>& lReturnMatrixTemp = *lReturnMatrix;
    
//     std::cout << "Harmonic count: " << cHarmonicCount << std::endl;
    
    // period
    double lT = 2 * PI / aFrequency;
    
//     std::cout << "Mass matrix: " << std::endl;
//     std::cout << mMassMatrix << std::endl;
//     std::cout << "Damping matrix: " << std::endl;
//     std::cout << mDampingMatrix << std::endl;
//     std::cout << "Stiffness matrix: " << std::endl;
//     std::cout << mStiffnessMatrix << std::endl;
    
    // row number in the physical domain
    for (int iDof = 0; iDof < mProblemParams.DofCountPhysical; iDof++)
    {
        // harmonic by which we differentiate
        // we only need to iterate over the harmonic index once, because we are differentiating a linear expression,
        // i.e. we will get an "identity"
        for (int iHarm = 0; iHarm < mProblemParams.HarmonicCount; iHarm++)
        {
//             std::cout << "Harm index: " << iHarm << " out of " << cHarmonicCount << std::endl;
            // dynamic matrix row index
            int lRowIndex = GetHBMDofIndex(iDof, iHarm, mProblemParams.HarmonicCount);
            
//             std::cout << "Row index: " << lRowIndex << std::endl;
            
            int lWaveNumber;
            if (iHarm == 0) lWaveNumber = 0;
            else if (iHarm % 2 == 1) lWaveNumber = (iHarm + 1) / 2; // cos wave
            else lWaveNumber = (iHarm + 1) / 2; // sin wave
            
            if (iHarm == 0)
            {
                // DC harmonic
                
                // col number in the physical domain
                for (int jDof = 0; jDof < mProblemParams.DofCountPhysical; jDof++)
                {
                    // index for stiffness contribution
                    int lColIndex = GetHBMDofIndex(jDof, iHarm, mProblemParams.HarmonicCount);
//                     std::cout << "DC col index: " << lColIndex << std::endl;
                    
                    lReturnMatrixTemp(lRowIndex, lColIndex) += lT * mStiffnessMatrix(iDof, jDof);
                }
            }
            else if (iHarm % 2 == 1)
            {
                // cos wave
                
                // col number in the physical domain
                for (int jDof = 0; jDof < mProblemParams.DofCountPhysical; jDof++)
                {
                    // index for mass and stiffness contribution
                    int lColIndex = GetHBMDofIndex(jDof, iHarm, mProblemParams.HarmonicCount);
                    // index for damping contribution
                    int lColIndex2 = GetHBMDofIndex(jDof, iHarm + 1, mProblemParams.HarmonicCount);
                    
//                     std::cout << "cos col index 1: " << lColIndex << std::endl;
//                     std::cout << "cos col index 2: " << lColIndex2 << std::endl;
                    
                    lReturnMatrixTemp(lRowIndex, lColIndex) += lT / 2 * mStiffnessMatrix(iDof, jDof);
                    lReturnMatrixTemp(lRowIndex, lColIndex) -= lT / 2 * aFrequency * aFrequency * lWaveNumber * lWaveNumber * mMassMatrix(iDof, jDof);
                    lReturnMatrixTemp(lRowIndex, lColIndex2) += lT / 2 * aFrequency * lWaveNumber * mDampingMatrix(iDof, jDof);
                }
            }
            else
            {
                // sin wave
                
                // col number in the physical domain
                for (int jDof = 0; jDof < mProblemParams.DofCountPhysical; jDof++)
                {
                    // index for mass and stiffness contribution
                    int lColIndex = GetHBMDofIndex(jDof, iHarm, mProblemParams.HarmonicCount);
                    // index for damping contribution
                    int lColIndex2 = GetHBMDofIndex(jDof, iHarm - 1, mProblemParams.HarmonicCount);
                    
//                     std::cout << "sin col index 1: " << lColIndex << std::endl;
//                     std::cout << "sin col index 2: " << lColIndex2 << std::endl;
                    
                    lReturnMatrixTemp(lRowIndex, lColIndex) += lT / 2 * mStiffnessMatrix(iDof, jDof);
                    lReturnMatrixTemp(lRowIndex, lColIndex) -= lT / 2 * aFrequency * aFrequency * lWaveNumber * lWaveNumber * mMassMatrix(iDof, jDof);
                    lReturnMatrixTemp(lRowIndex, lColIndex2) -= lT / 2 * aFrequency * lWaveNumber * mDampingMatrix(iDof, jDof);
                }
            }
        }
    }
    
    return lReturnMatrix;
}
