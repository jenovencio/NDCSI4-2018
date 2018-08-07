
#include "DSBuilderBlock.h"


void DSBuilderBlock::Init(const Config& aConfig, const ProblemParams& aParams)
{
    DSBuilderBase::Init(aConfig, aParams);
    
    mHarmonicWaveCount = aConfig.HarmonicWaveCount;
    
    if (aConfig.MassMatrices.size() != mHarmonicWaveCount - 1) throw "Exactly " + std::to_string(mHarmonicWaveCount - 1) + " mass matrices must be provided for this dynamic stiffness builder!";
    if (aConfig.DampingMatrices.size() != mHarmonicWaveCount - 1 && aConfig.DampingMatrices.size() != 0) throw "Either 0 or exactly " + std::to_string(mHarmonicWaveCount - 1) + " damping matrices must be provided for this dynamic stiffness builder!";
    if (aConfig.StiffnessMatrices.size() != mHarmonicWaveCount) throw "Exactly " + std::to_string(mHarmonicWaveCount) + " stiffness matrices must be provided for this dynamic stiffness builder!";
    
    // load the DC stiffness matrix
    mStiffnessMatrixDC = LoadSquareMatrix(aConfig.ConfigFilePath + "/" + aConfig.StiffnessMatrices[0].File, aConfig.StiffnessMatrices[0].Type);
    
    // load matrices for nonzero harmonics
    for (int iHarmWave = 1; iHarmWave < mHarmonicWaveCount; iHarmWave++)
    {
        NOX::LAPACK::Matrix<double> lMassMatrix = LoadSquareMatrix(aConfig.ConfigFilePath + "/" + aConfig.MassMatrices[iHarmWave - 1].File, aConfig.MassMatrices[iHarmWave - 1].Type);
        CheckMatrixSize(lMassMatrix, 2 * aParams.DofCountPhysical, "Mass");
        mMassMatricesDyn.push_back(lMassMatrix);
        
        // damping matrices are optional
        if (aConfig.DampingMatrices.size() > 0)
        {
            NOX::LAPACK::Matrix<double> lDampingMatrix = LoadSquareMatrix(aConfig.ConfigFilePath + "/" + aConfig.DampingMatrices[iHarmWave - 1].File, aConfig.DampingMatrices[iHarmWave - 1].Type);
            CheckMatrixSize(lDampingMatrix, 2 * aParams.DofCountPhysical, "Damping");
            mDampingMatricesDyn.push_back(lDampingMatrix);
        }
        
        NOX::LAPACK::Matrix<double> lStiffnessMatrix = LoadSquareMatrix(aConfig.ConfigFilePath + "/" + aConfig.StiffnessMatrices[iHarmWave].File, aConfig.StiffnessMatrices[iHarmWave].Type);
        CheckMatrixSize(lStiffnessMatrix, 2 * aParams.DofCountPhysical, "Stiffness");
        mStiffnessMatricesDyn.push_back(lStiffnessMatrix);
    }
    
    mDummyDamping = NOX::LAPACK::Matrix<double>(2 * aParams.DofCountPhysical, 2 * aParams.DofCountPhysical);
}

NOX::LAPACK::Matrix<double>* DSBuilderBlock::CreateDynamicStiffnessMatrixInner(const double& aFrequency)
{
    if (aFrequency <= 0) throw "Frequency must be a positive value! (attempted to set " + std::to_string(aFrequency) + ")";
    
    NOX::LAPACK::Matrix<double>* lReturnMatrix = new NOX::LAPACK::Matrix<double>(mProblemParams.DofCountHBM, mProblemParams.DofCountHBM);
    NOX::LAPACK::Matrix<double>& lReturnMatrixTemp = *lReturnMatrix;
    
//     std::cout << "Harmonic count: " << cHarmonicCount << std::endl;
    
    // period
    double lT = 2 * PI / aFrequency;
        
    // row number in the physical domain
    for (int iDof = 0; iDof < mProblemParams.DofCountPhysical; iDof++)
    {
        // 0th harmoninc
        for (int jDof = 0; jDof < mProblemParams.DofCountPhysical; jDof++)
        {
            int lRow = GetHBMDofIndex(iDof, 0, mProblemParams.HarmonicCount);
            int lCol = GetHBMDofIndex(jDof, 0, mProblemParams.HarmonicCount);
            
            lReturnMatrixTemp(lRow, lCol) += lT * mStiffnessMatrixDC(iDof, jDof);
        }
        
        // nonzero harmonics
        for (int iHarmWave = 1; iHarmWave < mHarmonicWaveCount; iHarmWave++)
        {            
            NOX::LAPACK::Matrix<double>& lMassMatrix = mMassMatricesDyn[iHarmWave - 1];
            NOX::LAPACK::Matrix<double>& lDampingMatrix = mDummyDamping;
            if (mDampingMatricesDyn.size() > 0) lDampingMatrix = mDampingMatricesDyn[iHarmWave - 1];
            NOX::LAPACK::Matrix<double>& lStiffnessMatrix = mStiffnessMatricesDyn[iHarmWave - 1];
                        
            std::function<double(double, double, double)> lFunc = [aFrequency, iHarmWave, lT](double aMass, double aDamp, double aStiff)
            {
                double lValue =
                    aStiff + 
                    aDamp * aFrequency * iHarmWave - 
                    aMass * aFrequency * aFrequency * iHarmWave * iHarmWave;
                    
                return lT / 2 * lValue;
            };
            
            // left top block
            for (int jDof = 0; jDof < mProblemParams.DofCountPhysical; jDof++)
            {
                int lRowBlock = iDof;
                int lColBlock = jDof;
                
                int lRowHBM = GetHBMDofIndex(iDof, 2 * iHarmWave - 1, mProblemParams.HarmonicCount);
                int lColHBM = GetHBMDofIndex(jDof, 2 * iHarmWave - 1, mProblemParams.HarmonicCount);
                
                lReturnMatrixTemp(lRowHBM, lColHBM) += 
                    lFunc(
                        lMassMatrix(lRowBlock, lColBlock),
                        lDampingMatrix(lRowBlock, lColBlock),
                        lStiffnessMatrix(lRowBlock, lColBlock));
            }
            
            // left bottom block
            for (int jDof = 0; jDof < mProblemParams.DofCountPhysical; jDof++)
            {
                int lRowBlock = iDof + mProblemParams.DofCountPhysical;
                int lColBlock = jDof;
                
                int lRowHBM = GetHBMDofIndex(iDof, 2 * iHarmWave, mProblemParams.HarmonicCount);
                int lColHBM = GetHBMDofIndex(jDof, 2 * iHarmWave - 1, mProblemParams.HarmonicCount);
                
                lReturnMatrixTemp(lRowHBM, lColHBM) += 
                    lFunc(
                        lMassMatrix(lRowBlock, lColBlock),
                        lDampingMatrix(lRowBlock, lColBlock),
                        lStiffnessMatrix(lRowBlock, lColBlock));
            }
            
            // right top block
            for (int jDof = 0; jDof < mProblemParams.DofCountPhysical; jDof++)
            {
                int lRowBlock = iDof;
                int lColBlock = jDof + mProblemParams.DofCountPhysical;
                
                int lRowHBM = GetHBMDofIndex(iDof, 2 * iHarmWave - 1, mProblemParams.HarmonicCount);
                int lColHBM = GetHBMDofIndex(jDof, 2 * iHarmWave, mProblemParams.HarmonicCount);
                
                lReturnMatrixTemp(lRowHBM, lColHBM) += 
                    lFunc(
                        lMassMatrix(lRowBlock, lColBlock),
                        lDampingMatrix(lRowBlock, lColBlock),
                        lStiffnessMatrix(lRowBlock, lColBlock));
            }
            
            // right bottom block
            for (int jDof = 0; jDof < mProblemParams.DofCountPhysical; jDof++)
            {
                int lRowBlock = iDof + mProblemParams.DofCountPhysical;
                int lColBlock = jDof + mProblemParams.DofCountPhysical;
                
                int lRowHBM = GetHBMDofIndex(iDof, 2 * iHarmWave, mProblemParams.HarmonicCount);
                int lColHBM = GetHBMDofIndex(jDof, 2 * iHarmWave, mProblemParams.HarmonicCount);
                
                lReturnMatrixTemp(lRowHBM, lColHBM) += 
                    lFunc(
                        lMassMatrix(lRowBlock, lColBlock),
                        lDampingMatrix(lRowBlock, lColBlock),
                        lStiffnessMatrix(lRowBlock, lColBlock));
            }
        }
    }
    
    return lReturnMatrix;
}
