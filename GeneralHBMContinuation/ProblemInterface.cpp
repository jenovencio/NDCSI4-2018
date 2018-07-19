
#include <cmath>
#include <functional>

#include "ProblemInterface.h"
#include "Functions.h"
#include "Misc.h"
#include "Nonlinearities/NonlinearitiesFactory.h"

const std::string ProblemInterface::cFrequencyName = "frequency";

ProblemInterface::ProblemInterface(const Config& aConfig, const std::vector<NonlinearBase*>& aNonlinearities, bool aSaveWholeSolutions)
 : cHarmonicCount(2 * aConfig.HarmonicWaveCount - 1), cNonlinearitiesOuter(aNonlinearities), cSaveWholeSolutions(aSaveWholeSolutions)
{
    int lDummy;
    
    mMassMatrix = LoadSquareMatrix(aConfig.ConfigFilePath + "/" + aConfig.MassMatrixFile, aConfig.MassMatrixType, mDOFCount);
    mDampingMatrix = LoadSquareMatrix(aConfig.ConfigFilePath + "/" + aConfig.DampingMatrixFile, aConfig.DampingMatrixType, lDummy);
    if (mDOFCount != lDummy) throw "Mass and damping matrices have different dimensions!";
    mStiffnessMatrix = LoadSquareMatrix(aConfig.ConfigFilePath + "/" + aConfig.StiffnessMatrixFile, aConfig.StiffnessMatrixType, lDummy);
    if (mDOFCount != lDummy) throw "Mass and stiffness matrices have different dimensions!";
    
    // ensure it's a negative value, so the function fills it with the value from the file
    mExcitationCoeffCount = -1;
    mExcitationAmp = LoadExcitationForce(aConfig.ConfigFilePath + "/" + aConfig.ExcitationForceFile, lDummy, mExcitationCoeffCount);
    if (mDOFCount != lDummy) throw "Mass matrix and excitation force have different dimensions!";
    if (mExcitationCoeffCount > cHarmonicCount) throw "Number of excitation force coeffs (" + std::to_string(mExcitationCoeffCount) + ") can not be larger than number of coeffs used for the HBM (" + std::to_string(cHarmonicCount) + ")!";
    
    mDOFCountHBM = cHarmonicCount * mDOFCount;
    
    mInitGuess = NOX::LAPACK::Vector(mDOFCountHBM);
        
    int lIntPointCount = aConfig.IntPointCount;
    if (lIntPointCount < 1) throw "Integration point count must be a positive number!";
    
    double lRelativeStep = 1.0 / lIntPointCount;
    
    mIntPointsRelative.reserve(lIntPointCount);
    
    double lCurrentPos = lRelativeStep / 2;
    
    for (int i = 0; i < lIntPointCount; i++)
    {
        mIntPointsRelative.push_back(lCurrentPos);
        lCurrentPos += lRelativeStep;
    }
    
    double* lTempArray1 = new double[lIntPointCount * cHarmonicCount];
    
    std::function<double(double)> lB;
    
    for (int i = 0; i < cHarmonicCount; i++)
    {
        int lWaveNumber = (i + 1) / 2;
        
        if (i == 0)             lB = [](double aX) { return 1.0; };
        else if (i % 2 == 1)    lB = [lWaveNumber](double aX) { return std::cos(2 * PI * lWaveNumber * aX); };
        else                    lB = [lWaveNumber](double aX) { return std::sin(2 * PI * lWaveNumber * aX); };
        
        for (int iIntPoint = 0; iIntPoint < lIntPointCount; iIntPoint++)
        {
            int lIndex = GetBValuesIndex(i, iIntPoint, cHarmonicCount, lIntPointCount);
            double lIntPointPos = mIntPointsRelative[iIntPoint];
            
            double lValue = lB(lIntPointPos);
            
            lTempArray1[lIndex] = lValue;
        }
    }
    mBValues.assign(lTempArray1, lTempArray1 + lIntPointCount * cHarmonicCount);
    
    double* lTempArray2 = new double[lIntPointCount * cHarmonicCount * cHarmonicCount];
    
    std::function<double(double)> lB1;
    std::function<double(double)> lB2;
    
    for (int i = 0; i < cHarmonicCount; i++)
    {
        int lWaveNumber = (i + 1) / 2;
        
        if (i == 0)             lB1 = [](double aX) { return 1.0; };
        else if (i % 2 == 1)    lB1 = [lWaveNumber](double aX) { return std::cos(2 * PI * lWaveNumber * aX); };
        else                    lB1 = [lWaveNumber](double aX) { return std::sin(2 * PI * lWaveNumber * aX); };
        
        for (int j = 0; j < cHarmonicCount; j++)
        {
            int lWaveNumber2 = (j + 1) / 2;
            
            if (j == 0)             lB2 = [](double aX) { return 1.0; };
            else if (j % 2 == 1)    lB2 = [lWaveNumber2](double aX) { return std::cos(2 * PI * lWaveNumber2 * aX); };
            else                    lB2 = [lWaveNumber2](double aX) { return std::sin(2 * PI * lWaveNumber2 * aX); };
            
            for (int iIntPoint = 0; iIntPoint < lIntPointCount; iIntPoint++)
            {
                int lIndex = GetBProductIndex(i, j, iIntPoint, cHarmonicCount, lIntPointCount);
                double lIntPointPos = mIntPointsRelative[iIntPoint];
                
                double lProduct = lB1(lIntPointPos) * lB2(lIntPointPos);
                lTempArray2[lIndex] = lProduct;
            }
        }
    }
    mBProducts.assign(lTempArray2, lTempArray2 + lIntPointCount * cHarmonicCount * cHarmonicCount);
    
    // initialise the outer nonlinearities with data from this class
    // also finalise them
    for (int i = 0; i < cNonlinearitiesOuter.size(); i++)
    {
        cNonlinearitiesOuter[i]->Init(mIntPointsRelative, mBValues, mBProducts, cHarmonicCount, mDOFCount);
        cNonlinearitiesOuter[i]->Finalise();
        
        // add the outer nonlinearity into the "all nonlinearities" vector
        mNonlinearities.push_back(cNonlinearitiesOuter[i]);
    }
    
    // load and init the nonlinearities from files (specified in the config)
    
    std::vector<NonlinearityDefinition> lNonlinDef = LoadNonlinearitiesDefinitions(aConfig.ConfigFilePath + "/" + aConfig.NonlinearitiesFile);
    
    for (int i = 0; i < lNonlinDef.size(); i++)
    {
        std::string lFilePath = lNonlinDef[i].File;
        std::string lType = lNonlinDef[i].Type;
        
        NonlinearBase* lNewNonlin = C_NonlinearitiesFactory.at(lType)();
        
        lNewNonlin->Init(mIntPointsRelative, mBValues, mBProducts, cHarmonicCount, mDOFCount);
        lNewNonlin->LoadFromFile(lFilePath);
        lNewNonlin->Finalise();
        
        mNonlinearitiesInner.push_back(lNewNonlin);
        
        // add the inner nonlinearity into the "all nonlinearities" vector
        mNonlinearities.push_back(lNewNonlin);
    }
    
    std::cout << "Problem interface successfully initialised" << std::endl;
    std::cout << BORDER << std::endl;
    std::cout << "Problem: " << std::endl;
    std::cout << "Number of physical DOFs: " << mDOFCount << std::endl;
    std::cout << "Total number of DOFs: " << mDOFCountHBM << std::endl;
    std::cout << "Number of nonlinearities loaded from files: " << lNonlinDef.size() << std::endl;
    std::cout << "Number of nonlinearities added from external code: " << cNonlinearitiesOuter.size() << std::endl;
    std::cout << "Total number of nonlinearities: " << mNonlinearities.size() << std::endl;
    std::cout << BORDER << std::endl;
}

ProblemInterface::~ProblemInterface()
{
    if (mDynamicStiffnessMatrix != nullptr)
        delete mDynamicStiffnessMatrix;
    if (mExcitationRHS != nullptr)
        delete mExcitationRHS;
    
    // delete the inner nonlinearities (loaded from files)
    for (int i = 0; i < mNonlinearitiesInner.size(); i++)
        delete mNonlinearitiesInner[i];
}


const NOX::LAPACK::Vector& ProblemInterface::getInitialGuess()
{
    return mInitGuess;
}
bool ProblemInterface::computeF(NOX::LAPACK::Vector& aRhs, const NOX::LAPACK::Vector& aX)
{
    if (mRecomputeDynamicStiffness)
    {
        delete mDynamicStiffnessMatrix;
        mDynamicStiffnessMatrix = nullptr;
        mDynamicStiffnessMatrix = CreateDynamicStiffnessMatrix(mFrequency);
        mRecomputeDynamicStiffness = false;
    }
    
    if (mRecomputeExcitationRHS)
    {
        delete mExcitationRHS;
        mExcitationRHS = nullptr;
        mExcitationRHS = CreateExcitationRHS(mFrequency);
        mRecomputeExcitationRHS = false;
    }
    
    // linear part
    aRhs = NOX::LAPACK::Vector(*mExcitationRHS);
        
    aRhs = aRhs.scale(-1.0);
    
    // dynamic stiffness matrix multiplication
    for (int i = 0; i < mDynamicStiffnessMatrix->numRows(); i++)
        for (int j = 0; j < mDynamicStiffnessMatrix->numCols(); j++)
            aRhs(i) += (*mDynamicStiffnessMatrix)(i, j) * aX(j);
        
    // nonlinear part
        
    for (int iNonlin = 0; iNonlin < mNonlinearities.size(); iNonlin++)
    {
        NOX::LAPACK::Vector lNonlinContrib = mNonlinearities[iNonlin]->ComputeF(aX, mFrequency);
        
        for (int i = 0; i < aRhs.length(); i++)
            aRhs(i) += lNonlinContrib(i);
    }
    
    return true;
}
bool ProblemInterface::computeJacobian(NOX::LAPACK::Matrix<double>& aJ, const NOX::LAPACK::Vector& aX)
{    
    if (mRecomputeDynamicStiffness)
    {
        delete mDynamicStiffnessMatrix;
        mDynamicStiffnessMatrix = nullptr;
        mDynamicStiffnessMatrix = CreateDynamicStiffnessMatrix(mFrequency);
        mRecomputeDynamicStiffness = false;
    }
    
    // linear part
    aJ = NOX::LAPACK::Matrix<double>(*mDynamicStiffnessMatrix);
    
    // nonlinear part
    
    for (int iNonlin = 0; iNonlin < mNonlinearities.size(); iNonlin++)
    {
        NOX::LAPACK::Matrix<double> lNonlinContrib = mNonlinearities[iNonlin]->ComputeJacobian(aX, mFrequency);
        
        for (int j = 0; j < aJ.numCols(); j++)
            for (int i = 0; i < aJ.numRows(); i++)
                aJ(i, j) += lNonlinContrib(i, j);
    }
    
    return true;
}
void ProblemInterface::setParams(const LOCA::ParameterVector& aParams)
{
    double lNewFrequency = aParams.getValue(cFrequencyName);
    if (lNewFrequency != mFrequency)
    {
        mFrequency = lNewFrequency;
        mRecomputeDynamicStiffness = true;
        mRecomputeExcitationRHS = true;
    }
}
void ProblemInterface::printSolution(const NOX::LAPACK::Vector& aX, const double aConParam)
{
    double lNorm = aX.norm();
    
    mSolutionFrequencies.push_back(aConParam);
    mSolutionNorms.push_back(lNorm);
    
    if (cSaveWholeSolutions) 
    {
        mSolutions.push_back(aX);
        mHasWholeSolutions = true;
    }
}
void ProblemInterface::ClearSolutions()
{
    mSolutionFrequencies.clear();
    mSolutionNorms.clear();
    mSolutions.clear();
    
    mHasWholeSolutions = false;
}
void ProblemInterface::WriteSolutionNorms(std::ostream& aStream) const
{
    for (int i = 0; i < mSolutionFrequencies.size(); i++)
    {
        aStream << mSolutionFrequencies[i] << "; " << mSolutionNorms[i] << std::endl;
    }
}
void ProblemInterface::WriteWholeSolutions(std::ostream& aStream) const
{
    if (!cSaveWholeSolutions) throw "Whole solutions were not saved, so they can not be outputted!";
    
    for (int i = 0; i < mSolutions.size(); i++)
    {
        aStream << mSolutionFrequencies[i] << "; " << mSolutions[i] << std::endl;
    }
}
bool ProblemInterface::HasWholeSolutions() const
{
    return mHasWholeSolutions;
}

NOX::LAPACK::Matrix<double>* ProblemInterface::CreateDynamicStiffnessMatrix(double aFrequency)
{
    if (aFrequency <= 0) throw "Frequency must be a positive value! (attempted to set " + std::to_string(aFrequency) + ")";
    
    NOX::LAPACK::Matrix<double>* lReturnMatrix = new NOX::LAPACK::Matrix<double>(mDOFCountHBM, mDOFCountHBM);
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
    for (int iDof = 0; iDof < mDOFCount; iDof++)
    {
        // harmonic by which we differentiate
        // we only need to iterate over the harmonic index once, because we are differentiating a linear expression,
        // i.e. we will get an "identity"
        for (int iHarm = 0; iHarm < cHarmonicCount; iHarm++)
        {
//             std::cout << "Harm index: " << iHarm << " out of " << cHarmonicCount << std::endl;
            // dynamic matrix row index
            int lRowIndex = GetHBMDofIndex(iDof, iHarm, cHarmonicCount);
            
//             std::cout << "Row index: " << lRowIndex << std::endl;
            
            if (iHarm == 0)
            {
                // DC harmonic
                
                // col number in the physical domain
                for (int jDof = 0; jDof < mDOFCount; jDof++)
                {
                    // index for stiffness contribution
                    int lColIndex = GetHBMDofIndex(jDof, iHarm, cHarmonicCount);
//                     std::cout << "DC col index: " << lColIndex << std::endl;
                    
                    lReturnMatrixTemp(lRowIndex, lColIndex) += lT * mStiffnessMatrix(iDof, jDof);
                }
            }
            else if (iHarm % 2 == 1)
            {
                // cos wave
                int lWaveNumber = (iHarm + 1) / 2;
                
                // col number in the physical domain
                for (int jDof = 0; jDof < mDOFCount; jDof++)
                {
                    // index for mass and stiffness contribution
                    int lColIndex = GetHBMDofIndex(jDof, iHarm, cHarmonicCount);
                    // index for damping contribution
                    int lColIndex2 = GetHBMDofIndex(jDof, iHarm + 1, cHarmonicCount);
                    
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
                int lWaveNumber = iHarm / 2;
                
                // col number in the physical domain
                for (int jDof = 0; jDof < mDOFCount; jDof++)
                {
                    // index for mass and stiffness contribution
                    int lColIndex = GetHBMDofIndex(jDof, iHarm, cHarmonicCount);
                    // index for damping contribution
                    int lColIndex2 = GetHBMDofIndex(jDof, iHarm - 1, cHarmonicCount);
                    
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
NOX::LAPACK::Vector* ProblemInterface::CreateExcitationRHS(double aFrequency)
{
    if (aFrequency <= 0) throw "Frequency must be a positive value! (attempted to set " + std::to_string(aFrequency) + ")";
    
    NOX::LAPACK::Vector* lReturnVector = new NOX::LAPACK::Vector(mDOFCountHBM);
    NOX::LAPACK::Vector& lReturnVectorTemp = *lReturnVector;
    
    // period
    double lT = 2 * PI / aFrequency;
    
    for (int iDof = 0; iDof < mDOFCount; iDof++)
    {
        // The mExcitationCoeffCount is equal or less than cHarmonicCount
        for (int iHarm = 0; iHarm < mExcitationCoeffCount; iHarm++)
        {
            int lExcitationIndex = GetHBMDofIndex(iDof, iHarm, mExcitationCoeffCount);
            int lSystemIndex = GetHBMDofIndex(iDof, iHarm, cHarmonicCount);
            
            double lForceAmp = mExcitationAmp[lExcitationIndex];
            
            if (iHarm == 0) lReturnVectorTemp(lSystemIndex) = lForceAmp * lT;
            else lReturnVectorTemp(lSystemIndex) = lForceAmp * lT / 2;
        }
    }
    
    return lReturnVector;
}
