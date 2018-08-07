
#include <cmath>
#include <functional>
#include <limits>

#include "matplotlibcpp.h"

#include "ProblemInterface.h"
#include "Functions.h"
#include "Misc.h"
#include "Nonlinearities/NonlinearitiesFactory.h"
#include "Aft/AftFactory.h"
#include "DSBuilder/DSBuilderFactory.h"

const std::string ProblemInterface::cContParameterName = "generic_parameter";

ProblemInterface::ProblemInterface(const Config& aConfig, const std::vector<NonlinearBase*>& aNonlinearities)
 : cConfig(aConfig), cProblemParams(aConfig.DofCount, 2 * aConfig.HarmonicWaveCount - 1), cNonlinearitiesOuter(aNonlinearities)
{
//     int lDummy;
//     
//     CheckMatrixCount("mass", aConfig.MassMatrices.size(), aConfig.HarmonicWaveCount);
//     CheckMatrixCount("damping", aConfig.DampingMatrices.size(), aConfig.HarmonicWaveCount);
//     CheckMatrixCount("stiffness", aConfig.StiffnessMatrices.size(), aConfig.HarmonicWaveCount);
//     
//     for (int i = 0; i < aConfig.MassMatrices.size(); i++)
//     {
//         NOX::LAPACK::Matrix<double> lMatrix = LoadSquareMatrix(aConfig.ConfigFilePath + "/" + aConfig.MassMatrices[i].File, aConfig.MassMatrices[i].Type);
//         
//         CheckMatrixSize(lMatrix, "Mass matrix");
//         
//         mMassMatrices.push_back(lMatrix);
//     }
//     for (int i = 0; i < aConfig.DampingMatrices.size(); i++)
//     {
//         NOX::LAPACK::Matrix<double> lMatrix = LoadSquareMatrix(aConfig.ConfigFilePath + "/" + aConfig.DampingMatrices[i].File, aConfig.DampingMatrices[i].Type);
//         
//         CheckMatrixSize(lMatrix, "Damping matrix");
//         
//         mDampingMatrices.push_back(lMatrix);
//     }
//     for (int i = 0; i < aConfig.StiffnessMatrices.size(); i++)
//     {
//         NOX::LAPACK::Matrix<double> lMatrix = LoadSquareMatrix(aConfig.ConfigFilePath + "/" + aConfig.StiffnessMatrices[i].File, aConfig.StiffnessMatrices[i].Type);
//         
//         CheckMatrixSize(lMatrix, "Stiffness matrix");
//         
//         mStiffnessMatrices.push_back(lMatrix);
//     }
    
    // ensure it's a negative value, so the function fills it with the value from the file
    mExcitationCoeffCount = -1;
    mExcitationAmp = LoadExcitationForce(aConfig.ConfigFilePath + "/" + aConfig.ExcitationForceFile, mExcitationCoeffCount);
    if (cConfig.DofCount != mExcitationAmp.size() / mExcitationCoeffCount) throw "Excitation force's number of dofs is different from number of dofs specified in the configuration!";
    if (mExcitationCoeffCount > cProblemParams.HarmonicCount) throw "Number of excitation force coeffs (" + std::to_string(mExcitationCoeffCount) + ") can not be larger than number of coeffs used for the HBM (" + std::to_string(cProblemParams.HarmonicCount) + ")!";
        
    mInitGuess = NOX::LAPACK::Vector(cProblemParams.DofCountHBM);
    
    int lIntPointCount = aConfig.IntPointCount;
    if (lIntPointCount < 1) throw "Integration point count must be a positive number!";
    
    // create the aft object
    mAft = C_AftFactory.at(aConfig.AftType)(lIntPointCount, cProblemParams);
    mDSBuilder = C_DSBuilderFactory.at(aConfig.DSBuilderType)();
    mDSBuilder->Init(aConfig, cProblemParams);
    
    // initialise the outer nonlinearities with data from this class
    // also finalise them
    for (int i = 0; i < cNonlinearitiesOuter.size(); i++)
    {
        cNonlinearitiesOuter[i]->Init(mAft, cProblemParams);
        cNonlinearitiesOuter[i]->Finalise();
        
        // add the outer nonlinearity into the "all nonlinearities" vector
        mNonlinearities.push_back(cNonlinearitiesOuter[i]);
    }
    
    // load and init the nonlinearities from files (specified in the config)
    
    for (int i = 0; i < aConfig.Nonlinearities.size(); i++)
    {
        std::string lFilePath = aConfig.Nonlinearities[i].File;
        std::string lType = aConfig.Nonlinearities[i].Type;
        
        NonlinearBase* lNewNonlin = C_NonlinearitiesFactory.at(lType)();
        
        lNewNonlin->Init(mAft, cProblemParams);
        lNewNonlin->LoadFromFile(aConfig.ConfigFilePath + "/" + lFilePath);
        lNewNonlin->Finalise();
        
        mNonlinearitiesInner.push_back(lNewNonlin);
        
        // add the inner nonlinearity into the "all nonlinearities" vector
        mNonlinearities.push_back(lNewNonlin);
    }
    
    std::cout << "Problem interface successfully initialised" << std::endl;
    std::cout << BORDER << std::endl;
    std::cout << "Problem: " << std::endl;
    std::cout << "Number of physical DOFs: " << cProblemParams.DofCountPhysical << std::endl;
    std::cout << "Total number of DOFs: " << cProblemParams.DofCountHBM << std::endl;
    std::cout << "Number of nonlinearities loaded from files: " << aConfig.Nonlinearities.size() << std::endl;
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
    
    if (mAft != nullptr)
        delete mAft;
    
    if (mDSBuilder != nullptr)
        delete mDSBuilder;
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
    
    const NOX::LAPACK::Matrix<double>& lDSMTemp = *mDynamicStiffnessMatrix;
    
    // dynamic stiffness matrix multiplication
    for (int i = 0; i < lDSMTemp.numRows(); i++)
        for (int j = 0; j < lDSMTemp.numCols(); j++)
            aRhs(i) += lDSMTemp(i, j) * aX(j);
        
    // nonlinear part
        
    for (int iNonlin = 0; iNonlin < mNonlinearities.size(); iNonlin++)
    {
        const NOX::LAPACK::Vector& lNonlinContrib = mNonlinearities[iNonlin]->ComputeF(aX, mFrequency);
        
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
        const NOX::LAPACK::Matrix<double>& lNonlinContrib = mNonlinearities[iNonlin]->ComputeJacobian(aX, mFrequency);
        
        for (int j = 0; j < aJ.numCols(); j++)
            for (int i = 0; i < aJ.numRows(); i++)
                aJ(i, j) += lNonlinContrib(i, j);
    }
    
    return true;
}
void ProblemInterface::setParams(const LOCA::ParameterVector& aParams)
{
    double lNewContParam = aParams.getValue(cContParameterName);
    if (lNewContParam != mCurrentContParam)
    {
        mCurrentContParam = lNewContParam;
        mFrequency = mCurrentContParam * (cConfig.FrequencyEnd - cConfig.FrequencyStart) + cConfig.FrequencyStart;
        
//         std::cout << "New continuation parameter set: " << mCurrentContParam << std::endl;
//         std::cout << "Corresponding frequency value : " << mFrequency << std::endl;
        
        mRecomputeDynamicStiffness = true;
        mRecomputeExcitationRHS = true;
    }
}
void ProblemInterface::printSolution(const NOX::LAPACK::Vector& aX, const double aConParam)
{
    double lFreq = aConParam * cConfig.FrequencyEnd + (1.0 - aConParam) * cConfig.FrequencyStart;
    
    double lNorm = aX.norm();
    
    mSolutionFrequencies.push_back(lFreq);
    mSolutionNorms.push_back(lNorm);
    
    if (cConfig.SaveWholeSolutions) 
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
    aStream.precision(std::numeric_limits<double>::max_digits10);
    
    for (int i = 0; i < mSolutionFrequencies.size(); i++)
    {
        aStream << mSolutionFrequencies[i] << "; " << mSolutionNorms[i] << std::endl;
    }
}
void ProblemInterface::WriteWholeSolutions(std::ostream& aStream) const
{
    if (!cConfig.SaveWholeSolutions) throw "Whole solutions were not saved, so they can not be outputted!";
    
    aStream.precision(std::numeric_limits<double>::max_digits10);
    
    for (int i = 0; i < mSolutions.size(); i++)
    {
        aStream << mSolutionFrequencies[i] << "; " << mSolutions[i] << std::endl;
    }
}
void ProblemInterface::PlotSolutionNorms() const
{
    matplotlibcpp::plot(mSolutionFrequencies, mSolutionNorms);
    matplotlibcpp::show();
}
bool ProblemInterface::HasWholeSolutions() const
{
    return mHasWholeSolutions;
}

NOX::LAPACK::Matrix<double>* ProblemInterface::CreateDynamicStiffnessMatrix(double aFrequency)
{
    if (aFrequency <= 0) throw "Frequency must be a positive value! (attempted to set " + std::to_string(aFrequency) + ")";
    NOX::LAPACK::Matrix<double>* lReturnMatrix = mDSBuilder->CreateDynamicStiffnessMatrix(aFrequency);
    return lReturnMatrix;
}
NOX::LAPACK::Vector* ProblemInterface::CreateExcitationRHS(double aFrequency)
{
    if (aFrequency <= 0) throw "Frequency must be a positive value! (attempted to set " + std::to_string(aFrequency) + ")";
    
    NOX::LAPACK::Vector* lReturnVector = new NOX::LAPACK::Vector(cProblemParams.DofCountHBM);
    NOX::LAPACK::Vector& lReturnVectorTemp = *lReturnVector;
    
    // period
    double lT = 2 * PI / aFrequency;
    
    for (int iDof = 0; iDof < cProblemParams.DofCountPhysical; iDof++)
    {
        // The mExcitationCoeffCount is equal or less than cHarmonicCount
        for (int iHarm = 0; iHarm < mExcitationCoeffCount; iHarm++)
        {
            int lExcitationIndex = GetHBMDofIndex(iDof, iHarm, mExcitationCoeffCount);
            int lSystemIndex = GetHBMDofIndex(iDof, iHarm, cProblemParams.HarmonicCount);
            
            double lForceAmp = mExcitationAmp[lExcitationIndex];
            
            if (iHarm == 0) lReturnVectorTemp(lSystemIndex) = lForceAmp * lT;
            else lReturnVectorTemp(lSystemIndex) = lForceAmp * lT / 2;
        }
    }
    
    return lReturnVector;
}

void ProblemInterface::CheckMatrixCount(const std::string& aMatrixName, const int& aMatrixCount, const int& aHarmonicWaveCount)
{
    if (aMatrixCount != 1 && aMatrixCount != aHarmonicWaveCount)
        throw "Number of \"" + aMatrixName + "\" matrices must be either 1 or equal to number of harmonic waves (\"" + std::to_string(aHarmonicWaveCount) + "\"). Instead, \"" + std::to_string(aMatrixCount) + "\" matrices were provided!";
}
void ProblemInterface::CheckMatrixSize(const NOX::LAPACK::Matrix<double>& aMatrix, const std::string& aMatrixName)
{
    if (aMatrix.numCols() != cConfig.DofCount)
        throw "Number of columns of \"" + aMatrixName + "\" does not equal to the number of dofs specified in the configuration!";
    
    if (aMatrix.numRows() != cConfig.DofCount)
        throw "Number of rows of \"" + aMatrixName + "\" does not equal to the number of dofs specified in the configuration!";
}
