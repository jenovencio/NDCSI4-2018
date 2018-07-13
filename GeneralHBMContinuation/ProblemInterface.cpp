
#include <cmath>
#include <functional>

#include "ProblemInterface.h"
#include "Functions.h"
#include "Misc.h"

ProblemInterface::ProblemInterface(const Config& aConfig)
 : cDOFCount(aConfig.DOFCount), cHarmonicCount(2 * aConfig.HarmonicWaveCount - 1), cDOFCountHBM(aConfig.DOFCount * (2 * aConfig.HarmonicWaveCount - 1))
{
    mMassMatrix = LoadSquareMatrix(aConfig.MassMatrixFilePath);
    mDampingMatrix = LoadSquareMatrix(aConfig.DampingMatrixFile);
    mStiffnessMatrix = LoadSquareMatrix(aConfig.StiffnessMatrixFile);
    
    mInitGuess = NOX::LAPACK::Vector(cDOFCountHBM);
        
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
    
    mBProducts = new double[lIntPointCount * lIntPointCount * cHarmonicCount];
    
    std::function<double(double)> lB1;
    std::function<double(double)> lB2;
    
    for (int i = 0; i < cHarmonicCount; i++)
    {
        if (i == 0)             lB1 = [](double aX) { return 2 * PI * std::cos(aX); };
        else if (i % 2 == 1)    lB1 = [](double aX) { return 2 * PI * std::cos(aX); };
        else                    lB1 = [](double aX) { return 2 * PI * std::sin(aX); };
        
        for (int j = i; j < cHarmonicCount; j++)
        {
            if (j == 0)             lB2 = [](double aX) { return 2 * PI * std::cos(aX); };
            else if (j % 2 == 1)    lB2 = [](double aX) { return 2 * PI * std::cos(aX); };
            else                    lB2 = [](double aX) { return 2 * PI * std::sin(aX); };
        
            for (int iIntPoint = 0; iIntPoint < lIntPointCount; iIntPoint++)
            {
                int lIndex = GetBProductIndex(i, j, iIntPoint, cHarmonicCount, lIntPointCount);
                double lIntPointPos = mIntPointsRelative[iIntPoint];
                
                double lProduct = lB1(lIntPointPos) * lB2(lIntPointPos);
                mBProducts[lIndex] = lProduct;
            }
        }
    }
    
    std::cout << "Problem interface successfully initialised" << std::endl;
}

ProblemInterface::~ProblemInterface()
{
    if (mBProducts != nullptr)
        delete[] mBProducts;
}


const NOX::LAPACK::Vector& ProblemInterface::getInitialGuess()
{    
    return mInitGuess;
}
bool ProblemInterface::computeF(NOX::LAPACK::Vector& aRhs, const NOX::LAPACK::Vector& aX)
{
}
bool ProblemInterface::computeJacobian(NOX::LAPACK::Matrix<double>& aJ, const NOX::LAPACK::Vector& aX)
{
}
void ProblemInterface::setParams(const LOCA::ParameterVector& aParams)
{
    double lNewFrequency = aParams.getValue("frequency");
    if (lNewFrequency != mFrequency)
    {
        mFrequency = lNewFrequency;
        mFrequencyChanged = true;
    }    
}

NOX::LAPACK::Matrix<double> ProblemInterface::CreateDynamicStiffnessMatrix(double aFrequency)
{
    NOX::LAPACK::Matrix<double> lReturnMatrix(cDOFCountHBM, cDOFCountHBM);
    
    
    
    return lReturnMatrix;
}

