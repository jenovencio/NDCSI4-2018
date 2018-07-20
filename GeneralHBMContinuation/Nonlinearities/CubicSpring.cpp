
#include <set>

#include "CubicSpring.h"
#include "../Functions.h"
#include "../Misc.h"

FResult CubicSpring::ComputeFTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const
{
    FResult lReturnResult;
    lReturnResult.FValues = NOX::LAPACK::Vector(aX.length());
    
    for (int iSpring = 0; iSpring < mSprings.size(); iSpring++)
    {
        int lDof1 = mSprings[iSpring].Dof1Index;
        int lDof2 = mSprings[iSpring].Dof2Index;
        
        double lStiffCoeff = mSprings[iSpring].StiffnessCoeff;
        double lDisp1;
        double lDisp2;
        
        if (lDof1 >= 0) lDisp1 = aX(lDof1);
        else lDisp1 = 0;
        
        if (lDof2 >= 0) lDisp2 = aX(lDof2);
        else lDisp2 = 0;
        
        double lForce1 = lStiffCoeff * (lDisp1 - lDisp2) * (lDisp1 - lDisp2) * (lDisp1 - lDisp2);
        double lForce2 = -lForce1;
        
        if (lDof1 >= 0) lReturnResult.FValues(lDof1) += lForce1;
        if (lDof2 >= 0) lReturnResult.FValues(lDof2) += lForce2;
    }
    
    return lReturnResult;
}

NOX::LAPACK::Matrix<double> CubicSpring::ComputeJacobianTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const
{
    NOX::LAPACK::Matrix<double> lReturnMatrix(aX.length(), aX.length());
    
    for (int iSpring = 0; iSpring < mSprings.size(); iSpring++)
    {
        int lDof1 = mSprings[iSpring].Dof1Index;
        int lDof2 = mSprings[iSpring].Dof2Index;
        
        double lStiffCoeff = mSprings[iSpring].StiffnessCoeff;
        
        double lDisp1;
        double lDisp2;
        
        if (lDof1 >= 0) lDisp1 = aX(lDof1);
        else lDisp1 = 0;
        
        if (lDof2 >= 0) lDisp2 = aX(lDof2);
        else lDisp2 = 0;
        
        double lDiff11 = 3 * lStiffCoeff * (lDisp1 - lDisp2) * (lDisp1 - lDisp2);
        double lDiff12 = -lDiff11;
//         double lDiff22 = 3 * lStiffCoeff * (lDisp2 - lDisp1) * (lDisp2 - lDisp1);
        double lDiff22 = lDiff11;
        double lDiff21 = -lDiff22;
        
        if (lDof1 >= 0) lReturnMatrix(lDof1, lDof1) += lDiff11;
        if (lDof1 >= 0 && lDof2 >= 0) lReturnMatrix(lDof1, lDof2) += lDiff12;
        if (lDof1 >= 0 && lDof2 >= 0) lReturnMatrix(lDof2, lDof1) += lDiff21;
        if (lDof2 >= 0) lReturnMatrix(lDof2, lDof2) += lDiff22;
    }
    
    return lReturnMatrix;
}

int CubicSpring::NumberOfPrepLoops() const
{
    return 0;
}

std::vector<int> CubicSpring::NonzeroFPositions() const
{
    std::set<int> lUnique;
    std::vector<int> lReturnVector;
    
    for (int iSpring = 0; iSpring < mSprings.size(); iSpring++)
    {
        int lDof1 = mSprings[iSpring].Dof1Index;
        int lDof2 = mSprings[iSpring].Dof2Index;
        
        if (lDof1 >= 0)
        {
            const bool lIsIn1 = lUnique.find(lDof1) != lUnique.end();
            if (!lIsIn1)
            {
                lReturnVector.push_back(lDof1);
                lUnique.insert(lDof1);
            }
        }
        
        if (lDof2 >= 0)
        {
            const bool lIsIn2 = lUnique.find(lDof2) != lUnique.end();
            if (!lIsIn2)
            {
                lReturnVector.push_back(lDof2);
                lUnique.insert(lDof2);
            }
        }
    }
    
    return lReturnVector;
}
bool CubicSpring::IsCorrectingX() const
{
    return false;
}

void CubicSpring::AddCubicSpring(const CubicSpringDef& aDef)
{
    if (IsFinalised()) throw "The nonlinearity is finalised, no modifications are allowed at this point!";
    
    if (aDef.Dof1Index < 0 && aDef.Dof2Index < 0) throw "At least one of the spring dof indices must be non negative!";
    if (aDef.Dof1Index == aDef.Dof2Index) throw "Dof indices are the same!";
    
    if (aDef.StiffnessCoeff < 0) throw "Cubic stiffness can not be negative!";
    if (aDef.StiffnessCoeff == 0) return;
    
    mSprings.push_back(aDef);
}

void CubicSpring::AddCubicSpring(const int& aDof1Index, const int& aDof2Index, const double& aStiffnessCoeff)
{    
    CubicSpringDef lNewSpring;
    lNewSpring.Dof1Index = aDof1Index;
    lNewSpring.Dof2Index = aDof2Index;
    lNewSpring.StiffnessCoeff = aStiffnessCoeff;
    
    AddCubicSpring(lNewSpring);
}
void CubicSpring::ClearSprings()
{
    if (IsFinalised()) throw "The nonlinearity is finalised, no modifications are allowed at this point!";
    mSprings.clear();
}
void CubicSpring::LoadFromFile(const std::string& aFilePath)
{
    
    std::ifstream lInputFile(aFilePath);
    if (!lInputFile.is_open()) throw "Unable to open file \"" + aFilePath + "\"!";
    
    while (!lInputFile.eof())
    {
        int lDof1;
        int lDof2;
        double lStiffCoeff;
        
        lInputFile >> lDof1 >> lDof2 >> lStiffCoeff;
        
        AddCubicSpring(lDof1, lDof2, lStiffCoeff);
    }
    
    lInputFile.close();
}
std::string CubicSpring::ClassName() const
{
    return "Cubic spring";
}
