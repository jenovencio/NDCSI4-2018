
#include <set>

#include "CubicSpring.h"
#include "../Functions.h"
#include "../Misc.h"

NOX::LAPACK::Vector CubicSpring::ComputeFTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const
{
    NOX::LAPACK::Vector lReturnVector(aX.length());
    
    for (int iSpring = 0; iSpring < mSprings.size(); iSpring++)
    {
        int lDof = mSprings[iSpring].DofIndex;
        double lStiffCoeff = mSprings[iSpring].StiffnessCoeff;
        double lDisp = aX(lDof) * aX(lDof) * aX(lDof);
        
        lReturnVector(lDof) += lStiffCoeff * lDisp;
    }
    
    return lReturnVector;
}

NOX::LAPACK::Matrix<double> CubicSpring::ComputeJacobianTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const
{
    NOX::LAPACK::Matrix<double> lReturnMatrix(aX.length(), aX.length());
    
    for (int iSpring = 0; iSpring < mSprings.size(); iSpring++)
    {
        int lDof = mSprings[iSpring].DofIndex;
        double lStiffCoeff = mSprings[iSpring].StiffnessCoeff;
        double lDisp = aX(lDof) * aX(lDof) * aX(lDof);
        
        lReturnMatrix(lDof, lDof) += 3 * lStiffCoeff * lDisp * lDisp;
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
        int lDof = mSprings[iSpring].DofIndex;
        const bool lIsIn = lUnique.find(lDof) != lUnique.end();
        
        if (!lIsIn)
        {
            lReturnVector.push_back(lDof);
            lUnique.insert(lDof);
        }
    }
    
    return lReturnVector;
}

void CubicSpring::AddCubicSpring(const CubicSpringDef& aDef)
{
    if (IsFinalised()) throw "The nonlinearity is finalised, no modifications are allowed at this point!";
    
    if (aDef.DofIndex < 0) throw "Dof index can not be negative!";
    if (aDef.StiffnessCoeff < 0) throw "Cubic stiffness can not be negative!";
    if (aDef.StiffnessCoeff == 0) return;
    
    mSprings.push_back(aDef);
}

void CubicSpring::AddCubicSpring(const int& aDofIndex, const double& aStiffnessCoeff)
{
    if (IsFinalised()) throw "The nonlinearity is finalised, no modifications are allowed at this point!";
    CubicSpringDef lNewSpring;
    lNewSpring.DofIndex = aDofIndex;
    lNewSpring.StiffnessCoeff = aStiffnessCoeff;
    
    AddCubicSpring(lNewSpring);
}
void CubicSpring::ClearSprings()
{
    if (IsFinalised()) throw "The nonlinearity is finalised, no modifications are allowed at this point!";
    mSprings.clear();
}
