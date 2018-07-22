
#include <algorithm>

#include "PenaltyWall.h"
#include "../Misc.h"

FResult PenaltyWall::ComputeFTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const
{
    FResult lReturnResult;
    lReturnResult.FValues = NOX::LAPACK::Vector(aX.length());
    
    for (int iWall = 0; iWall < mWalls.size(); iWall++)
    {
        int lDof = mWalls[iWall].DofIndex;
        double lStiffCoeff = mWalls[iWall].DirectionalStiffness;
        
        double lDiff = aX(lDof) - mWalls[iWall].Offset;
        
        if (lDiff * lStiffCoeff >= 0) continue;
        
//         std::cout << "Penetration amount: " << std::abs(lDiff) << std::endl;
//         std::cout << "Result force      : " << std::abs(lDiff * lDiff) * lStiffCoeff << std::endl;
//         
//         STOP
        
        lReturnResult.FValues(lDof) += std::abs(lDiff) * lStiffCoeff;
    }
    
    return lReturnResult;
}
NOX::LAPACK::Matrix<double> PenaltyWall::ComputeJacobianTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const
{   
    NOX::LAPACK::Matrix<double> lReturnMatrix(aX.length(), aX.length());
    
    for (int iWall = 0; iWall < mWalls.size(); iWall++)
    {
        int lDof = mWalls[iWall].DofIndex;
        double lStiffCoeff = mWalls[iWall].DirectionalStiffness;
        
        double lDiff = aX(lDof) - mWalls[iWall].Offset;
        
        if (lDiff * lStiffCoeff >= 0) continue;
        
        lReturnMatrix(lDof, lDof) += lStiffCoeff;
    }
    
    return lReturnMatrix;
}

int PenaltyWall::NumberOfPrepLoops() const
{
    return 0;
}
std::vector<int> PenaltyWall::NonzeroFPositions() const
{
    std::set<int> lUnique;
    std::vector<int> lReturnVector;
    
    for (int iSpring = 0; iSpring < mWalls.size(); iSpring++)
    {
        int lDof = mWalls[iSpring].DofIndex;
        
        if (lDof >= 0)
        {
            const bool lIsIn = lUnique.find(lDof) != lUnique.end();
            if (!lIsIn)
            {
                lReturnVector.push_back(lDof);
                lUnique.insert(lDof);
            }
        }
    }
    
    return lReturnVector;
}
bool PenaltyWall::IsCorrectingX() const
{
    return false;
}

void PenaltyWall::LoadFromFile(const std::string& aFilePath)
{
    std::ifstream lInputFile(aFilePath);
    if (!lInputFile.is_open()) throw "Unable to open file \"" + aFilePath + "\"!";
    
    while (!lInputFile.eof())
    {
        int lDof;
        double lStiffCoeff;
        double lOffset;
        
        lInputFile >> lDof >> lStiffCoeff >> lOffset;
        
        AddWall(lDof, lStiffCoeff, lOffset);
    }
    
    lInputFile.close();
}
std::string PenaltyWall::ClassName() const
{
    return "Penalty Wall";
}
void PenaltyWall::AddWall(const PenaltyWallDef& aWall)
{
    if (IsFinalised()) throw "The object is already finalised, can not add another wall definition!";
    if (aWall.DofIndex < 0) throw "Dof index must be non-negative!";
    
    mWalls.push_back(aWall);
}
void PenaltyWall::AddWall(const int& aDofIndex, const double& aDirStiffness, const double& aOffset)
{
    PenaltyWallDef lDef;
    lDef.DofIndex = aDofIndex;
    lDef.DirectionalStiffness = aDirStiffness;
    lDef.Offset = aOffset;
    
    AddWall(lDef);
}
