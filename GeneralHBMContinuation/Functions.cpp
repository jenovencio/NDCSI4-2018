
#include "Functions.h"
#include "Misc.h"

Config LoadConfig(const std::string& aConfigFilePath)
{
    std::ifstream lInputFile(aConfigFilePath);
    if (!lInputFile.is_open()) throw std::string("Unable to open file \"" + aConfigFilePath + "\"!");
    
    std::string lTempString;
    Config lReturnConfig;
    
    int lSeparator = aConfigFilePath.size() - 1;
    while (lSeparator > 0 && aConfigFilePath[lSeparator] != '/') lSeparator--;
    
    if (lSeparator == 0)
    {
        lReturnConfig.ConfigFilePath = ".";
        lReturnConfig.ConfigFileName = aConfigFilePath;
    }
    else
    {
        lReturnConfig.ConfigFilePath = aConfigFilePath.substr(0, lSeparator + 1);
        lReturnConfig.ConfigFileName = aConfigFilePath.substr(lSeparator + 1, aConfigFilePath.size() - lSeparator - 1);
    }
    
    lReturnConfig.MassMatrixFilePath = GetNextValidLine(lInputFile);
    lReturnConfig.DampingMatrixFile = GetNextValidLine(lInputFile);
    lReturnConfig.StiffnessMatrixFile = GetNextValidLine(lInputFile);
    lReturnConfig.ExcitationForceFile = GetNextValidLine(lInputFile);
    
    lTempString = GetNextValidLine(lInputFile);
    lReturnConfig.HarmonicWaveCount = std::stoi(lTempString);
    
    lTempString = GetNextValidLine(lInputFile);
    lReturnConfig.FrequencyMin = std::stod(lTempString);
    
    lTempString = GetNextValidLine(lInputFile);
    lReturnConfig.FrequencyMax = std::stod(lTempString);
    
    lTempString = GetNextValidLine(lInputFile);
    lReturnConfig.IntPointCount = std::stoi(lTempString);
    
    return lReturnConfig;
}
void PrintConfig(const Config& aConfig)
{   
    std::cout << BORDER << std::endl;
    std::cout << "Config: " << std::endl;
    std::cout << "Config path: " << aConfig.ConfigFilePath << std::endl;
    std::cout << "Config file name: " << aConfig.ConfigFileName << std::endl;
    std::cout << "Harmonic waves count: " << aConfig.HarmonicWaveCount << std::endl;
    std::cout << "Frequency range: <" << aConfig.FrequencyMin << ", " << aConfig.FrequencyMax << ">" << std::endl;
    std::cout << "Number of time integration points: " << aConfig.IntPointCount << std::endl;
    std::cout << "Mass matrix file: " << aConfig.MassMatrixFilePath << std::endl;
    std::cout << "Damping matrix file: " << aConfig.DampingMatrixFile << std::endl;
    std::cout << "Stiffness matrix file: " << aConfig.StiffnessMatrixFile << std::endl;
    std::cout << "Excitation force file: " << aConfig.ExcitationForceFile << std::endl;
    std::cout << BORDER << std::endl;
}
NOX::LAPACK::Matrix<double> LoadSquareMatrix(const std::string& aFilePath, int& aDim)
{
    std::ifstream lInputFile(aFilePath);
    if (!lInputFile.is_open()) throw "Unable to open file \"" + aFilePath + "\"!";
    
    int lDim;
    
    lInputFile >> lDim;
    
    aDim = lDim;
    
    NOX::LAPACK::Matrix<double> lReturnMatrix(lDim, lDim);
    
    int lRow = 0;
    int lCol = 0;
    for (int i = 0; i < lDim * lDim; i++)
    {
        lInputFile >> lReturnMatrix(lRow, lCol);
        
        lRow++;
        if (lRow == lDim)
        {
            lCol++;
            lRow = 0;
        }
    }
    
    return lReturnMatrix;
}
std::vector<double> LoadExcitationForce(const std::string& aFilePath, int& aDim, int& aHarmCoeffCount)
{
    std::ifstream lInputFile(aFilePath);
    if (!lInputFile.is_open()) throw "Unable to open file \"" + aFilePath + "\"!";
    
    int lDim;
    lInputFile >> lDim;
    
    aDim = lDim;
        
    int lHarmCount;
    lInputFile >> lHarmCount;
    
    aHarmCoeffCount = lHarmCount;
    
    std::vector<double> lReturnVector(lDim * lHarmCount);
    
    for (int iDof = 0; iDof < lDim; iDof++)
    {
        for (int iHarm = 0; iHarm < lHarmCount; iHarm++)
        {
            int lIndex = GetHBMDofIndex(iDof, iHarm, lHarmCount);
            lInputFile >> lReturnVector[lIndex];
        }
    }
    
    return lReturnVector;
}
std::string GetNextValidLine(std::ifstream& aFile)
{
    std::string lReturnValue = "";
    
    bool lIsComment = true;
    
    while (lIsComment)
    {
        lIsComment = false;
        std::getline(aFile, lReturnValue);
        if (lReturnValue.size() == 0 || lReturnValue[0] == '#') lIsComment = true;
    }
    
    return lReturnValue;
}
