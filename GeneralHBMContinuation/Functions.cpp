
#include "Functions.h"
#include <fstream>

Config LoadConfig(const std::string& aConfigFilePath)
{
    std::ifstream lInputFile(aConfigFilePath);
    if (!lInputFile.is_open()) throw std::string("Unable to open file \"" + aConfigFilePath + "\"!");
    
    std::string lTempString;
    Config lReturnConfig;
    
    std::getline(lInputFile, lTempString);
    lReturnConfig.DOFCount = std::stoi(lTempString);
    
    std::getline(lInputFile, lReturnConfig.MassMatrixFilePath);
    std::getline(lInputFile, lReturnConfig.DampingMatrixFile);
    std::getline(lInputFile, lReturnConfig.StiffnessMatrixFile);
    
    std::getline(lInputFile, lTempString);
    lReturnConfig.HarmonicWaveCount = std::stoi(lTempString);
    
    std::getline(lInputFile, lTempString);
    lReturnConfig.FrequencyMin = std::stoi(lTempString);
    
    std::getline(lInputFile, lTempString);
    lReturnConfig.FrequencyMax = std::stoi(lTempString);
    
    std::getline(lInputFile, lTempString);
    lReturnConfig.IntPointCount = std::stoi(lTempString);
    
    return lReturnConfig;
}
void PrintConfig(const Config& aConfig)
{
    std::string lBorder = "-----------------------------------------------------------------------------------------";
    std::cout << lBorder << std::endl;
    std::cout << "Config: " << std::endl;
    std::cout << "Physical dof count: " << aConfig.DOFCount << std::endl;
    std::cout << "Harmonic waves count: " << aConfig.HarmonicWaveCount << std::endl;
    std::cout << "Frequency range: <" << aConfig.FrequencyMin << ", " << aConfig.FrequencyMax << ">" << std::endl;
    std::cout << "Number of time integration points: " << aConfig.IntPointCount << std::endl;
    std::cout << "Mass matrix file: " << aConfig.MassMatrixFilePath << std::endl;
    std::cout << "Damping matrix file: " << aConfig.DampingMatrixFile << std::endl;
    std::cout << "Stiffness matrix file: " << aConfig.StiffnessMatrixFile << std::endl;
    
    std::cout << std::endl;
    
    std::cout << "Total HBM dof count: " << aConfig.DOFCount * (2 * aConfig.HarmonicWaveCount - 1) << std::endl;
    std::cout << lBorder << std::endl;
}
NOX::LAPACK::Matrix<double> LoadSquareMatrix(const std::string& aFilePath)
{
    std::ifstream lInputFile(aFilePath);
    if (!lInputFile.is_open()) throw "Unable to open file \"" + aFilePath + "\"!";
    
    int lDim;
    
    lInputFile >> lDim;
    
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
