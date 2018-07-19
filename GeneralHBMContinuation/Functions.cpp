
#include <algorithm>

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
    
    lReturnConfig.MassMatrixFile = GetNextValidLine(lInputFile);
    lReturnConfig.MassMatrixType = GetNextValidLine(lInputFile);
    CheckMatrixType(lReturnConfig.MassMatrixType);
    
    lReturnConfig.DampingMatrixFile = GetNextValidLine(lInputFile);
    lReturnConfig.DampingMatrixType = GetNextValidLine(lInputFile);
    CheckMatrixType(lReturnConfig.DampingMatrixType);
    
    lReturnConfig.StiffnessMatrixFile = GetNextValidLine(lInputFile);
    lReturnConfig.StiffnessMatrixType = GetNextValidLine(lInputFile);
    CheckMatrixType(lReturnConfig.StiffnessMatrixType);
    
    lReturnConfig.ExcitationForceFile = GetNextValidLine(lInputFile);
    lReturnConfig.ContinuationSettingsFile = GetNextValidLine(lInputFile);
    
    lTempString = GetNextValidLine(lInputFile);
    lReturnConfig.HarmonicWaveCount = std::stoi(lTempString);
    
    lTempString = GetNextValidLine(lInputFile);
    lReturnConfig.FrequencyMin = std::stod(lTempString);
    
    lTempString = GetNextValidLine(lInputFile);
    lReturnConfig.FrequencyMax = std::stod(lTempString);
    
    lTempString = GetNextValidLine(lInputFile);
    lReturnConfig.IntPointCount = std::stoi(lTempString);
    
    lTempString = GetNextValidLine(lInputFile);
    int lSaveWholeInt = std::stoi(lTempString);
    
    if (lSaveWholeInt == 1) lReturnConfig.SaveWholeSolutions = true;
    else if (lSaveWholeInt == 0) lReturnConfig.SaveWholeSolutions = false;
    else throw "Invalid bool flag for \"SaveWholeSolutions\", only 0 or 1 are accepted!";
    
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
    std::cout << "Mass      matrix file: " << aConfig.MassMatrixFile << " (" << aConfig.MassMatrixType << ")" << std::endl;
    std::cout << "Damping   matrix file: " << aConfig.DampingMatrixFile << " (" << aConfig.DampingMatrixType << ")" << std::endl;
    std::cout << "Stiffness matrix file: " << aConfig.StiffnessMatrixFile << " (" << aConfig.StiffnessMatrixType << ")" << std::endl;
    std::cout << "Excitation force file: " << aConfig.ExcitationForceFile << std::endl;
    std::cout << "Save whole solutions: " << (aConfig.SaveWholeSolutions ? "true" : "false") << std::endl;
    std::cout << BORDER << std::endl;
}
ContinuationSettings LoadContinuationSettings(const std::string& aContSetFilePath)
{
    std::ifstream lInputFile(aContSetFilePath);
    if (!lInputFile.is_open()) throw std::string("Unable to open file \"" + aContSetFilePath + "\"!");
    
    std::string lTempString;
    ContinuationSettings lReturnSettings;
    
    int MaxStepsContinuation;
    int MaxStepsNewton;
    double StepSizeInitial;
    double StepSizeMin;
    double StepSizeMax;
    
    std::string PredictorType; // secant, constant, tangent
    
    
    lTempString = GetNextValidLine(lInputFile);
    lReturnSettings.MaxStepsContinuation = std::stoi(lTempString);
    
    lTempString = GetNextValidLine(lInputFile);
    lReturnSettings.MaxStepsNewton = std::stoi(lTempString);
    
    lTempString = GetNextValidLine(lInputFile);
    lReturnSettings.StepSizeInitial = std::stod(lTempString);
    
    lTempString = GetNextValidLine(lInputFile);
    lReturnSettings.StepSizeMin = std::stod(lTempString);
    
    lTempString = GetNextValidLine(lInputFile);
    lReturnSettings.StepSizeMax = std::stod(lTempString);
        
    lTempString = GetNextValidLine(lInputFile);
    
    if (! (lTempString == "Constant" || lTempString == "Tangent" || lTempString == "Secant")) 
        throw "Invalid predictor option! Valid strings (case sensitive): Constant, Tangent, Secant";
    
    lReturnSettings.PredictorType = lTempString;
    
    return lReturnSettings;
}

void PrintContinuationSettings(const ContinuationSettings& aContSettings)
{
    std::cout << BORDER << std::endl;
    std::cout << "Continuation settings: " << std::endl;
    std::cout << "Max continuation steps: " << aContSettings.MaxStepsContinuation << std::endl;
    std::cout << "Max newton steps: " << aContSettings.MaxStepsNewton << std::endl;
    std::cout << "Init step: " << aContSettings.StepSizeInitial << std::endl;
    std::cout << "Min step: " << aContSettings.StepSizeMin << std::endl;
    std::cout << "Max step: " << aContSettings.StepSizeMax << std::endl;
    std::cout << "Predictor type: " << aContSettings.PredictorType << std::endl;
    std::cout << BORDER << std::endl;
}

NOX::LAPACK::Matrix<double> LoadSquareMatrixFull(const std::string& aFilePath, int& aDim)
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
NOX::LAPACK::Matrix<double> LoadSquareMatrixSparse(const std::string& aFilePath, int& aDim)
{
    std::ifstream lInputFile(aFilePath);
    if (!lInputFile.is_open()) throw "Unable to open file \"" + aFilePath + "\"!";
    
    int lDim;
    lInputFile >> lDim;
    
    aDim = lDim;
        
    NOX::LAPACK::Matrix<double> lReturnMatrix(lDim, lDim);
    
    int lRow;
    int lCol;
    double lValue;
    
    int lEntryNumber = 0;
    
    do
    {
        lInputFile >> lRow;
        lInputFile >> lCol;
        
        if (lRow < 0 || lRow >= lDim) throw "Invalid row value! (entry number " + std::to_string(lEntryNumber) + ")";
        if (lCol < 0 || lCol >= lDim) throw "Invalid col value! (entry number " + std::to_string(lEntryNumber) + ")";
        
        lInputFile >> lValue;
        
        lReturnMatrix(lRow, lCol) = lValue;
        
        lEntryNumber++;
    }
    while (!lInputFile.eof());
    
    return lReturnMatrix;
}
NOX::LAPACK::Matrix<double> LoadSquareMatrix(const std::string& aFilePath, const std::string& aMatrixType, int& aDim)
{
    if (aMatrixType == FULL_STRING)
        return LoadSquareMatrixFull(aFilePath, aDim);
    else if (aMatrixType == SPARSE_STRING)
        return LoadSquareMatrixSparse(aFilePath, aDim);
    else throw "Unknown matrix type: " + aMatrixType + " for file \"" + aFilePath + "\"!";
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

void CheckMatrixType(const std::string& aType)
{
    if (aType != FULL_STRING && aType != SPARSE_STRING)
        throw "String \"" + aType + "\" is not a valid matrix type! Valid types are: " + FULL_STRING + ", " + SPARSE_STRING;
}
