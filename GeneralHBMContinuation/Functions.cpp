
#include <algorithm>
#include <string>

#include "Functions.h"
#include "Misc.h"
#include "Nonlinearities/NonlinearitiesFactory.h"

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
    CheckString(lReturnConfig.MassMatrixType, C_MatrixTypes, "Matrix type");
    
    lReturnConfig.DampingMatrixFile = GetNextValidLine(lInputFile);
    lReturnConfig.DampingMatrixType = GetNextValidLine(lInputFile);
    CheckString(lReturnConfig.DampingMatrixType, C_MatrixTypes, "Matrix type");
    
    lReturnConfig.StiffnessMatrixFile = GetNextValidLine(lInputFile);
    lReturnConfig.StiffnessMatrixType = GetNextValidLine(lInputFile);
    CheckString(lReturnConfig.StiffnessMatrixType, C_MatrixTypes, "Matrix type");
    
    lReturnConfig.ExcitationForceFile = GetNextValidLine(lInputFile);
    lReturnConfig.ContinuationSettingsFile = GetNextValidLine(lInputFile);
    lReturnConfig.NonlinearitiesFile = GetNextValidLine(lInputFile);
    
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
    std::cout << "Nonlinearities file: " << aConfig.NonlinearitiesFile << std::endl;
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
    
    lInputFile.close();
    
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
    
    while (!lInputFile.eof())
    {
        lInputFile >> lRow;
        lInputFile >> lCol;
        
        if (lRow < 0 || lRow >= lDim) throw "Invalid row value! (entry number " + std::to_string(lEntryNumber) + ")";
        if (lCol < 0 || lCol >= lDim) throw "Invalid col value! (entry number " + std::to_string(lEntryNumber) + ")";
        
        lInputFile >> lValue;
        
        lReturnMatrix(lRow, lCol) = lValue;
        
        lEntryNumber++;
    }
    
    lInputFile.close();
    
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
    
    lInputFile.close();
    
    return lReturnVector;
}
std::vector<NonlinearityDefinition> LoadNonlinearitiesDefinitions(const std::string& aFilePath)
{
    std::ifstream lInputFile(aFilePath);
    if (!lInputFile.is_open()) throw "Unable to open file \"" + aFilePath + "\"!";
    
    std::cout << "Opened file: " << aFilePath << std::endl;
    
    std::vector<NonlinearityDefinition> lReturnVector;
    
    while (!lInputFile.eof())
    {
        std::string lFile = GetNextValidLine(lInputFile);
        
        // hack to deal with an empty file
        if (lFile == "" && lInputFile.eof()) break;
        
        std::string lType = GetNextValidLine(lInputFile);
        
        CheckString(lType, C_NonlinearitiesFactory, "Nonlinearity type");
        
        NonlinearityDefinition lDef;
        lDef.File = lFile;
        lDef.Type = lType;
        
        lReturnVector.push_back(lDef);
    }
    
    lInputFile.close();
    std::cout << "Closed file: " << aFilePath << std::endl;
    
    return lReturnVector;
}
std::string GetNextValidLine(std::ifstream& aFile)
{
    std::string lReturnValue = "";
    
    bool lSkipLine = true;
    
    while (lSkipLine && !aFile.eof())
    {
        lSkipLine = false;
        std::getline(aFile, lReturnValue);
        if (lReturnValue.size() == 0 || lReturnValue[0] == '#') lSkipLine = true;
    }
    
    return lReturnValue;
}

void CheckString(const std::string& aString, const std::vector<std::string>& aPossibilities, const std::string& aGroupName)
{
    bool lIsIn = std::find(aPossibilities.begin(), aPossibilities.end(), aString) != aPossibilities.end();
    
    if (!lIsIn)
    {
        std::stringstream lStringBuilder;
        lStringBuilder << "Value \"" + aString + "\" is not a valid option for \"" + aGroupName + "\"! Valid options are: " << std::endl;
        
        if (aPossibilities.size() == 0) lStringBuilder << " none (you are fucked)";
        
        for (int i = 0; i < aPossibilities.size(); i++)
        {
            lStringBuilder << aPossibilities[i];
            if (i < aPossibilities.size() - 1)
                lStringBuilder << ", ";
        }
        
        lStringBuilder << std::endl;
        
        throw lStringBuilder.str();
    }
}

template <class T>
void CheckString(const std::string& aString, const std::map<std::string, T>& aPossibilities, const std::string& aGroupName)
{
    bool lIsIn = aPossibilities.find(aString) != aPossibilities.end();
    
    if (!lIsIn)
    {
        std::stringstream lStringBuilder;
        lStringBuilder << "Value \"" + aString + "\" is not a valid option for \"" + aGroupName + "\"! Valid options are: " << std::endl;
        
        if (aPossibilities.size() == 0) lStringBuilder << " none (you are fucked)";
        
        int lCount = 0;
        for (auto nIt = aPossibilities.begin(); nIt != aPossibilities.end(); nIt++)
        {
            lStringBuilder << nIt->first;
            if (lCount < aPossibilities.size() - 1)
                lStringBuilder << ", ";
            
            lCount++;
        }
        
        lStringBuilder << std::endl;
        
        throw lStringBuilder.str();
    }
}
