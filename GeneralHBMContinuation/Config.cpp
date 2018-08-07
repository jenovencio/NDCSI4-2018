
#include "Config.h"
#include "Misc.h"
#include "Functions.h"
#include "Functions.hpp"
#include "Nonlinearities/NonlinearitiesFactory.h"
#include "Aft/AftFactory.h"
#include "DSBuilder/DSBuilderFactory.h"

const std::string Config::KEY_MASS      = "mass";
const std::string Config::KEY_DAMP      = "damp";
const std::string Config::KEY_STIFF     = "stiff";
const std::string Config::KEY_EXC       = "exc";
const std::string Config::KEY_DOF       = "dof";
const std::string Config::KEY_HARM      = "harm";
const std::string Config::KEY_NONLIN    = "nonlin";
const std::string Config::KEY_START     = "start";
const std::string Config::KEY_END       = "end";
const std::string Config::KEY_INT       = "int";
const std::string Config::KEY_WHOLE     = "whole";
const std::string Config::KEY_STEPSC    = "stepscont";
const std::string Config::KEY_STEPSN    = "stepsnewt";
const std::string Config::KEY_NORM      = "normnewt";
const std::string Config::KEY_STEPINI   = "stepinit";
const std::string Config::KEY_STEPMIN   = "stepmin";
const std::string Config::KEY_STEPMAX   = "stepmax";
const std::string Config::KEY_PRED      = "pred";
const std::string Config::KEY_SCALEARC  = "scalearc";
const std::string Config::KEY_AFT       = "aft";
const std::string Config::KEY_DSBUILD   = "dsbuild";

const std::map<std::string, std::function<std::vector<std::string>(const std::string&)>> Config::C_ConfigDefinition = 
{
    { Config::KEY_MASS,     THROW_NO_DEFAULT        },
    { Config::KEY_DAMP,     THROW_NO_DEFAULT        },
    { Config::KEY_STIFF,    THROW_NO_DEFAULT        },
    { Config::KEY_EXC,      THROW_NO_DEFAULT        },
    { Config::KEY_HARM,     THROW_NO_DEFAULT        },
    { Config::KEY_DOF,      THROW_NO_DEFAULT        },
    { Config::KEY_NONLIN,   RET_EMPTY_STR_VEC       },
    { Config::KEY_START,    THROW_NO_DEFAULT        },
    { Config::KEY_END,      THROW_NO_DEFAULT        },
    { Config::KEY_INT,      THROW_NO_DEFAULT        },
    { Config::KEY_WHOLE,    RET_ONE_STR("0")        },
    { Config::KEY_STEPSC,   RET_ONE_STR("1000")     },
    { Config::KEY_STEPSN,   RET_ONE_STR("10")       },
    { Config::KEY_NORM,     RET_ONE_STR("1e-6")     },
    { Config::KEY_STEPINI,  RET_ONE_STR("1e-4")     },
    { Config::KEY_STEPMIN,  RET_ONE_STR("1e-6")     },
    { Config::KEY_STEPMAX,  RET_ONE_STR("1e-1")     },
    { Config::KEY_PRED,     RET_ONE_STR("Tangent")  },
    { Config::KEY_SCALEARC, RET_ONE_STR("1")        },
    { Config::KEY_AFT,      RET_ONE_STR("Simple")   },
    { Config::KEY_DSBUILD,      RET_ONE_STR("Simple")   },
};

Config Config::LoadConfig(const std::string& aConfigFilePath)
{
    std::map<std::string, std::vector<std::string>> lMap = ParseGeneralConfigFile(aConfigFilePath);
    
    std::map<std::string, std::vector<std::string>> lConfigMap;
    
    for (auto nKeyValue : C_ConfigDefinition)
    {
        std::string lKey = nKeyValue.first;
        
        bool lIsInLoadedMap = lMap.find(lKey) != lMap.end();
        
        if (lIsInLoadedMap) lConfigMap[lKey] = lMap[lKey];
        else 
        {
            std::cout << "Key \"" + lKey + "\" not specified in \"" + aConfigFilePath + "\", attempting to use default" << std::endl;
            // get the default value (might raise an exception if the default value does not exist)
            lConfigMap[lKey] = nKeyValue.second(lKey);
            std::cout << "Default value used for the key \"" + lKey + "\": " << std::endl;
            for (int i = 0; i < lConfigMap[lKey].size(); i++)
            {
                std::cout << lConfigMap[lKey][i] << std::endl;
            }
            std::cout << std::endl;
        }
    }
    
    // list loaded keys that will be unused
    // this is just to inform the user that they configuration file contains some unused definitions
    std::cout << "Unused keys in file \"" + aConfigFilePath + "\": " << std::endl;
    bool lAnyUnused = false;
    for (auto nKeyValue : lMap)
    {
        std::string lKey = nKeyValue.first;
        
        bool lIsInUsedMap = lConfigMap.find(lKey) != lConfigMap.end();
        
        if (!lIsInUsedMap) 
        {
            std::cout << lKey << std::endl;
            lAnyUnused = true;
        }
    }
    if (!lAnyUnused)
        std::cout << "--- all keys will be used ---" << std::endl;
    
    Config lReturnConfig;
    
    // get the config path and file name
    int lSeparatorPos = aConfigFilePath.size() - 1;
    while (lSeparatorPos > 0 && aConfigFilePath[lSeparatorPos] != '/') lSeparatorPos--;
    
    if (lSeparatorPos == 0)
    {
        lReturnConfig.ConfigFilePath = ".";
        lReturnConfig.ConfigFileName = aConfigFilePath;
    }
    else
    {
        lReturnConfig.ConfigFilePath = aConfigFilePath.substr(0, lSeparatorPos + 1);
        lReturnConfig.ConfigFileName = aConfigFilePath.substr(lSeparatorPos + 1, aConfigFilePath.size() - lSeparatorPos - 1);
    }
    
    // parse loaded key-value pairs (value is a vector of strings here
    
    std::vector<std::string> lTempLines;
    int lTempLineInd = 0;
    
    // matrices
    lTempLines = lConfigMap[KEY_MASS];
    lReturnConfig.MassMatrices = ParseMatrixDefinition(lTempLines, "Mass");
    lTempLines = lConfigMap[KEY_DAMP];
    lReturnConfig.DampingMatrices = ParseMatrixDefinition(lTempLines, "Damping");
    lTempLines = lConfigMap[KEY_STIFF];
    lReturnConfig.StiffnessMatrices = ParseMatrixDefinition(lTempLines, "Stiffness");
    
    // excitation force
    lTempLines = lConfigMap[KEY_EXC];
    if (lTempLines.size() != 1) throw "Excitation force definition must have exactly one line!";
    lReturnConfig.ExcitationForceFile = lTempLines[0];
    
    // nonlinearities files
    lTempLines = lConfigMap[KEY_NONLIN];
    if (lTempLines.size() % 2 != 0) throw "Nonlinearities defininition must have an even number of lines!";
    lTempLineInd = 0;
    while (lTempLineInd < lTempLines.size())
    {
        NonlinearityDefinition lNewDef;
        lNewDef.File = lTempLines[lTempLineInd++];
        lNewDef.Type = lTempLines[lTempLineInd++];
        CheckString(lNewDef.Type, C_NonlinearitiesFactory, "nonlinearity type");
        lReturnConfig.Nonlinearities.push_back(lNewDef);
    }
    
    // number of physical dofs
    lReturnConfig.DofCount = ParseOneLineInt(lConfigMap[KEY_DOF], "Dof count");
    
    // number of harmonics
    lReturnConfig.HarmonicWaveCount = ParseOneLineInt(lConfigMap[KEY_HARM], "Harmonic count");
    
    // start frequency
    lReturnConfig.FrequencyStart = ParseOneLineDouble(lConfigMap[KEY_START], "Start frequency");
    
    // end frequency
    lReturnConfig.FrequencyEnd = ParseOneLineDouble(lConfigMap[KEY_END], "End frequency");
    
    // integration point count
    lReturnConfig.IntPointCount = ParseOneLineInt(lConfigMap[KEY_INT], "Integration point count");
    
    // save whole solutions
    lReturnConfig.SaveWholeSolutions = ParseOneLineBool(lConfigMap[KEY_WHOLE], "Save whole solutions");
    
    // continuation settings
    
    // max newton steps
    lReturnConfig.MaxStepsNewton = ParseOneLineInt(lConfigMap[KEY_STEPSN], "Newton step count");
    
    // max continuation steps
    lReturnConfig.MaxStepsContinuation = ParseOneLineInt(lConfigMap[KEY_STEPSC], "Continuation step count");
    
    // newton norm
    lReturnConfig.NewtonNormF = ParseOneLineDouble(lConfigMap[KEY_NORM], "Newton tolerance norm");
    
    // init step
    lReturnConfig.StepSizeInitial = ParseOneLineDouble(lConfigMap[KEY_STEPINI], "Initial conntinuation step");
    
    // min step
    lReturnConfig.StepSizeMin = ParseOneLineDouble(lConfigMap[KEY_STEPMIN], "Min conntinuation step");
    
    // max step
    lReturnConfig.StepSizeMax = ParseOneLineDouble(lConfigMap[KEY_STEPMAX], "Max conntinuation step");
    
    // predictor
    lTempLines = lConfigMap[KEY_PRED];
    if (lTempLines.size() != 1) throw "Predictor definition must have exactly one line!";
    lReturnConfig.PredictorType = lTempLines[0];
    CheckString(lReturnConfig.PredictorType, C_PredictorTypes, KEY_PRED);
    
    // arc length scaling
    lReturnConfig.EnableArcLengthScaling = ParseOneLineBool(lConfigMap[KEY_SCALEARC], "Arc length scaling");
    
    // aft type
    lTempLines = lConfigMap[KEY_AFT];
    if (lTempLines.size() != 1) throw "Aft type definition must have exactly one line!";
    lReturnConfig.AftType = lTempLines[0];
    CheckString(lReturnConfig.AftType, C_AftFactory, KEY_AFT);
    
    // dynamics stiffness matrix builder type
    lTempLines = lConfigMap[KEY_DSBUILD];
    if (lTempLines.size() != 1) throw "DS builder type definition must have exactly one line!";
    lReturnConfig.DSBuilderType = lTempLines[0];
    CheckString(lReturnConfig.DSBuilderType, C_DSBuilderFactory, KEY_DSBUILD);
    
    return lReturnConfig;
}

void Config::Print() const
{
    std::cout << BORDER << std::endl;
    std::cout << "Config: " << std::endl;
    std::cout << "Config path: " << ConfigFilePath << std::endl;
    std::cout << "Config file name: " << ConfigFileName << std::endl;
    std::cout << "Physical dof count: " << DofCount << std::endl;
    std::cout << "Harmonic waves count: " << HarmonicWaveCount << std::endl;
    std::cout << "Start frequency: " << FrequencyStart << std::endl;
    std::cout << "End frequency: " << FrequencyEnd << std::endl;
    std::cout << "Number of time integration points: " << IntPointCount << std::endl;
    std::cout << "AFT: " << AftType << std::endl;
    std::cout << "Mass matrix files: " << std::endl;
    PrintMatrixDefinitions(MassMatrices);
    std::cout << "Damping matrix files: " << std::endl;
    PrintMatrixDefinitions(DampingMatrices);
    std::cout << "Stiffness matrix files: " << std::endl;
    PrintMatrixDefinitions(StiffnessMatrices);
    std::cout << "DS builder: " << DSBuilderType << std::endl;
    std::cout << "Excitation force file: " << ExcitationForceFile << std::endl;
    std::cout << "Nonlinearities files: " << std::endl;
    PrintNonlinearityDefinitions(Nonlinearities);
    std::cout << "Save whole solutions: " << (SaveWholeSolutions ? "true" : "false") << std::endl;
    
    std::cout << "Continuation settings: " << std::endl;
    std::cout << "Max continuation steps: " << MaxStepsContinuation << std::endl;
    std::cout << "Max newton steps: " << MaxStepsNewton << std::endl;
    std::cout << "Newton norm F: " << NewtonNormF << std::endl;
    std::cout << "Init step: " << StepSizeInitial << std::endl;
    std::cout << "Min step: " << StepSizeMin << std::endl;
    std::cout << "Max step: " << StepSizeMax << std::endl;
    std::cout << "Predictor type: " << PredictorType << std::endl;
    std::cout << "Arc length scaling: " << (EnableArcLengthScaling ? "true" : "false") << std::endl;
    std::cout << BORDER << std::endl;
}

int Config::ParseOneLineInt(const std::vector<std::string>& aLines, const std::string& aPropetyName)
{
    if (aLines.size() != 1) throw "\"" + aPropetyName + "\" definition must have exactly one line!";
    return std::stoi(aLines[0]);
}
double Config::ParseOneLineDouble(const std::vector<std::string>& aLines, const std::string& aPropetyName)
{
    if (aLines.size() != 1) throw "\"" + aPropetyName + "\" definition must have exactly one line!";
    return std::stod(aLines[0]);
}
bool Config::ParseOneLineBool(const std::vector<std::string>& aLines, const std::string& aPropetyName)
{
    if (aLines.size() != 1) throw "\"" + aPropetyName + "\" definition must have exactly one line!";
    int lBoolInt = std::stoi(aLines[0]);
    if (lBoolInt == 0) return false;
    else if (lBoolInt == 1) return true;
    else throw "Invalid bool flag for \"" + aPropetyName + "\", only 0 or 1 are accepted!";
}
std::vector<MatrixDefinition> Config::ParseMatrixDefinition(const std::vector<std::string>& aLines, const std::string& aMatrixName)
{
    std::vector<MatrixDefinition> lReturnVector;
    
    int lLineInd = 0;
    
    if (aLines.size() % 2 != 0) throw aMatrixName + " matrix definition must have an even number of lines!";
    if (aLines.size() <= 0) throw aMatrixName + " matrix definition can't be empty!";
    
    while (lLineInd < aLines.size())
    {
        MatrixDefinition lNewDef;
        
        lNewDef.File = aLines[lLineInd++];
        lNewDef.Type = aLines[lLineInd++];
        
        CheckString(lNewDef.Type, C_MatrixTypes, "matrix type");
        
        lReturnVector.push_back(lNewDef);
    }
    
    return lReturnVector;
}
void Config::PrintMatrixDefinitions(const std::vector<MatrixDefinition>& aDef) const
{
    for (int i = 0; i < aDef.size(); i++)
    {
        std::cout << "\"" << aDef[i].File << "\" (" << aDef[i].Type << ")" << std::endl;
    }
}
void Config::PrintNonlinearityDefinitions(const std::vector<NonlinearityDefinition>& aDef) const
{
    for (int i = 0; i < aDef.size(); i++)
    {
        std::cout << "\"" << aDef[i].File << "\" (" << aDef[i].Type << ")" << std::endl;
    }
}
