
#include "Config.h"
#include "Misc.h"
#include "Functions.h"
#include "Functions.hpp"
#include "Nonlinearities/NonlinearitiesFactory.h"

const std::string Config::KEY_MASS    = "mass";
const std::string Config::KEY_DAMP    = "damp";
const std::string Config::KEY_STIFF   = "stiff";
const std::string Config::KEY_EXC     = "exc";
const std::string Config::KEY_HARM    = "harm";
const std::string Config::KEY_NONLIN  = "nonlin";
const std::string Config::KEY_START   = "start";
const std::string Config::KEY_END     = "end";
const std::string Config::KEY_INT     = "int";
const std::string Config::KEY_WHOLE   = "whole";
const std::string Config::KEY_STEPSC  = "stepscont";
const std::string Config::KEY_STEPSN  = "stepsnewt";
const std::string Config::KEY_NORM    = "normnewt";
const std::string Config::KEY_STEPINI = "stepinit";
const std::string Config::KEY_STEPMIN = "stepmin";
const std::string Config::KEY_STEPMAX = "stepmax";
const std::string Config::KEY_PRED    = "pred";

const std::map<std::string, std::function<std::vector<std::string>(const std::string&)>> Config::C_ConfigDefinition = 
{
    { Config::KEY_MASS,     THROW_NO_DEFAULT        },
    { Config::KEY_DAMP,     THROW_NO_DEFAULT        },
    { Config::KEY_STIFF,    THROW_NO_DEFAULT        },
    { Config::KEY_EXC,      THROW_NO_DEFAULT        },
    { Config::KEY_HARM,     THROW_NO_DEFAULT        },
    { Config::KEY_NONLIN,   RET_EMPTY_STR_VEC       },
    { Config::KEY_START,    THROW_NO_DEFAULT        },
    { Config::KEY_END,      THROW_NO_DEFAULT        },
    { Config::KEY_INT,      THROW_NO_DEFAULT        },
    { Config::KEY_WHOLE,    THROW_NO_DEFAULT        },
    { Config::KEY_STEPSC,   RET_ONE_STR("1000")     },
    { Config::KEY_STEPSN,   RET_ONE_STR("10")       },
    { Config::KEY_NORM,     RET_ONE_STR("1e-6")     },
    { Config::KEY_STEPINI,  RET_ONE_STR("1e-3")     },
    { Config::KEY_STEPMIN,  RET_ONE_STR("1e-5")     },
    { Config::KEY_STEPMAX,  RET_ONE_STR("1e-1")     },
    { Config::KEY_PRED,     RET_ONE_STR("Tangent")  },
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
    lReturnConfig.MassMatrices = ParseMatrixDefinition("Mass", lTempLines);
    lTempLines = lConfigMap[KEY_DAMP];
    lReturnConfig.DampingMatrices = ParseMatrixDefinition("Damping", lTempLines);
    lTempLines = lConfigMap[KEY_STIFF];
    lReturnConfig.StiffnessMatrices = ParseMatrixDefinition("Stiffness", lTempLines);
    
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
        CheckString(lNewDef.Type, C_NonlinearitiesFactory);
        lReturnConfig.Nonlinearities.push_back(lNewDef);
    }
    
    // number of harmonics
    lTempLines = lConfigMap[KEY_HARM];
    if (lTempLines.size() != 1) throw "Harmonic count definition must have exactly one line!";
    lReturnConfig.HarmonicWaveCount = std::stoi(lTempLines[0]);
    
    // start frequency
    lTempLines = lConfigMap[KEY_START];
    if (lTempLines.size() != 1) throw "Start frequency definition must have exactly one line!";
    lReturnConfig.FrequencyStart = std::stod(lTempLines[0]);
    
    // end frequency
    lTempLines = lConfigMap[KEY_END];
    if (lTempLines.size() != 1) throw "End frequency definition must have exactly one line!";
    lReturnConfig.FrequencyEnd = std::stod(lTempLines[0]);
    
    // integration point count
    lTempLines = lConfigMap[KEY_INT];
    if (lTempLines.size() != 1) throw "Integration point count definition must have exactly one line!";
    lReturnConfig.IntPointCount = std::stoi(lTempLines[0]);
    
    // save whole solutions
    lTempLines = lConfigMap[KEY_WHOLE];
    if (lTempLines.size() != 1) throw "Save whole solutions definition must have exactly one line!";
    int lBoolInt = std::stoi(lTempLines[0]);
    if (lBoolInt == 0) lReturnConfig.SaveWholeSolutions = false;
    else if (lBoolInt == 1) lReturnConfig.SaveWholeSolutions = true;
    else throw "Invalid bool flag for \"SaveWholeSolutions\", only 0 or 1 are accepted!";
    
    // continuation settings
    
    // max newton steps
    lTempLines = lConfigMap[KEY_STEPSN];
    if (lTempLines.size() != 1) throw "Newton step count definition must have exactly one line!";
    lReturnConfig.MaxStepsNewton = std::stoi(lTempLines[0]);
    
    // max continuation steps
    lTempLines = lConfigMap[KEY_STEPSC];
    if (lTempLines.size() != 1) throw "Continuation step count definition must have exactly one line!";
    lReturnConfig.MaxStepsContinuation = std::stoi(lTempLines[0]);
    
    // newton norm
    lTempLines = lConfigMap[KEY_NORM];
    if (lTempLines.size() != 1) throw "Newton tolerance norm definition must have exactly one line!";
    lReturnConfig.NewtonNormF = std::stod(lTempLines[0]);
    
    // init step
    lTempLines = lConfigMap[KEY_STEPINI];
    if (lTempLines.size() != 1) throw "Initial continuation step definition must have exactly one line!";
    lReturnConfig.StepSizeInitial = std::stod(lTempLines[0]);
    
    // min step
    lTempLines = lConfigMap[KEY_STEPMIN];
    if (lTempLines.size() != 1) throw "Min continuation step definition must have exactly one line!";
    lReturnConfig.StepSizeMin = std::stod(lTempLines[0]);
    
    // max step
    lTempLines = lConfigMap[KEY_STEPMAX];
    if (lTempLines.size() != 1) throw "Max continuation step definition must have exactly one line!";
    lReturnConfig.StepSizeMax = std::stod(lTempLines[0]);
    
    // predictor
    lTempLines = lConfigMap[KEY_PRED];
    if (lTempLines.size() != 1) throw "Predictor definition must have exactly one line!";
    lReturnConfig.PredictorType = lTempLines[0];
    CheckString(lReturnConfig.PredictorType, C_PredictorTypes);
        
    return lReturnConfig;
}

void Config::Print() const
{
    std::cout << BORDER << std::endl;
    std::cout << "Config: " << std::endl;
    std::cout << "Config path: " << ConfigFilePath << std::endl;
    std::cout << "Config file name: " << ConfigFileName << std::endl;
    std::cout << "Harmonic waves count: " << HarmonicWaveCount << std::endl;
    std::cout << "Start frequency: " << FrequencyStart << std::endl;
    std::cout << "End frequency: " << FrequencyEnd << std::endl;
    std::cout << "Number of time integration points: " << IntPointCount << std::endl;
    std::cout << "Mass matrix files: " << std::endl;
    PrintMatrixDefinitions(MassMatrices);
    std::cout << "Damping matrix files: " << std::endl;
    PrintMatrixDefinitions(DampingMatrices);
    std::cout << "Stiffness matrix files: " << std::endl;
    PrintMatrixDefinitions(StiffnessMatrices);
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
    std::cout << BORDER << std::endl;
}
