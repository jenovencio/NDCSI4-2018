
class Config;

#pragma once
#include <string>
#include <vector>
#include <map>
#include <functional>

#include "MatrixDefinition.h"
#include "Nonlinearities/NonlinearityDefinition.h"
#include "ProblemParams.h"

// the return after the throw is there just so the compiler knows what the function returns
// it should never get to that point because of the throw
#define THROW_NO_DEFAULT { [](const std::string& aStr) { throw std::string("\"") + aStr + "\" has no default value!"; return std::vector<std::string>(); } }
// this macro only works if you pass one string literal into it
#define RET_ONE_STR(STR_VAL) [](const std::string& aStr) { return std::vector<std::string>( { STR_VAL } ); }
#define RET_EMPTY_STR_VEC [](const std::string& aStr) { return std::vector<std::string>(); }

class Config
{
// static constants
public:
    // config keywords definition
    static const std::string KEY_MASS;      // mass matrix def
    static const std::string KEY_DAMP;      // damping matrix def
    static const std::string KEY_STIFF;     // stiffness matrix def
    static const std::string KEY_EXC;       // excitation force def
    static const std::string KEY_DOF;      // number of harmonice waves
    static const std::string KEY_HARM;      // number of harmonice waves
    static const std::string KEY_NONLIN;    // nonlinearities def
    static const std::string KEY_START;     // start continuation frequency
    static const std::string KEY_END;       // end continuation frequency
    static const std::string KEY_INT;       // number of integration points
    static const std::string KEY_WHOLE;     // save whole solutions (bool)
    static const std::string KEY_STEPSC;    // number of continuation steps
    static const std::string KEY_STEPSN;    // number of newton steps
    static const std::string KEY_NORM;      // norm criteria for newton convergence
    static const std::string KEY_STEPINI;   // initial continuation step
    static const std::string KEY_STEPMIN;   // minimum continuation step
    static const std::string KEY_STEPMAX;   // maximum continuation step
    static const std::string KEY_PRED;      // continuation predictor
    static const std::string KEY_SCALEARC;  // scaleing for arc length continuation
    static const std::string KEY_AFT;       // aft type
    
    // list of parameter keys and functions that return default values for those  parameters
    // the argument passed into the function should be the parameter name 
    // so it gets passed into the exception if there is no default so the user knows
    // which parameter they need to provide
    static const std::map<std::string, std::function<std::vector<std::string>(const std::string&)>> C_ConfigDefinition;

public:
    // Load function
    static Config LoadConfig(const std::string& aConfigFilePath);
    // print
    void Print() const;
    
public:
    
    std::string ConfigFilePath;
    std::string ConfigFileName;
    // relative paths with respect to the config path
    std::vector<MatrixDefinition> MassMatrices;
    std::vector<MatrixDefinition> DampingMatrices;
    std::vector<MatrixDefinition> StiffnessMatrices;
    std::string ExcitationForceFile;
    std::string ContinuationSettingsFile;
    std::vector<NonlinearityDefinition> Nonlinearities;
    
    // number of physical degrees of freedom (i.e. in time domain)
    int DofCount;
    // number of harmonic waves (first harmonic is the DC component)
    // that means that the number of harmonic dofs will be DOFCount * (2 * HarmonicCount - 1)
    int HarmonicWaveCount;
    // start continuation frequency
    double FrequencyStart;
    // end continuation frequency
    double FrequencyEnd;
    // integration point count
    int IntPointCount;
    // should the interface store the whole solutions? (all harmonic coefficients)
    bool SaveWholeSolutions;
    std::string AftType;
    
    // continuation settings
    int MaxStepsContinuation;
    int MaxStepsNewton;
    double NewtonNormF;
    double StepSizeInitial;
    double StepSizeMin;
    double StepSizeMax;
    std::string PredictorType; // secant, constant, tangent
    bool EnableArcLengthScaling;
    
private:
    static int ParseOneLineInt(const std::vector<std::string>& aLines, const std::string& aPropetyName = "");
    static double ParseOneLineDouble(const std::vector<std::string>& aLines, const std::string& aPropetyName = "");
    static bool ParseOneLineBool(const std::vector<std::string>& aLines, const std::string& aPropetyName = "");
    static std::vector<MatrixDefinition> ParseMatrixDefinition(const std::vector<std::string>& aLines, const std::string& aMatrixName);
    
    void PrintMatrixDefinitions(const std::vector<MatrixDefinition>& aDef) const;
    void PrintNonlinearityDefinitions(const std::vector<NonlinearityDefinition>& aDef) const;
};

