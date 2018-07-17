
struct ContinuationSettings;

#pragma once

struct ContinuationSettings
{
    int MaxStepsContinuation;
    int MaxStepsNewton;
    double StepSizeInitial;
    double StepSizeMin;
    double StepSizeMax;
    
    std::string PredictorType; // secant, constant, tangent
};
