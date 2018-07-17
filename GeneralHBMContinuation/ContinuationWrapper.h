
class ContinuationWrapper;

#pragma once
#include "LOCA_Stepper.H"

#include "ProblemInterface.h"
#include "Nonlinearities/NonlinearBase.h"

class ContinuationWrapper
{
private:
    LOCA::Stepper* mStepper = nullptr;
    // bound to the current stepper
    Teuchos::RCP<LOCA::GlobalData> mGlobalData;
    ProblemInterface* mInterface = nullptr;
    bool mCurrentStepperRan = false;
    
public:
    ~ContinuationWrapper();
    void Init(const std::string& aConfigPath, const std::vector<NonlinearBase*> aNonlinearities);
    LOCA::Abstract::Iterator::IteratorStatus RunContinuation();
    const ProblemInterface* const GetInterface();
    
private:
};
