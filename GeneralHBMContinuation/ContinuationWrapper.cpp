
#include "NOX_Utils.H"
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_NormF.H"
#include "NOX_StatusTest_MaxIters.H"
#include "NOX_StatusTest_Combo.H"

#include "LOCA_LAPACK_Factory.H"
#include "LOCA_GlobalData.H"
#include "LOCA_LAPACK_Group.H"

#include "ContinuationWrapper.h"
#include "Config.h"
#include "Functions.h"
#include "ProblemInterface.h"

ContinuationWrapper::~ContinuationWrapper()
{
    if (mStepper != nullptr)
    {
        delete mStepper;
        mStepper = nullptr;
        
        LOCA::destroyGlobalData(mGlobalData);
    }
    if (mInterface != nullptr)
    {
        delete mInterface;
        mInterface = nullptr;
    }
}
void ContinuationWrapper::Init(const std::string& aConfigPath, const std::vector<NonlinearBase*> aNonlinearities)
{
    Config lConfig = Config::LoadConfig(aConfigPath);
    lConfig.Print();
        
    // Create parameter list
    Teuchos::RCP<Teuchos::ParameterList> paramList =
        Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    stepperList.set("Continuation Method", "Arc Length");// Default
//         stepperList.set("Continuation Method", "Natural");
    stepperList.set("Continuation Parameter", ProblemInterface::cContParameterName);  // Must set
    stepperList.set("Initial Value", 0.0); // Must set
    stepperList.set("Min Value", 0.0);             // Must set
    stepperList.set("Max Value", 1.0);             // Must set
    stepperList.set("Max Steps", lConfig.MaxStepsContinuation);                    // Should set
    stepperList.set("Max Nonlinear Iterations", lConfig.MaxStepsNewton); // Should set
    stepperList.set("Compute Eigenvalues", false);        // Default
    stepperList.set("Enable Tangent Factor Step Size Scaling", true);        // Default
    stepperList.set("Enable Arc Length Scaling", lConfig.EnableArcLengthScaling);        // Default
    stepperList.set("Skip Parameter Derivative", false);        
//     stepperList.set("Min Tangent Factor", 0.1);        // Default

    // Create predictor sublist
    Teuchos::ParameterList& predictorList =
    locaParamsList.sublist("Predictor");
//         predictorList.set("Method", "Secant");               // Default
//         predictorList.set("Method", "Constant");
    predictorList.set("Method", lConfig.PredictorType);

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Method", "Adaptive");             // Default
    stepSizeList.set("Initial Step Size", lConfig.StepSizeInitial);   // Should set
    stepSizeList.set("Min Step Size", lConfig.StepSizeMin);    // Should set
    stepSizeList.set("Max Step Size", lConfig.StepSizeMax);      // Should set

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("Output Information",
//             NOX::Utils::Details +
//             NOX::Utils::OuterIteration +
//             NOX::Utils::InnerIteration +
            NOX::Utils::Warning +
            NOX::Utils::StepperIteration +
//             NOX::Utils::StepperDetails +
            NOX::Utils::StepperParameters);  // Should set
    
	Teuchos::ParameterList& lSearchParams = nlParams.sublist("Line Search");
	lSearchParams.set("Method", "Polynomial");


    // Create LAPACK Factory
    Teuchos::RCP<LOCA::LAPACK::Factory> lapackFactory =
    Teuchos::rcp(new LOCA::LAPACK::Factory);

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
    LOCA::createGlobalData(paramList, lapackFactory);

    // Set up the problem interface
    
    ProblemInterface* lInterface = new ProblemInterface(lConfig, aNonlinearities, lConfig.SaveWholeSolutions);
    
    LOCA::ParameterVector p;
    p.addParameter(ProblemInterface::cContParameterName, 0.0);

    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> grp =
    Teuchos::rcp(new LOCA::LAPACK::Group(globalData, *lInterface));

    grp->setParams(p);

    // Set up the status tests
    Teuchos::RCP<NOX::StatusTest::NormF> normF =
    Teuchos::rcp(new NOX::StatusTest::NormF(lConfig.NewtonNormF));
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxIters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(lConfig.MaxStepsNewton));
    Teuchos::RCP<NOX::StatusTest::Generic> comboOR =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
                        normF,
                        maxIters));

    std::cout << "Creating the stepper ... " << std::endl;
    // Create the stepper
    LOCA::Stepper* lStepper = new LOCA::Stepper(globalData, grp, comboOR, paramList);
    
    if (mStepper != nullptr)
    {
        delete mStepper;
        mStepper = nullptr;
        
        LOCA::destroyGlobalData(mGlobalData);
    }
    if (mInterface != nullptr)
    {
        delete mInterface;
        mInterface = nullptr;
    }
    
    mStepper = lStepper;
    mGlobalData = globalData;
    mInterface = lInterface;
    mCurrentStepperRan = false;
    
    std::cout << "Creating the stepper ... done" << std::endl;
}

LOCA::Abstract::Iterator::IteratorStatus ContinuationWrapper::RunContinuation()
{
    if (mStepper == nullptr) throw "Stepper not initialised! Call Init() first.";
    if (mCurrentStepperRan) throw "Stepper can not be run twice without re-initialisation. Call Init() before running again!";
    
    mInterface->ClearSolutions();
    
    LOCA::Abstract::Iterator::IteratorStatus lReturnStatus = mStepper->run();
    
    mCurrentStepperRan = true;
    return lReturnStatus;
}

const ProblemInterface* const ContinuationWrapper::GetInterface() const
{
    return mInterface;
}
