
#include <iostream>

#include "Teuchos_ParameterList.hpp"

#include "NOX_Utils.H"
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_NormF.H"
#include "NOX_StatusTest_MaxIters.H"
#include "NOX_StatusTest_Combo.H"

#include "LOCA_LAPACK_Factory.H"
#include "LOCA_LAPACK_Group.H"
#include "LOCA_GlobalData.H"
#include "LOCA_Stepper.H"

#include "Functions.h"
#include "ProblemInterface.h"
#include "Misc.h"

int main(int argc, char **argv)
{
    try
    {
        std::string lConfigPath = "config";
        if (argc >= 2) // first argument is always the name of the exe
        {
            lConfigPath = std::string(argv[1]);
            std::cout << "Config path specified: \"" << lConfigPath << "\"" << std::endl;
        }
        else
        {
            std::cout << "No config path specified, using the default: \"" << lConfigPath << "\"" << std::endl;
        }
        
        Config lConfig = LoadConfig(lConfigPath);
        PrintConfig(lConfig);
        
        // continuation parameters
        int lMaxNewtonIters = 20;
        
        // Create parameter list
        Teuchos::RCP<Teuchos::ParameterList> paramList =
        Teuchos::rcp(new Teuchos::ParameterList);

        // Create LOCA sublist
        Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

        // Create the stepper sublist and set the stepper parameters
        Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
//         stepperList.set("Continuation Method", "Arc Length");// Default
        stepperList.set("Continuation Method", "Natural");
        stepperList.set("Continuation Parameter", ProblemInterface::cFrequencyName);  // Must set
        stepperList.set("Initial Value", lConfig.FrequencyMin); // Must set
        stepperList.set("Min Value", lConfig.FrequencyMin);             // Must set
        stepperList.set("Max Value", lConfig.FrequencyMax);             // Must set
        stepperList.set("Max Steps", 400);                    // Should set
        stepperList.set("Max Nonlinear Iterations", lMaxNewtonIters); // Should set
        stepperList.set("Compute Eigenvalues", false);        // Default

        // Create predictor sublist
        Teuchos::ParameterList& predictorList =
        locaParamsList.sublist("Predictor");
//         predictorList.set("Method", "Secant");               // Default
        predictorList.set("Method", "Constant");
//         predictorList.set("Method", "Tangent");

        // Create step size sublist
        Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
        stepSizeList.set("Method", "Adaptive");             // Default
        stepSizeList.set("Initial Step Size", 1.0e-3);   // Should set
        stepSizeList.set("Min Step Size", 1.0e-4);    // Should set
        stepSizeList.set("Max Step Size", 4.0e-2);      // Should set

        // Create the "Solver" parameters sublist to be used with NOX Solvers
        Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
        Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
        nlPrintParams.set("Output Information",
                NOX::Utils::Details +
                NOX::Utils::OuterIteration +
                NOX::Utils::InnerIteration +
                NOX::Utils::Warning +
                NOX::Utils::StepperIteration +
                NOX::Utils::StepperDetails +
                NOX::Utils::StepperParameters);  // Should set

        // Create LAPACK Factory
        Teuchos::RCP<LOCA::LAPACK::Factory> lapackFactory =
        Teuchos::rcp(new LOCA::LAPACK::Factory);

        // Create global data object
        Teuchos::RCP<LOCA::GlobalData> globalData =
        LOCA::createGlobalData(paramList, lapackFactory);

        // Set up the problem interface
        
        ProblemInterface lInterface(lConfig);
        LOCA::ParameterVector p;
        p.addParameter(ProblemInterface::cFrequencyName, lConfig.FrequencyMin);

        // Create a group which uses that problem interface. The group will
        // be initialized to contain the default initial guess for the
        // specified problem.
        Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> grp =
        Teuchos::rcp(new LOCA::LAPACK::Group(globalData, lInterface));

        grp->setParams(p);

        // Set up the status tests
        Teuchos::RCP<NOX::StatusTest::NormF> normF =
        Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
        Teuchos::RCP<NOX::StatusTest::MaxIters> maxIters =
        Teuchos::rcp(new NOX::StatusTest::MaxIters(lMaxNewtonIters));
        Teuchos::RCP<NOX::StatusTest::Generic> comboOR =
        Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
                            normF,
                            maxIters));

        std::cout << "Creating the stepper ... " << std::endl;
        // Create the stepper
        LOCA::Stepper stepper(globalData, grp, comboOR, paramList);
        std::cout << "Creating the stepper ... done" << std::endl;
        
        std::cout << "Running the stepper ... " << std::endl;
        // Perform continuation run
        LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();
        
        std::cout << "Running the stepper ... done" << std::endl;
        
        // Check for convergence
        if (status == LOCA::Abstract::Iterator::Finished)
        std::cout << "All examples passed" << std::endl;
        else {
        if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
        globalData->locaUtils->out()
        << "Stepper failed to converge!" << std::endl;
        }

        // Get the final solution from the stepper
        Teuchos::RCP<const LOCA::LAPACK::Group> finalGroup =
        Teuchos::rcp_dynamic_cast<const LOCA::LAPACK::Group>(stepper.getSolutionGroup());
        const NOX::LAPACK::Vector& finalSolution =
        dynamic_cast<const NOX::LAPACK::Vector&>(finalGroup->getX());

        LOCA::destroyGlobalData(globalData);
        
        std::ofstream lOutputFile("continuation_output" + OUT_EXTENSION);
        lInterface.WriteSolutions(lOutputFile);
        lOutputFile.close();
    }
    catch (std::string aEx)
    {
        std::cout << "Exception: " << std::endl;
        std::cout << aEx << std::endl;
        return 1;
    }
    catch (char const* aEx)
    {
        std::cout << "Exception: " << std::endl;
        std::cout << aEx << std::endl;
        return 1;
    }
    
    return 0;
}
