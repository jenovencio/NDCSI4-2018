
#include <iostream>
#include <functional>

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

#include "ContinuationWrapper.h"
#include "Functions.h"
#include "ProblemInterface.h"
#include "Misc.h"
#include "Nonlinearities/CubicSpring.h"

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
        
        std::string lOutputPath = "./";
        if (argc >= 3) // output path
        {
            lOutputPath = std::string(argv[2]);
            std::cout << "Output path specified: \"" << lOutputPath << "\"" << std::endl;
        }
        else
        {
            std::cout << "No output path specified, using the default: \"" << lOutputPath << "\"" << std::endl;
        }
        
        std::vector<NonlinearBase*> lNonlinearities;
        CubicSpring lCS;
        lCS.AddCubicSpring(0, 1, 1);
        lNonlinearities.push_back(&lCS);
        
        ContinuationWrapper lCont;
        
        lCont.Init(lConfigPath);
        
        auto lStatus = lCont.RunContinuation();
        
        if (lStatus == LOCA::Abstract::Iterator::Finished)
        {
            std::cout << "Everything worked well!" << std::endl;
        }
        else
        {
            std::cout << "Something CRASHED!" << std::endl;
        }
        
        const ProblemInterface* const lInterface = lCont.GetInterface();
        
        std::ofstream lOutputFile(lOutputPath + "continuation_output_norms" + OUT_EXTENSION);
        lInterface->WriteSolutionNorms(lOutputFile);
        lOutputFile.close();
        
        if (lInterface->HasWholeSolutions())
        {
            std::ofstream lOutputFile2(lOutputPath + "continuation_output_raw" + OUT_EXTENSION);
            lInterface->WriteWholeSolutions(lOutputFile2);
            lOutputFile2.close();
        }
        
//         std::cout << "Solution norms: " << std::endl;
//         lInterface.WriteSolutionNorms(std::cout);
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
