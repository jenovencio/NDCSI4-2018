
#include <iostream>
#include <functional>

#include <fftw3.h>

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

#include "matplotlibcpp.h"

#include "ContinuationWrapper.h"
#include "Functions.h"
#include "ProblemInterface.h"
#include "Misc.h"
#include "Nonlinearities/CubicSpring.h"
#include "TestingFunctions.h"

void ProcessOutput(const ContinuationWrapper& aWrapper, const std::string& aOutputPath);

void TestWhatever();

int main(int argc, char **argv)
{
    ContinuationWrapper lCont;
    std::string lOutputPath = "./";
    
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
        
        if (argc >= 3) // output path
        {
            lOutputPath = std::string(argv[2]);
            std::cout << "Output path specified: \"" << lOutputPath << "\"" << std::endl;
        }
        else
        {
            std::cout << "No output path specified, using the default: \"" << lOutputPath << "\"" << std::endl;
        }
                
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
        
        ProcessOutput(lCont, lOutputPath);
        
//         std::cout << "Solution norms: " << std::endl;
//         lInterface.WriteSolutionNorms(std::cout);
    }
    catch (std::string aEx)
    {
        std::cout << "Exception: " << std::endl;
        std::cout << aEx << std::endl;
        
        ProcessOutput(lCont, lOutputPath);
        
        return 1;
    }
    catch (char const* aEx)
    {
        std::cout << "Exception: " << std::endl;
        std::cout << aEx << std::endl;
        
        ProcessOutput(lCont, lOutputPath);
        
        return 1;
    }
    
    return 0;
}

void ProcessOutput(const ContinuationWrapper& aWrapper, const std::string& aOutputPath)
{
    try
    {
        std::string lFile1Name = aOutputPath + "/" + "continuation_output_norms" + OUT_EXTENSION;
        
        const ProblemInterface* const lInterface = aWrapper.GetInterface();
        std::ofstream lOutputFile(lFile1Name);
        lInterface->WriteSolutionNorms(lOutputFile);
        lOutputFile.close();
        
        std::cout << "Solution norms written into: " << lFile1Name << std::endl;
        
        if (lInterface->HasWholeSolutions())
        {
            std::string lFile2Name = aOutputPath + "/" + "continuation_output_full" + OUT_EXTENSION;
            
            std::ofstream lOutputFile2(lFile2Name);
            lInterface->WriteWholeSolutions(lOutputFile2);
            lOutputFile2.close();
            
            std::cout << "Whole solutions written into: " << lFile2Name << std::endl;
        }
        
        lInterface->PlotSolutionNorms();
    }
    catch (std::string aEx)
    {
        std::cout << "Exception in processing outputs: " << std::endl;
        std::cout << aEx << std::endl;
    }
    catch (char const* aEx)
    {
        std::cout << "Exception in processing outputs: " << std::endl;
        std::cout << aEx << std::endl;
    }
}

void TestWhatever()
{
    try
    {
        std::string lPath = "/home/jiri/Desktop/config";
        
        Config lConfig = Config::LoadConfig(lPath);
        
        lConfig.Print();
    }
    catch (std::string aEx)
    {
        std::cout << "Exception in processing outputs: " << std::endl;
        std::cout << aEx << std::endl;
    }
    catch (char const* aEx)
    {
        std::cout << "Exception in processing outputs: " << std::endl;
        std::cout << aEx << std::endl;
    }
}

