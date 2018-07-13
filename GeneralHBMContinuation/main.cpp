#include <iostream>

#include "Functions.h"
#include "ProblemInterface.h"

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
        
        ProblemInterface lInterface(lConfig);
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
