#include <iostream>
#include <fstream>
#include <string>
#include "Web.h"

int main(int argc, char* argv[])
{
    if (argc != 2) 
    {
        std::cerr << "Usage: " << argv[0] << " <input_filepath>" << std::endl;
        return 1;
    }

    std::ifstream inputFile(argv[1]); 

    if (!inputFile.is_open()) 
    {
        std::cerr << "Error: Could not open file '" << argv[1] << "'" << std::endl;
        return 1;
    }

    Web web;

    try
    {
        web.ReadGraph(inputFile); 
        std::cout << "Maximal capacity: " << web.FindMaximalCapacity() << std::endl;
    }
    catch (const std::exception& err)
    {
        std::cerr << "An error occurred: " << err.what() << std::endl;
        if (inputFile.is_open()) {
            inputFile.close();
        }
        return 1;
    }

    if (inputFile.is_open()) 
    {
        inputFile.close();
    }

    return 0;
}