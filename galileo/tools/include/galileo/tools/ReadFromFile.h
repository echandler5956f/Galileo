
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <sstream>

namespace galileo
{
    namespace tools
    {

        std::map<std::string, std::string> readFromFile(std::string file_name);

        std::vector<std::string> readAsVector(std::string vector_as_string);

    }
}