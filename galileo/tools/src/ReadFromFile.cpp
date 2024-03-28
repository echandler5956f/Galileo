
#include "galileo/tools/ReadFromFile.h"

namespace galileo
{
    namespace tools
    {

        std::map<std::string, std::string> readFromFile(std::string file_name)
        {
            std::ifstream file(file_name);

            if (!file.is_open())
            {
                std::cerr << "Could not open file " << file_name << std::endl;
                return {};
            }

            std::map<std::string, std::string> dataMap;
            std::string line;
            while (std::getline(file, line))
            {
                std::istringstream iss(line);
                std::string key, value;
                if (std::getline(iss, key, ',') && std::getline(iss, value))
                {
                    dataMap[key] = value;
                }
            }

            return dataMap;
        }

        std::vector<std::string> readAsVector(std::string vector_as_string)
        {
            std::vector<std::string> vector_out;

            size_t firstParenthesis = vector_as_string.find_first_of("(");
            size_t lastParenthesis = vector_as_string.find_last_of(")");

            if (firstParenthesis == std::string::npos || lastParenthesis == std::string::npos)
            {
                std::cerr << "Could not find parenthesis in string; must be formatted as \"( value_1, value_2, value_3 ... value_l)\" " << std::endl;
                return {};
            }

            std::string vectorString = vector_as_string.substr(firstParenthesis + 1, lastParenthesis - firstParenthesis - 1);

            size_t comma = vectorString.find(",");
            while (comma != std::string::npos)
            {
                vector_out.push_back(vectorString.substr(0, comma));
                vectorString = vectorString.substr(comma + 1);
                comma = vectorString.find(",");
            }
            vector_out.push_back(vectorString);

            return vector_out;
        }
    }
}