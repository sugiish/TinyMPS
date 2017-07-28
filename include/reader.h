#ifndef MPS_READER_H_INCLUDED
#define MPS_READER_H_INCLUDED

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <regex>

class Reader
{
public:
	Reader(std::string path)
	{
		std::ifstream ifs(path);
		
		if(ifs.fail())
		{
			std::cerr << "Error: in Reader()" << std::endl;
			std::cerr << "Failed to read files: " << path << std::endl;
			return;
		}

		std::string tmp_str;
		std::regex re("\\(.*\\)"); // For removing (**)
		std::regex re2("-+\\w+-+");// For removing like --**--
		while(getline(ifs, tmp_str))
		{
			std::stringstream ss;
			ss.str(tmp_str);

			std::string item, value;
			ss >> item;
			{
				// Lines that begin with '#' are comments
				char first = item.at(0);
				if(first == '#') continue;
			}

			std::cout << item << std::endl;
			item = std::regex_replace(item, re, "");
			item = std::regex_replace(item, re2, "");
			std::cout << item << std::endl;
			
			ss >> value;
			data[item] = value;
		}
	}

	virtual ~Reader(){}

	inline int getValue(const std::string& item, int& value)
	{
		if(data.find(item) == data.end())
		{
			return 1;
		}
		stringstream ss;
		ss << data[item];
		ss >> value;
		return 0;
	}

	inline int getValue(const std::string& item, double& value)
	{
		if(data.find(item) == data.end())
		{
			return 1;
		}
		stringstream ss;
		ss << data[item];
		ss >> value;
		return 0;
	}

	int getValue(const std::string& item, bool& value)
	{
		if(data.find(item) == data.end())
		{
			return 1;
		}
		std::string tmp = data[item];
		std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
		if(tmp == "on" || tmp == "true")
		{
			value = true;
		}
		else
		{
			value = false;
		}
		return 0;
	}
	
	inline int getValue(const std::string& item, std::string& value)
	{
		if(data.find(item) == data.end())
		{
			return 1;
		}
		value = data[item];
		return 0;
	}

	inline void printValues()
	{
		for(std::unordered_map<std::string, std::string>::iterator itr = data.begin(); itr != data.end(); ++itr) {
        	std::cout << "key = " << itr->first << ", val = " << itr->second << "\n";
    	}
	}

private:
	std::unordered_map<std::string, std::string> data;
	
};

#endif //MPS_READER_H_INCLUDED