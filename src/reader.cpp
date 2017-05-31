#include "reader.h"
#include <string>

using namespace std;


int
Reader::readGridFile(const string &path)
{
    ifstream ifs(path);
    string tmp_str;

    if(ifs.fail())
    {
        cerr << "Error: in readGridFile()" << endl;
        cerr << "Failed to read files: " << path << endl;
    }
    while(getline(ifs, str))
    {
        cout << "[" << str << "]" << endl;
    }

    return 0;
}