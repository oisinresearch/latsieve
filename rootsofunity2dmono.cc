#include <stdint.h>	// int64_t
#include <cstdlib>	// atoi
#include <iostream> // cout
#include <math.h>	// sqrt
#include <fstream>	// file
#include <ctime>	// clock_t
#include <sstream>	// stringstream
#include <stack>	// stack
#include <string>
#include <vector>	// vector
#include <algorithm>	// sort

using std::cout;
using std::endl;
using std::flush;
using std::ifstream;
using std::string;
using std::sort;
using std::getline;
using std::vector;
using std::stringstream;

string commaSeparate(vector<string> v);
bool compareHexa(string x, string y);

int main (int argc, char** argv)
{
	if (argc == 1) {
		cout << endl << "Usage: ./rootsofunity2dmono relations" << endl;
		return 0;
	}

	ifstream rels(argv[1]);
	string separator1 = ":";
	string separator2 = ",";
	string line;
	while (rels >> line) {
		string side0 = line.substr(line.find(separator1) + 1);

		string Astr = line.substr(0, line.find(separator1));
		line.erase(0, Astr.length() + 1);
		int64_t a = strtoll(Astr.substr(0, Astr.find(separator2)).c_str(),NULL,10);  Astr.erase(0, Astr.find(separator2) + 1);
		int64_t b = strtoll(Astr.substr(0, Astr.find(separator1)).c_str(),NULL,10);  Astr.erase(0, Astr.find(separator1) + 1);

		vector<string> v1;
		string hexaprime;

		stringstream ss(side0);
		while (getline(ss, hexaprime, ',')) if (hexaprime != "") v1.push_back(hexaprime);

		ss.str("");
		ss.clear();

		sort(v1.begin(), v1.end(), compareHexa);

		side0 = commaSeparate(v1);

		if (a < 0) { a = -a; b = -b; }

		cout << a << "," << b << ":" << side0 << endl;
	}

	return 0;
}

string commaSeparate(vector<string> v)
{
	string output = v[0];
	for (int i = 1; i < v.size(); i++) output += "," + v[i];
	return output;
}


bool compareHexa(string x, string y)
{
	stringstream ss;
	int X, Y;
	ss << std::hex << x;
	ss >> X;
	ss.clear();
	ss << std::hex << y;
	ss >> Y;
	return (X < Y);
}

