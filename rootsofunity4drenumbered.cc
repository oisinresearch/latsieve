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
using std::hex;

string commaSeparate(vector<string> v);
bool compareHexa(string x, string y);

int main (int argc, char** argv)
{
	if (argc != 2) {
		cout << endl << "Usage: ./rootsofunity4drenumbered relations" << endl;
		return 0;
	}

	ifstream rels(argv[1]);
	string separator1 = ":";
	string separator2 = ",";
	string line;
	while (rels >> line) {
		string pid = line.substr(0, line.find(separator1));
		line.erase(0, pid.length() + 1);
		string Ahex = pid.substr(0, pid.find(separator2));
		string Bhex = pid.substr(pid.find(separator2)+1);
		int64_t A = stol(Ahex, NULL, 16);
		int64_t B = stol(Bhex, NULL, 16);
		
		int64_t a = (A>>24) - (1l<<23);
		int64_t b = (A&((1l<<24)-1)) - (1l<<23);
		int64_t c = (B>>24) - (1l<<23);
		int64_t d = (B&((1l<<24)-1)) - (1l<<23);

		// make a non-negative
		if (a == 0) {
			if (b < 0) {
				b = -b; c = -c; d = -d;
			}
		}
		if (a < 0) {
			a = -a; b = -b; c = -c; d = -d;
		}
		
		A = ((a+(1<<23))<<24)+(b+(1<<23));
		B = ((c+(1<<23))<<24)+(d+(1<<23));

		stringstream ss;
		ss.str("");
		ss << hex << A;
		string Astr = ss.str();
		ss.str("");
		ss << hex << B;
		string Bstr = ss.str();

		cout << Astr << "," << Bstr << ":" << line << endl;
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

