#include <stdint.h>	// int64_t
#include <iostream> // cout
#include <math.h>	// sqrt
#include <fstream>	// file
#include <ctime>	// clock_t
#include <sstream>	// stringstream
#include <stack>	// stack
#include <string>
#include <assert.h>

using std::cout;
using std::endl;
using std::flush;
using std::ifstream;
using std::string;
using std::to_string;
using std::stringstream;
auto npos = std::string::npos;

int main (int argc, char** argv)
{
	if (argc != 2 ) {	// argc must be 2
		cout << endl << "Usage: ./galoisexpand relations" << endl;
		return 0;
	}

	// read relations file
	ifstream rels(argv[1]);
	string separator1 = ":";
	string separator2 = ",";
	string line;
	while (rels >> line) {
		string line0 = line;
		bool isrel = true;

		string Astr = line.substr(0, line.find(separator1));
		line.erase(0, Astr.length() + 1);
		int64_t a = strtoll(Astr.substr(0, Astr.find(separator2)).c_str(), NULL, 10);
		Astr.erase(0, Astr.find(separator2) + 1);
		int64_t b = strtoll(Astr.substr(0, Astr.find(separator2)).c_str(), NULL, 10);
		Astr.erase(0, Astr.find(separator2) + 1);
		
		// print original relation
		cout << to_string(a) << "," << to_string(b) << ":" << line << endl;
		// print Galois conjugate
		cout << to_string(a) << "," << to_string(-b) << ":" << line << endl;
	}

	return 0;
}

