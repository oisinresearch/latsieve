#include <stdint.h>	// int64_t
#include <cstdlib>	// atoi
#include <iostream> // cout
#include <math.h>	// sqrt
#include <fstream>	// file
#include <ctime>	// clock_t
#include <sstream>	// stringstream
#include <stack>	// stack
#include <string>

using std::cout;
using std::endl;
using std::flush;
using std::ifstream;
using std::string;

int main (int argc, char** argv)
{
	if (argc == 1) {
		cout << endl << "Usage: ./encode3d relations" << endl << flush;
		return 0;
	}

	ifstream rels(argv[1]);
	string separator1 = ":";
	string separator2 = ",";
	string line;
	while (rels >> line) {
		string Astr = line.substr(0, line.find(separator1));
		line.erase(0, Astr.length() + 1);
		int64_t a = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		int64_t b = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		int64_t c = atoi(Astr.substr(0, Astr.find(separator2)).c_str());

        int ainfo = 0, binfo = 0, cinfo = 0;
        int abits = a; while (abits>>=1) ainfo++;    // a is always supposed to be non-negative
        int bbits = b; if (b < 0) { binfo++; bbits = -bbits; }; while (bbits>>=1) binfo++;
        int cbits = c; if (c < 0) { cinfo++; cbits = -cbits; }; while (cbits>>=1) cinfo++;

        if (ainfo <= 32 && binfo <= 32 && cinfo <= 32) {
            int64_t A = ((a+(1l<<31))<<16) + ((b+(1l<<31))>>16);
            uint64_t B = (b+(1l<<31))%(1l<<16)+((c+(1l<<31))<<16);

            cout << A << "," << B << ":" << line << endl;
        }
        else {
            cout << "#" << a << "," << b << "," << c << ":" << line << endl;
        }
	}

	return 0;
}

