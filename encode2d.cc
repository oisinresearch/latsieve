#include <stdint.h>	// int64_t
#include <cstdlib>	// atoi
#include <iomanip>	// hex
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
using std::hex;

int main (int argc, char** argv)
{
	if (argc == 1) {
		cout << endl << "Usage: ./encode2d relations" << endl;
		return 0;
	}

	ifstream rels(argv[1]);
	string separator1 = ":";
	string separator2 = ",";
	string line;
	while (rels >> line) {
		string Astr = line.substr(0, line.find(separator1));
		line.erase(0, Astr.length() + 1);
		int64_t a = strtoll(Astr.substr(0, Astr.find(separator2)).c_str(),NULL,10);  Astr.erase(0, Astr.find(separator2) + 1);
		int64_t b = strtoll(Astr.substr(0, Astr.find(separator2)).c_str(),NULL,10);  Astr.erase(0, Astr.find(separator2) + 1);

        int ainfo = 0, binfo = 0, cinfo = 0, dinfo = 0;
        int64_t abits = a; if (a < 0) { ainfo++; abits = -abits; }; while (abits>>=1) ainfo++;
        int64_t bbits = b; if (b < 0) { binfo++; bbits = -bbits; }; while (bbits>>=1) binfo++;

        if (ainfo <= 48 && binfo <= 48) {
            // 3d
			//int64_t A = ((a+(1l<<31))<<16) + ((b+(1l<<31))>>16);
            //uint64_t B = (b+(1l<<31))%(1l<<16)+((c+(1l<<31))<<16);

			// 4d
			uint64_t A = (a + (1l<<46));
			uint64_t B = (b + (1l<<46));

            cout << hex << A << "," << hex << B << ":" << line << endl;
        }
        else {
            cout << "#" << a << "," << b << ":" << line << endl;
        }
	}

	return 0;
}

