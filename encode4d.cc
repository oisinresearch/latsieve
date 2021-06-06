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
		cout << endl << "Usage: ./encode4d relations" << endl;
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
		int64_t c = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		int64_t d = atoi(Astr.substr(0, Astr.find(separator2)).c_str());

        int ainfo = 0, binfo = 0, cinfo = 0, dinfo = 0;
        int abits = a; if (a < 0) { ainfo++; abits = -abits; }; while (abits>>=1) ainfo++;
        int bbits = b; if (b < 0) { binfo++; bbits = -bbits; }; while (bbits>>=1) binfo++;
        int cbits = c; if (c < 0) { cinfo++; cbits = -cbits; }; while (cbits>>=1) cinfo++;
        int dbits = d; if (d < 0) { dinfo++; dbits = -dbits; }; while (dbits>>=1) dinfo++;

        if (ainfo <= 24 && binfo <= 24 && cinfo <= 24 && dinfo <= 24) {
            // 3d
			//int64_t A = ((a+(1l<<31))<<16) + ((b+(1l<<31))>>16);
            //uint64_t B = (b+(1l<<31))%(1l<<16)+((c+(1l<<31))<<16);

			// 4d
			uint64_t A = ((a + (1l<<23))<<24) + (b + (1l<<23));
			uint64_t B = ((c + (1l<<23))<<24) + (d + (1l<<23));

            cout << hex << A << "," << hex << B << ":" << line << endl;
        }
        else {
            cout << "#" << a << "," << b << "," << c << "," << d << ":" << line << endl;
        }
	}

	return 0;
}

