#include <stdint.h>	// int64_t
#include <iostream> // cout
#include <gmpxx.h>
#include "intpoly.h"
#include <math.h>	// sqrt
#include <fstream>	// file
#include <ctime>	// clock_t
#include "mpz_poly.h"
#include <sstream>	// stringstream
#include <stack>	// stack
#include <string>
#include "factorsmall.h"
#include <assert.h>
#include <vector>

using std::cout;
using std::endl;
using std::flush;
using std::ifstream;
using std::string;
using std::to_string;
using std::stringstream;
using std::vector;
using std::sort;
using std::hex;
auto npos = std::string::npos;

int main (int argc, char** argv)
{
	if (argc == 1) {
		cout << "Sort relations from smallest ideal to largest.  Outputs to stdout." << endl;
		cout << endl << "Usage: ./sortideals relations" << endl;
		return 0;
	}

	// read relations file
	mpz_t p_mpz; mpz_init(p_mpz); int BASE = 16;
	ifstream rels(argv[1]);
	string separator1 = ":";
	string separator2 = ",";
	int n = 0;
    string line;
	while (rels >> line) {
		n++;
		string ABstr = line.substr(0, line.find(separator1));
		line.erase(0, ABstr.length() + 1);

		vector<int> p;

		string idealstr = line;
		while (idealstr.length()) {
			string pstr = idealstr.substr(0, idealstr.find(separator2));
			if (pstr.length()) {
				mpz_set_str(p_mpz, pstr.c_str(), BASE);
				int pt = mpz_get_ui(p_mpz);
				p.push_back(pt);
				int pos = idealstr.find(separator2);
				if (pos == npos) pos = pstr.length() - 1;
				idealstr.erase(0, pos + 1);
			}
		}
		sort(p.begin(), p.end());
		
		// store relation
		cout << ABstr << hex << ":";
        int L = p.size();
        for (int i = 0; i < L-1; i++) cout << p[i] << ",";
        cout << p[L-1] << endl;
	}

	mpz_clear(p_mpz);

	return 0;
}


