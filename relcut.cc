#include <stdint.h>	// int64_t
#include <iostream> // cout
#include <math.h>	// sqrt
#include <fstream>	// file
#include <ctime>	// clock_t
#include <sstream>	// stringstream
#include <stack>	// stack
#include <string>
#include <assert.h>
#include <vector> // vector
#include <algorithm>  // sort

using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::to_string;
using std::stringstream;
using std::getline;
using std::hex;
using std::vector;
using std::unique;
using std::pair;
using std::sort;
using std::flush;
auto npos = std::string::npos;

string int2hex(int n);

int main (int argc, char** argv)
{
	if (argc < 8 ) {	// argc must be 6
		cout << endl << "Usage: ./relcut num_rels num_ideals keep relations_file "
			"output_file excess_file dlog" << endl;
		return 0;
	}

	bool verbose = false;

	string line;
	string separator0 = " ";
	string separator1 = ":";
	string separator2 = ",";

	int N = atoi(argv[1]);
	int T = atoi(argv[2]);
	int keep = atoi(argv[3]);

	// first pass, read relations and count unique ideals
	int* allideals = new int[T]();	// initialize to 0

	// read relations file
	int  num = 0; int mark = 1024; int maxideal = -1;
	ifstream rels(argv[4]);
	FILE* out = fopen(argv[5], "wt");
	FILE* excess = fopen(argv[6], "wt");
	bool dlog = false;
	if (argc == 9) dlog = true;
	cout << "Reading relations..." << endl;
	int  k = 0;
	bool done = false;
	while (std::getline(rels, line)) { //rels >> line) {
		string line0 = line;
		if (!done) {
			string Astr = line.substr(0, line.find(separator1));
			line.erase(0, Astr.length() + 1);

			stringstream ss(line);
			stringstream hh;
			vector<int> ideals0;
			vector<int> ideals;
			ideals.clear();

			while(ss.good()) {
				string str;
				getline(ss, str, ',');
				int ideal = stoi(str, NULL, 16);
				ideals0.push_back(ideal);
				ideals.push_back(ideal);
				if (ideal > maxideal) maxideal = ideal;
			}

			sort(ideals.begin(), ideals.end());
			ideals.erase(unique(ideals.begin(), ideals.end()), ideals.end());

			int s0 = ideals0.size();
			int s = ideals.size();
			for (int i = 0; i < s; i++) {
				if (allideals[ideals[i]] == 0) k++;
				allideals[ideals[i]] = 1;
			}
			string idealstr = to_string(ideals[0]);
			string idealrun = "," + int2hex(ideals[1]);
			int parity = 1;
			for (int i = 2; i < s0; i++) {
				if (ideals0[i] != ideals0[i-1]) {
					if (parity == 1) {
						idealstr = idealstr + idealrun;
					}
					idealrun = "," + int2hex(ideals0[i]);
					parity = 1;
				}
				else {
					idealrun = "," + int2hex(ideals0[i]);
					parity = 1 - parity;
				}
				if (i == s0 - 1 && parity == 1) {
					idealstr = idealstr + idealrun;
				}
			}
			string line1 = Astr + ":" + idealstr;
			fprintf(out, "%s\n", line1.c_str());
		}
		
		if (done)
			fprintf(excess, "%s\n", line0.c_str());

		if (num % mark == 0) {
			cout << num << " relations processed..." << endl;
			mark *= 2;
		}

		num++;
		if (num == k + keep) {
			done = true;
			cout << "fully determined for " << num << " relations and " << k
				<< " ideals!" << endl;
		}
	}
	cout << num << " relations processed." << endl;
	fclose(out);
	fclose(excess);

	delete[] allideals;

	return 0;
}

string int2hex(int n)
{
	string str1;
	stringstream stream;
	stream.str("");
	stream << hex << n;
	stream >> str1;
	return str1;
}

