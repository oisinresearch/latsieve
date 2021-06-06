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

using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::to_string;
using std::stringstream;
using std::getline;
using std::hex;
using std::vector;
auto npos = std::string::npos;

int main (int argc, char** argv)
{
	if (argc != 5 ) {	// argc must be 5
		cout << endl << "Usage: ./tnfsfilter2 num_rels taget_num num_sch_maps relations_file" << endl;
		return 0;
	}

	bool verbose = false;

	string line;
	string separator0 = " ";
	string separator1 = ":";
	string separator2 = ",";

	int N = atoi(argv[1]);
	int T = atoi(argv[2]);
	int numSch = atoi(argv[3]);

	// read relations file
	int nlast = -1; int  i = 0; int s = 0;
	ifstream rels(argv[4]);
	while (std::getline(rels, line)) { //rels >> line) {
		string line0 = line;

		string nstr = line.substr(0, line.find(separator0));
		int n = stoi(nstr);
		line.erase(0, nstr.length() + 1);
		if (nlast == -1) nlast = n;

		//string Astr = line.substr(0, line.find(separator1));
		//line.erase(0, Astr.length() + 1);
		//int a = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		//int b = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		//int c = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		//int d = atoi(Astr.substr(0, Astr.find(separator2)).c_str());

		if (i >= n || N - s == T + numSch) {
			cout << line << endl;
		}
		else {
			s++;	// skipped
		}

		i++;

		if (n < nlast) {
			nlast = n;
			i = 0;
		}
	}

	return 0;
}

