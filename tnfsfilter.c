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
	if (argc != 3 ) {	// argc must be 3
		cout << endl << "Usage: ./tnfsfilter ideal_count_file relations_file" << endl;
		return 0;
	}

	bool verbose = false;

	string line;
	string separator1 = ":";
	string separator2 = ",";

	// read ideal count file
	ifstream rrfile(argv[1]);
	getline(rrfile, line);
	int L = stoi(line);
	vector<int> C;
	for (int i = 0; i < L; i++) {
		getline(rrfile, line);
		line = line.substr(line.find(separator2)+1, line.length());
		int c = stoi(line);
		C.push_back(c);
	}		

	// read relations file
	ifstream rels(argv[2]);
	while (std::getline(rels, line)) { //rels >> line) {
		string line0 = line;

		string Astr = line.substr(0, line.find(separator1));
		line.erase(0, Astr.length() + 1);
		int a = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		int b = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		int c = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		int d = atoi(Astr.substr(0, Astr.find(separator2)).c_str());

		// find minimum global ideal count for current relation
	 	stringstream s1(line);
		stringstream s2;
		int min = L;
		while(s1.good()) {
			string id;
			getline(s1, id, ',');
			int intid = 0;
			s2.clear();
			s2 << hex << id;
			s2 >> intid;
			if (C[intid] < min) min = C[intid];
		}

		cout << to_string(min) << " " << line0 << endl;
	}

	return 0;
}

