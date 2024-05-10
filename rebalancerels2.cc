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

int bisection_search(vector<pair<int,int>> list, int val);

int main (int argc, char** argv)
{
	if (argc != 8 ) {	// argc must be 8
		cout << endl << "Usage: ./rebalancerels2 num_fb_ideals num_rels target_excess gap "
			"relations_file output_file unused_rels" << endl;
		return 0;
	}

	bool verbose = false;

	string line;
	string separator0 = " ";
	string separator1 = ":";
	string separator2 = ",";

	int numideals = atoi(argv[1]);
	int N = atoi(argv[2]);
	int T = atoi(argv[3]);
	int gap = atoi(argv[4]);

	// first pass, read relations and count unique ideals
	int* allideals = new int[numideals]();

	// read relations file
	int nlast = -1; int  r = 0; int num = 0; int mark = 1024; int maxideal = -1; int ti = 0;
	ifstream rels(argv[5]);
	FILE* out = fopen(argv[6], "wt");
	FILE* unused = fopen(argv[7], "wt");
	cout << "Reading relations..." << endl;
	while (std::getline(rels, line)) { //rels >> line) {
		string line0 = line;

		string Astr = line.substr(0, line.find(separator1));
		line.erase(0, Astr.length() + 1);

		stringstream ss(line);
		stringstream hh;
		vector<int> ideals;

		while(ss.good()) {
			string str;
			getline(ss, str, ',');
			int ideal = stoi(str, NULL, 16);
			ideals.push_back(ideal);
			if (ideal > maxideal) maxideal = ideal;
		}

		sort(ideals.begin(), ideals.end());
		ideals.erase(unique(ideals.begin(), ideals.end()), ideals.end());

		int ni = 0;
		for (int i = 0; i < ideals.size(); i++)
			if (allideals[ideals[i]] == 0)
				ni++;
		
		if (ni + ti - r >= T) {
			// add relation
			r++;
			fprintf(out, "%s\n", line0.c_str());
			// add ideals
			ti += ni;
			for (int i = 0; i < ideals.size(); i++)
				allideals[ideals[i]] = 1;
		}
		else {
			// skip this relation
			fprintf(unused, "%s\n", line0.c_str());
		}
		
		num++;
		
		if (ti - r < T)
			break;
		
		if (num % mark == 0) {
			cout << num << " relations read..." << endl;
			mark *= 2;
		}
	}
	cout << num << " relations processed." << endl;
	cout << r << " relations in output." << endl;
	cout << ti << " total ideals in relation set." << endl;
	fclose(out);
	fclose(unused);

	return 0;
}

// assume list is sorted on first element of pair
int bisection_search(vector<pair<int,int>> list, int val)
{
	int index = -1;

	int start = 0; int stop = list.size();
	while (start != stop) {
		int mid = start + (stop - start)/2;
		if (val < list[mid].first) {
			stop = mid;
		}
		else if (val > list[mid].first) {
			start = mid;
		}
		else {
			index = mid;
			break;
		}
	}

	return index;
}


