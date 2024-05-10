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

bool comparefirst(const pair<int,string> &a,const pair<int,string> &b);
int bisection_search(vector<pair<int,int>> list, int val);

int main (int argc, char** argv)
{
	if (argc != 6 ) {	// argc must be 6
		cout << endl << "Usage: ./smreorder numrels nsm ideals_file smfile_in smfile_out"
			<< endl;
		return 0;
	}

	bool verbose = false;

	string line;
	string separator0 = " ";
	string separator1 = ":";
	string separator2 = ",";

	int numrels = atoi(argv[1]);
	int nsm = atoi(argv[2]);

	// first pass, read ideals, unsorted sm maps
	int nlast = -1; int  r = 0; int num = 0; int mark = 1024; int maxideal = 0;
	ifstream idealsin(argv[3]);
	cout << "Reading ideals..." << endl;
	std::getline(idealsin, line);	// read first line in ideals file
	vector<int> ideals;
	while (std::getline(idealsin, line)) {
		string line0 = line;

		string str = line.substr(0, line.find(separator0));
		line.erase(0, str.length() + 1);

		int ideal = stoi(line, NULL, 16);
		if (ideal > maxideal) maxideal = ideal;
		ideals.push_back(ideal);
		num++;
		
		if (num % mark == 0) {
			cout << num << " ideals read..." << endl;
			mark *= 2;
		}
	}
	cout << num << " ideals read." << endl;
	ifstream smin(argv[4]);
	cout << "Reading unsorted Schirokauer maps..." << endl;
	vector<string> sm0;
	string smin0;
	std::getline(smin, smin0);	// read first line in ideals file
	num = 0; mark = 1024;
	while (std::getline(smin, line)) {

		sm0.push_back(line);
		num++;
		
		if (num % mark == 0) {
			cout << num << " maps read..." << endl;
			mark *= 2;
		}
	}
	cout << num << " maps read." << endl;

	// pair ideals/maps
	vector<pair<int, string>> sm1;
	for (int i = 0; i < numrels-nsm; i++) {
		sm1.push_back(make_pair(ideals[i], sm0[i]));
	}

	// sort list of pairs
	sort(sm1.begin(), sm1.end(), comparefirst);

	// append sm_exp
	for (int i = 0; i < nsm; i++) {
		sm1.push_back(make_pair(maxideal+1+i, sm0[numrels-nsm+i]));
	}

	FILE* out = fopen(argv[5], "wt");
	fprintf(out, "%s\n", smin0.c_str());
	num = 0; mark = 1024;
	for (int i = 0; i < numrels; i++) {
		fprintf(out, "%s\n", sm1[i].second.c_str());
		num++;
		
		if (num % mark == 0) {
			cout << num << " maps written..." << endl;
			mark *= 2;
		}
	}
	cout << num << " maps written." << endl;
		
	fclose(out);

	return 0;
}

bool comparefirst(const pair<int,string> &a,const pair<int,string> &b)
{
	return a.first < b.first;
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


