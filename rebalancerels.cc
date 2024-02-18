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
	if (argc != 6 ) {	// argc must be 6
		cout << endl << "Usage: ./rebalancerels num_rels target_num relations_file "
			"output_file excess_file" << endl;
		return 0;
	}

	bool verbose = false;

	string line;
	string separator0 = " ";
	string separator1 = ":";
	string separator2 = ",";

	int N = atoi(argv[1]);
	int T = atoi(argv[2]);

	// first pass, read relations and count unique ideals
	vector<int> allideals;

	// read relations file
	int nlast = -1; int  num = 0; int mark = 1024; int maxideal = -1;
	ifstream rels(argv[3]);
	cout << "Reading relations..." << endl;
	while (std::getline(rels, line)) { //rels >> line) {
		string line0 = line;

		string Astr = line.substr(0, line.find(separator1));
		line.erase(0, Astr.length() + 1);
		//int a = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		//int b = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		//int c = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		//int d = atoi(Astr.substr(0, Astr.find(separator2)).c_str());

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

		allideals.insert(allideals.end(), ideals.begin(), ideals.end());

		num++;
		if (num % mark == 0) {
			cout << num << " relations read..." << endl;
			mark *= 2;
		}
	}
	cout << num << " relations read." << endl;
	cout << "all relations read.  Sorting ideals..." << flush;
	sort(allideals.begin(), allideals.end());
	cout << "done." << endl;

	int* list1 = new int[maxideal+1]();  // initialise to 0

	cout << "creating ideal count lists..." << flush;
	int id = -1; int index = -1;
	for (int i = 0; i < allideals.size(); i++) {
		list1[allideals[i]]++;
	}
	cout << "done." << endl;
	
	int numideals = 0;
	for (int i = 0; i < maxideal+1; i++) {
		if (list1[i] > 0) numideals++;
	}
	cout << "There are " << numideals << " unique ideals in total." << endl;
	
	bool maybeskip = true;

	// second pass, reduce number of ideals by T
	cout << "final pass..." << endl;
	num = 0; mark = 1024;
	int k = 0;
	rels.clear();
	rels.seekg(0);
	FILE* out = fopen(argv[4], "wt");
	FILE* excess = fopen(argv[5], "wt");
	while (std::getline(rels, line)) { //rels >> line) {
		string line0 = line;
		
		if (maybeskip) {
			string Astr = line.substr(0, line.find(separator1));
			line.erase(0, Astr.length() + 1);
			//int a = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
			//int b = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
			//int c = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
			//int d = atoi(Astr.substr(0, Astr.find(separator2)).c_str());

			stringstream ss(line);
			stringstream hh;
			vector<int> ideals;

			while(ss.good()) {
				string str;
				getline(ss, str, ',');
				int ideal = stoi(str, NULL, 16);
				ideals.push_back(ideal);
			}

			sort(ideals.begin(), ideals.end());
			ideals.erase(unique(ideals.begin(), ideals.end()), ideals.end());

			int nsingle = 0;
			for (int i = 0; i < ideals.size(); i++) {
				if (list1[ideals[i]] == 1) nsingle++;
			}
			if (nsingle > 1) {
				k += nsingle;
				// now deplete global ideal counts of this removed relation
				for (int i = 0; i < ideals.size(); i++) {
					list1[ideals[i]]--;
				}
				if (k >= T) {
					maybeskip = false;
				}
				num++;
				if (num % mark == 0) {
					cout << num << " relations processed..." << endl;
					mark *= 2;
				}
				fprintf(excess, "%s\n", line0.c_str());
				continue;	// skip printing this relation
			}
		}
		fprintf(out, "%s\n", line0.c_str());
		num++;
		if (num % mark == 0) {
			cout << num << " relations processed..." << endl;
			mark *= 2;
		}
	}
	cout << num << " relations processed." << endl;
	fclose(out);
	fclose(excess);

	delete[] list1;

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


