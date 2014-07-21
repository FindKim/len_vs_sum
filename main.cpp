//
//  Main.cpp
//  
//	Driver program for Sequence class
//	Reads file; parses file; prints name, description, and sequence
//
//  Created by Kim Ngo on July 16, 2014.
//
//

#include "ExtractMMSeq.h"
#include "ExtractSigOrf.h"
#include "RCClust.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <fstream>	// Output file
#include <dirent.h>	// Traverse directory
#include <algorithm> // Sort

using namespace std;

const string DIRECTORY = "/afs/crc.nd.edu/user/k/kngo/orig_fasta/";
const string SIGORF_FILE = "/afs/crc.nd.edu/user/k/kngo/leRCC/results/sigOrfs_unmasked/sigOrfs_p1e-05.txt";
const string OUTPUTFILE_SIG_LEN_SUM = "/afs/crc.nd.edu/user/k/kngo/len_vs_min_sum/results/sig_len_min_sum_p1e-05.txt";
const string OUTPUTFILE_NSIG_LEN_SUM = "/afs/crc.nd.edu/user/k/kngo/len_vs_min_sum/results/nsig_len_min_sum_p1e-05.txt";

// Checks if provided directory exists
bool directory_exists(const char* pzPath) {
	if(pzPath == NULL) return false;
	
	DIR *pDir;
	bool bExists = false;
	
	pDir = opendir(pzPath);
	
	if(pDir != NULL) {
		bExists = true;
		(void) closedir (pDir);
	}
	return bExists;
}


// Traverses directory and creates a vector of filenames
vector<string> traverse_directory(string& directory) {

	vector<string> mmfiles;
	DIR *dpdf;
	struct dirent *epdf;
	dpdf = opendir(directory.c_str());
	if (dpdf) {
		while(true) {
			epdf = readdir(dpdf);
			if (epdf == NULL) break;
			string file = directory;
			file.append(string(epdf->d_name));	// Absolute path
			mmfiles.push_back(file);
		}
		closedir(dpdf);
		sort(mmfiles.begin(), mmfiles.end());
	}
	return mmfiles;
}


// Initializes outputfile with two column format
void initialize_outputfile(const string& outputfile) {

//	cout << "Creating " << outputfile << "..." << endl;
	ofstream ofile;
	ofile.open (outputfile.c_str());
	
	if (ofile.is_open()) {
	
		string x_coord = "rcclust length";
		string y_coord = "rcclust %min sum";
		ofile << x_coord << "," << y_coord << endl;
		ofile.close();

	} else cout << "Unable to open " << outputfile << endl;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// -------------------------------------CALCULATES LENGTH,SUM COORD

int main() {

	string directory = DIRECTORY;
	string sigOrf_file = SIGORF_FILE;

	vector<string> mmfiles;
	
	ExtractSigOrf sigOrf(sigOrf_file);
	vector<string> sigOrf_v = sigOrf.get_sigOrf();
	vector<string>::iterator sigOrf_it = sigOrf_v.begin();

	initialize_outputfile(OUTPUTFILE_SIG_LEN_SUM);
	initialize_outputfile(OUTPUTFILE_NSIG_LEN_SUM);

	if (directory_exists(directory.c_str())) {
	
		// Traverses directory and returns vector of filenames
		mmfiles = traverse_directory(directory);
		
		cout << "Computing length and %min sum coordinates for " << SIGORF_FILE << endl;
		
		// Iterates through each file
		vector<string>::iterator file_it = mmfiles.begin();
		for (file_it; file_it != mmfiles.end(); ++file_it) {

			ExtractMMSeq mm(*file_it);

			// Valid file with extension ".fasta.mm.mm"--sorted .mm file
			if (mm.valid_file_extension(*file_it)) {
//				cout << *file_it << endl;

// CALCULATES COORD FOR LEN,SUM
				RCClust coord(sigOrf_it, sigOrf_v, mm.get_mm_orfeome());

				vector<pair<int,float> > s = coord.get_sig_len_sum_v();
				vector<pair<int,float> > ns = coord.get_nsig_len_sum_v();
				sigOrf_it = coord.get_it_pos();

				coord.output_to_file(OUTPUTFILE_SIG_LEN_SUM, s);
				coord.output_to_file(OUTPUTFILE_NSIG_LEN_SUM, ns);
			}
		}
	}
}


