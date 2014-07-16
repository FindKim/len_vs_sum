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
#include "RunLength.h"
#include "MinValue.h"
#include "AvgSeqLength.h"
#include "StdDevLength.h"
#include "RunSum.h"
#include <iostream>
#include <iomanip>	// setprecision for decmial places
#include <string>
#include <vector>
#include <math.h>		// pow()
#include <fstream>	// Output file
#include <dirent.h>	// Traverse directory
#include <algorithm> // Sort
#include <limits>		// numeric_limits to find max value of type
//#include <regex>		// Regex to check for valid file ext WHY DOESN'T THE UNIV SUPPORT C++11?!

const string DIRECTORY = "/afs/crc.nd.edu/user/k/kngo/orig_fasta/";
const string SIGORF_FILE = "/afs/crc.nd.edu/user/k/kngo/leRCC/results/sigOrfs_unmasked/sigOrfs_p1e-05.txt";
const string OUTPUTFILE_CLUSTER_LENGTH = "/afs/crc.nd.edu/user/k/kngo/leRCC/results/sigOrfs_masked/rcclust_length/sig_pruned_masked_p1e-05_rcclust_length.txt";
const string OUTPUTFILE_CLUSTER_SUM = "/afs/crc.nd.edu/user/k/kngo/leRCC/results/sigOrfs_masked/rcclust_sum/sig_pruned_masked_p1e-05_rcclust_sum.txt";
const string OUTPUTFILE_MINVALUE = "/afs/crc.nd.edu/user/k/kngo/leRCC/results/sigOrfs_masked/minvalue/sig_pruned_masked_p1e-05_minvalue.txt";
const string OUTPUTFILE_AVG_LENGTH_T_TEST = "/afs/crc.nd.edu/user/k/kngo/leRCC/results/avg_seq_length_t_test/sigOrfs_avg_seq_length_t_test_p1e-05.txt";
const string OUTPUTFILE_SEQ_LENGTH = "/afs/crc.nd.edu/user/k/kngo/leRCC/results/avg_seq_length_t_test/sigOrfs_seq_length_p1e-05.txt";

using namespace std;


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


// Prints total sig & non-sig runs
void print_total_run (const vector<int>& tot_sig_runs, const vector<int>& tot_non_sig_runs) {
/*
	vector<int>::const_iterator tot_sig_it = tot_sig_runs.begin();
	vector<int>::const_iterator tot_non_sig_it = tot_non_sig_runs.begin();
	for (int i = 0; tot_sig_it != tot_sig_runs.end(); i++, ++tot_sig_it, ++tot_non_sig_it) {
		cout << i << " sig: " << *tot_sig_it << " ";
		cout << "non-sig: " << *tot_non_sig_it << endl;
	}
*/
// Specify print values
	for (int i = 0; i < 25; i++) {
		cout << i << " sig: " << tot_sig_runs[i] << " ";
		cout << "non-sig: " << tot_non_sig_runs[i] << endl;
	}
}


// Resizes vector with 0 as values
void resize_vector(vector<int> &sum, const vector<int> &runs) {
	// Increases size if needed
	if (runs.size() != sum.size()) {
		int resize = runs.size();
		sum.resize(resize, 0);
//		cout << "MAIN ONE RESIZED TO " << resize << endl;
	}
}


// Creates a column-based output file separated by ','
// last argument is what I am calculating--length or abs min value
void create_outputfile(const string& filename, const vector<int>& sig_runs, const vector<int>& non_sig_runs, string& type) {

	cout << "Creating " << filename << "..." << endl;
	
	ofstream ofile;
	ofile.open (filename.c_str());
	
	if (ofile.is_open()) {
	
		ofile << type << ",Sig,Non-Sig" << endl;
		vector<int>::const_iterator sig_it = sig_runs.begin();
		vector<int>::const_iterator non_sig_it = non_sig_runs.begin();
		for (int i = 0; sig_it != sig_runs.end();
		i++, ++sig_it, ++non_sig_it) {
			
			ofile << i << "," << *sig_it << "," << *non_sig_it << endl;
		}
		ofile.close();
		cout << filename << " has been created." << endl;

	} else cout << "Unable to open " << filename << endl;
}


// Creates output file that writes # elements, avg, variance, std dev,
// t value, and determines if diff btwn avg are significantly different
void create_outputfile_t_test(const string& outputfile, const float& sig_tot_num_seq, const float& non_sig_tot_num_seq, const float& sig_avg, const float& non_sig_avg, const float& sig_var, const float& non_sig_var, const float& sig_std_dev, const float& non_sig_std_dev, const float& t_value) {

	cout << "Creating " << outputfile << "..." << endl;
	ofstream ofile;
	ofile.open (outputfile.c_str());
	
	if (ofile.is_open()) {
		ofile << "# t test comparing average sequence lengths of significant & non-sig seqs\n# p-Value 1e-05" << endl << endl;
		ofile << setw(12) << "Variable" << setw(12) << "Sig" << setw(12) << "Non-Sig" << endl;

		ofile << setw(12) << "n elements" << setw(12)
			<< sig_tot_num_seq << setw(12) << non_sig_tot_num_seq << endl;

		ofile.precision(2); // fixed format for two decimal places
		ofile << fixed;
		
		ofile << setw(12) << "Average" << setw(12)
			<< sig_avg << setw(12) << non_sig_avg << endl;

		ofile << setw(12) << "Variance" << setw(12)
			<< sig_var << setw(12) << non_sig_var << endl;

		ofile << setw(12) << "Std Dev" << setw(12)
			<< sig_std_dev << setw(12) << non_sig_std_dev << endl;
		
		ofile << endl;
		ofile << "t Value: " << t_value << endl;
		ofile << endl;

		float degree_of_freedom = (sig_tot_num_seq + non_sig_tot_num_seq)-2;
		if (degree_of_freedom > 120)
			if (t_value > 3.29)
				ofile << "A significant difference at probability = 0.001" << endl;
			else if (t_value > 2.58)
				ofile << "A significant difference at probability = 0.01." << endl;
			else if (t_value > 1.96)
				ofile << "A significant difference at probability = 0.05" << endl;
			else if (t_value > 1.65)
				ofile << "A significant difference at probability = 0.1" << endl;
			else
				ofile << "No siginificant difference." << endl;
		
		else
			ofile << "Select a smaller degree of freedom value from t table." << endl;

		ofile.close();
		cout << outputfile << " has been created." << endl;

	} else cout << "Unable to open " << outputfile << endl;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ---------------------------------------CALCULATES AVG SEQ LENGTH
// ---------------------------------------CALCULATES STD_DEV LENGTH
// -------------------------------------------------COMPUTES T TEST

int main() {

	long double sig_tot_sum_length = 0;
	long double non_sig_tot_sum_length = 0;
	float sig_tot_num_seq = 0;
	float non_sig_tot_num_seq = 0;
	
	vector<int> tot_sig_num_occ(300,0);
	vector<int> tot_non_sig_num_occ(300,0);

	string directory = DIRECTORY;
	string sigOrf_file = SIGORF_FILE;
	string outputfile = OUTPUTFILE_AVG_LENGTH_T_TEST;

	vector<string> mmfiles;
	
	ExtractSigOrf sigOrf(sigOrf_file);
	vector<string> sigOrf_v = sigOrf.get_sigOrf();
	vector<string>::iterator sigOrf_it = sigOrf_v.begin();

	if (directory_exists(directory.c_str())) {
	
		// Traverses directory and returns vector of filenames
		mmfiles = traverse_directory(directory);
		
		cout << "Computing average for " << SIGORF_FILE << endl;
		
		// Iterates through each file
		vector<string>::iterator file_it = mmfiles.begin();
		for (file_it; file_it != mmfiles.end(); ++file_it) {

			ExtractMMSeq mm(*file_it);

			// Valid file with extension ".fasta.mm.mm"--sorted .mm file
			if (mm.valid_file_extension(*file_it)) {
//				cout << *file_it << endl;
// CALCULATES AVG SEQ LENGTH
				AvgSeqLength avg(sigOrf_it, sigOrf_v, mm.get_mm_orfeome());

				vector<int> s = avg.get_sig_num_occ();
				vector<int> ns = avg.get_non_sig_num_occ();

				if (tot_sig_num_occ.size() < s.size() || tot_non_sig_num_occ.size() < ns.size()) {
					resize_vector (tot_sig_num_occ, s);
					resize_vector (tot_non_sig_num_occ, ns);
//					print_total_run(tot_sig_num_occ, tot_non_sig_num_occ);
				}
				avg.add_num_occ_v (tot_sig_num_occ, s);
				avg.add_num_occ_v (tot_non_sig_num_occ, ns);
				sigOrf_it = avg.get_it_pos();

			}
		}
	}
	string type = "Seq Length";
	create_outputfile(OUTPUTFILE_SEQ_LENGTH, tot_sig_num_occ, tot_non_sig_num_occ, type);

}


