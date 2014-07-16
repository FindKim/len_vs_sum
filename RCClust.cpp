//
//  RRClust.h
//  
//
//  Created by Kim Ngo July 16, 2014.

//	Calculates rare codon cluster %min value sum
//	Calculates rare codon cluster length
//  Compares sum & length for sig sequences & non-sig seqs
//	SigOrf files are from Aaron's pruning results
//	mm_orfeome files are from minmax results
//	Outputs csv file with length as x and sum as y coords

#ifndef RCCLUST_H
#define RCCLUST_H

#include <string>
#include <vector>
#include <utility>	// pair
#include <iostream>
#include <fstream>	// output to file

using namespace std;

const int NUM_CODON_WINDOW = 17;
const int NUM_MASK_AA = 50 - NUM_CODON_WINDOW;

RCClust :: RCClust(const string& filename, vector<string>::iterator sigOrf_it, const vector<string>& sigOrf_v, const vector< pair< string, string> >& id_seq) {

	set_sig_len_sum_v(300);
	set_nsig_len_sum_v(300);
	calc_len_sum(sigOrf_it, sigOrf_v, id_seq);
}


// Instantiates len_sum_v with n # of 0's
void RCClust :: set_sig_len_sum_v(int n) {
	vector<int> temp(n, 0);
	sig_len_sum_v = temp;
}

void RCClust :: set_nsig_len_sum_v(int n) {
	vector<int> temp(n, 0);
	nsig_len_sum_v = temp;
}

// Returns coord vector
vector< pair<int, int> > RCClust :: get_sig_len_sum_v() {
	return sig_len_sum_v;
}

vector< pair<int, int> > RCClust :: get_nsig_len_sum_v() {
	return nsig_len_sum_v;
}

// Returns the position the significant organism iterator was last at
vector<string>::iterator RCClust :: get_it_pos() {
	return sigOrf_it;
}

// Sets where the last position the iterator pointed to
void RCClust :: set_sigOrf_it_pos(vector<string>::iterator it) {
	sigOrf_it = it;
}


// Parses minmax sequence string into vector of floats
vector<float> RCClust :: parse_mm_seq(const string& mm_seq) {

	float i;
	vector<float> mm_number_v;
	
	stringstream ss(mm_seq);
	while (ss >> i) {
		mm_number_v.push_back(i);
		if (ss.peek() == ',')
			ss.ignore();
	}
	return mm_number_v;
}


// Appends sig or nsig coord to output file
// Filename, vector of len, sum coord pair
void RCClust :: output_to_file(const string& filename, const vector< pair<int, int> >& len_sum) {

	fstream file;
	file.open (filename, fstream::app);
		// Appends to file
	
	if (file.is_open) {
		vector<pair<int, int> >::iterator coord_it;
		coord_it = len_sum.begin();
		for (coord_it; coord_it != len_sum.end(); ++coord_it) {
		
			file << coord_it->first << "," << coord_it->second << endl;
		}
	
		file.close();

	} else cout << "Unable to open " << filename << endl;
}


// Prints coordinate pnts for sig/nsig to terminal
void RCClust :: print_coord(const vector< pair<int, int> >& coord) {

	vector<pair<int, int> >::iterator coord_it;
	coord_it = coord.begin();
	for (coord_it; coord_it != coord.end(); ++coord_it) {
	
		cout << coord_it->first << "\t" << coord_it->second << endl;
	}
}


// 
void calc_len_sum(vector<string>::iterator sigOrf_it, const vector<string>& sigOrf, const vector< pair<string, string> >& id_seq) {

	
}


