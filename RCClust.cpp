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

#include "RCClust.h"
#include <string>
#include <vector>
#include <utility>	// pair
#include <iostream>
#include <fstream>	// output to file
#include <sstream> //stringstream peak, ignore

using namespace std;

const int NUM_CODON_WINDOW = 17;
const int NUM_MASK_AA = 50 - NUM_CODON_WINDOW;

RCClust :: RCClust(vector<string>::iterator sigOrf_it, const vector<string>& sigOrf_v, const vector< pair< string, string> >& id_seq) {

	calc_len_sum(sigOrf_it, sigOrf_v, id_seq);
}

// Returns coord vector
vector< pair<int, float> > RCClust :: get_sig_len_sum_v() {
	return sig_len_sum_v;
}

vector< pair<int, float> > RCClust :: get_nsig_len_sum_v() {
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
void RCClust :: output_to_file(const string& filename, const vector< pair<int,float> >& len_sum) {

	fstream file;
	file.open (filename.c_str(), ios::out | ios::app );
		// Appends to file
	
	if (file.is_open()) {
	
		vector< pair<int,float> >::const_iterator coord_it = len_sum.begin();
		for (coord_it; coord_it != len_sum.end(); ++coord_it) {
			file << coord_it->first << "," << coord_it->second << endl;
		}
		file.close();

	} else cout << "Unable to open " << filename << endl;
}


// Prints coordinate pnts for sig/nsig to terminal
void RCClust :: print_coord(const vector< pair<int,float> >& coord) {

	vector< pair<int,float> >::const_iterator coord_it = coord.begin();
	for (coord_it; coord_it != coord.end(); ++coord_it) {
	
		cout << coord_it->first << "\t" << coord_it->second << endl;
	}
}


// 
void RCClust :: calc_len_sum(vector<string>::iterator sigOrf_it, const vector<string>& sigOrf, const vector< pair<string, string> >& id_seq) {

	int rcclust_len_count = 0;
	float rcclust_sum_count = 0;

	// Iterates through pairs in the vector id_seq
	vector< pair< string, string> >::const_iterator id_seq_it = id_seq.begin();	
	for (id_seq_it; id_seq_it != id_seq.end(); ++id_seq_it) {

		// Parse the minmax values string by ','
		vector<float> mm_number_v;
		mm_number_v = parse_mm_seq(id_seq_it->second);

		// Masks the first 50 amino acids considering 17 codon windows
		// Iterates through min max values and counts for min runs
		if (mm_number_v.size() > NUM_MASK_AA) {
			
			//cout << id_seq_it->first << endl;
			vector<float>::iterator it = mm_number_v.begin();
			for (it; it != mm_number_v.end(); ++it) {
				rcclust_len_count = 0;	// Reset count every max value break
				rcclust_sum_count = 0;

				while (*it < 0) {				// Continuous min values
					rcclust_len_count++;
					rcclust_sum_count += *it;

					// Prevents from running off vec
					if (it+1 != mm_number_v.end())
						++it;
					else break;					
				}
				if (rcclust_len_count > 0) {
				
					if (*sigOrf_it == id_seq_it->first) {
						pair<int,float> coord = make_pair(rcclust_len_count, rcclust_sum_count);
						sig_len_sum_v.push_back(coord);
//						cout << *sigOrf_it << " " << coord.first << "," << coord.second << endl;
					
					} else {
						pair<int,float> coord = make_pair(rcclust_len_count, rcclust_sum_count);
						nsig_len_sum_v.push_back(coord);
					}
				}
			}
		}
		if (*sigOrf_it == id_seq_it->first && sigOrf_it+1 != sigOrf.end())
			++sigOrf_it;
	}
	set_sigOrf_it_pos(sigOrf_it);
}


