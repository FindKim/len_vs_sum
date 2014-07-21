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
#include <utility>

using namespace std;

class RCClust {

	public:
		RCClust(vector<string>::iterator, const vector<string>&, const vector< pair< string, string> >&);
			// filename, sigOrf iterator, sigOrf_id vec, pair of id & seq
		vector<string>::iterator get_it_pos();
			// Returns position iterator
		vector< pair<int,float> > get_sig_len_sum_v();
			// Returns vector of len,sum coordinates for sig
		vector< pair<int,float> > get_nsig_len_sum_v();
			// Returns vector of len,sum coord for non-sig
		void output_to_file(const string&, const vector< pair<int,float> >&);
			// Appends sig or nsig coord to output file
			// Filename, vector of len, sum coord pair
		void print_coord(const vector< pair<int,float> >&);
			// Prints coordinate pnts for sig/nsig to terminal
	
	private:
		vector< pair<int,float> > sig_len_sum_v;
		vector< pair<int,float> > nsig_len_sum_v;
			// vector of len, sum coord pair
		vector<string>::iterator sigOrf_it;
			// Position of sigOrf iterator
		void set_sigOrf_it_pos(vector<string>::iterator);
			// Sets where the last position the iterator
		vector<float> parse_mm_seq(const string&);
			// Parses minmax sequence string into vector of floats, delim=,
		void calc_len_sum(vector<string>::iterator, const vector<string>&, const vector< pair<string, string> >&);

};

#endif
