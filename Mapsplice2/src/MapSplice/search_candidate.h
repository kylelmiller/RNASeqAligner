#ifndef SEARCH_CANDIDATE_H
#define SEARCH_CANDIDATE_H

#include <string.h>


class Search_Candidate
{
public:
	Search_Candidate()
	{}

	Search_Candidate(string _chrom, string _strand, int _pos1, int _pos2, int _start_seg1, int _end_seg1,
	 int _start_seg2, int _end_seg2, string _left_seq, string _right_seq, string _middle_seq, size_t _bwt_index1, size_t _bwt_index2):
	 	chrom1(_chrom), strand1(_strand), pos1(_pos1), pos2(_pos2), start_seg_no1(_start_seg1), end_seg_no1(_end_seg1), start_seg_no2(_start_seg2),
	 	end_seg_no2(_end_seg2), left_sequence(_left_seq), right_sequence(_right_seq), middle_sequence(_middle_seq), bwt_index1(_bwt_index1), 
	 	bwt_index2(_bwt_index2)
	{
	}

	~Search_Candidate()
	{}

	string chrom1;
	string chrom2;
	string strand1;
	string strand2;
	int pos1;
	int pos2;
	int start_seg_no1;
	int end_seg_no1;
	int start_seg_no2;
	int end_seg_no2;
	string left_sequence;
	string right_sequence;
	string middle_sequence;
	size_t bwt_index1;
	size_t bwt_index2;
	bool fixed;
};


#endif