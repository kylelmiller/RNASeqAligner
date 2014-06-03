#ifndef __WORD_PAIR_H_INCLUDED__
#define __WORD_PAIR_H_INCLUDED__

#include <string>

#include "constantDefinitions.h"

using namespace std;

struct WordPair {
	/*
	 * Member Variables
	 */
	size_t upper;
	size_t lower;
	size_t bads;

	WordPair(const string &s);
	WordPair() : upper(0), lower(0), bads(0) {}

	size_t score(const WordPair &other, size_t mask) const;

	void left_shift(const size_t i);
	void ps_combine(const size_t prefix_mask, const size_t suffix_mask, const size_t big_buff_mask, const WordPair& suffix_wp, WordPair &wp);

	void duplicate_self(const size_t leftshift, WordPair &wp)
	{
		wp.upper = ((upper << leftshift)  + upper) ;
		wp.lower = ((lower << leftshift)  + lower);

		//why?
		wp.bads =  ((bads << leftshift) + bads);
	}

	static size_t base2int(char c)
	{
		switch(c) {
			case 'a' :
			case 'A' : return 0;
			case 'c' :
			case 'C' : return 1;
			case 'g' :
			case 'G' : return 2;
			case 't' :
			case 'T' : return 3;
			default	 : return 4;
		}
	}

	static bool isValid(char c) { return (base2int(c) != 4); }
	static inline size_t get_upper(const size_t i) {return i > 1;}
	static inline size_t get_lower(const size_t i) {return (i % 2);}
	static inline size_t get_bads(char c) {return !isValid(c);}
};

//return mismatches between two word pairs.
inline size_t
WordPair::score(const WordPair &other, size_t mask) const {
	register size_t bits = ((other.upper ^ upper) |
		(other.lower ^ lower) | other.bads | bads) & mask;

	bits = ((bits & 0xAAAAAAAAAAAAAAAA) >> 1)  + (bits & 0x5555555555555555);
	//  cerr << "bits " << endl << bits2string(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xCCCCCCCCCCCCCCCC) >> 2)  + (bits & 0x3333333333333333);
	//cerr << "bits " << endl << bits2string(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xF0F0F0F0F0F0F0F0) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0F);
	//cerr << "bits " << endl << bits2string(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xFF00FF00FF00FF00) >> 8)  + (bits & 0x00FF00FF00FF00FF);
	//cerr << "bits " << endl << bits2string(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xFFFF0000FFFF0000) >> 16) + (bits & 0x0000FFFF0000FFFF);
	//cerr << "bits " << endl << bits2string(ALL_BITS_ON, bits) << endl;
	// do you want to watch it making the sums?  This is the hypercube summation alg, wow, the reult from the score is not right, debugging
	// at this point right here you would have the sum of the top 32 bits and the bottom 32 bits, each sitting in their half of the 64 bit word
	return ((bits & 0xFFFFFFFF00000000) >> 32) + (bits & 0x00000000FFFFFFFF);
}

WordPair::WordPair(const string &s) : upper(0), lower(0), bads(0) {
	string::const_iterator i = s.begin();
	const string::const_iterator limit = s.end();
	while (i != limit) {
		const char c = base2int(*i) & static_cast<size_t>(3);
		upper = ((upper << 1) + get_upper(c));
		lower = ((lower << 1) + get_lower(c));
		bads  = ((bads << 1) + get_bads(*i));
		++i;
	}
}

//right most bits corresponds to the combined word

void WordPair::left_shift(const size_t i)
{
	upper = upper << i; //only the left seed_width + mid_buff_size is useful though
	lower = lower << i;
	bads =  bads << i;
}

inline
void WordPair::ps_combine(const size_t prefix_mask, const size_t suffix_mask, const size_t big_buff_mask, const WordPair& suffix_wp, WordPair &wp)
{
	wp.upper = ((upper &  prefix_mask) + (suffix_wp.upper & suffix_mask) ) & big_buff_mask;
	wp.lower = ((lower &  prefix_mask) + (suffix_wp.lower & suffix_mask) ) & big_buff_mask;
	wp.bads = ((bads &  prefix_mask) + (suffix_wp.bads & suffix_mask) ) & big_buff_mask;
}

#endif
