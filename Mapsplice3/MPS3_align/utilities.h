#ifndef __UTILITIES_H_INCLUDED__
#define __UTILITIES_H_INCLUDED__

#include <string>

#include "wordPair.h"
#include "constantDefinitions.h"

using namespace std;

class Utilities
{
public:
	Utilities() {}

	static string toUpper(string value)
	{
		for(int i=0; i<value.length();i++)
			value[i] = toupper(value[i]);

		return value;
	}

	static string int_to_str(int numerical)
	{
			char c[100];
			sprintf(c,"%d",numerical);
			string str(c);
			return str;
	}


	static bool scoreString(const string& s1, const string& s2,
		size_t max_mismatch, size_t& num_mismatch, size_t& comb_bits)
	{
		if (s1.length() != s2.length())
		{
			//cout << "different length of two strings in score_string from Double_anchored_score.h"<<endl;// <<s1 <<endl << s2 <<endl ;
			return false;
		}

		size_t mask = ALL_BITS_ON >> (BITS_SUPPORTED - s2.length());

		WordPair w1(s1), w2(s2);

		register size_t bits = ((w1.upper ^ w2.upper) |
			(w1.lower ^ w2.lower) | w1.bads | w2.bads) & mask;

		comb_bits = bits;
		bits = ((bits & 0xAAAAAAAAAAAAAAAA) >> 1)  + (bits & 0x5555555555555555);
		bits = ((bits & 0xCCCCCCCCCCCCCCCC) >> 2)  + (bits & 0x3333333333333333);
		bits = ((bits & 0xF0F0F0F0F0F0F0F0) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0F);
		bits = ((bits & 0xFF00FF00FF00FF00) >> 8)  + (bits & 0x00FF00FF00FF00FF);
		bits = ((bits & 0xFFFF0000FFFF0000) >> 16) + (bits & 0x0000FFFF0000FFFF);
		// do you want to watch it making the sums?  This is the hypercube summation alg, wow, the reult from the score is not right, debugging
		// at this point right here you would have the sum of the top 32 bits and the bottom 32 bits, each sitting in their half of the 64 bit word

		num_mismatch = ((bits & 0xFFFFFFFF00000000) >> 32) + (bits & 0x00000000FFFFFFFF);
		if(num_mismatch <= max_mismatch)
			return true;
		else
			return false;
	}

};

#endif
