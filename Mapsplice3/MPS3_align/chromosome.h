#ifndef __CHROMOSOME_H_INCLUDED__
#define __CHROMOSOME_H_INCLUDED__

/*****************
 * This class is used to store all
 * data related to the reference which
 * includes the Suffix Array, LCP array,
 * child tables and the reference
 ****************/

#include "index_info.h"
#include "constantDefinitions.h"
#include "read.h"

class Chromosome
{
private:

	// Constants
	static const short CHARACTER_OFFSET = 4;

	const int baseChar2intArray[NUMBER_OF_LETTERS_IN_THE_ALPHABET] = {0, 100, 1, 100, 100, 100, 2,
				100, 100, 100, 100, 100, 100, 100,
				100, 100, 100, 100, 100, 3,
				100, 100, 100, 100, 100, 100};


	// Member Variables
	unsigned int* _suffixArray;
	BYTE* _lcp;
	unsigned int* _child;
	BYTE* _verifyChild;
	char* _chrom;

	Index_Info* _indexInfo;
	int* _preIndexMappedLengthArray;
	unsigned int* _preIndexIntervalStartArray;
	unsigned int* _preIndexIntervalEndArray;

	unsigned int getChildDownValue(unsigned int index)
	{
		if(_verifyChild[index] == 4)
			return _child[_child[index]-1];
		else if (_verifyChild[index] == 2)
			return _child[index];
		else
			return 0;
	}
public:
	Chromosome(unsigned int* sa, BYTE* lcp, unsigned int* child, BYTE* verifyChild, char* chrom,
			Index_Info* indexInfo, int* preIndexMappedLengthArray,
			unsigned int* preIndexIntervalStartArray, unsigned int* preIndexIntervalEndArray)
	{
		_suffixArray = sa;
		_lcp = lcp;
		_child = child;
		_verifyChild = verifyChild;
		_chrom = chrom;

		_indexInfo = indexInfo;
		_preIndexMappedLengthArray = preIndexMappedLengthArray;
		_preIndexIntervalStartArray = preIndexIntervalStartArray;
		_preIndexIntervalEndArray = preIndexIntervalEndArray;
	}

	~Chromosome()
	{
		delete _suffixArray;
		delete _lcp;
		delete _child;
		delete _verifyChild;
		delete _chrom;

		delete _indexInfo;
		delete _preIndexMappedLengthArray;
		delete _preIndexIntervalStartArray;
		delete _preIndexIntervalEndArray;
	}

	/*
	 * Returns the suffix array
	 */
	unsigned int* getSuffixArray()
	{
		return _suffixArray;
	}

	/*
	 * Returns the reference genome
	 */
	char* getReference()
	{
		return _chrom;
	}

	/*
	 * Returns the index information
	 */
	Index_Info* getIndexInfo()
	{
		return _indexInfo;
	}

	unsigned int getLcp(unsigned int start, unsigned int end)
	{
		unsigned int index = (_verifyChild[end]==1) * _child[end];
		if(start >= index || end < index)
			index = getChildDownValue(start);

		return _lcp[index];
	}

	/*
	 *
	 */
	bool couldContainSegment(unsigned int length)
	{
		return length > FIRST_LEVEL_INDEX_KMER_LENGTH;
	}

	/*
	 * Gets the location of a given read based on INDEX_KMER_LENGTH
	 * nucleotides of the read starting at stop_loc_overall
	 */
	bool getMappedLocation(Read* read, unsigned int startLocation,
		int* mappedLength, unsigned int* indexIntervalStart,
		unsigned int* indexIntervalEnd)
	{
		string stringToSearch = read->getSequence().substr(startLocation, FIRST_LEVEL_INDEX_KMER_LENGTH);

		// If there is an N in our substring then return false since we can't find that
		if (stringToSearch.find("N") != stringToSearch.npos)
			return false;

		unsigned int preIndexNO = 0;
		int baseForCount = 1;

		// Build a hash based on INDEX_KMER_LENGTH nucleotides
		for(int i=FIRST_LEVEL_INDEX_KMER_LENGTH - 1; i>=0; i--)
		{
			preIndexNO = preIndexNO + baseChar2intArray[stringToSearch.at(i) - 'A'] * baseForCount;
			baseForCount = baseForCount * CHARACTER_OFFSET;
		}

		// Get the information about the read based on the generated hash
		*mappedLength = _preIndexMappedLengthArray[preIndexNO];

		// Given the length of the substring this is the first index in the SA that it occurs
		*indexIntervalStart = _preIndexIntervalStartArray[preIndexNO];

		// Given the length of the substring this is the last index in the SA that is occurs
		*indexIntervalEnd = _preIndexIntervalEndArray[preIndexNO];

		// Something was found, return true
		return true;
	}

	void getFirstInterval(char ch, unsigned int *interval_begin, unsigned int *interval_end)
	{
		switch(ch)
		{
			case 'C':
				*interval_begin = (_verifyChild[0]>2) * _child[0];
				*interval_end = (_verifyChild[*interval_begin]>2) * _child[*interval_begin] - 1;
				break;
			case 'A':
				*interval_begin = 0;
				*interval_end = (_verifyChild[0]>2) * _child[0] - 1;
				break;
			case 'G':
			{
				unsigned int g_child_next_tmp = (_verifyChild[0]>2) * _child[0];
				*interval_begin = (_verifyChild[g_child_next_tmp]>2) * _child[g_child_next_tmp];
				*interval_end = (_verifyChild[*interval_begin]>2) * _child[*interval_begin] - 1;
			}
				break;
			case 'T':
			default:
				unsigned int t_child_next_tmp = (_verifyChild[0]>2) * _child[0];
				unsigned int t_child_next_tmp2 = (_verifyChild[t_child_next_tmp]>2) * _child[t_child_next_tmp];
				*interval_begin = (_verifyChild[t_child_next_tmp2]>2) * _child[t_child_next_tmp2];
				*interval_end = (_verifyChild[*interval_begin]>2) * _child[*interval_begin] - 1;
				break;
		}
	}
};

#endif
