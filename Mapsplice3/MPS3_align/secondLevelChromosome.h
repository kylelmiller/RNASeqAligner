/*****************
 * This class is used to store all
 * data related to the second level of
 * the reference which
 * includes the second level of the: Suffix Array,
 * LCP array, child tables and the reference
 ****************/
#ifndef __SECONDLEVELCHROMOSOME_H_INCLUDED__
#define __SECONDLEVELCHROMOSOME_H_INCLUDED__

#include <vector>

#include "index_info.h"

using namespace std;

class SecondLevelChromosome
{
private:
	char* _secondLevelChrom;
	unsigned int* _secondLevelSa;
	BYTE* _secondLevelLcp;
	unsigned int* _secondLevelChildTab;
	BYTE* _secondLevelDetChild;
	static Index_Info* _indexInfo;

	unsigned int getChildDownValue(unsigned int index)
	{
		if(_secondLevelDetChild[index] == 4)
			return _secondLevelChildTab[_secondLevelChildTab[index]-1];
		else if (_secondLevelDetChild[index] == 2)
			return _secondLevelChildTab[index];
		else
			return 0;
	}

public:
	SecondLevelChromosome(
		char* secondLevelChrom,
		unsigned int* secondLevelSa,
		BYTE* secondLevelLcp,
		unsigned int* secondLevelChildTab,
		BYTE* secondLevelDetChild,
		Index_Info* indexInfo)
	{
		_secondLevelChrom = secondLevelChrom;
		_secondLevelSa = secondLevelSa;
		_secondLevelLcp = secondLevelLcp;
		_secondLevelChildTab = secondLevelChildTab;
		_secondLevelDetChild = secondLevelDetChild;
		_indexInfo = indexInfo;
	}

	/*
	 *
	 */
	bool couldContainSegment(unsigned int length)
	{
		return length > FIRST_LEVEL_INDEX_KMER_LENGTH;
	}

	Index_Info* getIndexInfo()
	{
		return _indexInfo;
	}

	/*
	 * Returns the second level suffix array
	 */
	unsigned int* getSuffixArray()
	{
		return _secondLevelSa;
	}

	/*
	 * Returns the second level reference genome
	 */
	char* getReference()
	{
		return _secondLevelChrom;
	}

	/*
	 * returns the length of the second level lcp inside a given interval
	 */
	unsigned int getLcpLength(unsigned int start, unsigned int end)
	{
		return start < _secondLevelChildTab[end + 1] && end >= _secondLevelChildTab[end +1]
			? _secondLevelLcp[_secondLevelChildTab[end+1]]
			: _secondLevelLcp[_secondLevelDetChild[start]];
	}

	void getFirstInterval(char ch, unsigned int *interval_begin, unsigned int *interval_end)
	{
		switch(ch)
		{
			case 'C':
				*interval_begin = _secondLevelDetChild[0];
				*interval_end = _secondLevelDetChild[*interval_begin] - 1;
				break;
			case 'A':
				*interval_begin = 0;
				*interval_end = _secondLevelDetChild[0] - 1;
				break;
			case 'G':
				*interval_begin = _secondLevelDetChild[_secondLevelDetChild[0]];
				*interval_end = _secondLevelDetChild[*interval_begin] - 1;
				break;
			case 'T':
			default:
				*interval_begin = _secondLevelDetChild[_secondLevelDetChild[_secondLevelDetChild[0]]];
				*interval_end = _secondLevelDetChild[*interval_begin] - 1;
				break;
		}
	}

	void getInterval(unsigned int start, unsigned int end, unsigned int position,
		char ch, unsigned int *interval_begin, unsigned int *interval_end)
	{
		unsigned int index_begin;
		unsigned int pos;
		*interval_end = 0;

		unsigned int child_up_value_tmp = _secondLevelDetChild[end]==1 * _secondLevelChildTab[end];
		*interval_begin = (start < child_up_value_tmp) && (end >= child_up_value_tmp)
			? child_up_value_tmp
			: getChildDownValue(start);

		if(ch == 'C')
		{
			if(_secondLevelChrom[_secondLevelSa[index_begin]+position] == 'C')
			{
				*interval_begin = index_begin;
				*interval_end = (_secondLevelDetChild[index_begin]>2) * _secondLevelChildTab[index_begin] - 1;
				if ((*interval_end < start) || (*interval_end > end))
					*interval_end = end;
			}
			else if(_secondLevelChrom[_secondLevelSa[start]+position] == 'C')
			{
				*interval_begin = start;
				*interval_end = index_begin - 1;
			}
		}
		else if(ch == 'A') // ch == 'A'
		{
			if(_secondLevelChrom[_secondLevelSa[start]+position] == 'A')
			{
				*interval_begin = start;
				*interval_end = index_begin - 1;
			}
		}
		else if(ch == 'G')
		{
			pos = (_secondLevelDetChild[index_begin]>2) * _secondLevelChildTab[index_begin];
			if(_secondLevelChrom[_secondLevelSa[pos]+position] == 'G')
			{
				*interval_begin = pos;
				*interval_end = (_secondLevelDetChild[pos]>2) * _secondLevelChildTab[pos] - 1;
			}
			else if(_secondLevelChrom[_secondLevelSa[start]+position] == 'G')
			{
				*interval_begin = start;
				*interval_end = index_begin - 1;
			}
			else if(_secondLevelChrom[_secondLevelSa[index_begin]+position] == 'G')
			{
				*interval_begin = index_begin;
				*interval_end = (_secondLevelDetChild[index_begin]>2) * _secondLevelChildTab[index_begin] - 1;
			}
		}
		else // if(ch == 'T')
		{
			if(_secondLevelChrom[_secondLevelSa[start]+position] == 'T')
			{
				*interval_begin = start;
				*interval_end = index_begin - 1;
			}
			else if(_secondLevelChrom[_secondLevelSa[index_begin]+position] == 'T')
			{
				*interval_begin = index_begin;
				*interval_end = ((_secondLevelDetChild[index_begin]>2)*_secondLevelChildTab[index_begin]) - 1;
			}
			else
			{
				pos = (_secondLevelDetChild[index_begin]>2) * _secondLevelChildTab[index_begin];
				if(_secondLevelChrom[_secondLevelSa[pos]+position] == 'T')
				{
					*interval_begin = pos;
					*interval_end = ((_secondLevelDetChild[pos]>2)*_secondLevelChildTab[pos]) - 1;
				}
				else
				{
					pos = (_secondLevelDetChild[pos]>2) * _secondLevelChildTab[pos];
					if(_secondLevelChrom[_secondLevelSa[pos]+position] == 'T')
					{
						*interval_begin = pos;
						*interval_end = ((_secondLevelDetChild[pos]>2)*_secondLevelChildTab[pos]) - 1;
					}
				}
			}
		}

		if ((*interval_end < start) || (*interval_end > end))
			*interval_end = end;
	}
};
Index_Info* SecondLevelChromosome::_indexInfo = NULL;

class SecondLevelChromosomeList
{
private:
	vector<SecondLevelChromosome*> _chroms;
	Index_Info* _indexInfo;

public:
	SecondLevelChromosomeList(vector<SecondLevelChromosome*> chroms, Index_Info* indexInfo)
	{
		_chroms = chroms;
		_indexInfo = indexInfo;
	}

	/*
	 * Gets the second level index information based on a given alignment
	 */
	SecondLevelChromosome* getSecondLevelChromosome(string chromosomeName, int chromosomePosition)
	{
		int chrNameInt = _indexInfo->convertStringToInt(chromosomeName);

		// Xinan: need to debug
		int secondLevelIndexNum = _indexInfo->getSecondLevelIndexFromChrAndPos(chrNameInt, chromosomePosition);
		int chrPosStartIn2ndLevelIndex = _indexInfo->getChrPosFromSecondLevelIndexPos(chrNameInt, secondLevelIndexNum, 1);

		if(_indexInfo->invalidSecondLevelIndexNOset.find(secondLevelIndexNum)
			!= _indexInfo->invalidSecondLevelIndexNOset.end())
			return NULL;

		return _chroms[secondLevelIndexNum];
	}
};

#endif
