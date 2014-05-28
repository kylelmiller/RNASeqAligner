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

#include "align_info.h"

using namespace std;

class SecondLevelChromosome
{
private:
	char* _secondLevelChrom;
	unsigned int* _secondLevelSa;
	BYTE* _secondLevelLcp;
	unsigned int* _secondLevelChildTab;
	BYTE* _secondLevelDetChild;
	Index_Info* _indexInfo;

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

	Index_Info* getIndexInfo()
	{
		return _indexInfo;
	}
};

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
	SecondLevelChromosome* getSecondLevelChromosome(Alignment_Info* alignmentInfo)
	{
		int chrNameInt = _indexInfo->convertStringToInt(alignmentInfo->alignChromName);

		// Xinan: need to debug
		int secondLevelIndexNum = _indexInfo->getSecondLevelIndexFromChrAndPos(chrNameInt, alignmentInfo->alignChromPos);
		int chrPosStartIn2ndLevelIndex = _indexInfo->getChrPosFromSecondLevelIndexPos(chrNameInt, secondLevelIndexNum, 1);

		if(_indexInfo->invalidSecondLevelIndexNOset.find(secondLevelIndexNum)
			!= _indexInfo->invalidSecondLevelIndexNOset.end())
			return NULL;

		return _chroms[secondLevelIndexNum];
	}
};
#endif
