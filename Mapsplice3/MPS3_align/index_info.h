#ifndef __INDEX_INFO_H_INCLUDED__
#define __INDEX_INFO_H_INCLUDED__

#include <stdlib.h>
#include <string>
#include <string.h>
#include <map>

#include "constantDefinitions.h"
#include "utilities.h"

using namespace std;

class Index_Info
{
private:

	/*
	 * Member Variables
	 */
	int baseCharCount2intArray[FIRST_LEVEL_INDEX_KMER_LENGTH][NUMBER_OF_LETTERS_IN_THE_ALPHABET] = { { 0 } };


public:
	string chromString;
	unsigned int genomeLength;

	int chromNum;
	vector<string> chrNameStr; // size = chromNum
	vector<int> chromLength; // size = chromNum
	vector<string> chromStr;
	vector<unsigned int> chrEndPosInGenome;

	map<string, int> chrNameMap;
	//map<string, int>::iterator chrNameMapIter;

	int secondLevelIndexNormalSize;// = 3000000;
	vector<int> secondLevelIndexPartsNum;
	int secondLevelIndexPartsNumSum;

	set<int> invalidSecondLevelIndexNOset;
	//vector<int> secondLevelIndexLengthVec;

	unsigned int null_num; // 2654911540 for mm9_noRandom genome
	unsigned int indexSize; //2654911539  //sequence length + 1, the length of sa-lcp-down-next 

	string getInvalidSecondLevelIndexNOstr()
	{
		string tmpStr = "invalidSecondLevelIndexNO: \n";
		for(set<int>::iterator setIter = invalidSecondLevelIndexNOset.begin(); setIter != invalidSecondLevelIndexNOset.end(); setIter ++)
		{
			tmpStr += Utilities::int_to_str(*setIter);
			tmpStr += ",";
		} 
		tmpStr += "\n";
		return tmpStr;
	}

	Index_Info()
	{}

	Index_Info(ifstream& inputIndexInfoFile)
	{
		for(int tmp1 = 0; tmp1 < FIRST_LEVEL_INDEX_KMER_LENGTH; tmp1++)
		{
			for(int tmp2 = 0; tmp2 < NUMBER_OF_LETTERS_IN_THE_ALPHABET; tmp2++)
			{
				baseCharCount2intArray[tmp1][tmp2] = 100;
			}
		}
		int tmpBaseCount = 1;
		for(int tmp3 = 0; tmp3 < FIRST_LEVEL_INDEX_KMER_LENGTH; tmp3++)
		{
			baseCharCount2intArray[tmp3][0] = 0*tmpBaseCount;
			baseCharCount2intArray[tmp3][2] = 1*tmpBaseCount;
			baseCharCount2intArray[tmp3][6] = 2*tmpBaseCount;
			baseCharCount2intArray[tmp3][19] = 3*tmpBaseCount;
			tmpBaseCount = 4*tmpBaseCount;
		}


		string s;
		string chromNumLine;
		string chromNameLine;
		string chromEndPosInGenomeLine;
		string secondLevelIndexSizeLine;
		string chrom2ndLevelIndexNumLine;
		getline(inputIndexInfoFile, s);
		getline(inputIndexInfoFile, chromNumLine);
		getline(inputIndexInfoFile, s);
		getline(inputIndexInfoFile, chromNameLine);
		getline(inputIndexInfoFile, s);
		getline(inputIndexInfoFile, chromEndPosInGenomeLine);
		getline(inputIndexInfoFile, s);
		getline(inputIndexInfoFile, secondLevelIndexSizeLine);
		getline(inputIndexInfoFile, s);
		getline(inputIndexInfoFile, chrom2ndLevelIndexNumLine);

		chromNum = atoi( (chromNumLine.substr(0, chromNumLine.length())).c_str() );
		//cout << "chromNum: " << chromNum << endl;

		int startSearchPos = 0;
		int foundSearchPos;
		string tmpChromNameStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromNameLine.find(",", startSearchPos);
			tmpChromNameStr = chromNameLine.substr(startSearchPos+1, foundSearchPos - 2 - startSearchPos - 1 + 1);
			chrNameStr.push_back(tmpChromNameStr);
			//cout << tmp+1 << " tmpChromNameStr: " << tmpChromNameStr << " strLen: " << tmpChromNameStr.length() << endl;
			startSearchPos = foundSearchPos + 1;
		}

		startSearchPos = 0;
		string tmpChromEndPosStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromEndPosInGenomeLine.find(",", startSearchPos);
			tmpChromEndPosStr = chromEndPosInGenomeLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			unsigned int tmpChromEndPos = strtoul(tmpChromEndPosStr.c_str(), NULL, 10);
			chrEndPosInGenome.push_back(tmpChromEndPos);
			//cout << tmp+1 << " tmpChromEndPos: " << tmpChromEndPos << endl;
			startSearchPos = foundSearchPos + 1;
		}		

		secondLevelIndexNormalSize = atoi( (secondLevelIndexSizeLine.substr(0, secondLevelIndexSizeLine.length())).c_str() );
		//cout << "secondLevelIndexNormalSize: " << secondLevelIndexNormalSize << endl;

		startSearchPos = 0;
		string tmpChrom2ndLevelIndexNumStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chrom2ndLevelIndexNumLine.find(",", startSearchPos);
			tmpChrom2ndLevelIndexNumStr = chrom2ndLevelIndexNumLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			int tmpChrom2ndLevelIndexNum = atoi(tmpChrom2ndLevelIndexNumStr.c_str());
			secondLevelIndexPartsNum.push_back(tmpChrom2ndLevelIndexNum);
			//cout << tmp+1 << "tmp2ndLevelIndexNum: " << tmpChrom2ndLevelIndexNum << endl;
			startSearchPos = foundSearchPos + 1;
		}

		indexSize = chrEndPosInGenome[chrEndPosInGenome.size()-1] + 2;
		//cout << "MAX: " << indexSize << endl;
		null_num = indexSize + 1;
		//cout << "NULL_NUM: " << null_num << endl;
		this->buildChrNameMap();
	}

	void buildChrNameMap()
	{
		for(int tmp = 0; tmp < chromNum; tmp++)
			chrNameMap.insert(pair <string, int> (chrNameStr[tmp], tmp));
	}

	int getSecondLevelIndexFromChrAndPos(int chrNameInt, int chrMapPos)
	{
		int tmpTimes = chrMapPos/secondLevelIndexNormalSize;
		int partsTimeBase = 0;
		for(int tmp = 0; tmp < chrNameInt; tmp++)
			partsTimeBase += secondLevelIndexPartsNum[tmp];

		return partsTimeBase + tmpTimes + 1;
	}

	int getChrPosFromSecondLevelIndexPos(int chrNameInt, int secondLevelIndexNum, int secondLevelIndexPos)
	{
		int partsTimeBase = 0;
		for(int tmp = 0; tmp < chrNameInt; tmp++)
			partsTimeBase += secondLevelIndexPartsNum[tmp];

		int tmpSecondLevelIndexNO = secondLevelIndexNum - partsTimeBase;
		return (tmpSecondLevelIndexNO - 1) * secondLevelIndexNormalSize + secondLevelIndexPos;
	}

	void getChrLocation(unsigned int locationInWholeGenome, unsigned int *chr_name_int, unsigned int *chr_local_location)
	{
		if(locationInWholeGenome <= chrEndPosInGenome[0])
		{
			(*chr_name_int) = 0;
			*chr_local_location = locationInWholeGenome;
		}
		else
		{
			for(int tmp = 1; tmp < chrEndPosInGenome.size(); tmp++)
				if( (locationInWholeGenome >= chrEndPosInGenome[tmp-1] + 2) 
					&& (locationInWholeGenome <= chrEndPosInGenome[tmp]) )
				{
					*chr_name_int = tmp;
					*chr_local_location = locationInWholeGenome - chrEndPosInGenome[tmp-1] - 2;
				}
		}
	}

	unsigned int getWholeGenomeLocation(unsigned int chr_name_int, unsigned int locationInWholeGenome)
	{
		unsigned int chr_local_location = 0;
		if(chr_name_int == 0)
		{
			chr_local_location = locationInWholeGenome;
		}
		else if(chr_name_int < chromNum)
		{
			chr_local_location = locationInWholeGenome 
				+ chrEndPosInGenome[chr_name_int-1] + 2;
		}
		else
		{
			cout << "chr_name_int error: " << chr_name_int << endl;
		}
		return chr_local_location;
	}

	int convertStringToInt(const string& chrName)
	{
		
		//omp_init_lock(&lock);
		map<string, int>::iterator chrNameMapIter;
		//omp_set_lock(&lock);
		int chrNameInt = 1000;
		//cout << "chrName = " << chrName << endl;
		chrNameMapIter = chrNameMap.find(chrName);
		if(chrNameMapIter != chrNameMap.end())
		{
			chrNameInt = chrNameMapIter->second;
		}
		else
		{
			cout << "...... chrom name error! ...... " << endl;
		}
		//cout << "chrNameInt: " << chrNameInt;
		//omp_unset_lock(&lock);
		//omp_destroy_lock(&lock);
		return chrNameInt;
	}
};

#endif
