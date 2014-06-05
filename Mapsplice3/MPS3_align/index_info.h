#ifndef __INDEX_INFO_H_INCLUDED__
#define __INDEX_INFO_H_INCLUDED__

#include <stdlib.h>
#include <stdexcept>
#include <string>
#include <string.h>
#include <map>

#include "chromosome.h"
#include "constantDefinitions.h"
#include "utilities.h"

using namespace std;

class Index_Info
{
private:

	/*
	 * Member Variables
	 */
	int _baseCharCount2intArray[FIRST_LEVEL_INDEX_KMER_LENGTH][NUMBER_OF_LETTERS_IN_THE_ALPHABET] = { { 0 } };

	string _reference;
	vector<Chromosome*> _chromosomes;

	int _secondLevelIndexNormalSize;
	set<int> _invalidSecondLevelIndexNOset;

	/*
	 * Private Methods
	 */


public:

	/*
	 * Constructor
	 */
	Index_Info(ifstream& inputIndexInfoFile, ifstream& chromosomeBitFile)
	{
		for(int i=0; i<FIRST_LEVEL_INDEX_KMER_LENGTH; i++)
			for(int j=0; j<NUMBER_OF_LETTERS_IN_THE_ALPHABET; j++)
				_baseCharCount2intArray[i][j] = 100;

		int baseCount = 1;
		for(int i=0; i<FIRST_LEVEL_INDEX_KMER_LENGTH; i++)
		{
			_baseCharCount2intArray[i][0] = 0 * baseCount;
			_baseCharCount2intArray[i][2] = 1 * baseCount;
			_baseCharCount2intArray[i][6] = 2 * baseCount;
			_baseCharCount2intArray[i][19] = 3 * baseCount;
			baseCount = 4 * baseCount;
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

		int numberOfChromosomes = (atoi((chromNumLine.substr(0, chromNumLine.length())).c_str()));

		int nameSearchPos = 0;
		int endPositionSearchPos = 0;
		int secondLevelIndexSearchPos = 0;

		int foundNameSearchPos;
		int foundEndPositionSearchPos;
		int secondLevelIndexFoundPos;

		for(int i=0; i<numberOfChromosomes; i++)
		{
			foundNameSearchPos = chromNameLine.find(",", nameSearchPos);
			foundEndPositionSearchPos = chromEndPosInGenomeLine.find(",", endPositionSearchPos);
			secondLevelIndexFoundPos = chrom2ndLevelIndexNumLine.find(",", secondLevelIndexSearchPos);

			string endPosition = chromEndPosInGenomeLine.substr(
				endPositionSearchPos,
				foundEndPositionSearchPos - endPositionSearchPos);

			string secondLevelIndex = chrom2ndLevelIndexNumLine.substr(
				secondLevelIndexSearchPos,
				secondLevelIndexFoundPos - secondLevelIndexSearchPos);

			_chromosomes.push_back(new Chromosome(
				chromNameLine.substr(nameSearchPos+1, foundNameSearchPos - 2 - nameSearchPos),
				(int)strtoul(endPosition.c_str(), NULL, 10),
				(int)atoi(secondLevelIndex.c_str())));

			nameSearchPos = foundNameSearchPos + 1;
			endPositionSearchPos = foundEndPositionSearchPos + 1;
			secondLevelIndexSearchPos = secondLevelIndexFoundPos + 1;
		}

		_secondLevelIndexNormalSize = atoi((secondLevelIndexSizeLine.substr(0, secondLevelIndexSizeLine.length())).c_str());

		char *chrom = (char*)malloc(getSize() * sizeof(char));
		chromosomeBitFile.read((char*)chrom, getSize() * sizeof(char));
		_reference = chrom;

		int previousStart = 0;
		for (vector<Chromosome*>::iterator it = _chromosomes.begin(); it != _chromosomes.end(); ++it)
		{
			int endPosition = (*it)->getEndPosition() + 1;
			(*it)->setSequence(_reference.substr(previousStart, endPosition));
			previousStart = endPosition;
		}

		/*
		(indexInfo->chromStr).push_back((indexInfo->chromString).substr(0, (indexInfo->chrEndPosInGenome)[0]+1));
		(indexInfo->chromLength).push_back(((indexInfo->chrEndPosInGenome)[0]+1));
		for(int tmp = 1; tmp < indexInfo->chromNum; tmp++)
		{
			(indexInfo->chromStr).push_back((indexInfo->chromString).substr((indexInfo->chrEndPosInGenome)[tmp-1]+2,
				(indexInfo->chrEndPosInGenome)[tmp]-(indexInfo->chrEndPosInGenome)[tmp-1]-1));
			(indexInfo->chromLength).push_back(((indexInfo->chrEndPosInGenome)[tmp]-(indexInfo->chrEndPosInGenome)[tmp-1]-1));
		}
		*/
		delete chrom;
	}

	~Index_Info()
	{
		for (vector<Chromosome*>::iterator it = _chromosomes.begin(); it != _chromosomes.end(); ++it)
			delete *it;
		_chromosomes.clear();
	}

	/*
	 * FIXME - KLM 6/4/14
	 * THIS WILL BE DELETED
	 */
	string getChromSequence(int index)
	{
		return _chromosomes[index]->getSequence();
	}

	int getChromLength(int index)
	{
		return _chromosomes[index]->getLength();
	}

	string getChromName(int index)
	{
		return _chromosomes[index]->getName();
	}

	const char* getReference()
	{
		return _reference.c_str();
	}

	/*
	 * End OF DELETION
	 */

	unsigned int getSize()
	{
		return _chromosomes.size() == 0
			? 0
			: _chromosomes.back()->getEndPosition() + 2;
	}

	int getNumberOfChromosomes()
	{
		return _chromosomes.size();
	}

	int getSecondLevelIndexNormalSize()
	{
		return _secondLevelIndexNormalSize;
	}

	set<int> getInvalidSecondLevelIndexNOset()
	{
		return _invalidSecondLevelIndexNOset;
	}

	vector<Chromosome*> getChromosomes()
	{
		return _chromosomes;
	}

	// FIXME - THIS SHOULD GO
	int convertStringToInt(const string& value)
	{
		for(int i=0; i <_chromosomes.size(); i++)
			if(_chromosomes[i]->getName() == value)
				return i;

		throw invalid_argument("Invalid Chromosome Name!");
		return -1;
	}

	// FIXME - THIS IS WRONG
	int getSecondLevelIndexFromChrAndPos(int chromosomeIndex, int positionInChromosome)
	{
		int tmpTimes = positionInChromosome / _secondLevelIndexNormalSize;
		int previousChromTotal = 0;
		for (vector<Chromosome*>::iterator it = _chromosomes.begin(); it != _chromosomes.end(); ++it)
			previousChromTotal += (*it)->getPartNumber();

		return previousChromTotal + tmpTimes + 1;
	}

	// FIXME - THIS IS WRONG
	int getChrPosFromSecondLevelIndexPos(int chrNameInt, int secondLevelIndexNum, int secondLevelIndexPos)
	{
		int previousChromTotal = 0;
		for (vector<Chromosome*>::iterator it = _chromosomes.begin(); it != _chromosomes.end(); ++it)
			previousChromTotal += (*it)->getPartNumber();

		int tmpSecondLevelIndexNO = secondLevelIndexNum - previousChromTotal;
		return (tmpSecondLevelIndexNO - 1) * _secondLevelIndexNormalSize + secondLevelIndexPos;
	}

	void getChromosomeLocation(unsigned int locationInWholeGenome, unsigned int *chr_name_int,
		unsigned int *chr_local_location)
	{
		*chr_name_int = 0;
		*chr_local_location = locationInWholeGenome;

		for(int i=1; i<_chromosomes.size(); i++)
			if(locationInWholeGenome >= _chromosomes[i-1]->getEndPosition() + 2
				&& locationInWholeGenome <= _chromosomes[i]->getEndPosition())
				{
					*chr_name_int = i;
					*chr_local_location = locationInWholeGenome - _chromosomes[i-1]->getEndPosition() - 2;
					break;
				}
	}

	/*


	int getSecondLevelIndexFromChrAndPos(int chromosomeIndex, int positionInChromosome)
	{
		int tmpTimes = positionInChromosome / _secondLevelIndexNormalSize;
		int previousChromTotal = 0;
		for(int i=0; i<chromosomeIndex; i++)
			previousChromTotal += _secondLevelIndexPartsNum[i];

		return previousChromTotal + tmpTimes + 1;
	}

	int getChrPosFromSecondLevelIndexPos(int chrNameInt, int secondLevelIndexNum, int secondLevelIndexPos)
	{
		int previousChromTotal = 0;
		for(int i=0; i<chrNameInt; i++)
			previousChromTotal += _secondLevelIndexPartsNum[i];

		int tmpSecondLevelIndexNO = secondLevelIndexNum - previousChromTotal;
		return (tmpSecondLevelIndexNO - 1) * _secondLevelIndexNormalSize + secondLevelIndexPos;
	}

	void getChromosomeLocation(unsigned int locationInWholeGenome, unsigned int *chr_name_int,
		unsigned int *chr_local_location)
	{

		*chr_name_int = 0;
		*chr_local_location = locationInWholeGenome;

		for(int i=1; i<getNumberOfChromosomes(); i++)
			if(locationInWholeGenome >= _chrEndPosInGenome[i-1] + 2
				&& locationInWholeGenome <= _chrEndPosInGenome[i])
				{
					*chr_name_int = i;
					*chr_local_location = locationInWholeGenome - _chrEndPosInGenome[i-1] - 2;
					break;
				}
	}

	int convertStringToInt(const string& chrName)
	{
		map<string, int>::iterator chrNameMapIter;
		int chrNameInt = 1000;
		chrNameMapIter = chrNameMap.find(chrName);
		if(chrNameMapIter != chrNameMap.end())
			chrNameInt = chrNameMapIter->second;
		else
			cout << "...... chrom name error! ...... " << endl;

		return chrNameInt;
	}*/


	string getInvalidSecondLevelIndexNOstr()
	{
		string tmpStr = "invalidSecondLevelIndexNO: \n";
		for(set<int>::iterator setIter = _invalidSecondLevelIndexNOset.begin(); setIter != _invalidSecondLevelIndexNOset.end(); setIter ++)
		{
			tmpStr += Utilities::int_to_str(*setIter);
			tmpStr += ",";
		}
		tmpStr += "\n";
		return tmpStr;
	}
};

#endif
