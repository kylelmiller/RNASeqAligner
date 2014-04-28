#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <map>
#include <hash_map>

#define SPLICE_MIN_LENGTH 50
using namespace std;


typedef map<int, int> SpliceWeightMap; 
typedef map <int, SpliceWeightMap > SpliceJunctionHash; 
typedef SpliceJunctionHash::iterator SpliceJunctionHashIter;
typedef SpliceWeightMap::iterator SpliceWeightMapIter;

SpliceJunctionHash spliceJunction[22];
//SpliceJunctionHash spliceJunctionAcceptor[22]; //<spliceStartPos, spliceEndPos>
SpliceJunctionHashIter iter;
SpliceWeightMapIter weightMapIter;
int groundTruthJunctionNum = 0;
int groundTruthSpliceNum = 0;
int groundTruthDeletionNum = 0;
//int deletionNum = 0;
int tofindJunctionNum = 0;
int tofindSpliceNum = 0;
int tofindDeletionNum = 0;
int covertStringToInt(string chrName)
{
	if(chrName == "chr1")
		return 0;
	else if(chrName == "chr2")
		return 1;
	else if(chrName == "chr3")
		return 2;
	else if(chrName == "chr4")
		return 3;
	else if(chrName == "chr5")
		return 4;
	else if(chrName == "chr6")
		return 5;
	else if(chrName == "chr7")
		return 6;
	else if(chrName == "chr8")
		return 7;
	else if(chrName == "chr9")
		return 8;
	else if(chrName == "chr10")
		return 9;
	else if(chrName == "chr11")
		return 10;
	else if(chrName == "chr12")
		return 11;
	else if(chrName == "chr13")
		return 12;
	else if(chrName == "chr14")
		return 13;
	else if(chrName == "chr15")
		return 14;
	else if(chrName == "chr16")
		return 15;
	else if(chrName == "chr17")
		return 16;
	else if(chrName == "chr18")
		return 17;
	else if(chrName == "chr19")
		return 18;
	else if(chrName == "chrX")
		return 19;
	else if(chrName == "chrY")
		return 20;
	else if(chrName == "chrM")
		return 21;
	else
		return 22;
}

bool insertSpliceJunction2Hash(SpliceJunctionHash* spliceJunction, //SpliceJunctionHash* spliceJunctionAcceptor,
	int chrInt, int spliceStartPos, int spliceEndPos) 
{
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////insert to spliceJunction///////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	groundTruthJunctionNum ++;
	if((spliceEndPos - spliceStartPos) < SPLICE_MIN_LENGTH)
	{
		groundTruthDeletionNum ++;
		return true;
	}
	groundTruthSpliceNum ++;
	iter = spliceJunction[chrInt].find(spliceStartPos); // to find spliceStartPos in Hash
	if(iter == spliceJunction[chrInt].end()) 
	{
		SpliceWeightMap* newSpliceWeightMap = new SpliceWeightMap;
		(*newSpliceWeightMap).insert(pair<unsigned int, unsigned int> (spliceEndPos, 1));  // insert spliceEndPos and weight(1) into Map
		spliceJunction[chrInt].insert(pair<unsigned int, SpliceWeightMap> (spliceStartPos, (*newSpliceWeightMap)));
	}
	else
	{
		weightMapIter = (iter->second).find(spliceEndPos);
		if(weightMapIter == (iter->second).end())
		{
			(iter->second).insert(pair<unsigned int, unsigned int> (spliceEndPos, 1));
		}
		else
		{
			(weightMapIter->second) ++;
		}
	}
	return true;
}

bool findSpliceJunctionInHash(SpliceJunctionHash* spliceJunction, //SpliceJunctionHash* spliceJunctionAcceptor,
	int chrInt, int spliceStartPos, int spliceEndPos) 
{
	//tofindSpliceNum ++;
	if(chrInt > 21)
		return false;
	iter = spliceJunction[chrInt].find(spliceStartPos);
	if(iter == spliceJunction[chrInt].end()) 
	{
		return false;
	}
	else
	{
		weightMapIter = (iter->second).find(spliceEndPos);
		if(weightMapIter == (iter->second).end())
		{
			return false;
		}
		else
		{
			return true;	
		}
	}
}

int main(int argc, char**argv)
{
    if(argc < 4)
	{
		cout << "Executable <InputSpliceJunction1> <InputSpliceJunction2> <Output(Not)FoundSpliceJunctionPrefix>"<< endl;
		exit(0);
	}

	char* Input1 = argv[1];
	char* Input2 = argv[2];
	string outputFile = argv[3];

	string outputFoundSpliceJunction = outputFile + ".found";
	ofstream outputFoundSpliceJunction_ofs(outputFoundSpliceJunction.c_str());
	string outputNotFoundSpliceJunction = outputFile + ".notFound";
	ofstream outputNotFoundSpliceJunction_ofs(outputNotFoundSpliceJunction.c_str());

	FILE *fp_in1 = fopen(Input1, "r");
	string entryString;
	int tabLocation1;
	int tabLocation2;
	int tabLocation3;
	char entry[500];
	int chrInt;
	int spliceStartPos;
	int spliceEndPos;
	string chrIntString;
	string spliceStartPosString;
	string spliceEndPosString;

	int found = 0, notFound = 0;
	//int groundTruthSpliceNum = 0;
	//int toFindSpliceNum = 0;

	fgets(entry, sizeof(entry), fp_in1);
	while(!feof(fp_in1))
	{
		fgets(entry, sizeof(entry), fp_in1);
		//groundTruthSpliceNum ++;
		entryString = entry;
		tabLocation1 = entryString.find('\t', 0);
		tabLocation2 = entryString.find('\t', tabLocation1+1);
		tabLocation3 = entryString.find('\t', tabLocation2+1);
		chrIntString = entryString.substr(0, tabLocation1);
		spliceStartPosString = entryString.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
		spliceEndPosString = entryString.substr(tabLocation2+1, tabLocation3-tabLocation2-1);
		chrInt = covertStringToInt(chrIntString);
		spliceStartPos = atoi(spliceStartPosString.c_str());
		spliceEndPos = atoi(spliceEndPosString.c_str());
		bool insertion = insertSpliceJunction2Hash(spliceJunction, //SpliceJunctionHash* spliceJunctionAcceptor,
			chrInt, spliceStartPos, spliceEndPos);
		if(!insertion)
		{
			cout << "insertion failed " << endl;
			cout << chrInt << " " << spliceStartPos << " " << spliceEndPos << endl;
		}
	}
	FILE *fp_in2 = fopen(Input2, "r");
	fgets(entry, sizeof(entry), fp_in2);
	while(!feof(fp_in2))
	{
		fgets(entry, sizeof(entry), fp_in2);
		tofindJunctionNum++;
		entryString = entry;
		tabLocation1 = entryString.find('\t', 0);
		tabLocation2 = entryString.find('\t', tabLocation1+1);
		tabLocation3 = entryString.find('\t', tabLocation2+1);
		chrIntString = entryString.substr(0, tabLocation1);
		spliceStartPosString = entryString.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
		spliceEndPosString = entryString.substr(tabLocation2+1, tabLocation3-tabLocation2-1);
		chrInt = covertStringToInt(chrIntString);
		spliceStartPos = atoi(spliceStartPosString.c_str());
		spliceEndPos = atoi(spliceEndPosString.c_str());
		if((spliceEndPos-spliceStartPos) < SPLICE_MIN_LENGTH)
		{
			tofindDeletionNum ++;
			continue;
		}
		tofindSpliceNum ++;
		bool findSplice = false;
		for(int tmp1 = -5; tmp1 <= 5; tmp1++)
		{
			if(findSplice)
				break;
			for(int tmp2 = -5; tmp2 <= 5; tmp2++)
			{
				findSplice = findSpliceJunctionInHash(spliceJunction, //SpliceJunctionHash* spliceJunctionAcceptor,
					chrInt, spliceStartPos+tmp1, spliceEndPos+tmp2);
				if(findSplice)
				{
					break;			
				}
			}
		}
		//findSplice = findSpliceJunctionInHash(spliceJunction, //SpliceJunctionHash* spliceJunctionAcceptor,
		//	chrInt, spliceStartPos, spliceEndPos);
		if(findSplice)
		{
			found++;
			outputFoundSpliceJunction_ofs << entryString ;
		}

		else
		{
			notFound++;
			outputNotFoundSpliceJunction_ofs << entryString;
		}

	}
	outputFoundSpliceJunction_ofs << endl;
	outputNotFoundSpliceJunction_ofs << endl;
	cout << "GroundTruthJunctionNum = " << groundTruthJunctionNum << endl;
	cout << "GroundTruthSpliceNum = " << groundTruthSpliceNum << endl;
	cout << "GroundTruthDeletionNum = " << groundTruthDeletionNum << endl << endl;
	cout << "tofindJunctionNum = " << tofindJunctionNum << endl;
	cout << "tofindSpliceNum = " << tofindSpliceNum << endl;
	cout << "tofindDeletionNum = " << tofindDeletionNum << endl;
	cout << "spliceFound = " << found << endl << endl 
	<< "sensitivity: " << 100*(float)(found)/(float)(groundTruthSpliceNum) << endl	
	//<< "sensitivity2: " << 100*(float)(found)/(float)groundTruthSpliceNum << endl
	<< "specificity: " << 100*(float)(found)/(float)tofindSpliceNum << endl << endl; 
	cout << "notFound = " << notFound << " " << (float)(notFound)/(float)tofindSpliceNum << endl;

	return 0;
}