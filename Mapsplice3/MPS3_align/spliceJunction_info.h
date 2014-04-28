#include <string>
#include <string.h>
#include <map>
#include <set>

using namespace std;

typedef map<string, vector<int> > SpliceEndStrHash; // SJ hash (in 2nd hash, key = anchor string)
typedef map<int, SpliceEndStrHash> SplicePosHash;
typedef SplicePosHash::iterator SplicePosHashIter; 
typedef SpliceEndStrHash::iterator SpliceEndStrHashIter;

typedef map<int, set<int> > SJintHash;  // SJ hash (only 1 hash, value is all possible splice junction positions in the other side)
typedef SJintHash::iterator SJintHashIter;

typedef map<int, set<int> > SJareaHash; //( areaNO = pos/100 ) intermediate hash to directly get all possible splice junctions near a certain position
typedef SJareaHash::iterator SJareaHashIter;


class SJhash_Info
{
public: 
	vector<SplicePosHash> spliceJunctionNormal; // size = chromNum in index_info file
	vector<SplicePosHash> spliceJunctionReverse; // size = chromNum in index_info file

	vector<SJintHash> SJintHashNormal;
	vector<SJintHash> SJintHashReverse;

	int areaSize;

	vector<SJareaHash> SJstartPosAreaHash;
	vector<SJareaHash> SJendPosAreaHash;

	void initiateSpliceJunctionStringHash(int chromNum)
	{
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			SplicePosHash newSplicePosHash;
			spliceJunctionNormal.push_back(newSplicePosHash);
		}
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			SplicePosHash newSplicePosHash;
			spliceJunctionReverse.push_back(newSplicePosHash);
		}
	}

	void initiateSJintHash(int chromNum)
	{
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			SJintHash newSJintHash;
			SJintHashNormal.push_back(newSJintHash);	
		}
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			SJintHash newSJintHash;
			SJintHashReverse.push_back(newSJintHash);	
		}		
	}

	void insert2SJintHash(
		int chrInt, int spliceStartPos, int spliceEndPos)
	{

		SJintHashIter foundIntHashIter;
		// insert to SJintHashNormal
		foundIntHashIter = SJintHashNormal[chrInt].find(spliceStartPos);
		if(foundIntHashIter == SJintHashNormal[chrInt].end())
		{
			set<int> newPosSet; 
			newPosSet.insert(spliceEndPos);
			SJintHashNormal[chrInt].insert(pair<int, set<int> > (spliceStartPos, newPosSet));
		}
		else
		{
			if((foundIntHashIter->second).find(spliceEndPos) == (foundIntHashIter->second).end())
			{
				(foundIntHashIter->second).insert(spliceEndPos);
			}
			else
			{}
		}

		//insert to SJintHashReverse
		foundIntHashIter = SJintHashReverse[chrInt].find(spliceEndPos);
		if(foundIntHashIter == SJintHashReverse[chrInt].end())
		{
			set<int> newPosSet; 
			newPosSet.insert(spliceStartPos);
			SJintHashReverse[chrInt].insert(pair<int, set<int> > (spliceEndPos, newPosSet));
		}
		else
		{
			if((foundIntHashIter->second).find(spliceStartPos) == (foundIntHashIter->second).end())
			{
				(foundIntHashIter->second).insert(spliceStartPos);
			}
			else
			{}
		}		
	}

	void insert2SJintHashAsClass(Splice_Junction& newSJ)
	{
		//int chrInt =
		//int spliceStartPos = 
		//int spliceEndPos = 
		//this->insert2SJintHash(chrInt, spliceStartPos, spliceEndPos)
	}

	void insert2SJareaHash(
		int chrInt, int spliceStartPos, int spliceEndPos)
	{
		int spliceStartPosArea = spliceStartPos/areaSize;
		int spliceEndPosArea = spliceEndPos/areaSize;

		SJareaHashIter foundAreaHashIter;
		
		// insert to SJstartPosAreaHash;
		foundAreaHashIter = SJstartPosAreaHash[chrInt].find(spliceStartPosArea);
		if(foundAreaHashIter == SJstartPosAreaHash[chrInt].end())
		{
			set<int> newPosSet;
			newPosSet.insert(spliceStartPos);
			SJstartPosAreaHash[chrInt].insert(pair<int, set<int> > (spliceStartPosArea, newPosSet));
		}
		else
		{
			if( (foundAreaHashIter->second).find(spliceStartPos) == (foundAreaHashIter->second).end() )
			{
				(foundAreaHashIter->second).insert(spliceStartPos);
			}
			else
			{}
		}

		
		// insert to SJendPosAreaHash;
		foundAreaHashIter = SJendPosAreaHash[chrInt].find(spliceEndPosArea);
		if(foundAreaHashIter == SJendPosAreaHash[chrInt].end())
		{
			set<int> newPosSet;
			newPosSet.insert(spliceEndPos);
			SJendPosAreaHash[chrInt].insert(pair<int, set<int> > (spliceEndPosArea, newPosSet));
		}
		else
		{
			if( (foundAreaHashIter->second).find(spliceStartPos) == (foundAreaHashIter->second).end() )
			{
				(foundAreaHashIter->second).insert(spliceStartPos);
			}
			else
			{}
		}
	}

	void insert2SJareaHashAsClass(Splice_Junction& newSJ)
	{
		//int chrInt =
		//int spliceStartPos = 
		//int spliceEndPos = 
		//this->insert2SJintHash(chrInt, spliceStartPos, spliceEndPos)
	}

	
};