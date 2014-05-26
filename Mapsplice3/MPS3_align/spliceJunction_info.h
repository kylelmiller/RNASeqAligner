#include <string>
#include <string.h>
#include <map>
#include <set>

using namespace std;

typedef map<string, set<int> > SpliceEndStrHash; // SJ hash (in 2nd hash, key = anchor string)
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

	int anchorStringLength;

	vector<SJintHash> SJintHashNormal;
	vector<SJintHash> SJintHashReverse;

	int areaSize;

	vector<SJareaHash> SJstartPosAreaHash;
	vector<SJareaHash> SJendPosAreaHash;

	SJhash_Info()
	{
		areaSize = 1000;
		anchorStringLength = 3;
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

};
