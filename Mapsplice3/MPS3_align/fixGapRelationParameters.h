
/*#ifdef DOREMAPPING_ANNOTATION_RELEASE
bool DO_REMAPPING_DuringFirstMapping= false;
bool BUILD_SPLICE_HASH_FROM_FILE = true;
bool BUILD_SPLICE_HASH_FROM_FIRST_MAPPING = false;
bool PRINT_JUNC = false;
#endif

#ifdef DOREMAPPING_NO_ANNOTATION_RELEASE
bool DO_REMAPPING_DuringFirstMapping = false;
bool BUILD_SPLICE_HASH_FROM_FILE = false;
bool BUILD_SPLICE_HASH_FROM_FIRST_MAPPING = true;
bool PRINT_JUNC = true;
#endif

#ifdef NOT_DOREMAPPING_RELEASE
bool DO_REMAPPING_DuringFirstMapping = false;
bool BUILD_SPLICE_HASH_FROM_FILE = false;
bool BUILD_SPLICE_HASH_FROM_FIRST_MAPPING = false;
bool PRINT_JUNC = false;
#endif

#ifdef DEBUG
bool DO_REMAPPING_DuringFirstMapping= false;
bool BUILD_SPLICE_HASH_FROM_FILE = true;
bool BUILD_SPLICE_HASH_FROM_FIRST_MAPPING = false;
bool PRINT_JUNC = false;
#endif*/

bool DO_REMAPPING_DuringFirstMapping= false;
bool BUILD_SPLICE_HASH_FROM_FILE = true;
bool BUILD_SPLICE_HASH_FROM_FIRST_MAPPING = false;
bool PRINT_JUNC = false;

typedef map<unsigned int, unsigned int> SpliceWeightMapForPrint; 
typedef map <unsigned int, SpliceWeightMapForPrint > SpliceJunctionHashForPrint; 
typedef SpliceJunctionHashForPrint::iterator SpliceJunctionHashIterForPrint;
typedef SpliceWeightMapForPrint::iterator SpliceWeightMapIterForPrint;

SpliceJunctionHashForPrint spliceJunctionForPrint[22];
SpliceJunctionHashIterForPrint spliceJuncHashIterForPrint;
SpliceWeightMapIterForPrint weightMapIterForPrint;


bool insertSpliceJunction2HashForPrint(
	unsigned int chrInt, unsigned int spliceStartPos, unsigned int spliceEndPos) 
{
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////insert to spliceJunction For Print///////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//groundTruthJunctionNum ++;
	if((spliceEndPos - spliceStartPos) < SPLICE_MIN_LENGTH)
	{
		//groundTruthDeletionNum ++;
		return true;
	}
	//groundTruthSpliceNum ++;
	SpliceJunctionHashIterForPrint iter = spliceJunctionForPrint[chrInt].find(spliceStartPos); // to find spliceStartPos in Hash
	if(iter == spliceJunctionForPrint[chrInt].end()) 
	{
		SpliceWeightMapForPrint* newSpliceWeightMap = new SpliceWeightMapForPrint;
		(*newSpliceWeightMap).insert(pair<unsigned int, unsigned int> (spliceEndPos, 1));  // insert spliceEndPos and weight(1) into Map
		spliceJunctionForPrint[chrInt].insert(pair<unsigned int, SpliceWeightMapForPrint> (spliceStartPos, (*newSpliceWeightMap)));
	}
	else
	{
		weightMapIterForPrint = (iter->second).find(spliceEndPos);
		if(weightMapIterForPrint == (iter->second).end())
		{
			(iter->second).insert(pair<unsigned int, unsigned int> (spliceEndPos, 1));
		}
		else
		{
			(weightMapIterForPrint->second) ++;
		}
	}
	return true;
}

bool DO_NONCANONICAL = false;

bool DO_NONCANONICAL_SHORT_ANCHOR = false;

const int FIX_NO_RELATIONSHIP = 0;
const int FIX_TOO_CLOSE = 1;
const int FIX_DOUBLE_ANCHOR = 2;
const int FIX_TOO_FAR = 3;
const int FIX_INSERTION_GAP = 4;
const int FIX_INSERTION_NEIGHBOUR = 5;
const int FIX_DELETION_GAP = 6;
const int FIX_DELETION_NEIGHBOUR = 7;
const int FIX_SPLICE_GAP = 8;
const int FIX_SPLICE_NEIGHBOUR = 9;
const int FIX_MATCH = 10;
const int FIX_SOFTCLIPPING = 11;
const int FIX_REMAPPING_SHORT_HEAD = 12;
const int FIX_REMAPPING_SHORT_TAIL = 13;

