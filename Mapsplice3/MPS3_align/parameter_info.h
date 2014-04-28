#include <string>
#include <string.h>

using namespace std;

class Parameter_Info
{
public:
	int SEGMENTNUM;
	int POSSIBLE_MAP_CASES_MAX;
	int MAX_SPLICE_LENGTH;

	int START_READ_NUM;
	int SPLICEDISTANCE;

	unsigned int NULL_NUM;
	int MAX_INSERTION_LENGTH;
	int MAX_DELETION_LENGTH;
	int MAPCASES_MAX;

	int ANCHOR_LENGTH;
	int MIN_LENGTH_TO_STITCH;
	int FIX_SPLICE_BUFFER;

	int SHORT_EXON_BUFFER;
	int FIX_MULTI_SPLICE_BUFFER;
	int SHORT_EXON_MIN;

	int SPLICE_MIN_LENGTH; 

	int CANDALILOC;

	int mindValSegLength;

	int SPLICE_MIN_LENGTH;

	int PE_MAP_DISTANCE;

};