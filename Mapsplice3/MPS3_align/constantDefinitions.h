#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>

using namespace std;

#define SEGMENTNUM 20
#define POSSIBLE_MAP_CASES_MAX 40
#define MAX_SPLICE_LENGTH 300000
#define READ_ALIGN_AREA_LENGTH 500000
#define START_READ_NUM 0
#define SPLICEDISTANCE 300000
#define PAIR_READ_DISTANCE_MAX 500000
//#define SEGMENTNUM 20
//#define CANDALILOC 100
//#define NULL_NUM 2654911540
  
//#define READ_LENGTH 100
#define MAX_INSERTION_LENGTH 10
//#define MAX_SPLICE_LENGTH 300000
#define MAX_DELETION_LENGTH 10
#define MAPCASES_MAX 200
//#define POSSIBLE_MAP_CASES_MAX 40
#define ANCHOR_LENGTH 10
#define MIN_LENGTH_TO_STITCH 20
#define FIX_SPLICE_BUFFER 5

#define SHORT_EXON_BUFFER 2
#define FIX_MULTI_SPLICE_BUFFER 3
#define SHORT_EXON_MIN 4

#define SPLICE_MIN_LENGTH 50
//#define FIX_SPLICE_BUFFER 5

//#define READ_LENGTH 100
#define CANDALILOC 40

#define minValSegLength 20

typedef unsigned char BYTE;

#define min_anchor_length 8

#define SPLICE_MIN_LENGTH 50

/*
typedef map<unsigned int, unsigned int> SpliceWeightMapForPrint; 
typedef map <unsigned int, SpliceWeightMapForPrint > SpliceJunctionHashForPrint; 
typedef SpliceJunctionHashForPrint::iterator SpliceJunctionHashIterForPrint;
typedef SpliceWeightMapForPrint::iterator SpliceWeightMapIterForPrint;

SpliceJunctionHashForPrint spliceJunctionForPrint[22];
//SpliceJunctionHash spliceJunctionReverse[22];

SpliceJunctionHashIterForPrint spliceJuncHashIterForPrint;
SpliceWeightMapIterForPrint weightMapIterForPrint;*/

/*
string chromStr[22];
string readString;
string chromString;

bool shortHeadFixedBool = false;
bool shortTailFixedBool = false;
*/



//#define CHROM_NUM 22
//#define MAX 2654911539  //sequence length + 1, the length of sa-lcp-down-next
//#define MAX 16300 // for chrM.fa
/*
char chr_name[CHROM_NUM][10];
string chrNameStr[CHROM_NUM];
unsigned int chrEndPosInGenome[CHROM_NUM];
*/
/*
void getChrLocation(unsigned int locationInWholeGenome, unsigned int *chr_name_int, unsigned int *chr_local_location);

int getChr(unsigned int locationInWholeGenome);

unsigned int max(unsigned int a, unsigned int b);
*/

#define PE_MAP_DISTANCE 300000