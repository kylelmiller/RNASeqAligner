
#ifndef __CONSTANT_DEFINITIONS_H_INCLUDED__
#define __CONSTANT_DEFINITIONS_H_INCLUDED__


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

#define PE_MAP_DISTANCE 300000


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

const int NUMBER_OF_LETTERS_IN_THE_ALPHABET = 26;
const int FIRST_LEVEL_INDEX_KMER_LENGTH = 14;

const size_t BITS_SUPPORTED = 64;
const size_t LAST_TWO_BIT = 3;
const size_t LAST_THIRD_FOUTH = LAST_TWO_BIT << 2;

const size_t ALL_BITS_ON = static_cast<size_t>(-1);
const size_t LEAST_SIG_BIT = static_cast<size_t>(1);
const size_t MOST_SIG_BIT = static_cast<size_t>(0x8000000000000000);

#endif
