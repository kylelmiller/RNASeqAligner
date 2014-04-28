#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <map>

#include "constantDefinitions.h"
//#include "defineCompilePattern.h"
//#include "segmentMapping.h"
//#include "read_block_test.h"
//#include "bwtmap_info.h"
//#include "DoubleAnchorScore.h"
//#include "sbndm.h"
//#include "splice_info.h"

using namespace std;  

//#define CHROM_NUM 22
//#define MAX 2654911539  //sequence length + 1, the length of sa-lcp-down-next
//#define MAX 16300 // for chrM.fa
/*
const char chr_name[CHROM_NUM][10] = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
	 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
	 "chr16", "chr17", "chr18", "chr19", "chrX", "chrY", "chrM"};

const string chrNameStr[CHROM_NUM] = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
	 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
	 "chr16", "chr17", "chr18", "chr19", "chrX", "chrY", "chrM"};	 

const int secondLevelIndexPartsNum[CHROM_NUM] = {66, 61, 54, 52, 51, 50, 51, 44, 42, 44,
	 41, 41, 41, 42, 35, 33, 32, 31, 21, 56, 6, 1};
*/	 
/*const int secondLevelIndexPartsNumSum = 66 + 61 + 54 + 52 + 51 + 50 + 51 + 44 + 42 + 44 +
	 41 + 41 + 41 + 42 + 35 + 33 + 32 + 31 + 21 + 56 + 6 + 1;*/
/*
const int secondLevelIndexNormalSize = 3000000;	 

int getSecondLevelIndexFromChrAndPos(int chrNameInt, int chrMapPos)
{
	int tmpTimes = chrMapPos/secondLevelIndexNormalSize;
	int partsTimeBase = 0;
	for(int tmp = 0; tmp < chrNameInt; tmp++)
	{
		partsTimeBase += secondLevelIndexPartsNum[tmp];
	}

	return	(partsTimeBase + tmpTimes + 1); 
}

int getChrPosFromSecondLevelIndexPos(int chrNameInt, int secondLevelIndexNum, int secondLevelIndexPos)
{
	int partsTimeBase = 0;
	for(int tmp = 0; tmp < chrNameInt; tmp++)
	{
		partsTimeBase += secondLevelIndexPartsNum[tmp];
	}

	int tmpSecondLevelIndexNO = secondLevelIndexNum - partsTimeBase;
	return ( (tmpSecondLevelIndexNO-1) * secondLevelIndexNormalSize + secondLevelIndexPos);	
}*/
/*
const unsigned int chrEndPosInGenome[CHROM_NUM] = {197195431, 378943519, 538543303, 694173424,
	846710684, 996227722, 1148752276, 1280491148, 1404567321, 1534560577, 1656404434,
	1777661965, 1897946278, 2023141143, 2126636118, 2224955269, 2320227921, 2410999953,
	2472342384, 2638992681, 2654895237, 2654911538};	*/
/*
void getChrLocation(unsigned int locationInWholeGenome, unsigned int *chr_name_int, unsigned int *chr_local_location)
{
	if ((0 <= locationInWholeGenome) && (locationInWholeGenome <= 197195431))
	{
		*chr_name_int = 0;
		*chr_local_location = locationInWholeGenome;
	}
	else if ((197195433 <= locationInWholeGenome) && (locationInWholeGenome <= 378943519))
	{
		*chr_name_int = 1;
		*chr_local_location = locationInWholeGenome-197195433;
	}
	else if ((378943521 <= locationInWholeGenome) && (locationInWholeGenome <= 538543303))
	{
		*chr_name_int = 2;
		*chr_local_location = locationInWholeGenome-378943521;
	}
	else if ((538543305 <= locationInWholeGenome) && (locationInWholeGenome <= 694173424))
	{
		*chr_name_int = 3;
		*chr_local_location = locationInWholeGenome-538543305;
	}
	else if ((694173426 <= locationInWholeGenome) && (locationInWholeGenome <= 846710684))
	{
		*chr_name_int = 4;
		*chr_local_location = locationInWholeGenome-694173426;
	}
	else if ((846710686 <= locationInWholeGenome) && (locationInWholeGenome <= 996227722))
	{
		*chr_name_int = 5;
		*chr_local_location = locationInWholeGenome-846710686;
	}
	else if ((996227724 <= locationInWholeGenome) && (locationInWholeGenome <= 1148752276))
	{
		*chr_name_int = 6;
		*chr_local_location = locationInWholeGenome-996227724;
	}
	else if ((1148752278 <= locationInWholeGenome) && (locationInWholeGenome <= 1280491148))
	{
		*chr_name_int = 7;
		*chr_local_location = locationInWholeGenome-1148752278;
	}
	else if ((1280491150 <= locationInWholeGenome) && (locationInWholeGenome <= 1404567321))
	{
		*chr_name_int = 8;
		*chr_local_location = locationInWholeGenome-1280491150;
	}
	else if ((1404567323 <= locationInWholeGenome) && (locationInWholeGenome <= 1534560577))
	{
		*chr_name_int = 9;
		*chr_local_location = locationInWholeGenome-1404567323;
	}
	else if ((1534560579 <= locationInWholeGenome) && (locationInWholeGenome <= 1656404434))
	{
		*chr_name_int = 10;
		*chr_local_location = locationInWholeGenome-1534560579;
	}
	else if ((1656404436 <= locationInWholeGenome) && (locationInWholeGenome <= 1777661965))
	{
		*chr_name_int = 11;
		*chr_local_location = locationInWholeGenome-1656404436;
	}
	else if ((1777661967 <= locationInWholeGenome) && (locationInWholeGenome <= 1897946278))
	{
		*chr_name_int = 12;
		*chr_local_location = locationInWholeGenome-1777661967;
	}
	else if ((1897946280 <= locationInWholeGenome) && (locationInWholeGenome <= 2023141143))
	{
		*chr_name_int = 13;
		*chr_local_location = locationInWholeGenome-1897946280;
	}
	else if ((2023141145 <= locationInWholeGenome) && (locationInWholeGenome <= 2126636118))
	{
		*chr_name_int = 14;
		*chr_local_location = locationInWholeGenome-2023141145;
	}
	else if ((2126636120 <= locationInWholeGenome) && (locationInWholeGenome <= 2224955269))
	{
		*chr_name_int = 15;
		*chr_local_location = locationInWholeGenome-2126636120;
	}
	else if ((2224955271 <= locationInWholeGenome) && (locationInWholeGenome <= 2320227921))
	{
		*chr_name_int = 16;
		*chr_local_location = locationInWholeGenome-2224955271;
	}
	else if ((2320227923 <= locationInWholeGenome) && (locationInWholeGenome <= 2410999953))
	{
		*chr_name_int = 17;
		*chr_local_location = locationInWholeGenome-2320227923;
	}		
	else if ((2410999955 <= locationInWholeGenome) && (locationInWholeGenome <= 2472342384))
	{
		*chr_name_int = 18;
		*chr_local_location = locationInWholeGenome-2410999955;
	}
	else if ((2472342386 <= locationInWholeGenome) && (locationInWholeGenome <= 2638992681))
	{
		*chr_name_int = 19;
		*chr_local_location = locationInWholeGenome-2472342386;
	}
	else if ((2638992683 <= locationInWholeGenome) && (locationInWholeGenome <= 2654895237))
	{
		*chr_name_int = 20;
		*chr_local_location = locationInWholeGenome-2638992683;
	}
	else if ((2654895239 <= locationInWholeGenome) && (locationInWholeGenome <= 2654911538))
	{
		*chr_name_int = 21;
		*chr_local_location = locationInWholeGenome-2654895239;
	}
}

unsigned int getWholeGenomeLocation(unsigned int chr_name_int, unsigned int locationInWholeGenome)
{
	//the names of locationInWholeGenome and chr_local_location should be exchanged. 
	//
	unsigned int chr_local_location;
	if (chr_name_int = 0)
	{
		chr_local_location = locationInWholeGenome;
	}
	else if (chr_name_int = 1)
	{
		chr_local_location = locationInWholeGenome + 197195433;
	}
	else if (chr_name_int = 2)
	{
		chr_local_location = locationInWholeGenome+378943521;
	}
	else if (chr_name_int = 3)
	{
		chr_local_location = locationInWholeGenome+538543305;
	}
	else if (chr_name_int = 4)
	{
		chr_local_location = locationInWholeGenome+694173426;
	}
	else if (chr_name_int = 5)
	{
		chr_local_location = locationInWholeGenome+846710686;
	}
	else if (chr_name_int = 6)
	{
		chr_local_location = locationInWholeGenome+996227724;
	}
	else if (chr_name_int = 7)
	{
		chr_local_location = locationInWholeGenome+1148752278;
	}
	else if (chr_name_int = 8)
	{
		chr_local_location = locationInWholeGenome+1280491150;
	}
	else if (chr_name_int = 9)
	{
		chr_local_location = locationInWholeGenome+1404567323;
	}
	else if (chr_name_int = 10)
	{
		chr_local_location = locationInWholeGenome+1534560579;
	}
	else if (chr_name_int = 11)
	{
		chr_local_location = locationInWholeGenome+1656404436;
	}
	else if (chr_name_int = 12)
	{
		chr_local_location = locationInWholeGenome+1777661967;
	}
	else if (chr_name_int = 13)
	{
		chr_local_location = locationInWholeGenome+1897946280;
	}
	else if (chr_name_int = 14)
	{
		chr_local_location = locationInWholeGenome+2023141145;
	}
	else if (chr_name_int = 15)
	{
		chr_local_location = locationInWholeGenome+2126636120;
	}
	else if (chr_name_int = 16)
	{
		chr_local_location = locationInWholeGenome+2224955271;
	}
	else if (chr_name_int = 17)
	{
		chr_local_location = locationInWholeGenome+2320227923;
	}		
	else if (chr_name_int = 18)
	{
		chr_local_location = locationInWholeGenome+2410999955;
	}
	else if (chr_name_int = 19)
	{
		chr_local_location = locationInWholeGenome+2472342386;
	}
	else if (chr_name_int = 20)
	{
		chr_local_location = locationInWholeGenome+2638992683;
	}
	else if (chr_name_int = 21)
	{
		chr_local_location = locationInWholeGenome+2654895239;
	}
	else
	{
		cout << "error: no that chromosomes" << endl;
	}
	return chr_local_location;
}

int getChr(unsigned int locationInWholeGenome)
{
	if ((0 <= locationInWholeGenome) && (locationInWholeGenome <= 197195431))
	{
		return 0;
	}
	else if ((197195433 <= locationInWholeGenome) && (locationInWholeGenome <= 378943519))
	{
		return 1;
	}
	else if ((378943521 <= locationInWholeGenome) && (locationInWholeGenome <= 538543303))
	{
		return 2;
	}
	else if ((538543305 <= locationInWholeGenome) && (locationInWholeGenome <= 694173424))
	{
		return 3;
	}
	else if ((694173426 <= locationInWholeGenome) && (locationInWholeGenome <= 846710684))
	{
		return 4;
	}
	else if ((846710686 <= locationInWholeGenome) && (locationInWholeGenome <= 996227722))
	{
		return 5;
	}
	else if ((996227724 <= locationInWholeGenome) && (locationInWholeGenome <= 1148752276))
	{
		return 6;
	}
	else if ((1148752278 <= locationInWholeGenome) && (locationInWholeGenome <= 1280491148))
	{
		return 7;
	}
	else if ((1280491150 <= locationInWholeGenome) && (locationInWholeGenome <= 1404567321))
	{
		return 8;
	}
	else if ((1404567323 <= locationInWholeGenome) && (locationInWholeGenome <= 1534560577))
	{
		return 9;
	}
	else if ((1534560579 <= locationInWholeGenome) && (locationInWholeGenome <= 1656404434))
	{
		return 10;
	}
	else if ((1656404436 <= locationInWholeGenome) && (locationInWholeGenome <= 1777661965))
	{
		return 11;
	}
	else if ((1777661967 <= locationInWholeGenome) && (locationInWholeGenome <= 1897946278))
	{
		return 12;
		//*chr_local_location = locationInWholeGenome-1777661967;
	}
	else if ((1897946280 <= locationInWholeGenome) && (locationInWholeGenome <= 2023141143))
	{
		return 13;
		//*chr_local_location = locationInWholeGenome-1897946280;
	}
	else if ((2023141145 <= locationInWholeGenome) && (locationInWholeGenome <= 2126636118))
	{
		return 14;
		//*chr_local_location = locationInWholeGenome-2023141145;
	}
	else if ((2126636120 <= locationInWholeGenome) && (locationInWholeGenome <= 2224955269))
	{
		return 15;
		//*chr_local_location = locationInWholeGenome-2126636120;
	}
	else if ((2224955271 <= locationInWholeGenome) && (locationInWholeGenome <= 2320227921))
	{
		return 16;
		//*chr_local_location = locationInWholeGenome-2224955271;
	}
	else if ((2320227923 <= locationInWholeGenome) && (locationInWholeGenome <= 2410999953))
	{
		return 17;
		//*chr_local_location = locationInWholeGenome-2320227923;
	}		
	else if ((2410999955 <= locationInWholeGenome) && (locationInWholeGenome <= 2472342384))
	{
		return 18;
		//*chr_local_location = locationInWholeGenome-2410999955;
	}
	else if ((2472342386 <= locationInWholeGenome) && (locationInWholeGenome <= 2638992681))
	{
		return 19;
		//*chr_local_location = locationInWholeGenome-2472342386;
	}
	else if ((2638992683 <= locationInWholeGenome) && (locationInWholeGenome <= 2654895237))
	{
		return 20;
		//*chr_local_location = locationInWholeGenome-2638992683;
	}
	else if ((2654895239 <= locationInWholeGenome) && (locationInWholeGenome <= 2654911538))
	{
		return 21;
		//*chr_local_location = locationInWholeGenome-2654895239;
	}
}

int convertStringToInt(const string& chrName)
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
}*/