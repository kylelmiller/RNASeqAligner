
#ifndef __FIX_PHASE_ONE_H_INCLUDED__
#define __FIX_PHASE_ONE_H_INCLUDED__

#include <string>
#include <string.h>

#include "pairedEndRead.h"
#include "mappedRead.h"
#include "path.h"
#include "gap.h"
#include "chromosome.h"

using namespace std;

class FixPhase1Info
{
public:

	MappedRead* firstNormalMappedRead;
	MappedRead* firstReverseCompMappedRead;
	MappedRead* secondNormalMappedRead;
	MappedRead* secondReverseCompMappedRead;

	Path* pathInfo_Nor1;
	Path* pathInfo_Rcm1;
	Path* pathInfo_Nor2;
	Path* pathInfo_Rcm2;	

	Gap* gapInfo_Nor1;
	Gap* gapInfo_Rcm1;
	Gap* gapInfo_Nor2;
	Gap* gapInfo_Rcm2;

	FixPhase1Info(PairedEndRead* pairedEndRead, Chromosome* chrom)
	{
	   	firstNormalMappedRead = new MappedRead(
			pairedEndRead->getFirstRead(),
			chrom);
	    firstReverseCompMappedRead = new MappedRead(
				pairedEndRead->getFirstReadReverseComplement(),
				chrom);

		secondNormalMappedRead = new MappedRead(
				pairedEndRead->getSecondRead(),
				chrom);
		secondReverseCompMappedRead = new MappedRead(
				pairedEndRead->getSecondReadReverseComplement(),
				chrom);

		pathInfo_Nor1 = new Path();
		pathInfo_Rcm1 = new Path();
		pathInfo_Nor2 = new Path();
		pathInfo_Rcm2 = new Path();

		gapInfo_Nor1 = new Gap();
		gapInfo_Rcm1 = new Gap();			
		gapInfo_Nor2 = new Gap();
		gapInfo_Rcm2 = new Gap();		
	}

	~FixPhase1Info()
	{
		delete gapInfo_Nor1;
		delete gapInfo_Rcm1;
		delete gapInfo_Nor2;
		delete gapInfo_Rcm2;

		delete pathInfo_Nor1;
		delete pathInfo_Rcm1;
		delete pathInfo_Nor2;
		delete pathInfo_Rcm2;

		delete firstNormalMappedRead;
		delete secondNormalMappedRead;
		delete firstReverseCompMappedRead;
		delete secondReverseCompMappedRead;
	}

	void fixPhase1_pathInfo()
	{
		pathInfo_Nor1->getPossiPathFromSeg(firstNormalMappedRead);
		pathInfo_Rcm1->getPossiPathFromSeg(firstReverseCompMappedRead);

		pathInfo_Nor2->getPossiPathFromSeg(secondNormalMappedRead);
		pathInfo_Rcm2->getPossiPathFromSeg(secondReverseCompMappedRead);
	}

	void fixPhase1_gapInfo(Index_Info* indexInfo)
	{
		gapInfo_Nor1->fixGapInPath(
			pathInfo_Nor1,
			firstNormalMappedRead,
			indexInfo);

		gapInfo_Rcm1->fixGapInPath(
			pathInfo_Rcm1,
			firstReverseCompMappedRead,
			indexInfo);

		gapInfo_Nor2->fixGapInPath(
			pathInfo_Nor2,
			secondNormalMappedRead,
			indexInfo);

		gapInfo_Rcm2->fixGapInPath(
			pathInfo_Rcm2,
			secondReverseCompMappedRead,
			indexInfo);
	}
};

#endif
