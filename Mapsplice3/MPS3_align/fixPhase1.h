#include <string>
#include <string.h>
#include "splice_info.h"

using namespace std;

class FixPhase1Info
{
public:

	Seg_Info* segInfo_Nor1;
	Seg_Info* segInfo_Rcm1;
	Seg_Info* segInfo_Nor2;
	Seg_Info* segInfo_Rcm2;

	bool normalMapMain;
	bool rcmMapMain;
	bool normalMapMain_PE;
	bool rcmMapMain_PE;

	Path_Info* pathInfo_Nor1;
	Path_Info* pathInfo_Rcm1;
	Path_Info* pathInfo_Nor2;
	Path_Info* pathInfo_Rcm2;	

	Gap_Info* gapInfo_Nor1;
	Gap_Info* gapInfo_Rcm1;
	Gap_Info* gapInfo_Nor2;
	Gap_Info* gapInfo_Rcm2;

	unsigned int norValLength;// = 0;
	unsigned int norValLength_PE;// = 0;
	unsigned int rcmValLength;// = 0;
	unsigned int rcmValLength_PE;// = 0;
	unsigned int minValLengthToStitch;// = MIN_LENGTH_TO_STITCH;
	
	FixPhase1Info()
	{
		norValLength = 0;
		norValLength_PE = 0;
		rcmValLength = 0;
		rcmValLength_PE = 0;
		minValLengthToStitch = MIN_LENGTH_TO_STITCH;

	   	segInfo_Nor1 = new Seg_Info();		
	    segInfo_Rcm1 = new Seg_Info();
		segInfo_Nor2 = new Seg_Info();
		segInfo_Rcm2 = new Seg_Info();	

		pathInfo_Nor1 = new Path_Info();
		pathInfo_Rcm1 = new Path_Info();
		pathInfo_Nor2 = new Path_Info();
		pathInfo_Rcm2 = new Path_Info();

		gapInfo_Nor1 = new Gap_Info();
		gapInfo_Rcm1 = new Gap_Info();			
		gapInfo_Nor2 = new Gap_Info();
		gapInfo_Rcm2 = new Gap_Info();		
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

		delete segInfo_Nor1;
		delete segInfo_Nor2;
		delete segInfo_Rcm1;
		delete segInfo_Rcm2;
	}

	/*
	 * Attempt to map the paired-end read to the reference
	 */
	void fixPhase1_segInfo(PairedEndRead* pairedEndRead, Chromosome* chrom)
	{
		// Attempt to map the first read and it's reverse complement
		normalMapMain = segInfo_Nor1->mapMain_SegInfo_preIndex(
			pairedEndRead->getFirstRead(),
			chrom);
		rcmMapMain = segInfo_Rcm1->mapMain_SegInfo_preIndex(
			pairedEndRead->getFirstReadReverseComplement(),
			chrom);

		// Attempt to map the second read and it's reverse complement
		normalMapMain_PE = segInfo_Nor2->mapMain_SegInfo_preIndex(
			pairedEndRead->getSecondRead(),
			chrom);
		rcmMapMain_PE = segInfo_Rcm2->mapMain_SegInfo_preIndex(
			pairedEndRead->getSecondReadReverseComplement(),
			chrom);
	}

	void fixPhase1_pathInfo()
	{
		if(normalMapMain)
			pathInfo_Nor1->getPossiPathFromSeg(segInfo_Nor1);
		if(rcmMapMain)
			pathInfo_Rcm1->getPossiPathFromSeg(segInfo_Rcm1);
		if(normalMapMain_PE)
			pathInfo_Nor2->getPossiPathFromSeg(segInfo_Nor2);
		if(rcmMapMain_PE)
			pathInfo_Rcm2->getPossiPathFromSeg(segInfo_Rcm2);		
	}

	void fixPhase1_gapInfo(PairedEndRead* pairedEndRead, Index_Info* indexInfo)
	{
		gapInfo_Nor1->fixGapInPath(
			pathInfo_Nor1,
			segInfo_Nor1,
			indexInfo,
			pairedEndRead->getFirstRead());

		gapInfo_Rcm1->fixGapInPath(
			pathInfo_Rcm1,
			segInfo_Rcm1,
			indexInfo,
			pairedEndRead->getFirstReadReverseComplement());

		gapInfo_Nor2->fixGapInPath(
			pathInfo_Nor2,
			segInfo_Nor2,
			indexInfo,
			pairedEndRead->getSecondRead());

		gapInfo_Rcm2->fixGapInPath(
			pathInfo_Rcm2,
			segInfo_Rcm2,
			indexInfo,
			pairedEndRead->getSecondReadReverseComplement());
	}
};
