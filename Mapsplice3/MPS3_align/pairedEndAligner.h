#ifndef __PAIRED_END_ALIGNER_H_INCLUDED__
#define __PAIRED_END_ALIGNER_H_INCLUDED__

#include <string>
#include <string.h>

#include "pairedEndRead.h"
#include "mappedRead.h"
#include "path.h"
#include "gap.h"
#include "referenceGenome.h"

using namespace std;

class PairedEndAligner
{
private:

	/*
	 * Member Variables
	 */
	MappedRead* _firstNormalMappedRead;
	MappedRead* _firstReverseCompMappedRead;
	MappedRead* _secondNormalMappedRead;
	MappedRead* _secondReverseCompMappedRead;

	ReferenceGenome* _reference;

	Gap* gapInfo_Nor1;
	Gap* gapInfo_Rcm1;
	Gap* gapInfo_Nor2;
	Gap* gapInfo_Rcm2;

public:

	Path* _pathInfo_Nor1;
	Path* pathInfo_Rcm1;
	Path* pathInfo_Nor2;
	Path* pathInfo_Rcm2;

	/*
	 * Constructor
	 */
	PairedEndAligner(PairedEndRead* pairedEndRead, ReferenceGenome* reference)
	{
		_reference = reference;

	   	_firstNormalMappedRead = new MappedRead(
			pairedEndRead->getFirstRead(),
			reference);
	    _firstReverseCompMappedRead = new MappedRead(
				pairedEndRead->getFirstReadReverseComplement(),
				reference);

		_secondNormalMappedRead = new MappedRead(
				pairedEndRead->getSecondRead(),
				reference);
		_secondReverseCompMappedRead = new MappedRead(
				pairedEndRead->getSecondReadReverseComplement(),
				reference);

		_pathInfo_Nor1 = new Path();
		pathInfo_Rcm1 = new Path();
		pathInfo_Nor2 = new Path();
		pathInfo_Rcm2 = new Path();

		gapInfo_Nor1 = new Gap();
		gapInfo_Rcm1 = new Gap();			
		gapInfo_Nor2 = new Gap();
		gapInfo_Rcm2 = new Gap();		
	}

	/*
	 * Deconstructor
	 */
	~PairedEndAligner()
	{
		delete gapInfo_Nor1;
		delete gapInfo_Rcm1;
		delete gapInfo_Nor2;
		delete gapInfo_Rcm2;

		delete _pathInfo_Nor1;
		delete pathInfo_Rcm1;
		delete pathInfo_Nor2;
		delete pathInfo_Rcm2;

		delete _firstNormalMappedRead;
		delete _secondNormalMappedRead;
		delete _firstReverseCompMappedRead;
		delete _secondReverseCompMappedRead;
	}



	void fixPhase1_pathInfo()
	{
		_pathInfo_Nor1->getPossiPathFromSeg(_firstNormalMappedRead);
		pathInfo_Rcm1->getPossiPathFromSeg(_firstReverseCompMappedRead);

		pathInfo_Nor2->getPossiPathFromSeg(_secondNormalMappedRead);
		pathInfo_Rcm2->getPossiPathFromSeg(_secondReverseCompMappedRead);
	}

	void fixPhase1_gapInfo()
	{
		gapInfo_Nor1->fixGapInPath(
			_pathInfo_Nor1,
			_firstNormalMappedRead,
			_reference->getIndexInfo());

		gapInfo_Rcm1->fixGapInPath(
			pathInfo_Rcm1,
			_firstReverseCompMappedRead,
			_reference->getIndexInfo());

		gapInfo_Nor2->fixGapInPath(
			pathInfo_Nor2,
			_secondNormalMappedRead,
			_reference->getIndexInfo());

		gapInfo_Rcm2->fixGapInPath(
			pathInfo_Rcm2,
			_secondReverseCompMappedRead,
			_reference->getIndexInfo());
	}
};

#endif
