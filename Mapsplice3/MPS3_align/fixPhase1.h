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

	void memoryFree()
	{
			delete(gapInfo_Nor1); delete(gapInfo_Rcm1); delete(gapInfo_Nor2); delete(gapInfo_Rcm2);
			pathInfo_Nor1->memoryFree(); delete(pathInfo_Nor1); pathInfo_Rcm1->memoryFree(); delete(pathInfo_Rcm1);
			pathInfo_Nor2->memoryFree(); delete(pathInfo_Nor2); pathInfo_Rcm2->memoryFree(); delete(pathInfo_Rcm2);
			delete(segInfo_Nor1); delete(segInfo_Nor2); delete(segInfo_Rcm1); delete(segInfo_Rcm2);		
	}

	void fixPhase1_segInfo(char* read, char* read_RC, char* read_PE, char* read_RC_PE, unsigned int* sa, BYTE* lcpCompress, 
		unsigned int* childTab, char* chrom, 
		//unsigned int* valLength, 
		BYTE* verifyChild, //int readLength, 
		Index_Info* indexInfo,
		int* preIndexMapLengthArray, unsigned int* preIndexIntervalStartArray,
		unsigned int* preIndexIntervalEndArray, PE_Read_Info* readInfo)
	{
		//cout << endl << "## do segment mapping for Nor_1 ##" << endl;
		//cout << "start to fix segInfo_Nor1 ..." << endl << endl;
		//bool 
		normalMapMain = segInfo_Nor1->mapMain_SegInfo_preIndex(read, sa, lcpCompress, childTab, 
			chrom, &norValLength, verifyChild, (readInfo->readInfo_pe1).readSeqLength, indexInfo, preIndexMapLengthArray,
			preIndexIntervalStartArray, preIndexIntervalEndArray);
		
		//cout << (segInfo_Nor1)->segInfoStr(indexInfo) << endl << endl;

		//cout << endl << "## do segment mapping for Rcm_1 ##" << endl;
		//cout << "start to fix segInfo_Rcm1 ..." << endl << endl;
	   	//bool 
	   	rcmMapMain = segInfo_Rcm1->mapMain_SegInfo_preIndex(read_RC, sa, lcpCompress, childTab,
			chrom, &rcmValLength, verifyChild, (readInfo->readInfo_pe1).readSeqLength, indexInfo, preIndexMapLengthArray,
			preIndexIntervalStartArray, preIndexIntervalEndArray);		
	   
	   	//cout << (segInfo_Rcm1)->segInfoStr(indexInfo) << endl << endl;

		//cout << endl << "## do segment mapping for Nor_2 ##" << endl;
	   	//cout << "start to fix segInfo_Nor2 ..." << endl << endl;
		//bool 
		normalMapMain_PE = segInfo_Nor2->mapMain_SegInfo_preIndex(read_PE, sa, lcpCompress, childTab, 
			chrom, &norValLength_PE, verifyChild, (readInfo->readInfo_pe2).readSeqLength, indexInfo, preIndexMapLengthArray,
			preIndexIntervalStartArray, preIndexIntervalEndArray);
		
		//cout << (segInfo_Nor2)->segInfoStr(indexInfo) << endl << endl;		
		
		//cout << endl << "## do segment mapping for Rcm_2 ##" << endl;
		//cout << "start to fix segInfo_Rcm2 ..." << endl << endl;
		//bool 
		rcmMapMain_PE = segInfo_Rcm2->mapMain_SegInfo_preIndex(read_RC_PE, sa, lcpCompress, childTab, 
			chrom, &rcmValLength_PE, verifyChild, (readInfo->readInfo_pe2).readSeqLength, indexInfo, preIndexMapLengthArray,
			preIndexIntervalStartArray, preIndexIntervalEndArray);			
		
		//cout << (segInfo_Rcm2)->segInfoStr(indexInfo) << endl << endl;

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

	void fixPhase1_gapInfo(PE_Read_Info* readInfo, Index_Info* indexInfo)
	{
		//cout << "start to fix gapInfo_Nor1 ..." << endl;
		gapInfo_Nor1->fixGapInPath(pathInfo_Nor1, segInfo_Nor1, 
			indexInfo, (readInfo->readInfo_pe1).readSeq, (readInfo->readInfo_pe1).readSeqLength);
		//cout << "start to fix gapInfo_Rcm1 ..." << endl;
		gapInfo_Rcm1->fixGapInPath(pathInfo_Rcm1, segInfo_Rcm1, 
			indexInfo, (readInfo->readInfo_pe1).rcmReadSeq, (readInfo->readInfo_pe1).readSeqLength);
		//cout << "start to fix gapInfo_Nor2 ..." << endl;	
		gapInfo_Nor2->fixGapInPath(pathInfo_Nor2, segInfo_Nor2, 
			indexInfo, (readInfo->readInfo_pe2).readSeq, (readInfo->readInfo_pe2).readSeqLength);
		//cout << "start to fix gapInfo_Rcm2 ..." << endl;
		gapInfo_Rcm2->fixGapInPath(pathInfo_Rcm2, segInfo_Rcm2, 
			indexInfo, (readInfo->readInfo_pe2).rcmReadSeq, (readInfo->readInfo_pe2).readSeqLength);			
	}

	void coutDebugInfo(PE_Read_Info* readInfo, Index_Info* indexInfo)
	{
				cout << endl << "##### readName_1: " << (readInfo->readInfo_pe1).readName << " #####" << endl;
				cout << "##### readName_2: " << (readInfo->readInfo_pe2).readName << " #####"<< endl;

				cout << endl << "## do segment mapping for Nor_1 ##" << endl;
				cout << segInfo_Nor1->segInfoStr(indexInfo) << endl;
				cout << pathInfo_Nor1->possiPathStr() << endl;
				cout << pathInfo_Nor1->fixedPathVecStr(indexInfo, segInfo_Nor1) << endl;
				cout << pathInfo_Nor1->finalFixedPathStr(indexInfo) << endl;

		    	cout << endl << "## do segment mapping for Rcm_1 ##" << endl;
		    	cout << segInfo_Rcm1->segInfoStr(indexInfo) << endl;
				cout << pathInfo_Rcm1->possiPathStr() << endl;
				cout << pathInfo_Rcm1->fixedPathVecStr(indexInfo, segInfo_Rcm1) << endl;
				cout << pathInfo_Rcm1->finalFixedPathStr(indexInfo) << endl;

				cout << endl << "## do segment mapping for Nor_2 ##" << endl;
				cout << segInfo_Nor2->segInfoStr(indexInfo) << endl;
				cout << pathInfo_Nor2->possiPathStr() << endl;
				cout << pathInfo_Nor2->fixedPathVecStr(indexInfo, segInfo_Nor2) << endl;
				cout << pathInfo_Nor2->finalFixedPathStr(indexInfo) << endl;

				cout << endl << "## do segment mapping for Rcm_2 ##" << endl;
				cout << segInfo_Rcm2->segInfoStr(indexInfo) << endl;
				cout << pathInfo_Rcm2->possiPathStr() << endl;
				cout << pathInfo_Rcm2->fixedPathVecStr(indexInfo, segInfo_Rcm2) << endl;
				cout << pathInfo_Rcm2->finalFixedPathStr(indexInfo) << endl;
	}
	/*void fixPhase1()
	{
			gapInfo_Nor1->fixGapInPath(pathInfo_Nor1, segInfo_Nor1, 
				indexInfo, (readInfo->readInfo_pe1).readSeq, (readInfo->readInfo_pe1).readSeqLength);

			gapInfo_Rcm1->fixGapInPath(pathInfo_Rcm1, segInfo_Rcm1, 
				indexInfo, (readInfo->readInfo_pe1).rcmReadSeq, (readInfo->readInfo_pe1).readSeqLength);
			
			gapInfo_Nor2->fixGapInPath(pathInfo_Nor2, segInfo_Nor2, 
				indexInfo, (readInfo->readInfo_pe2).readSeq, (readInfo->readInfo_pe2).readSeqLength);

			gapInfo_Rcm2->fixGapInPath(pathInfo_Rcm2, segInfo_Rcm2, 
				indexInfo, (readInfo->readInfo_pe2).rcmReadSeq, (readInfo->readInfo_pe2).readSeqLength);	
	}*/

};