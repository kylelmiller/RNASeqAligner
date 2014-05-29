#include <string>
#include <string.h>
#include "splice_info.h"
#include "secondLevelChromosome.h"

using namespace std;

class FixOneEndUnmappedInfo
{
public:
	FixOneEndUnmappedInfo()
	{}

	void fixOneEndReadForward(Read read, Read incompleteEndRead, PE_Read_Alignment_Info* peAlignInfo,
		vector<Alignment_Info*> alignmentVector, SecondLevelChromosomeList* secondLevelChromosomeList)
	{
		string chrNameStr;
		int chrMapPos_start;
		int chrMapPos_end;

		int mapPosIntervalStart;
		int mapPosIntervalEnd;
		int secondLevelIndexNum; // start from 1

		int chrPosStartIn2ndLevelIndex;

		for(int i = 0; i < alignmentVector.size();i++)
		{
			Alignment_Info* alignmentInfo = alignmentVector[i];

			mapPosIntervalStart = alignmentInfo->alignChromPos;
			mapPosIntervalEnd = alignmentInfo->getEndMatchedPosInChr() + READ_ALIGN_AREA_LENGTH;

			SecondLevelChromosome* secondLevelChromosome =
				secondLevelChromosomeList->getSecondLevelChromosome(
						alignmentInfo->alignChromName,
						alignmentInfo->alignChromPos);

			if(secondLevelChromosome == NULL)
				continue;

			Seg2ndOri_Info* seg2ndOriInfo = new Seg2ndOri_Info();
			if(!seg2ndOriInfo->mapMainSecondLevel_compressedIndex(read, secondLevelChromosome))
			{
				delete seg2ndOriInfo;
				continue;
			}

			Seg_Info* segInfo = new Seg_Info(seg2ndOriInfo, mapPosIntervalStart,
				mapPosIntervalEnd, chrPosStartIn2ndLevelIndex,
				secondLevelChromosome->getIndexInfo(), chrNameStr);

			Path_Info* pathInfo = new Path_Info();
			pathInfo->getPossiPathFromSeg(segInfo);

			int pathValidNum = pathInfo->pathValidNumInt();
			if(pathValidNum > 10)
			{
				delete(pathInfo);
				delete(segInfo);
				delete(seg2ndOriInfo);
				continue;
			}

			Gap_Info* gapInfo = new Gap_Info();
			gapInfo->fixGapInPath(pathInfo, segInfo,
					secondLevelChromosome->getIndexInfo(), incompleteEndRead);


			/* FIX ME - FIX THIS LATER 5/28/14 KLM
			if(pathInfo->finalPathVec.size() == 1)
				peAlignInfo->pushBackPathInfo2PeAlignInfo(pathInfo, End1OrEnd2, NorOrRcm, indexInfo);

			peAlignInfo->pushBackPathInfo2PeAlignInfo(pathInfo, End1OrEnd2, NorOrRcm, indexInfo);
			 */
			delete(gapInfo);
			delete(pathInfo);
			delete(segInfo);
			delete(seg2ndOriInfo);
		}
	}

	void fixOneEndUnmappedReadBackwards(Read read, Read incompleteEndRead,
		PE_Read_Alignment_Info* peAlignInfo, vector<Alignment_Info*> alignmentVector,
		SecondLevelChromosomeList* secondLevelChromosomeList)
	{

	}

	void fixOneEndUnmapped(PairedEndRead* peReadInfo, PE_Read_Alignment_Info* peAlignInfo,
			SecondLevelChromosomeList* secondLevelChromosomeList)
	{
		// Fix the first read's ends
		fixOneEndReadForward(peReadInfo->getFirstRead(), peReadInfo->getSecondReadReverseComplement(),
			peAlignInfo, peAlignInfo->norAlignmentInfo_PE_1, secondLevelChromosomeList); // End1OrEnd2=true, NorOrRcm=true
		fixOneEndUnmappedReadBackwards(peReadInfo->getFirstReadReverseComplement(), peReadInfo->getSecondRead(),
			peAlignInfo, peAlignInfo->rcmAlignmentInfo_PE_1, secondLevelChromosomeList); // End1OrEnd2=true, NorOrRcm=false

		// Fix the second read's ends
		fixOneEndUnmappedReadBackwards(peReadInfo->getSecondRead(), peReadInfo->getFirstReadReverseComplement(),
			peAlignInfo, peAlignInfo->norAlignmentInfo_PE_2, secondLevelChromosomeList); // End1OrEnd2=false, NorOrRcm=true
		fixOneEndReadForward(peReadInfo->getSecondReadReverseComplement(), peReadInfo->getFirstRead(),
			peAlignInfo, peAlignInfo->rcmAlignmentInfo_PE_2, secondLevelChromosomeList); // End1OrEnd2=false, NorOrRcm=false
	}
};
