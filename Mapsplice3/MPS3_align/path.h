
#ifndef __PATH_H_INCLUDED__
#define __PATH_H_INCLUDED__

#include <stdlib.h>
#include <stdio.h>

#include "secondLevelChromosome.h"
#include "chromosome.h"
#include "segment.h"
#include "splice_info.h"
#include "mappedRead.h"

using namespace std;

class Path
{
public:

	vector< vector< pair<int, int> > > PathVec_seg;
	vector< bool > PathValidBoolVec;

	vector < pair<int, pair<int, int> > > validPathVec_toPair;
	vector <int> validPathVec;

	vector< pair<int, Splice_Info*> > fixedPathVec;
	vector< int > fixedPathMismatchVec;

	vector< bool > PathFixedBoolVec;

	vector < pair< pair<int, int>, Splice_Info*> > finalPathVec;

	Path()
	{}

	~Path()
	{
		for(int tmp = 0; tmp < fixedPathVec.size(); tmp++)
			delete fixedPathVec[tmp].second;

		for(int tmp = 0; tmp < finalPathVec.size(); tmp++)
			delete finalPathVec[tmp].second;
	}

	int pathValidNumInt()
	{
		int pathValidNum = 0;
		for(int tmp = 0; tmp < PathValidBoolVec.size(); tmp++)
			if(PathValidBoolVec[tmp])
				pathValidNum++;

		return pathValidNum;
	}

	// FIXME - THIS METHOD NEEDS TO BE CLEANED UP KLM 5/29/14
	void getFinalPath_extend2HeadTail(Index_Info* indexInfo, MappedRead* mappedRead)
	{
		if(mappedRead->getNumberOfSegments() < 1)
			return;

		for(int tmpPath = 0; tmpPath < fixedPathVec.size(); tmpPath++)
		{
			int mismatchNumToAdd = 0;
			int fixedPathNO = fixedPathVec[tmpPath].first;

			int tmpPath1stSegGroupNO = (PathVec_seg[fixedPathNO])[0].first;
			int tmpPath1stSegCandiNO = (PathVec_seg[fixedPathNO])[0].second;

			Segment* currentSegment = mappedRead->getSegment(tmpPath1stSegGroupNO);

			int tmpPathElementSize = (PathVec_seg[fixedPathNO]).size();
			int tmpPathLastSegGroupNO = (PathVec_seg[fixedPathNO])[tmpPathElementSize - 1].first;

			unsigned int PathMapPos = currentSegment->getAlignmentLocation(tmpPath1stSegCandiNO);
			unsigned int tmpChrNameInt, tmpChrPosInt;
			indexInfo->getChrLocation(PathMapPos, &tmpChrNameInt, &tmpChrPosInt);

			int tmpPathFinalMapPos = tmpChrPosInt;
			int tmpPathFinalMapChr = tmpChrNameInt;

			Splice_Info* tmpSpliceInfo = new Splice_Info();
			tmpSpliceInfo->jump_code.clear();

			int tmpUnfixedHeadLength = currentSegment->getLocationInRead();
			string readSubSeqInProcess_head = mappedRead->getRead()->getSequence().substr(0, tmpUnfixedHeadLength);

			bool scoreStringBool_head;
			string chromSubSeqInProcess_head;

			size_t max_mismatch_head = (tmpUnfixedHeadLength)/8 + 1;
			size_t mismatch_bits_head = 0;
			size_t comb_bits_head = 0;

			if(tmpPathFinalMapPos - tmpUnfixedHeadLength - 1 < 0)
			{
				scoreStringBool_head = false;
			}
			else
			{
				chromSubSeqInProcess_head = (indexInfo->chromStr[tmpChrNameInt]).substr(tmpPathFinalMapPos - tmpUnfixedHeadLength - 1, tmpUnfixedHeadLength);

				// So, this is important
				// it gives a score to the alignment, i think
				scoreStringBool_head = Utilities::scoreString(readSubSeqInProcess_head, chromSubSeqInProcess_head, max_mismatch_head, mismatch_bits_head, comb_bits_head);
			}

			//cout << "scoreStringBool_head " << scoreStringBool_head << endl;

			if(tmpUnfixedHeadLength > 0)
			{
				if(scoreStringBool_head)
				{
					Jump_Code tmpHeadJumpCode(tmpUnfixedHeadLength, "M");
					tmpSpliceInfo->jump_code.push_back(tmpHeadJumpCode);
					mismatchNumToAdd = mismatchNumToAdd + mismatch_bits_head;
				}
				else
				{
					Jump_Code tmpHeadJumpCode(tmpUnfixedHeadLength, "S");
					tmpSpliceInfo->jump_code.push_back(tmpHeadJumpCode);
				}
			}

			tmpSpliceInfo->appendJumpCode(fixedPathVec[tmpPath].second);

			//////////////////////////// add last jump code /////////////////////////////////////////////////////
			int endMappedBaseMapPos = (fixedPathVec[tmpPath].second)->getEndBaseMapPos_jump_code(tmpPathFinalMapPos);
			int tmpUnfixedTailLength =
				mappedRead->getRead()->getLength() -
				(currentSegment->getLocationInRead() + currentSegment->getLength());

			string readSubSeqInProcess_tail
				= mappedRead->getRead()->getSequence().substr(mappedRead->getRead()->getLength() - tmpUnfixedTailLength, tmpUnfixedTailLength);

			bool scoreStringBool_tail;
			string chromSubSeqInProcess_tail;

			size_t max_mismatch_tail = (tmpUnfixedTailLength)/8 + 1;
			size_t mismatch_bits_tail = 0;
			size_t comb_bits_tail = 0;

			if(endMappedBaseMapPos + tmpUnfixedTailLength >= indexInfo->chromLength[tmpChrNameInt])
			{
				scoreStringBool_tail = false;
			}
			else
			{
				chromSubSeqInProcess_tail = (indexInfo->chromStr[tmpChrNameInt]).substr(endMappedBaseMapPos, tmpUnfixedTailLength);

				// this is important
				// it gives a score to the alignment, i think
				scoreStringBool_tail
					= Utilities::scoreString(readSubSeqInProcess_tail, chromSubSeqInProcess_tail,
						max_mismatch_tail, mismatch_bits_tail, comb_bits_tail);
			}

			if(tmpUnfixedTailLength > 0)
			{
				if(scoreStringBool_tail)
				{
					Jump_Code tmpTailJumpCode(tmpUnfixedTailLength, "M");
					tmpSpliceInfo->jump_code.push_back(tmpTailJumpCode);
					mismatchNumToAdd = mismatchNumToAdd + mismatch_bits_tail;
				}
				else
				{
					Jump_Code tmpTailJumpCode(tmpUnfixedTailLength, "S");
					tmpSpliceInfo->jump_code.push_back(tmpTailJumpCode);
				}
			}

			tmpSpliceInfo->getFinalJumpCode();

			if(tmpUnfixedHeadLength > 0 && scoreStringBool_head)
				tmpPathFinalMapPos -= tmpUnfixedHeadLength;

			int oldMismatchNum = fixedPathMismatchVec[tmpPath];
			int newMismatchNum = oldMismatchNum + mismatchNumToAdd;


			fixedPathMismatchVec[tmpPath] = newMismatchNum;

			finalPathVec.push_back(pair< pair<int, int>, Splice_Info*> (pair<int,int> (tmpPathFinalMapChr, tmpPathFinalMapPos), tmpSpliceInfo) );
		}
	}

	void copyOldPath_AddSegCandi_AddNewPath(int pathElementNO, int segGroupNO, int segCandiNO)
	{
		//copy path;
		vector< pair<int,int> > newPath;
		for(int tmpElementNO = 0; tmpElementNO < PathVec_seg[pathElementNO].size(); tmpElementNO++)
		{
			newPath.push_back((PathVec_seg[pathElementNO])[tmpElementNO]);
		}
		newPath.push_back(pair<int, int> (segGroupNO, segCandiNO));
		PathVec_seg.push_back(newPath);
		PathValidBoolVec.push_back(true);
		vector< pair<int,int> >().swap(newPath);
	}

	// FIXME - THIS NEEDS TO GET REWORKED KLM 5/30/14
	int minDistanceWithCurrentPath(Segment* segment, int segCandiNO, MappedRead* mappedRead,  int currentPathNum, int* minSegNumGap)
	{
		int minDistance_value = 900000;
		int minDistance_path = 1000;
		int tmpMinSegNumGap = 100;
		for(int tmpPathElementNO = 0; tmpPathElementNO < currentPathNum; tmpPathElementNO++)
		{
			int tmpPathElementVecSize = (PathVec_seg[tmpPathElementNO]).size();
			int tmpPathLastSegGroupNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).first;
			int tmpPathLastSegCandiNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).second;
			int tmpSegDistance = mappedRead->distanceBetweenSegment(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO,
					segment, segCandiNO);
			int tmpSegNumGap = mappedRead->getGroupNumber(segment) - tmpPathLastSegGroupNO;

			if(tmpSegNumGap <= 0)
				continue;

			if(tmpSegDistance < SPLICEDISTANCE)
			{
				if(tmpSegNumGap < tmpMinSegNumGap)
				{
					tmpMinSegNumGap = tmpSegNumGap;
					minDistance_value = tmpSegDistance;
					minDistance_path = tmpPathElementNO;
				}
				else if(tmpSegNumGap == tmpMinSegNumGap && tmpSegDistance < minDistance_value)
				{
					minDistance_value = tmpSegDistance;
					minDistance_path = tmpPathElementNO;
				}
			}

		}
		*minSegNumGap = tmpMinSegNumGap;
		return minDistance_value;
	}

	// FIXME - THIS NEEDS TO GET REWORKED KLM 5/30/14
	void addNewSegCandiToCurrentPathInfoVec_matchIndelUniqueSpliceMultiPath(Segment* segment,
			int segCandiNO, MappedRead* mappedRead, bool longSegBool, int currentPathNum)
	{
		bool relatedToSomePath = false;
		int minSegNumGap = 0;

		int minDistance = minDistanceWithCurrentPath(segment, segCandiNO, mappedRead, currentPathNum, &minSegNumGap);
		int compareDistance = minDistance <= MAX_DELETION_LENGTH
			? MAX_DELETION_LENGTH
			: SPLICEDISTANCE;

		for(int tmpPathElementNO = 0; tmpPathElementNO < currentPathNum; tmpPathElementNO++)
		{
			int tmpPathElementVecSize = (PathVec_seg[tmpPathElementNO]).size();
			int tmpPathLastSegGroupNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).first;
			int tmpPathLastSegCandiNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).second;
			int tmpSegDistance = mappedRead->distanceBetweenSegment(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO,
					segment, segCandiNO);
			int tmpSegNumGap = mappedRead->getGroupNumber(segment) - tmpPathLastSegGroupNO;
			if(tmpSegDistance <= compareDistance && tmpSegNumGap <= minSegNumGap)
			{
				PathValidBoolVec[tmpPathElementNO] = false;
				this->copyOldPath_AddSegCandi_AddNewPath(tmpPathElementNO, mappedRead->getGroupNumber(segment), segCandiNO);
				relatedToSomePath = true;
			}
		}

		if(!relatedToSomePath && longSegBool)
		{
			vector< pair<int,int> > tmpPathElementVec;
			tmpPathElementVec.push_back( pair<int,int> (mappedRead->getGroupNumber(segment), segCandiNO) );
			PathVec_seg.push_back(tmpPathElementVec);
			PathValidBoolVec.push_back(true);
			vector< pair<int,int> >().swap(tmpPathElementVec);
		}
	}

	// FIXME - KLM 6/2/14 WHAT IS THE POINT OF THIS?
	bool getPossiPathFromSeg(MappedRead* mappedRead)
	{
		int firstLongSegNO = mappedRead->getFirstLongSegNO();

		if(firstLongSegNO < 0)
			return false;

		Segment* currentSegment = mappedRead->getSegment(firstLongSegNO);

		for(int i=0; i<currentSegment->getAlignmentNumber(); i++)
		{
			vector< pair<int,int> > tmpPathElementVec;
			tmpPathElementVec.push_back( pair<int,int> (firstLongSegNO, i));
			PathVec_seg.push_back(tmpPathElementVec);
			PathValidBoolVec.push_back(true);
			vector< pair<int,int> > ().swap(tmpPathElementVec);
		}

		return true;
	}
};

#endif
