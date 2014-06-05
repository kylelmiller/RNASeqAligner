
#ifndef __PATH_H_INCLUDED__
#define __PATH_H_INCLUDED__

#include <stdlib.h>
#include <stdio.h>

#include "secondLevelReferenceGenome.h"
#include "referenceGenome.h"
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

		// for each potential path between segments
		// this will contain something like (1, 5), (1, 6), (2, 5), (2, 6)
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
			indexInfo->getChromosomeLocation(PathMapPos, &tmpChrNameInt, &tmpChrPosInt);

			int tmpPathFinalMapPos = tmpChrPosInt;
			int tmpPathFinalMapChr = tmpChrNameInt;

			Splice_Info* tmpSpliceInfo = new Splice_Info();
			tmpSpliceInfo->jump_code.clear();

			// So, this gets the start of our mapped segment which is say 25.
			// We then try and map locations 0 - 24
			int tmpUnfixedHeadLength = currentSegment->getLocationInRead();
			string readSubSeqInProcess_head = mappedRead->getRead()->getSequence().substr(0, tmpUnfixedHeadLength);

			bool scoreStringBool_head;
			string chromSubSeqInProcess_head;

			// We allow a single error every 8 locations, this needs to be configurable
			size_t max_mismatch_head = (tmpUnfixedHeadLength)/8 + 1;
			size_t mismatch_bits_head = 0;
			size_t comb_bits_head = 0;

			// If the head goes off the chromosome then we aren't going to map it
			if(tmpPathFinalMapPos - tmpUnfixedHeadLength - 1 < 0)
			{
				scoreStringBool_head = false;
			}
			else
			{
				// Gets the chromosome values for region in front of the current segment being
				// considered
				// The way this is written we may be doing this multiple times for the same segment
				// Why not just match the segment to begin with and say it starts at a different location
				chromSubSeqInProcess_head = (indexInfo->getChromSequence(tmpChrNameInt)).substr(tmpPathFinalMapPos - tmpUnfixedHeadLength - 1, tmpUnfixedHeadLength);

				// you then check the beginning of the read with the reference genome
				scoreStringBool_head = Utilities::stringsWithinMismatchLimit(
					readSubSeqInProcess_head,
					chromSubSeqInProcess_head,
					max_mismatch_head,
					mismatch_bits_head,
					comb_bits_head);
			}

			// if there was a segment at the beginning which wasn't matched
			if(tmpUnfixedHeadLength > 0)
			{
				// and we matched it
				if(scoreStringBool_head)
				{
					// add #M to the jump code
					// this will make that segment mapped
					// again, jump codes should not be here at all
					Jump_Code tmpHeadJumpCode(tmpUnfixedHeadLength, "M");
					tmpSpliceInfo->jump_code.push_back(tmpHeadJumpCode);

					// add the head mismatches to the total mismatches
					mismatchNumToAdd = mismatchNumToAdd + mismatch_bits_head;
				}
				else
				{
					// if we didn't find it then we consider the beginning "soft" clipped
					Jump_Code tmpHeadJumpCode(tmpUnfixedHeadLength, "S");
					tmpSpliceInfo->jump_code.push_back(tmpHeadJumpCode);
				}
			}

			tmpSpliceInfo->appendJumpCode(fixedPathVec[tmpPath].second);

			// So, this uses the jump codes to figure out the last mapped region based on
			// a given location
			int endMappedBaseMapPos = (fixedPathVec[tmpPath].second)->getEndBaseMapPos_jump_code(tmpPathFinalMapPos);

			// Gets the length of the end of the second segment which has not yet been mapped
			int tmpUnfixedTailLength =
				mappedRead->getRead()->getLength() -
				(currentSegment->getLocationInRead() + currentSegment->getLength());

			// Get the portion of the read that hasn't been mapped
			string readSubSeqInProcess_tail
				= mappedRead->getRead()->getSequence().substr(mappedRead->getRead()->getLength() - tmpUnfixedTailLength, tmpUnfixedTailLength);

			bool scoreStringBool_tail;
			string chromSubSeqInProcess_tail;

			size_t max_mismatch_tail = (tmpUnfixedTailLength)/8 + 1;
			size_t mismatch_bits_tail = 0;
			size_t comb_bits_tail = 0;

			// If the tail exits the chromosome
			if(endMappedBaseMapPos + tmpUnfixedTailLength >= indexInfo->getChromLength(tmpChrNameInt))
			{
				// no mapping the tail
				scoreStringBool_tail = false;
			}
			else
			{
				// get the tail piece of the reference genome
				chromSubSeqInProcess_tail = (indexInfo->getChromSequence(tmpChrNameInt)).substr(endMappedBaseMapPos, tmpUnfixedTailLength);

				// you then check the beginning of the read with the reference genome
				scoreStringBool_tail
					= Utilities::stringsWithinMismatchLimit(readSubSeqInProcess_tail, chromSubSeqInProcess_tail,
						max_mismatch_tail, mismatch_bits_tail, comb_bits_tail);
			}

			// if there was a tail to map
			if(tmpUnfixedTailLength > 0)
			{
				// if we did map the tail then add a jump code saying that
				// again, there shouldn't be any jump codes here at all
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

			// combine all of the jump codes
			tmpSpliceInfo->getFinalJumpCode();

			// if we mapped some head stuff then change the final map position
			// to be closer
			if(tmpUnfixedHeadLength > 0 && scoreStringBool_head)
				tmpPathFinalMapPos -= tmpUnfixedHeadLength;

			// add the new mismatches to the path
			fixedPathMismatchVec[tmpPath] += mismatchNumToAdd;

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

	// this is all about paths of mapped regions within a single read
	// FIXME - THIS NEEDS TO GET REWORKED KLM 5/30/14
	void addNewSegCandiToCurrentPathInfoVec_matchIndelUniqueSpliceMultiPath(Segment* segment,
			int segCandiNO, MappedRead* mappedRead, bool longSegBool, int currentPathNum)
	{
		// We have a current segment that we are looking at and we want to add potential paths
		// from the previous segments to this one
		bool relatedToSomePath = false;
		int minSegNumGap = 0;

		// I'm not quite sure what this is doing
		// I thin it's determining whether this segment should be identified as a
		// splice or a deletion and then matching each path that fits that criteria
		// and not making paths with the others
		int minDistance = minDistanceWithCurrentPath(segment, segCandiNO, mappedRead, currentPathNum, &minSegNumGap);
		int compareDistance = minDistance <= MAX_DELETION_LENGTH
			? MAX_DELETION_LENGTH
			: SPLICEDISTANCE;

		// This looks through all of the previous segments and tries
		// to find the one that is closest to the current segment
		// It then determines whether the gap is a deletion or a splice
		for(int tmpPathElementNO = 0; tmpPathElementNO < currentPathNum; tmpPathElementNO++)
		{
			int tmpPathElementVecSize = (PathVec_seg[tmpPathElementNO]).size();
			int tmpPathLastSegGroupNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).first;
			int tmpPathLastSegCandiNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).second;
			int tmpSegDistance = mappedRead->distanceBetweenSegment(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO,
					segment, segCandiNO);
			int tmpSegNumGap = mappedRead->getGroupNumber(segment) - tmpPathLastSegGroupNO;

			// if the distance between the two paths meets the criteria then add the path
			if(tmpSegDistance <= compareDistance && tmpSegNumGap <= minSegNumGap)
			{
				PathValidBoolVec[tmpPathElementNO] = false;
				this->copyOldPath_AddSegCandi_AddNewPath(tmpPathElementNO, mappedRead->getGroupNumber(segment), segCandiNO);
				relatedToSomePath = true;
			}
		}

		// If we created some new paths then add them to the existing paths
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
		Segment* currentSegment = mappedRead->getMostConfidentSegment();

		if(currentSegment == NULL)
			return false;

		/*
		// THIS DOES NOTHING. This information is already in the mapped read
		// object
		for(int i=0; i<currentSegment->getAlignmentNumber(); i++)
		{
			vector< pair<int,int> > tmpPathElementVec;
			tmpPathElementVec.push_back( pair<int,int> (firstLongSegNO, i));
			PathVec_seg.push_back(tmpPathElementVec);
			PathValidBoolVec.push_back(true);
			vector< pair<int,int> > ().swap(tmpPathElementVec);
		}*/

		/*
		 * This loops through the segments longer than 20 after the long segment
		 * and adds them to the path vector
		 * We are handling this already by not adding segments that are shorter than 20
		for(int tmpSegNO = firstLongSegNO + 1; tmpSegNO < segInfo->segmentNum; tmpSegNO ++)
		{
			if(//(!segInfo->checkSegLongOrNot(tmpSegNO)) &&
				(segInfo->norSegmentAlignNum[tmpSegNO]) > CANDALILOC)
				continue;
			this->addNewSegGroupToCurrentPathInfoVec(segInfo, tmpSegNO);
		}
		*/
		return true;
	}

	/* AS THE CODE CURRENTLY IS THIS IS NEVER HIT
	void addNewSegGroupToCurrentPathInfoVec(Seg_Info* segInfo, int segGroupNO)
	{
		bool longSegBool = segInfo->checkSegLongOrNot(segGroupNO);

		// So, for each segment in this group we send it off
		int currentPathNum = PathVec_seg.size();
		for(int tmpSegCandiLoc = 0; tmpSegCandiLoc < segInfo->norSegmentAlignNum[segGroupNO]; tmpSegCandiLoc++)
		{
			this->addNewSegCandiToCurrentPathInfoVec_matchIndelUniqueSpliceMultiPath(segInfo, segGroupNO, tmpSegCandiLoc, longSegBool, currentPathNum);
		}
	}
	*/
};

#endif
