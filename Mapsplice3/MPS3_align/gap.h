#ifndef __GAP_H_INCLUDED__
#define __GAP_H_INCLUDED__

/* 
To fix:
	1. get # of mismatches: score_string, ...
	2. in fix_insertion: if secondMatchLength == 0, then .......
	3. extend back 2nd segment to avoid pendingSeq.length() > 31;
*/

#include <stdlib.h>
#include <stdio.h>

#include "genomeScan.h"
#include "path.h"
#include "mappedRead.h"
#include "index_info.h"
#include "read.h"

using namespace std;

class Gap
{
public:
	int LengthOfSeqPerMismatchAllowed;
	int FixSpliceBuffer;

	Gap()
	{
		LengthOfSeqPerMismatchAllowed = 10;
		FixSpliceBuffer = 4;
	}

	/*
	 * FIXME - KLM 6/4/14 removing path
	 */
	bool fixGapInPath(MappedRead* mappedRead, Index_Info* indexInfo)
	{
		Segment* longSegment = mappedRead->getMostConfidentSegment();
		if(longSegment == NULL)
		{
			//pathInfo->getFinalPath_extend2HeadTail(indexInfo, segInfo, readSeq_inProcess);
			return true;
		}

		for(int i=0; i<longSegment->getAlignmentNumber();i++)
		{
			Splice_Info* newPathSpliceInfo = new Splice_Info();
			Jump_Code firstJumpCode(longSegment->getLength(), "M");
			newPathSpliceInfo->jump_code.push_back(firstJumpCode);

			// this code is never hit, need to figure out what is going
			// on with it
			//for(int j=0; j<longSegment->)

			newPathSpliceInfo->getFinalJumpCode();
			bool allJumpCodeValidBool = newPathSpliceInfo->allFinalJumpCodeValid();
			if(allJumpCodeValidBool)
			{
				//(pathInfo->PathFixedBoolVec).push_back(allJumpCodeValidBool);
				//(pathInfo->fixedPathVec).push_back(pair <int, Splice_Info*> (tmpPathNO, newPathSpliceInfo) );
				//(pathInfo->fixedPathMismatchVec).push_back(newPathMismatchNum);
			}
			else
			{
				//(pathInfo->PathFixedBoolVec).push_back(allJumpCodeValidBool);
			}
		}

		//pathInfo->getFinalPath_extend2HeadTail(indexInfo, segInfo, readSeq_inProcess);
		return true;
	}

	// at this point we have segments in the mapped read
	// and we have potential path information in pathInfo
	// which contains whether we think they are splices or deletions
	// we haven't done any mismatch comparisons or paired end read logic
	// we just have a list of vectors of number pairs which forms a bunch of
	// line segments, essentially
	bool fixGapInPath(Path* pathInfo, MappedRead* mappedRead, Index_Info* indexInfo)
	{
		// this is the list of segment group paths so [seg1group, seg2group, seg3group]
		// the group contains integers which correspond to indexes in other groups
		int pathVecSize = pathInfo->PathVec_seg.size();

		// for each group of path segments found (each group may contain many alignments)
		for(int tmpPathNO = 0; tmpPathNO < pathVecSize; tmpPathNO ++)
		{
			// There is also a vector of bools which states whether there
			// are paths from this segment group
			if(!(pathInfo->PathValidBoolVec[tmpPathNO]))
			{
				pathInfo->PathFixedBoolVec.push_back(false);				
				continue;
			}

			// The first segment group is going to be the first long segment from
			// earlier, so everything builds off of this and this will always be
			// present in your
			int firstSegGroupNO = (pathInfo->PathVec_seg[tmpPathNO])[0].first;
			int firstSegLength = mappedRead->getSegment(firstSegGroupNO)->getLength();

			// So, you have an alignment you call it "M". So, a length
			// of read which is aligned is called #M. If the read is length
			// 50 and the entire read is aligned the cigar string is 50M
			// "JumpCode" stuff should not be here. this is specific to sam/bam files
			Splice_Info* newPathSpliceInfo = new Splice_Info();
			Jump_Code firstJumpCode(firstSegLength, "M");
			newPathSpliceInfo->jump_code.push_back(firstJumpCode);
			// The weird thing is that each mapped region inside a group is the same length
			// even if there are perfectly mapped areas right next to it you
			// treat each segment inside a group as the same even if one has a mapped area of 20
			// and another is perfectly mapped for 70

			int newPathMismatchNum = 0;
			bool newPathFixed = true;
			int tmpPathSegSize = pathInfo->PathVec_seg[tmpPathNO].size();

			// Now for each segment inside of the group
			for(int tmpPathSegNO = 0; tmpPathSegNO < tmpPathSegSize-1; tmpPathSegNO++)
			{
				int tmpMismatchNum = 0;
				// I think this is the group of segments that were mapped
				int tmpSegGroupNO = (pathInfo->PathVec_seg[tmpPathNO])[tmpPathSegNO].first;
				int tmpSegCandiNO = (pathInfo->PathVec_seg[tmpPathNO])[tmpPathSegNO].second;

				// the next potential segment there is a path to
				int tmpSegGroupNO_next = (pathInfo->PathVec_seg[tmpPathNO])[tmpPathSegNO+1].first;
				int tmpSegCandiNO_next = (pathInfo->PathVec_seg[tmpPathNO])[tmpPathSegNO+1].second;

				// This checks to make sure that the segment groups are next to each other
				// then it checks to see if the segments line up perfectly, have an insertion or
				// have a deletion. We aren't going to have the perfectly matching because we
				// expand our segments as large as they can go immediately
				// It also determines if the segments are too far apart to be related
				// That should never happen if you are highly confident about the mappings
				// The algorithm really breaks down here. If a first segment you map is
				// length 20 but it shows up 20 times are wasting a lot of time using that as your
				// anchor. You may have a uniquely mapped 2nd segment which is of length 70 that
				// you should be using
				//int tmpRelation = mappedRead->checkSegRelation(tmpSegGroupNO, tmpSegCandiNO, tmpSegGroupNO_next, tmpSegCandiNO_next);
				int tmpRelation = 1;

				Segment* firstSegment = mappedRead->getSegment(tmpSegGroupNO);
				Segment* secondSegment = mappedRead->getSegment(tmpSegGroupNO_next);

				int tmpSegmentLocInRead_1 = firstSegment->getLocationInRead();
				int tmpSegmentLocInRead_2 = secondSegment->getLocationInRead();
				int tmpSegmentLength_1 = firstSegment->getLength();
				int tmpSegmentLength_2 = secondSegment->getLength();

				// ok, so this just gets the specific alignment within the alignment group
				// so there may be many mappings like 554, 789 and 1102 for example
				// and it gets the second segments group with may have more
				// alignments just like the first
				unsigned int tmpSegmentMapPosInWholeGenome_1 = firstSegment->getAlignmentLocation(tmpSegCandiNO);
				unsigned int tmpSegmentMapPosInWholeGenome_2 = secondSegment->getAlignmentLocation(tmpSegCandiNO_next);

				// so, this is also strange
				// it returns the indexes of the chromosome that the segments are mapped to
				// if they are in different chromosomes then there is no path
				// why not just ask index info if they are in the same one and work off
				// of a returned bool? No idea.
				unsigned int tmpChrNameInt, tmpChrPosInt;
				indexInfo->getChromosomeLocation(tmpSegmentMapPosInWholeGenome_1, &tmpChrNameInt, &tmpChrPosInt);
				string tmpChrNameStr_1 = indexInfo->getChromName(tmpChrNameInt);
				int tmpSegmentMapPos_1 = tmpChrPosInt;
				indexInfo->getChromosomeLocation(tmpSegmentMapPosInWholeGenome_2, &tmpChrNameInt, &tmpChrPosInt);
				string tmpChrNameStr_2 = indexInfo->getChromName(tmpChrNameInt);
				int tmpSegmentMapPos_2 = tmpChrPosInt;
				
				// If they are in the same chromosome then we set the name of the
				// chromosome
				string tmpChrNameStr;
				if(tmpChrNameStr_1 == tmpChrNameStr_2)
				{
					tmpChrNameStr = tmpChrNameStr_1;
				}
				else
				{
					// So, report this match as not possible and try and next one
					newPathFixed = false;
					pathInfo->PathFixedBoolVec.push_back(newPathFixed);
					break;
				}

				// So, now we have two potential segments which are within the same
				// chromosome
				// We are passing spliceInfo which is essentially the first segments
				// length + "M".
				// Relation which may even be invalid at this point. It also tells you
				// whether there is an insertion or deletion
				// The rest of this is basically segment and read information
				// it also tells you the number of mismatches if you were to put these
				// two segments together
				// And, for some reason the chromosome name is passed
				bool tmpDoubleAnchorFixed = fixDoubleAnchor_extendBack(
					newPathSpliceInfo,
					tmpRelation,
					tmpSegmentLocInRead_1,
					tmpSegmentLocInRead_2,
					tmpSegmentLength_1,
					tmpSegmentLength_2,
					tmpSegmentMapPos_1,
					tmpSegmentMapPos_2,
					mappedRead->getRead()->getSequence(),
					indexInfo,
					tmpChrNameStr,
					&tmpMismatchNum);

				// if this didn't work then put the path back in
				// the path list?
				if(!tmpDoubleAnchorFixed)
				{
					newPathFixed = false;
					pathInfo->PathFixedBoolVec.push_back(newPathFixed);
					break;					
				}

				// add the mismatch to the mismatch count
				newPathMismatchNum += tmpMismatchNum;
			}

			// if we get through all this with a working path
			if(newPathFixed)
			{
				// store the jump code
				newPathSpliceInfo->getFinalJumpCode();
				if(newPathSpliceInfo->allFinalJumpCodeValid())
				{
					(pathInfo->PathFixedBoolVec).push_back(true);
					(pathInfo->fixedPathVec).push_back(pair <int, Splice_Info*> (tmpPathNO, newPathSpliceInfo));
					(pathInfo->fixedPathMismatchVec).push_back(newPathMismatchNum);
				}
				else
				{
					(pathInfo->PathFixedBoolVec).push_back(false);
				}
			}
		}

		pathInfo->getFinalPath_extend2HeadTail(indexInfo, mappedRead);

		return true;
	}

	int extendBackInChromSeq(int readLoc, const string& readSeq, int chromLoc, const string& chromSeq, int extendBackLengthMax)
	{
		for (int i=0; i<extendBackLengthMax; i++)
			if(readSeq.at(readLoc - i) != chromSeq.at(chromLoc - i))
				return i;

		return extendBackLengthMax;
	}

	/*
	 * We are attempting to put two segments in the same read together
	 * We have the current cigar (not required), the two segments with information,
	 * the relation (insertion, deletion, neither or invalid), the entire read string,
	 * and the name of the chromosome and the index
	 */
	bool fixDoubleAnchor_extendBack(Splice_Info* cigarInfo, int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, Index_Info* indexInfo, const string& chromName, int* mismatchNum)
	{
		bool fixDoubleAnchorBool = false;

		int chrNameInt = indexInfo->convertStringToInt(chromName);
		// This determines the gap between the end of the first segment and the beginning of the
		// second segment
		// We have done this already with determining the relation stuff so I don't see the
		// point in doing it twice
		int extendBackNumMax = segmentLocInRead_2 - (segmentLocInRead_1 + segmentLength_1);

		// Why does the map position have anything to do with the gap between segment
		// one and segment two? The map position is the position in the chromosome
		if(extendBackNumMax > segmentMapPos_2 - 1)
		{
			extendBackNumMax = segmentMapPos_2 - 1; 
		}

		// Ok, so this extends the current segment back as far as it can while perfectly matching
		// We are currently doing this already in the mapping code
		// the only difference is that it won't extend back over the first segment
		int extendBackNum = extendBackInChromSeq(segmentLocInRead_2, readSeq_inProcess, 
			segmentMapPos_2, indexInfo->getChromSequence(chrNameInt), extendBackNumMax);

		// Now we modify the second segments length, location in the read
		// and location in the chromosome
		segmentLocInRead_2 = segmentLocInRead_2 - extendBackNum;
		segmentLength_2 = segmentLength_2 + extendBackNum;
		segmentMapPos_2 = segmentMapPos_2 - extendBackNum;

		// so, we all of the work above and then toss if based on
		// what relation is. That doesn't make sense either
		// all of the relation stuff is now old since we modified the
		// second segments data
		if((relation == FIX_TOO_CLOSE) || (relation == FIX_TOO_FAR) || (relation == FIX_NO_RELATIONSHIP))
		{
			//return false;
		}
		else if(relation == FIX_MATCH)
		{
			// Ok, so now all of the segment stuff is passed into this method
			// it's the same as before except for the modified second segment
			fixDoubleAnchorBool = fixDoubleAnchor_Match(cigarInfo, relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chromName, mismatchNum);
		}
		else if((relation == FIX_INSERTION_NEIGHBOUR) || (relation == FIX_INSERTION_GAP))
		{
			// So, we think that extra stuff was added to our read
			fixDoubleAnchorBool = fixDoubleAnchor_Insertion(cigarInfo, relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chromName, mismatchNum);
		}
		else if((relation == FIX_DELETION_NEIGHBOUR) || (relation == FIX_DELETION_GAP))
		{
			// We think that stuff was deleted from our read
			fixDoubleAnchorBool = fixDoubleAnchor_Deletion(cigarInfo, relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chromName, mismatchNum);
		}
		else if((relation == FIX_SPLICE_NEIGHBOUR) || (relation == FIX_SPLICE_GAP))
		{
			// there is a large gap so we think we have found a splice
			fixDoubleAnchorBool = fixDoubleAnchor_Splice(cigarInfo, relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chromName, mismatchNum);
		}
		else
		{
			cout << "error in fixDoubleAnchor ... " << endl;
		}
		return fixDoubleAnchorBool;
	}

	bool fixDoubleAnchor_Match(Splice_Info* cigarInfo, int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, //Index_Info* indexInfo, 
		Index_Info* indexInfo, const string& chromName, int* mismatchNum
		)
	{
		// This happens if the end of the first and start of the second match
		// so a lot of this code seems unnessisary
		bool fixDoubleAnchorBool_Match = false;

		// We get the new difference between segment one and segment two
		int subSeqLengthInProcess = segmentLocInRead_2 - (segmentLocInRead_1 + segmentLength_1);

		// more chrom name/int conversions, great
		int chrNameInt = indexInfo->convertStringToInt(chromName);

		// If the difference between the two is 0 or 1 we just count the different as a read
		// error and add that to the mismatch count
		// it can't be less than 0 because our modification of segment two earlier accounted for that
		// our code will allow for negative
		if(subSeqLengthInProcess < 2)
		{
			// So, we modify the jump code by adding the second segment to it
			// very important we include the mismatched character to the cigar string
			// so, we just say this whole segment is matched and gloss over the mismatch
			// value
			Jump_Code matchJumpCode(subSeqLengthInProcess + segmentLength_2, "M");
			cigarInfo->jump_code.push_back(matchJumpCode);
			fixDoubleAnchorBool_Match = true;
			*mismatchNum = subSeqLengthInProcess;
		}
		else
		{
			// Ok, so the gap is large, we snag the read string which is not covered by either
			// of these segments
			string readSubSeqInProcess = readSeq_inProcess.substr(segmentLocInRead_1 + segmentLength_1 - 1,
				subSeqLengthInProcess);
		
			// Convert the chrom name to an integer... again
			int chrNameInt = indexInfo->convertStringToInt(chromName);
		
			// Gets the reference at this location
			// I don't know why we are converting from the reference genome to the chromosome
			//  You don't seem to gain anything, just complexity
			string chromSubSeqInProcess = indexInfo->getChromSequence(chrNameInt).substr(segmentMapPos_1 + segmentLength_1 - 1,
				subSeqLengthInProcess);

			// Ok, another problem here:
			// We get the reference from the first segment, not the second segment.
			// So, if the mismatch is on the second segment size we are going to lose that
			// We should check the mismatch before the second segment and after the first
			// segment
			// I also don't like the + 2 on the mismatches. We allow 2 mismatches on the string of length
			// 2?
			// This also means that if there are many segments we will increase our
			// mismatch limit. The mismatches should be constant with the number of mapped regions
			// for example one mismatch for each 10 nucleotides rather than two+ for each segment
			// which may be as short as length 3 or 4.
			size_t max_mismatch = (subSeqLengthInProcess)/LengthOfSeqPerMismatchAllowed + 2;
			size_t mismatch_bits = 0;
			size_t comb_bits = 0;

			// Compares the two strings, returns true if they are within the mismatch limit
			bool scoreStringBool = Utilities::stringsWithinMismatchLimit(readSubSeqInProcess, chromSubSeqInProcess, max_mismatch, mismatch_bits, comb_bits);// need to debug
			
			if(scoreStringBool)
			{
				// If the strings check out then add a new jump code which will be merged
				// later and the mismatch amount associated with this jump code
				*mismatchNum = mismatch_bits;
				Jump_Code matchJumpCode(subSeqLengthInProcess + segmentLength_2, "M");
				cigarInfo->jump_code.push_back(matchJumpCode);
			}
			else // score string failed, insert sudo-match jump code
			{
				// I have no idea what this is
				// It takes the region that it can't map and creates a
				// new jump code with a lower case m
				// The second segment is then added with the upper case M
				Jump_Code midMatchJumpCode(subSeqLengthInProcess, "m");
				Jump_Code secondMatchJumpCode(segmentLength_2, "M");
				cigarInfo->jump_code.push_back(midMatchJumpCode);
				cigarInfo->jump_code.push_back(secondMatchJumpCode);			
			}
			fixDoubleAnchorBool_Match = scoreStringBool;
		}

		// If the area between the two segments were mapped return true, else return false
		return fixDoubleAnchorBool_Match;
	}

	bool fixDoubleAnchor_Insertion(Splice_Info* cigarInfo, int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, //Index_Info* indexInfo, 
		Index_Info* indexInfo, const string& chromName, int* mismatchNum
		)
	{
		bool fixDoubleAnchorBool_Insertion = false;
		int subSeqLengthInProcess = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1;
		int insertionLength = (segmentMapPos_1 - segmentLocInRead_1) - (segmentMapPos_2 - segmentLocInRead_2);
		int chrNameInt = indexInfo->convertStringToInt(chromName);

		if(subSeqLengthInProcess <= insertionLength)
		{
			//(*mismatchNum) = subSeqLengthInProcess;
			//Jump_Code firstMatchJumpCode(segmentLength_1, "M");
			//cigarInfo->jump_code.push_back(firstMatchJumpCode);

			int secondMatchLength = segmentLength_2 + subSeqLengthInProcess - insertionLength;
			if(secondMatchLength > 0)
			{
				Jump_Code midInsertionJumpCode(insertionLength, "I");	
				Jump_Code secondMatchJumpCode(segmentLength_2 + subSeqLengthInProcess - insertionLength, "M");							
				cigarInfo->jump_code.push_back(midInsertionJumpCode);
				cigarInfo->jump_code.push_back(secondMatchJumpCode);
				fixDoubleAnchorBool_Insertion = true;
				(*mismatchNum) = 0;
			}
			else
			{
				Jump_Code midInsertionJumpCode(insertionLength, "i");
				Jump_Code midMatchJumpCode(subSeqLengthInProcess, "m");
				Jump_Code secondMatchJumpCode(segmentLength_2, "M");
				cigarInfo->jump_code.push_back(midInsertionJumpCode);
				cigarInfo->jump_code.push_back(midMatchJumpCode);
				cigarInfo->jump_code.push_back(secondMatchJumpCode);
				fixDoubleAnchorBool_Insertion = false;
			}
			//fixDoubleAnchorBool_Insertion = true;	
		}	  
		else
		{
			// gets the read within the range
			string readSubSeqInProcess = readSeq_inProcess.substr(segmentLocInRead_1 + segmentLength_1 - 1, subSeqLengthInProcess);

			// get the reference genome within the range
			string chromSubSeqInProcess = indexInfo->getChromSequence(chrNameInt).substr(segmentMapPos_1 + segmentLength_1 - 1, segmentMapPos_2 - 1 - (segmentMapPos_1 + segmentLength_1) + 1);
			size_t prefix_length = 0;
			size_t mismatch_bits = 0; //?
			size_t max_ins_mismatch = (subSeqLengthInProcess)/LengthOfSeqPerMismatchAllowed + 2;
			size_t comb_bits_ins = 0;

			GenomeScan* genome_scan = new GenomeScan;
			bool insertion_fixed = (*genome_scan).doubleAnchorInsertionWithinMismatchLimit(
				readSubSeqInProcess,
				chromSubSeqInProcess,
				max_ins_mismatch,
				prefix_length,
				comb_bits_ins,
				mismatch_bits);

			if(insertion_fixed)
			{
				int firstMatchLength = //segmentLength_1 + 
						prefix_length;
				int secondMatchLength = segmentLength_2 + subSeqLengthInProcess - prefix_length - insertionLength;

				if(secondMatchLength > 0)
				{
					Jump_Code firstMatchJumpCode(firstMatchLength, "M");
					Jump_Code insertionJumpCode(insertionLength, "I");
					Jump_Code secondMatchJumpCode(secondMatchLength, "M");
					cigarInfo->jump_code.push_back(firstMatchJumpCode);
					cigarInfo->jump_code.push_back(insertionJumpCode);
					cigarInfo->jump_code.push_back(secondMatchJumpCode);
					(*mismatchNum) = mismatch_bits;
				}	
				else
				{
					insertion_fixed = false;
					//Jump_Code firstMatchJumpCode(segmentLength_1, "M");
					Jump_Code insertionJumpCode(insertionLength, "i");
					Jump_Code midMatchJumpCode(subSeqLengthInProcess - insertionLength, "m");
					Jump_Code secondMatchJumpCode(segmentLength_2, "M");
					//cigarInfo->jump_code.push_back(firstMatchJumpCode);
					cigarInfo->jump_code.push_back(insertionJumpCode);
					cigarInfo->jump_code.push_back(midMatchJumpCode);
					cigarInfo->jump_code.push_back(secondMatchJumpCode);
				}			
			}
			else
			{
				//Jump_Code firstMatchJumpCode(segmentLength_1, "M");
				Jump_Code insertionJumpCode(insertionLength, "i");
				Jump_Code midMatchJumpCode(subSeqLengthInProcess - insertionLength, "m");
				Jump_Code secondMatchJumpCode(segmentLength_2, "M");
				//cigarInfo->jump_code.push_back(firstMatchJumpCode);
				cigarInfo->jump_code.push_back(insertionJumpCode);
				cigarInfo->jump_code.push_back(midMatchJumpCode);
				cigarInfo->jump_code.push_back(secondMatchJumpCode);				
			}	
			delete genome_scan;
			fixDoubleAnchorBool_Insertion = insertion_fixed;
		}	
		return fixDoubleAnchorBool_Insertion;
	}

	bool fixDoubleAnchor_Deletion(Splice_Info* cigarInfo, int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, //Index_Info* indexInfo, 
		Index_Info* indexInfo, const string& chromName, int* mismatchNum
		)
	{
		//cout << "start to fix deletion ... " << endl;
		bool fixDoubleAnchorBool_Deletion = false;
		int subSeqLengthInProcess = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1;
		//cout << "subSeqLengthInProcess: " << subSeqLengthInProcess << endl;
		int deletionLength = (segmentMapPos_2 - segmentLocInRead_2) - (segmentMapPos_1 - segmentLocInRead_1);
		//cout << "deletionLength: " << deletionLength << endl;
		int chrNameInt = indexInfo->convertStringToInt(chromName);

		if(subSeqLengthInProcess < 2)
		{
			//cout << " subSeqLengthInProcess < 2 " << endl;
			//Jump_Code firstMatchJumpCode(segmentLength_1, "M");
			Jump_Code deletionJumpCode(deletionLength, "D");
			Jump_Code secondMatchJumpCode(segmentLength_2 + subSeqLengthInProcess, "M");
			//cigarInfo->jump_code.push_back(firstMatchJumpCode);
			cigarInfo->jump_code.push_back(deletionJumpCode);
			cigarInfo->jump_code.push_back(secondMatchJumpCode);
			fixDoubleAnchorBool_Deletion = true;
			(*mismatchNum) = subSeqLengthInProcess;
		}
		else
		{
			string readSubSeqInProcess = readSeq_inProcess.substr(segmentLocInRead_1 + segmentLength_1 - 1, subSeqLengthInProcess);
			int chromSubSeqLengthInProcess= subSeqLengthInProcess + 2;
			string left_chrom_seq = indexInfo->getChromSequence(chrNameInt).substr(segmentMapPos_1 + segmentLength_1 - 1, chromSubSeqLengthInProcess);
			string right_chrom_seq = indexInfo->getChromSequence(chrNameInt).substr(segmentMapPos_2 - 1 - chromSubSeqLengthInProcess, chromSubSeqLengthInProcess);
		
			bool small_deletion = true;
			size_t prefix_length = 0;
			size_t max_double_splice_mismatch = subSeqLengthInProcess/LengthOfSeqPerMismatchAllowed + 2;
			size_t mismatch_bits = 0;
			size_t comb_bits = 0;
			GenomeScan* genome_scan = new GenomeScan;
			bool deletion_fixed = (*genome_scan).doubleAnchoredLeastMisWithinMismatchLimit(readSubSeqInProcess, left_chrom_seq, right_chrom_seq, 
				prefix_length, max_double_splice_mismatch, comb_bits, small_deletion, mismatch_bits);
			
			if(deletion_fixed)
			{
				int firstMatchLength = //segmentLength_1 + 
						prefix_length;
				int secondMatchLength = segmentLength_2 + subSeqLengthInProcess - prefix_length;
				Jump_Code firstMatchJumpCode(firstMatchLength, "M");
				Jump_Code deletionJumpCode(deletionLength, "D");
				Jump_Code secondMatchJumpCode(secondMatchLength, "M");
				cigarInfo->jump_code.push_back(firstMatchJumpCode);
				cigarInfo->jump_code.push_back(deletionJumpCode);
				cigarInfo->jump_code.push_back(secondMatchJumpCode);
				(*mismatchNum) = mismatch_bits;
			}
			else
			{
				//int firstMatchLength = segmentLength_1;
				int midMatchLength = subSeqLengthInProcess;
				int secondMatchLength = segmentLength_2;
				//Jump_Code firstMatchJumpCode(firstMatchLength, "M");
				Jump_Code deletionJumpCode(deletionLength, "d");
				Jump_Code midMatchJumpCode(midMatchLength, "m");
				Jump_Code secondMatchJumpCode(secondMatchLength, "M");				
				//cigarInfo->jump_code.push_back(firstMatchJumpCode);
				cigarInfo->jump_code.push_back(deletionJumpCode);
				cigarInfo->jump_code.push_back(midMatchJumpCode);
				cigarInfo->jump_code.push_back(secondMatchJumpCode);				
			}
			delete genome_scan;
			fixDoubleAnchorBool_Deletion = deletion_fixed;
		}
		return fixDoubleAnchorBool_Deletion;
	}

	bool fixDoubleAnchor_Splice(Splice_Info* cigarInfo, int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, //Index_Info* indexInfo, 
		Index_Info* indexInfo, const string& chromName, int* mismatchNum
		)
	{
		//cout << "start to fix splice" << endl;
		bool fixDoubleAnchorBool_Splice = false;
		int chrNameInt = indexInfo->convertStringToInt(chromName);
		//cout << "chrNameInt: " << chrNameInt << endl;

		int tmpBuffer_left = FixSpliceBuffer;
		if(tmpBuffer_left > segmentLength_1 - 2) //anchor >= 2
		{
			tmpBuffer_left = segmentLength_1 - 2;
		}
		int tmpBuffer_right = FixSpliceBuffer;
		if(tmpBuffer_right > segmentLength_2 - 2)
		{
			tmpBuffer_right = segmentLength_2 - 2;
		}

		int subSeqLengthInProcess = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1 + tmpBuffer_left + tmpBuffer_right;
		int spliceJunctionLength = (segmentMapPos_2 - segmentLocInRead_2) - (segmentMapPos_1 - segmentLocInRead_1);
		int chromSubSeqLengthInProcess = subSeqLengthInProcess + 2;

		//cout << "subSeqLengthInProcess: " << subSeqLengthInProcess << " spliceJunctionLength: " << spliceJunctionLength << endl;
		//cout << "segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left - 1: " << segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left - 1 << endl;
		string readSubSeqInProcess = readSeq_inProcess.substr(segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left - 1, subSeqLengthInProcess);
		//cout << "segmentMapPos_1 + segmentLength_1 - tmpBuffer_left - 1: " << segmentMapPos_1 + segmentLength_1 - tmpBuffer_left - 1 << endl;
		string left_chrom_seq = (indexInfo->getChromSequence(chrNameInt)).substr(segmentMapPos_1 + segmentLength_1 - tmpBuffer_left - 1, chromSubSeqLengthInProcess);
		//cout << "segmentMapPos_2 - 1 - chromSubSeqLengthInProcess: " << segmentMapPos_2 - 1 - chromSubSeqLengthInProcess << endl;
		string right_chrom_seq = (indexInfo->getChromSequence(chrNameInt)).substr(segmentMapPos_2 + tmpBuffer_right - 1 - chromSubSeqLengthInProcess, chromSubSeqLengthInProcess);

		//cout << "readSubSeqInProcess: " << endl << readSubSeqInProcess << endl;
		//cout << "left_chrom_seq: " << endl << left_chrom_seq << endl;
		//cout << "right_chrom_seq: " << endl << right_chrom_seq << endl; 

		size_t prefix_length = 0;
		size_t max_double_splice_mismatch = subSeqLengthInProcess/LengthOfSeqPerMismatchAllowed + 2;
		size_t mismatch_bits = 0;
		size_t comb_bits = 0;
		//bool adjacent_segments = false;
		bool double_anchor_noncanonical = true;//false;//DO_NONCANONICAL; ////debug
		string flank_seq;
		GenomeScan* genome_scan = new GenomeScan;
		bool splice_fixed = (*genome_scan).doubleAnchoredWithinMismatchLimit(readSubSeqInProcess, left_chrom_seq, right_chrom_seq, prefix_length, 
			max_double_splice_mismatch, comb_bits,
		 	//(!adjacent_segments) && 
			double_anchor_noncanonical, flank_seq, mismatch_bits);
		//cout << "splice_fixed: " << endl;
		if(splice_fixed)
		{
			int firstMatchLength = //segmentLength_1 
				- tmpBuffer_left + prefix_length;
			int secondMatchLength = segmentLength_2 - tmpBuffer_right + subSeqLengthInProcess - prefix_length;
			//cout << "prefix_length: " << prefix_length << ", spliceJunctionLength: " << spliceJunctionLength << endl;

			Jump_Code firstMatchJumpCode(firstMatchLength, "M");
			Jump_Code spliceJumpCode(spliceJunctionLength, "N");
			Jump_Code secondMatchJumpCode(secondMatchLength, "M"); 
			cigarInfo->jump_code.push_back(firstMatchJumpCode);
			cigarInfo->jump_code.push_back(spliceJumpCode);
			cigarInfo->jump_code.push_back(secondMatchJumpCode);
			(*mismatchNum) = mismatch_bits;
		}
		else
		{
			//int firstMatchLength = segmentLength_1;
			int midMatchLength = subSeqLengthInProcess - tmpBuffer_right - tmpBuffer_left;
			int secondMatchLength = segmentLength_2;
			//Jump_Code firstMatchJumpCode(firstMatchLength, "M");
			Jump_Code spliceJumpCode(spliceJunctionLength, "n");
			Jump_Code midMatchJumpCode(midMatchLength, "m");
			Jump_Code secondMatchJumpCode(secondMatchLength, "M"); 		
			//cigarInfo->jump_code.push_back(firstMatchJumpCode);
			cigarInfo->jump_code.push_back(spliceJumpCode);
			cigarInfo->jump_code.push_back(midMatchJumpCode);
			cigarInfo->jump_code.push_back(secondMatchJumpCode);	
		}
		delete genome_scan;
		fixDoubleAnchorBool_Splice = splice_fixed;
		return fixDoubleAnchorBool_Splice;
	}
};

#endif
