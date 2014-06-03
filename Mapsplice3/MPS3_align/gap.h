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
	 * removing path
	 */
	bool fixGapInPath(MappedRead* mappedRead, Index_Info* indexInfo)
	{
		Segment* longSegment = mappedRead->getFirstLongSegment();
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

	bool fixGapInPath(Path* pathInfo, MappedRead* mappedRead, Index_Info* indexInfo)
	{
		int pathVecSize = pathInfo->PathVec_seg.size();

		for(int tmpPathNO = 0; tmpPathNO < pathVecSize; tmpPathNO ++)
		{
			if(!(pathInfo->PathValidBoolVec[tmpPathNO]))
			{
				pathInfo->PathFixedBoolVec.push_back(false);				
				continue;
			}

			int firstSegGroupNO = (pathInfo->PathVec_seg[tmpPathNO])[0].first;
			int firstSegLength = mappedRead->getSegment(firstSegGroupNO)->getLength();

			Splice_Info* newPathSpliceInfo = new Splice_Info();
			Jump_Code firstJumpCode(firstSegLength, "M");
			newPathSpliceInfo->jump_code.push_back(firstJumpCode);

			int newPathMismatchNum = 0;
			bool newPathFixed = true;
			int tmpPathSegSize = pathInfo->PathVec_seg[tmpPathNO].size();

			for(int tmpPathSegNO = 0; tmpPathSegNO < tmpPathSegSize-1; tmpPathSegNO++)
			{
				int tmpMismatchNum = 0;
				int tmpSegGroupNO = (pathInfo->PathVec_seg[tmpPathNO])[tmpPathSegNO].first;
				int tmpSegCandiNO = (pathInfo->PathVec_seg[tmpPathNO])[tmpPathSegNO].second;

				int tmpSegGroupNO_next = (pathInfo->PathVec_seg[tmpPathNO])[tmpPathSegNO+1].first;
				int tmpSegCandiNO_next = (pathInfo->PathVec_seg[tmpPathNO])[tmpPathSegNO+1].second;

				int tmpRelation = mappedRead->checkSegRelation(tmpSegGroupNO, tmpSegCandiNO, tmpSegGroupNO_next, tmpSegCandiNO_next);

				Segment* firstSegment = mappedRead->getSegment(tmpSegGroupNO);
				Segment* secondSegment = mappedRead->getSegment(tmpSegGroupNO_next);

				int tmpSegmentLocInRead_1 = firstSegment->getLocationInRead();
				int tmpSegmentLocInRead_2 = secondSegment->getLocationInRead();
				int tmpSegmentLength_1 = firstSegment->getLength();
				int tmpSegmentLength_2 = secondSegment->getLength();
				unsigned int tmpSegmentMapPosInWholeGenome_1 = firstSegment->getAlignmentLocation(tmpSegCandiNO);
				unsigned int tmpSegmentMapPosInWholeGenome_2 = secondSegment->getAlignmentLocation(tmpSegCandiNO_next);

				unsigned int tmpChrNameInt, tmpChrPosInt;
				indexInfo->getChrLocation(tmpSegmentMapPosInWholeGenome_1, &tmpChrNameInt, &tmpChrPosInt);
				string tmpChrNameStr_1 = indexInfo->chrNameStr[tmpChrNameInt];
				int tmpSegmentMapPos_1 = tmpChrPosInt;
				indexInfo->getChrLocation(tmpSegmentMapPosInWholeGenome_2, &tmpChrNameInt, &tmpChrPosInt);
				string tmpChrNameStr_2 = indexInfo->chrNameStr[tmpChrNameInt];
				int tmpSegmentMapPos_2 = tmpChrPosInt;
				
				string tmpChrNameStr;
				if(tmpChrNameStr_1 == tmpChrNameStr_2)
				{
					tmpChrNameStr = tmpChrNameStr_1;
				}
				else
				{
					newPathFixed = false;
					pathInfo->PathFixedBoolVec.push_back(newPathFixed);
					break;
				}

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

				newPathMismatchNum = newPathMismatchNum + tmpMismatchNum;
				if(!tmpDoubleAnchorFixed)
				{
					newPathFixed = false;
					pathInfo->PathFixedBoolVec.push_back(newPathFixed);
					break;					
				}
			}

			if(newPathFixed)
			{
				newPathSpliceInfo->getFinalJumpCode();
				bool allJumpCodeValidBool = newPathSpliceInfo->allFinalJumpCodeValid();
				if(allJumpCodeValidBool)
				{
					(pathInfo->PathFixedBoolVec).push_back(allJumpCodeValidBool);
					(pathInfo->fixedPathVec).push_back(pair <int, Splice_Info*> (tmpPathNO, newPathSpliceInfo) );
					(pathInfo->fixedPathMismatchVec).push_back(newPathMismatchNum);
				}
				else
				{
					(pathInfo->PathFixedBoolVec).push_back(allJumpCodeValidBool);
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

	bool fixDoubleAnchor_extendBack(Splice_Info* cigarInfo, int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, Index_Info* indexInfo, const string& chromName, int* mismatchNum)
	{
		//cout << "fixDoubleAnchor_extendBack starts!" << endl;
		bool fixDoubleAnchorBool = false;

		int chrNameInt = indexInfo->convertStringToInt(chromName);
		int extendBackNumMax = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1;
		if(extendBackNumMax > segmentMapPos_2 - 1)
		{
			extendBackNumMax = segmentMapPos_2 - 1; 
		}
		int extendBackNum = extendBackInChromSeq(segmentLocInRead_2, readSeq_inProcess, 
			segmentMapPos_2, indexInfo->chromStr[chrNameInt], extendBackNumMax);

		segmentLocInRead_2 = segmentLocInRead_2 - extendBackNum;
		segmentLength_2 = segmentLength_2 + extendBackNum;
		segmentMapPos_2 = segmentMapPos_2 - extendBackNum;

		if((relation == FIX_TOO_CLOSE) || (relation == FIX_TOO_FAR) || (relation == FIX_NO_RELATIONSHIP))
		{
			//return false;
		}
		else if(relation == FIX_MATCH)
		{
			fixDoubleAnchorBool = fixDoubleAnchor_Match(cigarInfo, relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chromName, mismatchNum);
		}
		else if((relation == FIX_INSERTION_NEIGHBOUR) || (relation == FIX_INSERTION_GAP))
		{
			fixDoubleAnchorBool = fixDoubleAnchor_Insertion(cigarInfo, relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chromName, mismatchNum);
		}
		else if((relation == FIX_DELETION_NEIGHBOUR) || (relation == FIX_DELETION_GAP))
		{
			fixDoubleAnchorBool = fixDoubleAnchor_Deletion(cigarInfo, relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chromName, mismatchNum);
		}
		else if((relation == FIX_SPLICE_NEIGHBOUR) || (relation == FIX_SPLICE_GAP))
		{
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
		//cout << " fixDoubleAnchor_Match starts ... " << endl;
		bool fixDoubleAnchorBool_Match = false;
		int subSeqLengthInProcess = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1;
		int chrNameInt = indexInfo->convertStringToInt(chromName);
		//cout << "subSeqLengthInProcess: " << subSeqLengthInProcess << endl;

		if(subSeqLengthInProcess < 2)
		{
			//(*mismatchNum) = subSeqLengthInProcess;
			
			Jump_Code matchJumpCode(//segmentLength_1 + 
				subSeqLengthInProcess + segmentLength_2, "M");
			cigarInfo->jump_code.push_back(matchJumpCode);
			fixDoubleAnchorBool_Match = true;
			(*mismatchNum) = subSeqLengthInProcess;
		}
		else
		{
			string readSubSeqInProcess = readSeq_inProcess.substr(segmentLocInRead_1 + segmentLength_1 - 1,
				subSeqLengthInProcess);
		
			int chrNameInt = indexInfo->convertStringToInt(chromName);
		
			string chromSubSeqInProcess = indexInfo->chromStr[chrNameInt].substr(segmentMapPos_1 + segmentLength_1 - 1,
				subSeqLengthInProcess);

			size_t max_mismatch = (subSeqLengthInProcess)/LengthOfSeqPerMismatchAllowed + 2;
			size_t mismatch_bits = 0;
			size_t comb_bits = 0;
			//cout << "readSeqInProcess: " << endl << readSubSeqInProcess << endl;
			//cout << "chromSeqInProcess: " << endl << chromSubSeqInProcess << endl;

			bool scoreStringBool = Utilities::scoreString(readSubSeqInProcess, chromSubSeqInProcess, max_mismatch, mismatch_bits, comb_bits);// need to debug
			
			//cout << "scoreStringBool: " << scoreStringBool << endl;
			if(scoreStringBool)
			{
				(*mismatchNum) = mismatch_bits;
				Jump_Code matchJumpCode(//segmentLength_1 + 
					subSeqLengthInProcess + segmentLength_2, "M");
				cigarInfo->jump_code.push_back(matchJumpCode);
			}
			else // score string failed, insert sudo-match jump code
			{
				//cout << " score_string failed !" << endl;
				//Jump_Code firstMatchJumpCode(segmentLength_1, "M");
				Jump_Code midMatchJumpCode(subSeqLengthInProcess, "m");
				Jump_Code secondMatchJumpCode(segmentLength_2, "M");	
				//cigarInfo->jump_code.push_back(firstMatchJumpCode);
				cigarInfo->jump_code.push_back(midMatchJumpCode);
				cigarInfo->jump_code.push_back(secondMatchJumpCode);			
			}
			fixDoubleAnchorBool_Match = scoreStringBool;
		}
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
			string readSubSeqInProcess = readSeq_inProcess.substr(segmentLocInRead_1 + segmentLength_1 - 1, subSeqLengthInProcess);
			string chromSubSeqInProcess = indexInfo->chromStr[chrNameInt].substr(segmentMapPos_1 + segmentLength_1 - 1, segmentMapPos_2 - 1 - (segmentMapPos_1 + segmentLength_1) + 1);
			size_t prefix_length = 0;
			size_t mismatch_bits = 0; //?
			size_t max_ins_mismatch = (subSeqLengthInProcess)/LengthOfSeqPerMismatchAllowed + 2;
			size_t comb_bits_ins = 0;

			GenomeScan* genome_scan = new GenomeScan;
			bool insertion_fixed = (*genome_scan).Double_anchored_score_ins(readSubSeqInProcess, chromSubSeqInProcess, max_ins_mismatch, prefix_length, comb_bits_ins, mismatch_bits); //X: fix insertion

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
			string left_chrom_seq = indexInfo->chromStr[chrNameInt].substr(segmentMapPos_1 + segmentLength_1 - 1, chromSubSeqLengthInProcess);
			string right_chrom_seq = indexInfo->chromStr[chrNameInt].substr(segmentMapPos_2 - 1 - chromSubSeqLengthInProcess, chromSubSeqLengthInProcess);
		
			bool small_deletion = true;
			size_t prefix_length = 0;
			size_t max_double_splice_mismatch = subSeqLengthInProcess/LengthOfSeqPerMismatchAllowed + 2;
			size_t mismatch_bits = 0;
			size_t comb_bits = 0;
			GenomeScan* genome_scan = new GenomeScan;
			bool deletion_fixed = (*genome_scan).Double_anchored_score_least_mis(readSubSeqInProcess, left_chrom_seq, right_chrom_seq, 
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
		string left_chrom_seq = (indexInfo->chromStr[chrNameInt]).substr(segmentMapPos_1 + segmentLength_1 - tmpBuffer_left - 1, chromSubSeqLengthInProcess);
		//cout << "segmentMapPos_2 - 1 - chromSubSeqLengthInProcess: " << segmentMapPos_2 - 1 - chromSubSeqLengthInProcess << endl;
		string right_chrom_seq = (indexInfo->chromStr[chrNameInt]).substr(segmentMapPos_2 + tmpBuffer_right - 1 - chromSubSeqLengthInProcess, chromSubSeqLengthInProcess); 

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
		bool splice_fixed = (*genome_scan).Double_anchored_score(readSubSeqInProcess, left_chrom_seq, right_chrom_seq, prefix_length, 
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
