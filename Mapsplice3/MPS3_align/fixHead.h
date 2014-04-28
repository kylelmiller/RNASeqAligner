#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>

bool fixInsertionGapHead(int anchorLength, unsigned int anchorMapPos, unsigned int mainMapPos, 
	int mainLocInRead, int extendLength, char* read, char* chrom, Splice_Info* headCigarInfo, 
	const string& readString)
	// mainLocInRead is the location in read after extension
{
	//cout << "test chromString start..." << endl;
	bool insertion_fixed = false;
	int headLength = mainLocInRead - 1 + extendLength;
	unsigned int readSeqLeftBound = anchorLength+1;
	unsigned int readSeqRightBound = mainLocInRead - 1;
	unsigned int readFixSeqLength = readSeqRightBound - readSeqLeftBound + 1;
	unsigned int chromSeqLeftBound = anchorMapPos + anchorLength;
	unsigned int chromSeqRightBound = mainMapPos + mainLocInRead - 2;
	unsigned int chromFixSeqLength = chromSeqRightBound - chromSeqLeftBound + 1; 

	////debug("readSeqLeftBound = "); //debug(readSeqLeftBound); //debug(" readSeqRightBound = "); 
	////debugln(readSeqRightBound); //debug("chromSeqLeftBound = "); //debug(chromSeqLeftBound); 
	////debug("chromSeqRightBound = "); //debugln(chromSeqRightBound); 
	////debug("readFixSeqLength = "); //debugln(readFixSeqLength);
	////debug("chromFixSeqLength = "); //debugln(chromFixSeqLength);

	if(chromSeqLeftBound >= chromSeqRightBound) ////debug
	{
		//int matchPrefixLength = chromSeqRightBound - anchorMapPos;
		int matchPrefixLength = anchorLength + 1 + chromSeqRightBound - chromSeqLeftBound;
		if(matchPrefixLength < 0)
		{
			////debug("matchPrefixLength<0 = ");
			////debugln(matchPrefixLength);
			return false;
		}
		int insertionLength = readSeqRightBound - matchPrefixLength;
		int matchSuffixLength = extendLength;
		Jump_Code matchPrefixJumpCode(matchPrefixLength, "M");
		Jump_Code insertionJumpCode(insertionLength, "I");
		Jump_Code matchSuffixJumpCode(matchSuffixLength, "M");
		headCigarInfo->jump_code.push_back(matchPrefixJumpCode);
		headCigarInfo->jump_code.push_back(insertionJumpCode);
		headCigarInfo->jump_code.push_back(matchSuffixJumpCode);		
		return true;
	}
	else
	{	
		string pending_seq_ins = readString.substr(readSeqLeftBound-1, readFixSeqLength); 
		string chrom_seq_ins = chromString.substr(chromSeqLeftBound-1, chromFixSeqLength);

		////debug("pending_seq_ins = "); //debugln(pending_seq_ins);
		/////debug("chrom_seq_ins = "); //debugln(chrom_seq_ins);

		size_t prefix_length = 0;
		size_t mismatch_bits = 0; //?
		size_t max_ins_mismatch = 2;

		////debugln("start Double_anchored_score_ins function ...");
		GenomeScan* genome_scan = new GenomeScan;
		insertion_fixed = (*genome_scan).Double_anchored_score_ins
		(pending_seq_ins, chrom_seq_ins, max_ins_mismatch, prefix_length, mismatch_bits); //X: fix insertion
		////debugln("finish Double_anchored_score_ins function ...");

		if(insertion_fixed)
		{
			////debugln("insertion_fixed = true");
			int matchPrefixLength = anchorLength + prefix_length;
			int insertionLength = anchorMapPos - mainMapPos;
			int matchSuffixLength = headLength - insertionLength - prefix_length - anchorLength;
			Jump_Code matchPrefixJumpCode(matchPrefixLength, "M");
			Jump_Code insertionJumpCode(insertionLength, "I");
			Jump_Code matchSuffixJumpCode(matchSuffixLength, "M");
			headCigarInfo->jump_code.push_back(matchPrefixJumpCode);
			headCigarInfo->jump_code.push_back(insertionJumpCode);
			headCigarInfo->jump_code.push_back(matchSuffixJumpCode);
		}
		//delete genome_scan;
		////debug ("pending_seq_ins = "); //debugln(pending_seq_ins); 
		return insertion_fixed;
	}
}

bool fixDeletionGapHead(int anchorLength, unsigned int anchorMapPos, unsigned int mainMapPos,
	int mainLocInRead, int extendLength, char* read, char* chrom, Splice_Info* headCigarInfo,
	const string& readString 
	//string chromString, 
	//string readString
	)
{
	////debugln("start fixDeletionGap function ...");
	int headLength = mainLocInRead -1 + extendLength;
	bool deletion_fixed = false;
	unsigned int readSeqLeftBound = anchorLength + 1;
	unsigned int readSeqRightBound = mainLocInRead - 1;
	unsigned int readFixSeqLength = readSeqRightBound - readSeqLeftBound + 1;
	unsigned int chromSeqLeftBound = anchorMapPos + anchorLength;
	unsigned int chromSeqRightBound = mainMapPos + mainLocInRead - 2;
	unsigned int chromFixSeqLength = readFixSeqLength + 2;

	string pending_seq = readString.substr(readSeqLeftBound-1, readFixSeqLength);
	string left_chrom_seq = chromString.substr(chromSeqLeftBound-1, chromFixSeqLength);
	string right_chrom_seq = chromString.substr(chromSeqRightBound-chromFixSeqLength, chromFixSeqLength);	

	bool small_deletion = true;
	size_t prefix_length = 0;
	size_t max_double_splice_mismatch = 2;
	size_t mismatch_bits = 0;
	GenomeScan* genome_scan = new GenomeScan;
	deletion_fixed = (*genome_scan).Double_anchored_score_least_mis
	(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_double_splice_mismatch, mismatch_bits, small_deletion);

	if(deletion_fixed)
	{
		int matchPrefixLength = anchorLength + prefix_length;
		int deletionLength = mainMapPos - anchorMapPos;
		int matchSuffixLength = headLength - prefix_length-anchorLength;
		Jump_Code matchPrefixJumpCode(matchPrefixLength, "M");
		Jump_Code deletionJumpCode(deletionLength, "D");
		Jump_Code matchSuffixJumpCode(matchSuffixLength, "M");
		headCigarInfo->jump_code.push_back(matchPrefixJumpCode);
		headCigarInfo->jump_code.push_back(deletionJumpCode);
		headCigarInfo->jump_code.push_back(matchSuffixJumpCode);
	}
	//delete genome_scan;
	return deletion_fixed;

}

bool fixSpliceNeighbourHead(int anchorLength, unsigned int anchorMapPos, unsigned int mainMapPos,
	int mainLocInRead, int extendLength, char* read, char* chrom, Splice_Info* headCigarInfo, 
	const string& readString, Index_Info* indexInfo)
{
	#ifdef debug
	cout << "start fixSpliceNeighbourHead function ......" << endl;
	#endif
	bool splice_fixed = false;
	int buffer = FIX_SPLICE_BUFFER;
	int headLength = anchorLength + extendLength;
	unsigned int readSeqLeftBound = anchorLength + 1 - buffer;
	unsigned int readSeqRightBound = mainLocInRead - 1 + buffer;
	unsigned int readFixSeqLength = readSeqRightBound - readSeqLeftBound + 1;
	unsigned int chromSeqLeftBound = anchorMapPos + anchorLength - buffer;
	unsigned int chromSeqRightBound = mainMapPos + mainLocInRead - 2 + buffer;
	unsigned int chromFixSeqLength = readFixSeqLength + 2;

	string pending_seq = readString.substr(readSeqLeftBound-1, readFixSeqLength);
	string left_chrom_seq = chromString.substr(chromSeqLeftBound-1, chromFixSeqLength);
	string right_chrom_seq = chromString.substr(chromSeqRightBound-chromFixSeqLength, chromFixSeqLength);	

	size_t prefix_length = 0;
	size_t max_double_splice_mismatch = 2;
	size_t mismatch_bits = 0;
	bool adjacent_segments = false;
	bool double_anchor_noncanonical = DO_NONCANONICAL;////debug
	string flank_seq;
	GenomeScan* genome_scan = new GenomeScan;
	splice_fixed = (*genome_scan).Double_anchored_score
	(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_double_splice_mismatch, 
	mismatch_bits, (!adjacent_segments) && double_anchor_noncanonical, flank_seq);

	if(splice_fixed)
	{
		int matchPrefixLength = anchorLength + prefix_length - buffer;
		int spliceLength = mainMapPos - anchorMapPos;
		int matchSuffixLength = headLength - matchPrefixLength;
		Jump_Code matchPrefixJumpCode(matchPrefixLength, "M");
		Jump_Code spliceJumpCode(spliceLength, "N");
		Jump_Code matchSuffixJumpCode(matchSuffixLength, "M");
		headCigarInfo->jump_code.push_back(matchPrefixJumpCode);
		headCigarInfo->jump_code.push_back(spliceJumpCode);
		headCigarInfo->jump_code.push_back(matchSuffixJumpCode);
		
	    if((!BUILD_SPLICE_HASH_FROM_FILE) && (BUILD_SPLICE_HASH_FROM_FIRST_MAPPING))
	    {
	    	unsigned int spliceStartPosInWholeGenome = anchorMapPos + matchPrefixLength - 1;
	    	unsigned int spliceJunctionStartPos, spliceJunctionEndPos; 
	    	unsigned int spliceJunctionChrInt;
	    	indexInfo->getChrLocation(spliceStartPosInWholeGenome, &spliceJunctionChrInt, &spliceJunctionStartPos);
	    	(spliceJunctionEndPos) = (spliceJunctionStartPos) + spliceLength + 1;
			bool insertion = insertSpliceJunction2ReverseHashForHead(spliceJunctionReverse, spliceJunctionChrInt, 
				spliceJunctionStartPos, spliceJunctionEndPos)
			&& insertSpliceJunction2NormalHashForTail(spliceJunctionNormal, spliceJunctionChrInt, 
				spliceJunctionStartPos, spliceJunctionEndPos);

			if(PRINT_JUNC)
			{
				insertSpliceJunction2HashForPrint(spliceJunctionChrInt, 
				spliceJunctionStartPos, spliceJunctionEndPos);
				//cout << (spliceJunctionChrInt) << " " << (spliceJunctionStartPos) << " " << (spliceJunctionEndPos) << " case 1"<< endl; 
			}

			/*
			if(!insertion)
			{
				cout << "insertion failed " << endl;
				cout << chrInt << " " << spliceStartPos << " " << spliceEndPos << endl;
			}*/
	    }
	}
	delete genome_scan;
	return splice_fixed;	
}

bool fixSpliceGapHead(int anchorLength, unsigned int anchorMapPos, unsigned int mainMapPos,
	int mainLocInRead, int extendLength, char* read, char* chrom, Splice_Info* headCigarInfo,
	const string& readString, Index_Info* indexInfo)
{
	bool splice_fixed = false;
	int buffer = FIX_SPLICE_BUFFER;
	int headLength = mainLocInRead - 1 + extendLength;
	unsigned int readSeqLeftBound = anchorLength + 1 - buffer;
	unsigned int readSeqRightBound = mainLocInRead - 1 + buffer;
	unsigned int readFixSeqLength = readSeqRightBound - readSeqLeftBound + 1;
	unsigned int chromSeqLeftBound = anchorMapPos + anchorLength - buffer;
	unsigned int chromSeqRightBound = mainMapPos + mainLocInRead - 2 + buffer;
	unsigned int chromFixSeqLength = readFixSeqLength + 2;

	string pending_seq = readString.substr(readSeqLeftBound-1, readFixSeqLength);
	string left_chrom_seq = chromString.substr(chromSeqLeftBound-1, chromFixSeqLength);
	string right_chrom_seq = chromString.substr(chromSeqRightBound-chromFixSeqLength, chromFixSeqLength);	

	size_t prefix_length = 0;
	size_t max_double_splice_mismatch = 2;
	size_t mismatch_bits = 0;
	bool adjacent_segments = false;
	bool double_anchor_noncanonical = false;////debug
	string flank_seq;
	GenomeScan* genome_scan = new GenomeScan;
	splice_fixed = (*genome_scan).Double_anchored_score
	(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_double_splice_mismatch, 
	mismatch_bits, (!adjacent_segments) && double_anchor_noncanonical, flank_seq);

	if(splice_fixed)
	{
		int matchPrefixLength = anchorLength + prefix_length - buffer;
		int spliceLength = mainMapPos - anchorMapPos;
		int matchSuffixLength = headLength - matchPrefixLength;
		Jump_Code matchPrefixJumpCode(matchPrefixLength, "M");
		Jump_Code spliceJumpCode(spliceLength, "N");
		Jump_Code matchSuffixJumpCode(matchSuffixLength, "M");
		headCigarInfo->jump_code.push_back(matchPrefixJumpCode);
		headCigarInfo->jump_code.push_back(spliceJumpCode);
		headCigarInfo->jump_code.push_back(matchSuffixJumpCode);

	    if((!BUILD_SPLICE_HASH_FROM_FILE) && (BUILD_SPLICE_HASH_FROM_FIRST_MAPPING))
	    {
	    	unsigned int spliceStartPosInWholeGenome = anchorMapPos + matchPrefixLength - 1;
	    	unsigned int spliceJunctionStartPos, spliceJunctionEndPos; 
	    	unsigned int spliceJunctionChrInt;
	    	indexInfo->getChrLocation(spliceStartPosInWholeGenome, &spliceJunctionChrInt, &spliceJunctionStartPos);
	    	(spliceJunctionEndPos) = (spliceJunctionStartPos) + spliceLength + 1;
			bool insertion = insertSpliceJunction2ReverseHashForHead(spliceJunctionReverse, spliceJunctionChrInt, 
				spliceJunctionStartPos, spliceJunctionEndPos)
			&& insertSpliceJunction2NormalHashForTail(spliceJunctionNormal, spliceJunctionChrInt, 
				spliceJunctionStartPos, spliceJunctionEndPos);

			if(PRINT_JUNC)
			{
				//cout << (spliceJunctionChrInt) << " " << (spliceJunctionStartPos) << " " << (spliceJunctionEndPos) << " case 2"<< endl; 
				insertSpliceJunction2HashForPrint(spliceJunctionChrInt, spliceJunctionStartPos, spliceJunctionEndPos);			
			}

			/*
			if(!insertion)
			{
				cout << "insertion failed " << endl;
				cout << chrInt << " " << spliceStartPos << " " << spliceEndPos << endl;
			}*/
	    }
		////debugln(" "); //debug("flankSeq = "); //debugln(flank_seq);
	}
	delete genome_scan;
	return splice_fixed;	
}

bool fixDoubleAnchorHead(int relation, Splice_Info* headCigarInfo, 
	unsigned int anchorMapPos, unsigned int mainMapPos, 
	unsigned int mainLocInRead, int anchorLength, 
	char* read, char* chrom, int extendLength, const string& readString, Index_Info* indexInfo)
	// mainLocInRead is the location in read after extension
{
	////debug("//debug -- start fixDoubleAnchorHead, relation = "); //debugln(relation);
	if((relation == FIX_TOO_CLOSE) || (relation == FIX_TOO_FAR) || (relation == FIX_NO_RELATIONSHIP))
	{
		////debugln(" No relation or too close or too far!!! ");
		return false;
	}
	else if(relation == FIX_REMAPPING_SHORT_HEAD)
	{
		int spliceDistance = anchorLength;
		int headMatchLength = mainLocInRead;
		Jump_Code headMatchJumpCode(headMatchLength, "M");
		Jump_Code spliceDistanceJumpCode(spliceDistance, "N");
		Jump_Code extendMatchJumpCode(extendLength, "M");
		headCigarInfo->jump_code.push_back(headMatchJumpCode);
		headCigarInfo->jump_code.push_back(spliceDistanceJumpCode);
		headCigarInfo->jump_code.push_back(extendMatchJumpCode);
		return true;
	}
	else if(relation == FIX_SOFTCLIPPING)
	{
		int headMatchLength = mainLocInRead - 1; //extendLength;
		Jump_Code headMatchJumpCode(headMatchLength, "S");
		Jump_Code extendMatchJumpCode(extendLength, "M");
		headCigarInfo->jump_code.push_back(headMatchJumpCode);
		headCigarInfo->jump_code.push_back(extendMatchJumpCode);
		return true;		
	}
	else if(relation == FIX_MATCH)
	{
		int headMatchLength = mainLocInRead - 1 + extendLength;
		Jump_Code headMatchJumpCode(headMatchLength, "M");
		headCigarInfo->jump_code.push_back(headMatchJumpCode);
		return true;
	}
	else if (relation == FIX_INSERTION_NEIGHBOUR)
	{
		int insertionLength = anchorMapPos - mainMapPos;
		Jump_Code matchPrefixJumpCode(anchorLength, "M");
		Jump_Code insertionJumpCode(insertionLength, "I");
		Jump_Code matchSuffixJumpCode(extendLength-insertionLength, "M");
		headCigarInfo->jump_code.push_back(matchPrefixJumpCode);
		headCigarInfo->jump_code.push_back(insertionJumpCode);
		headCigarInfo->jump_code.push_back(matchSuffixJumpCode);
		return true;
	}
	else if (relation == FIX_INSERTION_GAP)
	{
		return fixInsertionGapHead(anchorLength, anchorMapPos, mainMapPos, mainLocInRead, 
			extendLength, read, chrom, headCigarInfo, //chromString, 
			readString
			);		
	}
	else if (relation == FIX_DELETION_NEIGHBOUR)
	{
		int deletionLength = mainMapPos - anchorMapPos;
		Jump_Code matchPrefixJumpCode(anchorLength, "M");
		Jump_Code deletionJumpCode(deletionLength, "D");
		Jump_Code matchSuffixJumpCode(extendLength, "M");
		headCigarInfo->jump_code.push_back(matchPrefixJumpCode);
		headCigarInfo->jump_code.push_back(deletionJumpCode);
		headCigarInfo->jump_code.push_back(matchSuffixJumpCode);
		return true;
	}
	else if (relation == FIX_DELETION_GAP)
	{
		return fixDeletionGapHead(anchorLength, anchorMapPos, mainMapPos, mainLocInRead, 
			extendLength, read, chrom, headCigarInfo, readString );		
	}
	else if (relation == FIX_SPLICE_NEIGHBOUR)
	{
		/*
		int spliceLength = mainMapPos - anchorMapPos;
		Jump_Code matchPrefixJumpCode(anchorLength, "M");
		Jump_Code spliceJumpCode(spliceLength, "N");
		Jump_Code matchSuffixJumpCode(extendLength, "M");
		headCigarInfo->jump_code.push_back(matchPrefixJumpCode);
		headCigarInfo->jump_code.push_back(spliceJumpCode);
		headCigarInfo->jump_code.push_back(matchSuffixJumpCode);
		return true;*/
		return fixSpliceNeighbourHead(anchorLength, anchorMapPos, mainMapPos, mainLocInRead, 
			extendLength, read, chrom, headCigarInfo, readString, indexInfo);			
 
	}
	else if (relation == FIX_SPLICE_GAP)
	{
		return fixSpliceGapHead(anchorLength, anchorMapPos, mainMapPos, mainLocInRead, 
			extendLength, read, chrom, headCigarInfo, readString, indexInfo);			
	}
	else
	{
		////debugln("wrong relation");
		return FIX_NO_RELATIONSHIP;
	}		
}

int checkRelationHead(int mainLocInRead, unsigned int anchorMapLoc, unsigned int mainMapLoc, int anchorLength,
		Index_Info* indexInfo) 
	// mainLocInRead is the location in read after extension
{
	////debugln("//debug -- start checkRelationHead"); 
	if((anchorMapLoc>(indexInfo->indexSize))||(mainMapLoc>(indexInfo->indexSize)))
	{
		////debug("mapPos < 0 -- anchorMapLoc: "); //debug(anchorMapLoc); //debug(" mainMapLoc: "); //debugln(mainMapLoc);
		return FIX_NO_RELATIONSHIP;
	}
	if (anchorMapLoc == mainMapLoc)
		return FIX_MATCH;
	else if(anchorMapLoc > mainMapLoc)
	{
		if (anchorMapLoc - mainMapLoc > MAX_INSERTION_LENGTH)
			return FIX_TOO_CLOSE;
		else
		{
			if(mainLocInRead == anchorLength + 1)
				return FIX_INSERTION_NEIGHBOUR;
			else
				return FIX_INSERTION_GAP;
		}
	}
	else if(anchorMapLoc < mainMapLoc)
	{
		unsigned int mapDistance = mainMapLoc - anchorMapLoc;
		if(mapDistance > MAX_SPLICE_LENGTH)
		{
			return FIX_TOO_FAR;
		}
		else if(mapDistance > MAX_DELETION_LENGTH)
		{
			if(mainLocInRead == anchorLength + 1)
			{	
				return FIX_SPLICE_NEIGHBOUR;
			}
			else
			{
				return FIX_SPLICE_GAP;
			}
		}
		else 
		{
			if(mainLocInRead == anchorLength + 1)		
			{
				return FIX_DELETION_NEIGHBOUR;
			}
			else
			{
				return FIX_DELETION_GAP;
			}
		}
	}
	else 
	{
		////debugln("check head relation error!!!");
		return FIX_NO_RELATIONSHIP;
	}
	//debugln("//debug -- finish checkRelationHead");
}