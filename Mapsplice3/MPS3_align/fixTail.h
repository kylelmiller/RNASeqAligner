/*int checkRelation(unsigned int segStartNum2, unsigned int segEndNum1, unsigned int alignLoc1, unsigned int alignLoc2)
{
	if (segEndNum1 < segStartNum2)
	{
		if ((segStartNum2 - segEndNum1) == 1)
		{
			if (alignLoc1 > alignLoc2)
			{
				unsigned int gapInChr = alignLoc1 - alignLoc2;
				if (gapInChr > MAX_INSERTION_LENGTH)
					return FIX_TOO_CLOSE;
				else
					return FIX_INSERTION_NEIGHBOUR; 
			}
			else
			{
				unsigned int gapInChr = alignLoc2 - alignLoc1;
				if (gapInChr > MAX_SPLICE_LENGTH)
					return FIX_TOO_FAR;
				else if(gapInChr > MAX_DELETION_LENGTH)
					return FIX_SPLICE_NEIGHBOUR;
				else
					return FIX_DELETION_NEIGHBOUR;
			}
		}
		else
		{
			if (alignLoc1 > alignLoc2)
			{
				unsigned int gapInChr = alignLoc1 - alignLoc2;
				if (gapInChr > MAX_INSERTION_LENGTH)
					return FIX_TOO_CLOSE;
				else
					return FIX_INSERTION_GAP; 
			}
			else
			{
				unsigned int gapInChr = alignLoc2 - alignLoc1;
				if (gapInChr > MAX_SPLICE_LENGTH)
					return FIX_TOO_FAR;
				else if(gapInChr > MAX_DELETION_LENGTH)
					return FIX_SPLICE_GAP;
				else
					return FIX_DELETION_GAP;
			}
		}		
	}

	else
		return FIX_NO_RELATIONSHIP;
}*/

bool fixInsertionGapTail(Splice_Info* tailCigarInfo, 
	int anchorLength, unsigned int mainEndLocInRead, 
	unsigned int anchorMapPos, unsigned int mainMapPos, char* read, char* chrom, 
	const string& readString, int readLength )
	// extend function can be added when finding the anchor
{
	//int readLength = READ_LENGTH;
	bool insertion_fixed = false;
	int tailLength = readLength - mainEndLocInRead;
	unsigned int readSeqLeftBound = mainEndLocInRead + 1;
	unsigned int readSeqRightBound = readLength - anchorLength;
	unsigned int readFixSeqLength = readSeqRightBound - readSeqLeftBound + 1;
	unsigned int chromSeqLeftBound = mainMapPos + mainEndLocInRead;
	unsigned int chromSeqRightBound = anchorMapPos + readLength - anchorLength - 1;
	unsigned int chromFixSeqLength = chromSeqRightBound - chromSeqLeftBound + 1; 

	if(chromSeqLeftBound >= chromSeqRightBound)
	{
		//int matchPrefixLength = chromSeqRightBound - anchorMapPos;
		//if(matchPrefixLength < 0)
		//	return false;
		int insertionLength = readFixSeqLength + chromSeqLeftBound - chromSeqRightBound - 1;
		int matchSuffixLength = tailLength - insertionLength;
		if(matchSuffixLength < 0)
		{
			////debug("matchSuffixLength<0 = "); //debugln(matchSuffixLength);
			return false;
		}
		//Jump_Code matchPrefixJumpCode(matchPrefixLength, "M");
		Jump_Code insertionJumpCode(insertionLength, "I");
		Jump_Code matchSuffixJumpCode(matchSuffixLength, "M");
		//headCigarInfo->jump_code.push_back(matchPrefixJumpCode);
		tailCigarInfo->jump_code.push_back(insertionJumpCode);
		tailCigarInfo->jump_code.push_back(matchSuffixJumpCode);		
		return true;		
	}
	else 
	{
		string pending_seq_ins = readString.substr(readSeqLeftBound-1, readFixSeqLength); //cout << "readSequence = " << pending_seq_ins << endl;
		string chrom_seq_ins = chromString.substr(chromSeqLeftBound-1, chromFixSeqLength); //cout << "chromSequence = " << chrom_seq_ins << endl;

		size_t prefix_length = 0;
		size_t mismatch_bits = 0; 
		size_t max_ins_mismatch = 2;

		GenomeScan* genome_scan = new GenomeScan;
		insertion_fixed = (*genome_scan).Double_anchored_score_ins
		(pending_seq_ins, chrom_seq_ins, max_ins_mismatch, prefix_length, mismatch_bits); //X: fix insertion

		if(insertion_fixed)
		{
			int matchPrefixLength = prefix_length;
			int insertionLength = mainMapPos - anchorMapPos;
			int matchSuffixLength = tailLength - prefix_length - insertionLength;
			Jump_Code matchPrefixJumpCode(matchPrefixLength, "M");
			Jump_Code insertionJumpCode(insertionLength, "I");
			Jump_Code matchSuffixJumpCode(matchSuffixLength, "M");
			tailCigarInfo->jump_code.push_back(matchPrefixJumpCode);
			tailCigarInfo->jump_code.push_back(insertionJumpCode);
			tailCigarInfo->jump_code.push_back(matchSuffixJumpCode);
		}
		//delete genome_scan;
		return insertion_fixed;
	}
}

bool fixDeletionGapTail(Splice_Info* tailCigarInfo, 
	int anchorLength, unsigned int mainEndLocInRead, 
	unsigned int anchorMapPos, unsigned int mainMapPos, char* read, char* chrom,
	const string& readString, int readLength)
{
	//int readLength = READ_LENGTH;
	int tailLength = readLength - mainEndLocInRead;
	bool deletion_fixed = false;
	int readSeqLeftBound = mainEndLocInRead + 1;
	int readSeqRightBound = readLength - anchorLength;
	int readFixSeqLength = readSeqRightBound - readSeqLeftBound + 1;
	unsigned int chromSeqLeftBound = mainMapPos + mainEndLocInRead;
	unsigned int chromSeqRightBound = anchorMapPos + readLength - anchorLength - 1;
	int chromFixSeqLength = readFixSeqLength + 2;

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
		int matchPrefixLength = prefix_length;
		int deletionLength = anchorMapPos - mainMapPos;
		int matchSuffixLength = tailLength - prefix_length;
		Jump_Code matchPrefixJumpCode(matchPrefixLength, "M");
		Jump_Code deletionJumpCode(deletionLength, "D");
		Jump_Code matchSuffixJumpCode(matchSuffixLength, "M");
		tailCigarInfo->jump_code.push_back(matchPrefixJumpCode);
		tailCigarInfo->jump_code.push_back(deletionJumpCode);
		tailCigarInfo->jump_code.push_back(matchSuffixJumpCode);
	}
	//delete genome_scan;
	return deletion_fixed;
}

bool fixSpliceNeighbourTail(Splice_Info* tailCigarInfo, 
	int anchorLength, unsigned int mainEndLocInRead, 
	unsigned int anchorMapPos, unsigned int mainMapPos, char* read, char* chrom,
	const string& readString, int readLength, Index_Info* indexInfo
	)
{
	int buffer = FIX_SPLICE_BUFFER;
	//int readLength = READ_LENGTH;
	int tailLength = anchorLength;
	bool splice_fixed = false;
	int readSeqLeftBound = mainEndLocInRead + 1 - buffer;
	int readSeqRightBound = mainEndLocInRead + buffer;
	int readFixSeqLength = 2 * buffer;
	unsigned int chromSeqLeftBound = mainMapPos + mainEndLocInRead -buffer;
	unsigned int chromSeqRightBound = anchorMapPos + readLength - anchorLength - 1 + buffer;
	int chromFixSeqLength = readFixSeqLength + 2;

	string pending_seq = readString.substr(readSeqLeftBound-1, readFixSeqLength);
	string left_chrom_seq = chromString.substr(chromSeqLeftBound-1, chromFixSeqLength);
	string right_chrom_seq = chromString.substr(chromSeqRightBound-chromFixSeqLength, chromFixSeqLength);	

	size_t prefix_length = 0;
	size_t max_double_splice_mismatch = 2;
	size_t mismatch_bits = 0;
	bool adjacent_segments = false;
	bool double_anchor_noncanonical = DO_NONCANONICAL; ////debug 
	string flank_seq;
	GenomeScan* genome_scan = new GenomeScan;
	splice_fixed = (*genome_scan).Double_anchored_score
	(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_double_splice_mismatch, 
	mismatch_bits, (!adjacent_segments) && double_anchor_noncanonical, flank_seq);

	if(splice_fixed)
	{
		int matchPrefixLength = prefix_length - buffer;
		int spliceLength = anchorMapPos - mainMapPos;
		int matchSuffixLength = tailLength + buffer - prefix_length;
		Jump_Code matchPrefixJumpCode(matchPrefixLength, "M");
		Jump_Code spliceJumpCode(spliceLength, "N");
		Jump_Code matchSuffixJumpCode(matchSuffixLength, "M");
		tailCigarInfo->jump_code.push_back(matchPrefixJumpCode);
		tailCigarInfo->jump_code.push_back(spliceJumpCode);
		tailCigarInfo->jump_code.push_back(matchSuffixJumpCode);

	    if((!BUILD_SPLICE_HASH_FROM_FILE) && (BUILD_SPLICE_HASH_FROM_FIRST_MAPPING))
	    {
	    	unsigned int spliceEndPosInWholeGenome = anchorMapPos + readLength - matchSuffixLength;
	    	unsigned int spliceJunctionStartPos, spliceJunctionEndPos; 
	    	unsigned int spliceJunctionChrInt;
	    	indexInfo->getChrLocation(spliceEndPosInWholeGenome, &spliceJunctionChrInt, &spliceJunctionEndPos);
	    	(spliceJunctionStartPos) = (spliceJunctionEndPos) - spliceLength - 1;
			bool insertion = insertSpliceJunction2ReverseHashForHead(spliceJunctionReverse, spliceJunctionChrInt, 
				spliceJunctionStartPos, spliceJunctionEndPos)
			&& insertSpliceJunction2NormalHashForTail(spliceJunctionNormal, spliceJunctionChrInt, 
				spliceJunctionStartPos, spliceJunctionEndPos);

			if(PRINT_JUNC)
			{
				insertSpliceJunction2HashForPrint(spliceJunctionChrInt, 
				spliceJunctionStartPos, spliceJunctionEndPos);				
				//cout << (spliceJunctionChrInt) << " " << (spliceJunctionStartPos) << " " << (spliceJunctionEndPos) << " case 3"<<endl; 
			}
	    }

		//debugln(" "); //debug("flankSeq = "); //debugln(flank_seq);
	}
	delete genome_scan;
	return splice_fixed;
}

bool fixSpliceGapTail(Splice_Info* tailCigarInfo, 
	int anchorLength, unsigned int mainEndLocInRead, 
	unsigned int anchorMapPos, unsigned int mainMapPos, char* read, char* chrom, 
	const string& readString, int readLength, Index_Info* indexInfo
	)
{
	int buffer = FIX_SPLICE_BUFFER;
	//int readLength = READ_LENGTH;
	int tailLength = readLength - mainEndLocInRead;
	bool splice_fixed = false;
	int readSeqLeftBound = mainEndLocInRead + 1 -buffer;
	int readSeqRightBound = readLength - anchorLength + buffer;
	int readFixSeqLength = readSeqRightBound - readSeqLeftBound + 1;
	unsigned int chromSeqLeftBound = mainMapPos + mainEndLocInRead -buffer;
	unsigned int chromSeqRightBound = anchorMapPos + readLength - anchorLength - 1 + buffer;
	int chromFixSeqLength = readFixSeqLength + 2;

	string pending_seq = readString.substr(readSeqLeftBound-1, readFixSeqLength);
	////debugln("pending_seq = "); //debugln(pending_seq);
	string left_chrom_seq = chromString.substr(chromSeqLeftBound-1, chromFixSeqLength);
	////debugln("left_chrom = "); //debugln(left_chrom_seq);
	string right_chrom_seq = chromString.substr(chromSeqRightBound-chromFixSeqLength, chromFixSeqLength);	
	////debugln("right_chrom = "); //debugln(right_chrom_seq);
	size_t prefix_length = 0;
	size_t max_double_splice_mismatch = 2;
	size_t mismatch_bits = 0;
	bool adjacent_segments = false;
	bool double_anchor_noncanonical = DO_NONCANONICAL; ////debug 
	string flank_seq;
	GenomeScan* genome_scan = new GenomeScan;
	splice_fixed = (*genome_scan).Double_anchored_score
	(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_double_splice_mismatch, 
	mismatch_bits, (!adjacent_segments) && double_anchor_noncanonical, flank_seq);
	////debugln("prefix_length = "); //debugln((int)prefix_length);
	if(splice_fixed)
	{
		int matchPrefixLength = prefix_length - buffer;
		int spliceLength = anchorMapPos - mainMapPos;
		int matchSuffixLength = tailLength + buffer - prefix_length;
		Jump_Code matchPrefixJumpCode(matchPrefixLength, "M");
		Jump_Code spliceJumpCode(spliceLength, "N");
		Jump_Code matchSuffixJumpCode(matchSuffixLength, "M");
		tailCigarInfo->jump_code.push_back(matchPrefixJumpCode);
		tailCigarInfo->jump_code.push_back(spliceJumpCode);
		tailCigarInfo->jump_code.push_back(matchSuffixJumpCode);

	    if((!BUILD_SPLICE_HASH_FROM_FILE) && (BUILD_SPLICE_HASH_FROM_FIRST_MAPPING))
	    {
	    	unsigned int spliceEndPosInWholeGenome = anchorMapPos + readLength - matchSuffixLength;
	    	unsigned int spliceJunctionStartPos, spliceJunctionEndPos; 
	    	unsigned int spliceJunctionChrInt;
	    	indexInfo->getChrLocation(spliceEndPosInWholeGenome, &spliceJunctionChrInt, &spliceJunctionEndPos);
	    	(spliceJunctionStartPos) = (spliceJunctionEndPos) - spliceLength - 1;
			bool insertion = insertSpliceJunction2ReverseHashForHead(spliceJunctionReverse, spliceJunctionChrInt, 
				spliceJunctionStartPos, spliceJunctionEndPos)
			&& insertSpliceJunction2NormalHashForTail(spliceJunctionNormal, spliceJunctionChrInt, 
				spliceJunctionStartPos, spliceJunctionEndPos);

			if(PRINT_JUNC)
			{
				insertSpliceJunction2HashForPrint(spliceJunctionChrInt, 
				spliceJunctionStartPos, spliceJunctionEndPos);				
				//cout << (spliceJunctionChrInt) << " " << (spliceJunctionStartPos) << " " << (spliceJunctionEndPos) << " case 4"<< endl; 
			}
	    }
		////debugln(" "); //debug("flankSeq = "); //debugln(flank_seq);
	}
	delete genome_scan;
	return splice_fixed;	
}



bool fixDoubleAnchorTail(int relation, Splice_Info* tailCigarInfo,
	unsigned int anchorMapPos, unsigned int mainMapPos,
	unsigned int mainEndLocInRead, int anchorLength,
	char* read, char* chrom, const string& readString, int readLength, Index_Info* indexInfo)
{
	//int readLength = READ_LENGTH;
	if((relation == FIX_TOO_CLOSE) || (relation == FIX_TOO_FAR) || (relation == FIX_NO_RELATIONSHIP))
	{
		return false;
	}
	else if(relation == FIX_REMAPPING_SHORT_TAIL)
	{
		int spliceDistance = mainEndLocInRead;
		int middlePartExtendLength = anchorLength;
		int finalTailLength = anchorMapPos;
		Jump_Code midPartExtendJumpCode(middlePartExtendLength, "M");
		Jump_Code spliceDistanceJumpCode(spliceDistance, "N");
		Jump_Code tailMatchJumpCode(finalTailLength, "M");
		tailCigarInfo->jump_code.push_back(midPartExtendJumpCode);
		tailCigarInfo->jump_code.push_back(spliceDistanceJumpCode);
		tailCigarInfo->jump_code.push_back(tailMatchJumpCode);
		return true;
	}
	else if(relation == FIX_SOFTCLIPPING)
	{
		int tailMatchLength = readLength - mainEndLocInRead;
		Jump_Code tailMatchJumpCode(tailMatchLength, "S");
		tailCigarInfo->jump_code.push_back(tailMatchJumpCode);
		return true;		
	}
	else if(relation == FIX_MATCH)
	{
		int tailMatchLength = readLength - mainEndLocInRead;
		////debug("readLength = "); //debugln(readLength);
		////debug("mainEndLocInRead = "); //debugln(mainEndLocInRead);		
		////debug("tailMatchLength = "); //debugln(tailMatchLength);
		Jump_Code tailMatchJumpCode(tailMatchLength, "M");
		tailCigarInfo->jump_code.push_back(tailMatchJumpCode);
		return true;
	}
	else if(relation == FIX_INSERTION_NEIGHBOUR)
	{
		int insertionLength = mainMapPos - anchorMapPos;
		int tailLength = readLength - mainEndLocInRead;
		Jump_Code insertionJumpCode(insertionLength, "I");
		Jump_Code matchJumpCode(tailLength-insertionLength, "M");
		tailCigarInfo->jump_code.push_back(insertionJumpCode);
		tailCigarInfo->jump_code.push_back(matchJumpCode);
		return true;
	}
	else if(relation == FIX_INSERTION_GAP)
	{
		return fixInsertionGapTail(tailCigarInfo, 
			anchorLength, mainEndLocInRead, anchorMapPos, mainMapPos,
			read, chrom,//chromString, 
			readString, readLength
			);
	}
	else if(relation == FIX_DELETION_NEIGHBOUR)
	{
		int deletionLength = anchorMapPos- mainMapPos;
		int tailLength = readLength - mainEndLocInRead + 1;
		Jump_Code deletionJumpCode(deletionLength, "D");
		Jump_Code matchJumpCode(tailLength, "M");
		tailCigarInfo->jump_code.push_back(deletionJumpCode);
		tailCigarInfo->jump_code.push_back(matchJumpCode);
		return true;
	}
	else if(relation == FIX_DELETION_GAP)
	{
		return fixDeletionGapTail(tailCigarInfo, 
			anchorLength, mainEndLocInRead, anchorMapPos, mainMapPos,
			read, chrom, readString, readLength);
	}
	else if(relation == FIX_SPLICE_NEIGHBOUR)
	{
		/*
		unsigned int spliceLength = anchorMapPos - mainMapPos;
		int tailLength = readLength - mainEndLocInRead;
		Jump_Code spliceJumpCode(spliceLength, "N");
		Jump_Code matchJumpCode(tailLength, "M");
		tailCigarInfo->jump_code.push_back(spliceJumpCode);
		tailCigarInfo->jump_code.push_back(matchJumpCode);*/
		return fixSpliceNeighbourTail(tailCigarInfo, 
			anchorLength, mainEndLocInRead, anchorMapPos, mainMapPos,
			read, chrom, readString, readLength, indexInfo);
	}
	else if(relation == FIX_SPLICE_GAP)
	{
		return fixSpliceGapTail(tailCigarInfo, 
			anchorLength, mainEndLocInRead, anchorMapPos, mainMapPos,
			read, chrom, readString, readLength, indexInfo);
	}
	else
	{
		////debugln("fix double anchor tail error!!!");
		return false;
	}
}



int checkRelationTail(int mainEndLocInRead, unsigned int anchorMapLoc, unsigned int mainMapLoc, int anchorLength, int readLength, Index_Info* indexInfo)
{
	//int readLength = READ_LENGTH;
	
	if((anchorMapLoc>(indexInfo->indexSize))||(mainMapLoc>(indexInfo->indexSize)))
	{
		////debugln("mapPos < 0 or > (indexInfo->indexSize)-- anchorMapLoc: ");
		return FIX_NO_RELATIONSHIP;
	}	
	else if (anchorMapLoc == mainMapLoc)
		return FIX_MATCH;
	else if (mainMapLoc > anchorMapLoc)
	{
		if (mainMapLoc - anchorMapLoc > MAX_INSERTION_LENGTH)
			return FIX_TOO_CLOSE;
		else
		{
			if((mainEndLocInRead + anchorLength) == readLength)
				return FIX_INSERTION_NEIGHBOUR;
			else
				return FIX_INSERTION_GAP;
		}
	}
	else if (mainMapLoc < anchorMapLoc)
	{
		unsigned int mapDistance = anchorMapLoc - mainMapLoc;
		if(mapDistance > MAX_SPLICE_LENGTH)
			return FIX_TOO_FAR;
		else if(mapDistance > MAX_DELETION_LENGTH)
		{
			if((mainEndLocInRead + anchorLength) == readLength)
				return FIX_SPLICE_NEIGHBOUR;
			else
				return FIX_SPLICE_GAP;			
		}
		else
		{
			if((mainEndLocInRead + anchorLength) == readLength)
				return FIX_DELETION_NEIGHBOUR;
			else
				return FIX_DELETION_GAP;				
		}
	}
	else
	{
		//debugln("check tail relation error!!!");
		return false;
	}
}

