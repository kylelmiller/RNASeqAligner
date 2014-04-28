int checkRelation(unsigned int segStartNum2, unsigned int segEndNum1, unsigned int alignLoc1, unsigned int alignLoc2)
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
}

bool fixInsertionGap(int firstFragNo, int secondFragNo, unsigned int* segMapRangeStart, unsigned int* segMapRangeEnd, 
	unsigned int* segmentLocInRead, unsigned int* segMapLoc, char* read, char* chrom, int segmentNum, 
	Splice_Info* cigarInfo,
	//string chromString, 
	const string& readString, int readLength, Index_Info* indexInfo
	)
{
	bool insertion_fixed = false;
	int firstFragEndSegNo = *(segMapRangeEnd + (firstFragNo-1));
	int secondFragStartSegNo = *(segMapRangeStart + (secondFragNo-1));
	int firstFragStartSegNo = *(segMapRangeStart + (firstFragNo-1));

	unsigned int readSeqLeftBound = segmentLocInRead[firstFragEndSegNo] - 1;
	unsigned int readSeqRightBound = segmentLocInRead[secondFragStartSegNo-1] - 1;
	unsigned int readFixSeqLength = readSeqRightBound - readSeqLeftBound + 1;
	unsigned int chromSeqLeftBound = segMapLoc[firstFragNo-1] + segmentLocInRead[firstFragStartSegNo-1]
	 	+ getFragmentLength(firstFragNo, segmentNum, segmentLocInRead, segMapRangeStart, segMapRangeEnd, readLength) - 1;
	unsigned int chromSeqRightBound = segMapLoc[secondFragNo-1] + segmentLocInRead[secondFragStartSegNo-1] - 2;
	
	unsigned int chromFixSeqLength = chromSeqRightBound - chromSeqLeftBound + 1;
	if(chromSeqLeftBound >= chromSeqRightBound)
	{
		int insertionLength = readFixSeqLength + chromSeqLeftBound - chromSeqRightBound - 1;
		int matchSuffixLength = insertionLength - readFixSeqLength;
		if(matchSuffixLength < 0)
		{
			////debug("matchSuffixLength < 0 = "); //debugln(matchSuffixLength);
			//return false;
			#ifdef DEBUG
			cout << "match Suffix length < 0" << endl;
			#endif
		}
		//Jump_Code matchPrefixJumpCode(matchPrefixLength, "M");
		Jump_Code insertionJumpCode(insertionLength, "I");
		Jump_Code matchSuffixJumpCode(matchSuffixLength, "M");
		//headCigarInfo->jump_code.push_back(matchPrefixJumpCode);
		cigarInfo->jump_code.push_back(insertionJumpCode);
		cigarInfo->jump_code.push_back(matchSuffixJumpCode);		
		return true;				
	}
	else
	{
		if((readSeqLeftBound-1 > (indexInfo->indexSize)) //|| (readSeqRightBound - 1 > readLength)
		) // Xinan: two short in one side, need to be fixed later;
		{
			#ifdef DEBUG
			cout << "two short in one side" << readSeqLeftBound << " " << readSeqRightBound << endl;
			#endif
			return false;
		}
		string pending_seq_ins = readString.substr(readSeqLeftBound-1, readFixSeqLength); 
		string chrom_seq_ins = chromString.substr(chromSeqLeftBound-1, chromFixSeqLength); 

		size_t prefix_length = 0;
		size_t mismatch_bits = 0; //?
		size_t max_ins_mismatch = 2;

		GenomeScan* genome_scan = new GenomeScan;
		insertion_fixed = (*genome_scan).Double_anchored_score_ins
		(pending_seq_ins, chrom_seq_ins, max_ins_mismatch, prefix_length, mismatch_bits); //X: fix insertion

		if(insertion_fixed)
		{
			int frag2Length = getFragmentLength(secondFragNo, segmentNum, segmentLocInRead, segMapRangeStart, segMapRangeEnd, readLength);
			int insertionLength = segMapLoc[firstFragNo-1] - segMapLoc[secondFragNo-1];	
			Jump_Code firstMatchJumpCode(prefix_length, "M");
			Jump_Code insertionJumpCode(insertionLength, "I");
			Jump_Code secondMatchJumpCode(frag2Length + readFixSeqLength - prefix_length - insertionLength, "M");
			cigarInfo->jump_code.push_back(firstMatchJumpCode);
			cigarInfo->jump_code.push_back(insertionJumpCode);
			cigarInfo->jump_code.push_back(secondMatchJumpCode);
		}
		//delete genome_scan;
		return insertion_fixed;
	}
} 

bool fixDeletionGap(int firstFragNo, int secondFragNo, unsigned int* segMapRangeStart, 
	unsigned int* segMapRangeEnd, unsigned int* segmentLocInRead,
	unsigned int* segMapLoc, char* read, char* chrom, int segmentNum, Splice_Info* cigarInfo,
	const string& readString, int readLength, Index_Info* indexInfo)
{
	bool deletion_fixed = false;
	int firstFragEndSegNo = *(segMapRangeEnd + (firstFragNo-1));
	int secondFragStartSegNo = *(segMapRangeStart + (secondFragNo-1));
	int firstFragStartSegNo = *(segMapRangeStart + (firstFragNo-1));
	unsigned int readSeqLeftBound = segmentLocInRead[firstFragEndSegNo] - 1;
	unsigned int readSeqRightBound = segmentLocInRead[secondFragStartSegNo-1] - 1;
	unsigned int readFixSeqLength = readSeqRightBound - readSeqLeftBound + 1;

	if((readSeqLeftBound-1 > (indexInfo->indexSize)) //|| (readSeqRightBound - 1 > readLength)
	) // Xinan: two short in one side, need to be fixed later;
	{
		#ifdef DEBUG
		cout << "two short in one side" << readSeqLeftBound << " " << readSeqRightBound << endl;
		#endif
		return false;
	}

	unsigned int chromSeqLeftBound = segMapLoc[firstFragNo-1] + segmentLocInRead[firstFragStartSegNo-1]
	 	+ getFragmentLength(firstFragNo, segmentNum, segmentLocInRead, segMapRangeStart, segMapRangeEnd, readLength) - 1;
	unsigned int chromSeqRightBound = segMapLoc[secondFragNo-1] + segmentLocInRead[secondFragStartSegNo-1] - 2;	
	unsigned int chromFixSeqLength = readFixSeqLength+2;

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
		int frag2Length = getFragmentLength(secondFragNo, segmentNum, segmentLocInRead, segMapRangeStart, segMapRangeEnd, readLength);
		int deletionLength = segMapLoc[secondFragNo-1] - segMapLoc[firstFragNo-1];	
		int matchSuffixLength = frag2Length + readFixSeqLength - prefix_length;
		Jump_Code firstMatchJumpCode(prefix_length, "M");
		Jump_Code deletionJumpCode(deletionLength, "D");
		Jump_Code secondMatchJumpCode(matchSuffixLength, "M");
		cigarInfo->jump_code.push_back(firstMatchJumpCode);
		cigarInfo->jump_code.push_back(deletionJumpCode);
		cigarInfo->jump_code.push_back(secondMatchJumpCode);
	}
	//delete genome_scan;
	return deletion_fixed;
}

bool fixSpliceNeighbour(int firstFragNo, int secondFragNo, unsigned int* segMapRangeStart, 
	unsigned int* segMapRangeEnd, unsigned int* segmentLocInRead,
	unsigned int* segMapLoc, char* read, char* chrom, int segmentNum, Splice_Info* cigarInfo,
	const string& readString, int readLength, Index_Info* indexInfo)
{
	#ifdef DEBUG
	cout << "start fixSpliceNeighbour function ..." << endl;
	#endif
	bool splice_fixed = false;
	int buffer = FIX_SPLICE_BUFFER;
	int firstFragEndSegNo = *(segMapRangeEnd + (firstFragNo-1));
	int secondFragStartSegNo = firstFragEndSegNo + 1;
	int firstFragStartSegNo = *(segMapRangeStart + (firstFragNo-1));
	
	//int tmpReadSeqLeftBound = segmentLocInRead[firstFragEndSegNo] - 1 - buffer;
	//int tmpReadSeqRightBound = segmentLocInRead[secondFragStartSegNo-1] - 1 + buffer;
	//if(tmpReadSeqLeftBound )
	
	unsigned int readSeqLeftBound = segmentLocInRead[firstFragEndSegNo] - 1 - buffer;
	unsigned int readSeqRightBound = segmentLocInRead[secondFragStartSegNo-1] - 1 + buffer;
	unsigned int readFixSeqLength = buffer*2 + 1;
	
	if((readSeqLeftBound-1 > (indexInfo->indexSize)) //|| (readSeqRightBound - 1 > readLength)
	) // Xinan: two short in one side, need to be fixed later;
	{
		#ifdef DEBUG
		cout << "two short in one side" << readSeqLeftBound << " " << readSeqRightBound << endl;
		#endif
		return false;
	}

	int mid_seq_length = 1;
	unsigned int chromSeqLeftBound = segMapLoc[firstFragNo-1] + segmentLocInRead[firstFragStartSegNo-1]
	 	+ getFragmentLength(firstFragNo, segmentNum, segmentLocInRead, segMapRangeStart, segMapRangeEnd, readLength) 
	 	- 1 -buffer;
	unsigned int chromSeqRightBound = segMapLoc[secondFragNo-1] + segmentLocInRead[secondFragStartSegNo-1] - 2 
		+ buffer;	

	if((chromSeqLeftBound-1 > (indexInfo->indexSize)) //|| (chromSeqRightBound-1 > (indexInfo->indexSize))
	) // Xinan: two short in one side, need to be fixed later;
	{
		#ifdef DEBUG
		cout << "two short in one side " << chromSeqLeftBound << " " << chromSeqRightBound << endl;
		#endif
		return false;
	}


	
	unsigned int chromFixSeqLength = readFixSeqLength+2;
	#ifdef DEBUG 
	cout << "start to extract two read and chrom strings  " << endl;
	#endif
	string pending_seq = readString.substr(readSeqLeftBound-1, readFixSeqLength);
	string left_chrom_seq = chromString.substr(chromSeqLeftBound-1, chromFixSeqLength);
	string right_chrom_seq = chromString.substr(chromSeqRightBound-chromFixSeqLength, chromFixSeqLength);	

	#ifdef DEBUG
	cout << "start Double_anchored_score function ...... " << endl;
	#endif
	size_t prefix_length = 0;
	size_t max_double_splice_mismatch = 2;
	size_t mismatch_bits = 0;
	bool adjacent_segments = false;
	bool double_anchor_noncanonical = false;//DO_NONCANONICAL; ////debug
	string flank_seq;
	GenomeScan* genome_scan = new GenomeScan;
	splice_fixed = (*genome_scan).Double_anchored_score
	(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_double_splice_mismatch, 
	mismatch_bits, (!adjacent_segments) && double_anchor_noncanonical, flank_seq);
	/*
	if(splice_fixed && (flank_seq == "ATAC" || 
		flank_seq == "GTAT" || flank_seq == "CTGC" || 
		flank_seq == "GCAG") && ((int)(right_bound - left_bound - pending_seq.length() + 1) > max_semicanonical_intron_length))
					splice_fixed = false;
	*/

	if(splice_fixed)
	{
		#ifdef DEBUG
		cout << "splice_fixed for fixSpliceNeighbour function ... " << endl;
		#endif
		int frag2Length = getFragmentLength(secondFragNo, segmentNum, segmentLocInRead, segMapRangeStart, segMapRangeEnd, readLength);
		int spliceLength = segMapLoc[secondFragNo-1] - segMapLoc[firstFragNo-1];
		int matchSuffixLength = frag2Length + 1 + buffer - prefix_length ;	
		int matchPrefixLength = prefix_length - buffer;
		Jump_Code firstMatchJumpCode(matchPrefixLength, "M");
		Jump_Code spliceJumpCode(spliceLength, "N");
		Jump_Code secondMatchJumpCode(matchSuffixLength, "M");
		cigarInfo->jump_code.push_back(firstMatchJumpCode);
		cigarInfo->jump_code.push_back(spliceJumpCode);
		cigarInfo->jump_code.push_back(secondMatchJumpCode);

	    if((!BUILD_SPLICE_HASH_FROM_FILE) && (BUILD_SPLICE_HASH_FROM_FIRST_MAPPING))
	    {

	    	unsigned int spliceStartPosInWholeGenome = segMapLoc[firstFragNo-1] + segmentLocInRead[firstFragEndSegNo] - 2 + matchPrefixLength - 1; 
	    	//if(PRINT_JUNC)
	    	//{
	    	//	cout << "spliceStartPosInWholeGenome : " << spliceStartPosInWholeGenome << endl;
	    	//}
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
				//cout << (spliceJunctionChrInt) << " " << (spliceJunctionStartPos) << " " << (spliceJunctionEndPos) << " case 5"<< endl; 
			}
	    }
		////debugln(" "); //debug("flankSeq = "); //debugln(flank_seq);
	}
	delete genome_scan;
	return splice_fixed;	
}

bool fixSpliceGap(int firstFragNo, int secondFragNo, unsigned int* segMapRangeStart, 
	unsigned int* segMapRangeEnd, unsigned int* segmentLocInRead,
	unsigned int* segMapLoc, char* read, char* chrom, int segmentNum, Splice_Info* cigarInfo,
	const string& readString, int readLength, Index_Info* indexInfo//, string strand
	)
{
	#ifdef DEBUG
	cout << "start fixSpliceGap function ... " << endl; 
	#endif
	bool splice_fixed = false;
	int buffer = FIX_SPLICE_BUFFER;
	int firstFragEndSegNo = *(segMapRangeEnd + (firstFragNo-1));
	int secondFragStartSegNo = *(segMapRangeStart + (secondFragNo-1));
	int firstFragStartSegNo = *(segMapRangeStart + (firstFragNo-1));
	unsigned int readSeqLeftBound = segmentLocInRead[firstFragEndSegNo] - 1 - buffer;
	unsigned int readSeqRightBound = segmentLocInRead[secondFragStartSegNo-1] - 1 + buffer;
	unsigned int readFixSeqLength = readSeqRightBound - readSeqLeftBound + 1;
	int mid_seq_length = segmentLocInRead[secondFragStartSegNo-1] - segmentLocInRead[firstFragEndSegNo] + 1;
	unsigned int chromSeqLeftBound = segMapLoc[firstFragNo-1] + segmentLocInRead[firstFragStartSegNo-1]
	 	+ getFragmentLength(firstFragNo, segmentNum, segmentLocInRead, segMapRangeStart, segMapRangeEnd, readLength) - 1 -buffer;
	unsigned int chromSeqRightBound = segMapLoc[secondFragNo-1] + segmentLocInRead[secondFragStartSegNo-1] - 2 + buffer;	
	
	unsigned int chromFixSeqLength = readFixSeqLength+2;

	if((readSeqLeftBound-1 > (indexInfo->indexSize)) //|| (readSeqRightBound - 1 > readLength)
	) // Xinan: two short in one side, need to be fixed later;
	{
		#ifdef DEBUG
		cout << "too short in one side" << readSeqLeftBound << " " << readSeqRightBound << endl;
		#endif
		return false;
	}

	string pending_seq = readString.substr(readSeqLeftBound-1, readFixSeqLength);
	string left_chrom_seq = chromString.substr(chromSeqLeftBound-1, chromFixSeqLength);
	string right_chrom_seq = chromString.substr(chromSeqRightBound-chromFixSeqLength, chromFixSeqLength);	

	size_t prefix_length = 0;
	size_t max_double_splice_mismatch = 2;
	size_t mismatch_bits = 0;
	bool adjacent_segments = false;
	bool double_anchor_noncanonical = false;//DO_NONCANONICAL; ////debug
	string flank_seq;
	GenomeScan* genome_scan = new GenomeScan;
	splice_fixed = (*genome_scan).Double_anchored_score
	(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_double_splice_mismatch, 
	mismatch_bits, (!adjacent_segments) && double_anchor_noncanonical, flank_seq);
	/*
	if(splice_fixed && (flank_seq == "ATAC" || 
		flank_seq == "GTAT" || flank_seq == "CTGC" || 
		flank_seq == "GCAG") && ((int)(right_bound - left_bound - pending_seq.length() + 1) > max_semicanonical_intron_length))
					splice_fixed = false;
	*/
	int frag2Length = getFragmentLength(secondFragNo, segmentNum, segmentLocInRead, segMapRangeStart, segMapRangeEnd, readLength);
	if(splice_fixed)
	{
		//int frag2Length = getFragmentLength(secondFragNo, segmentNum, segmentLocInRead, segMapRangeStart, segMapRangeEnd);
		int spliceLength = segMapLoc[secondFragNo-1] - segMapLoc[firstFragNo-1];
		int matchSuffixLength = frag2Length + mid_seq_length + buffer - prefix_length ;	
		int matchPrefixLength = prefix_length - buffer;
		Jump_Code firstMatchJumpCode(matchPrefixLength, "M");
		Jump_Code spliceJumpCode(spliceLength, "N");
		Jump_Code secondMatchJumpCode(matchSuffixLength, "M");
		cigarInfo->jump_code.push_back(firstMatchJumpCode);
		cigarInfo->jump_code.push_back(spliceJumpCode);
		cigarInfo->jump_code.push_back(secondMatchJumpCode);

	    if((!BUILD_SPLICE_HASH_FROM_FILE) && (BUILD_SPLICE_HASH_FROM_FIRST_MAPPING))
	    {
	    	unsigned int spliceStartPosInWholeGenome = segMapLoc[firstFragNo-1] + segmentLocInRead[firstFragEndSegNo] - 2 + matchPrefixLength - 1; 
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
				//cout << (spliceJunctionChrInt) << " " << (spliceJunctionStartPos) << " " << (spliceJunctionEndPos) << " case 5"<< endl; 
			}				
			/*
			if(!insertion)
			{
				cout << "insertion failed " << endl;
				cout << chrInt << " " << spliceStartPos << " " << spliceEndPos << endl;
			}*/
	    }
	}
	else //fix multi-splice junction
	{
		/*#ifdef DEBUG
		cout << "fix spliceGap failed .... may try fixing multi-splice ....." << endl;
		#endif

		//cout << "readNO: " << read_num << endl;

		int firstFragLength = getFragmentLength(firstFragNo, segmentNum, segmentLocInRead, segMapRangeStart, segMapRangeEnd, readLength); 
		int secondFragLength = getFragmentLength(secondFragNo, segmentNum, segmentLocInRead, segMapRangeStart, segMapRangeEnd, readLength); 
		int firstFragEndPosInRead = segmentLocInRead[firstFragEndSegNo] - 2; // one base jump 
		int secondFragStartPosInRead = segmentLocInRead[secondFragStartSegNo-1];
		int firstFragEndMapPos = segMapLoc[firstFragNo-1] + firstFragEndPosInRead - 1;
		int secondFragStartMapPos = segMapLoc[secondFragNo-1] + secondFragStartPosInRead - 1;

		//cout << "firstFragLength: " << firstFragLength << endl;
		//cout << "secondFragLength: " << secondFragLength << endl;
		//cout << "firstFragEndPosInRead: " << firstFragEndPosInRead << endl;
		//cout << "secondFragStartPosInRead: " << secondFragStartPosInRead << endl;
		//cout << "firstFragEndMapPos: " << firstFragEndMapPos << endl;
		//cout << "secondFragStartMapPos: " << secondFragStartMapPos << endl;

		Small_Exon* smallExon = new Small_Exon(firstFragLength, secondFragLength, firstFragEndPosInRead, secondFragStartPosInRead,
							firstFragEndMapPos, secondFragStartMapPos, chromString, readString);
		
		unsigned int firstFragEndMapPosInChr, smallExonMapChromNameInt;
		indexInfo->getChrLocation(firstFragEndMapPos, &smallExonMapChromNameInt, &firstFragEndMapPosInChr);
		int secondFragStartMapPosInChr = secondFragStartMapPos + firstFragEndMapPosInChr - firstFragEndMapPos;

		string targetMappingExonWithFlankString;
		int tmpSmallExonStartPosInRead;
		int tmpSmallExonEndPosInRead;


		int PossibleSmallExonNum = 0;

		//cout << "start to search for GTAG SJ .... .... " << endl;

		for(int tmp = 0; tmp < (smallExon->possiSJposInReadGTAG).size(); tmp++)
		{

			if(PossibleSmallExonNum > 1)
			{
				break;
			}


			tmpSmallExonStartPosInRead = (smallExon->possiSJposInReadGTAG)[tmp].first;
			tmpSmallExonEndPosInRead = (smallExon->possiSJposInReadGTAG)[tmp].second;
			//cout << "GTAG small Exon from: " << tmpSmallExonStartPosInRead << " to: " << tmpSmallExonEndPosInRead << endl; 
			//cout << "tmpSmallExonStartPosInRead: " << tmpSmallExonStartPosInRead << endl;
			//cout << "tmpSmallExonEndPosInRead: " << tmpSmallExonEndPosInRead << endl;
			
			targetMappingExonWithFlankString = "AG";
			targetMappingExonWithFlankString += 
				readString.substr(tmpSmallExonStartPosInRead-1, tmpSmallExonEndPosInRead-tmpSmallExonStartPosInRead+1);
			targetMappingExonWithFlankString += "GT";


			////////////  search for target string, should be replaced by 2nd-level index mapping method ///////////////
			SBNDM_FORWARD* anchor_search = new SBNDM_FORWARD;
			string strand_anchor = "+";
			int max_anchor_hits = 2;
			string anchor_string = targetMappingExonWithFlankString;
			int left_bound_forward = firstFragEndMapPosInChr + tmpSmallExonStartPosInRead - firstFragEndPosInRead;
			int left_bound_reverse = secondFragStartMapPosInChr + tmpSmallExonStartPosInRead - secondFragStartPosInRead;
			int right_bound = left_bound_reverse; 
			int anchor_pos;
			bool mapped_fw = true;	
			
			(*anchor_search).set(anchor_string, chromStr[smallExonMapChromNameInt], strand_anchor == "+", 
				strand_anchor == "-", max_anchor_hits);

			while((*anchor_search).search_next_simple(chromStr[smallExonMapChromNameInt], left_bound_forward,
				left_bound_reverse, right_bound, anchor_pos, mapped_fw))
			{
				//cout << endl << "... small Exon found !!!" << endl << "anchor_pos = " << anchor_pos << endl;
				//cout << "left_bound_forward = " << left_bound_forward << endl;
				//cout << "left_bound_reverse = " << left_bound_reverse << endl;	
				//cout << "right_bound =        " << right_bound << endl;
				PossibleSmallExonNum ++;
				if(PossibleSmallExonNum > 1)
				{
					break;
				}
				anchor_pos ++;	
				(smallExon->finalSJVec).push_back(pair <int , pair <int, int> >( anchor_pos+2, pair <int, int > (tmpSmallExonStartPosInRead, tmpSmallExonEndPosInRead))); 
			}		

		}

		//cout << "start to search for CTAC SJ .... .... " << endl;

		for(int tmp = 0; tmp < (smallExon->possiSJposInReadCTAC).size(); tmp++)
		{
			if(PossibleSmallExonNum > 1)
			{
				break;
			}

			tmpSmallExonStartPosInRead = (smallExon->possiSJposInReadCTAC)[tmp].first;
			tmpSmallExonEndPosInRead = (smallExon->possiSJposInReadCTAC)[tmp].second;
			//cout << "CTAC small Exon from: " << tmpSmallExonStartPosInRead << " to: " << tmpSmallExonEndPosInRead << endl; 
			//cout << "tmpSmallExonStartPosInRead: " << tmpSmallExonStartPosInRead << endl;
			//cout << "tmpSmallExonEndPosInRead: " << tmpSmallExonEndPosInRead << endl;
			
			targetMappingExonWithFlankString = "AC";
			targetMappingExonWithFlankString += 
				readString.substr(tmpSmallExonStartPosInRead-1, tmpSmallExonEndPosInRead-tmpSmallExonStartPosInRead+1);
			targetMappingExonWithFlankString += "CT";
		
			SBNDM_FORWARD* anchor_search = new SBNDM_FORWARD;
			string strand_anchor = "+";
			int max_anchor_hits = 2;
			string anchor_string = targetMappingExonWithFlankString;
			int left_bound_forward = firstFragEndMapPosInChr + tmpSmallExonStartPosInRead - firstFragEndPosInRead;
			int left_bound_reverse = secondFragStartMapPosInChr + tmpSmallExonStartPosInRead - secondFragStartPosInRead;
			int right_bound = left_bound_reverse; 
			int anchor_pos;
			bool mapped_fw = true;	

			(*anchor_search).set(anchor_string, chromStr[smallExonMapChromNameInt], strand_anchor == "+", 
				strand_anchor == "-", max_anchor_hits);
			
			while((*anchor_search).search_next_simple(chromStr[smallExonMapChromNameInt], left_bound_forward,
				left_bound_reverse, right_bound, anchor_pos, mapped_fw))
			{
				//cout << endl << "... small Exon found !!!" << endl << "anchor_pos = " << anchor_pos << endl;
				//cout << "left_bound_forward = " << left_bound_forward << endl;
				//cout << "left_bound_reverse = " << left_bound_reverse << endl;	
				//cout << "right_bound =        " << right_bound << endl;
				PossibleSmallExonNum ++;
				if(PossibleSmallExonNum > 1)
				{
					break;
				}

				anchor_pos ++;	
				(smallExon->finalSJVec).push_back(pair <int , pair <int, int> >( anchor_pos+2, pair <int, int > (tmpSmallExonStartPosInRead, tmpSmallExonEndPosInRead))); 
			}		

		}

		if(PossibleSmallExonNum == 1)//((smallExon->finalSJVec).size() == 1)
		{
			//cout << "readNO: " << read_num << endl;
			int tmpSmallExonMapPos = (smallExon->finalSJVec)[0].first;
			int tmpSmallExonStartPosInRead = ((smallExon->finalSJVec)[0].second).first;
			int tmpSmallExonEndPosInRead = ((smallExon->finalSJVec)[0].second).second;

			int matchPrefixLength = tmpSmallExonStartPosInRead - firstFragEndPosInRead - 1;
			int firstSpliceLength = tmpSmallExonMapPos - (firstFragEndMapPosInChr + matchPrefixLength) - 1;
			int smallExonLength = tmpSmallExonEndPosInRead - tmpSmallExonStartPosInRead + 1;
			int matchSuffixLength = secondFragLength + secondFragStartPosInRead - 1 - tmpSmallExonEndPosInRead;
			int secondSpliceLength = secondFragStartMapPosInChr - (matchSuffixLength-secondFragLength) - (tmpSmallExonMapPos + smallExonLength - 1) - 1;

			cout << "add JumpCode: " << matchPrefixLength << "M" << firstSpliceLength << "N" << smallExonLength << "M" 
				<< secondSpliceLength << "N" << matchSuffixLength << "M" << endl;
			Jump_Code firstMatchJumpCode(matchPrefixLength, "M");
			Jump_Code firstSpliceJumpCode(firstSpliceLength, "N");
			Jump_Code smallExonJumpCode(smallExonLength, "M");
			Jump_Code secondSpliceJumpCode(secondSpliceLength, "N");
			Jump_Code secondMatchJumpCode(matchSuffixLength, "M");
			

			cigarInfo->jump_code.push_back(firstMatchJumpCode);
			cigarInfo->jump_code.push_back(firstSpliceJumpCode);
			cigarInfo->jump_code.push_back(smallExonJumpCode);
			cigarInfo->jump_code.push_back(secondSpliceJumpCode);
			cigarInfo->jump_code.push_back(secondMatchJumpCode);

			splice_fixed = true;
 		}*/

	}
	delete genome_scan;
	return splice_fixed;
}

bool fixDoubleAnchor(int relation, Splice_Info* cigarInfo, int firstFragNo, int secondFragNo, 
	unsigned int* segMapRangeStart, unsigned int* segMapRangeEnd, unsigned int* segmentLocInRead,
	unsigned int* segMapLoc, char* read, char* chrom, int segmentNum, const string& readString, int readLength,
	Index_Info* indexInfo//, string strand
	)
{
	////debug("//debug -- start fixDoubleAnchor, relation = "); //debugln(relation);
	if((relation == FIX_TOO_CLOSE) || (relation == FIX_TOO_FAR) || (relation == FIX_NO_RELATIONSHIP))
	{
		////debugln(" No relation or too close or too far!!! ");
		return false;
	}
	else if (relation == FIX_INSERTION_NEIGHBOUR)//need //debug, for the base between those two segments
	{
		int frag2Length = 
		getFragmentLength(secondFragNo, segmentNum, segmentLocInRead, segMapRangeStart, segMapRangeEnd, readLength);
		int insertionLength = segMapLoc[firstFragNo-1] - segMapLoc[secondFragNo-1];
		Jump_Code insertionJumpCode(insertionLength, "I");
		Jump_Code matchJumpCode(frag2Length+1-insertionLength, "M");
		cigarInfo->jump_code.push_back(insertionJumpCode);
		cigarInfo->jump_code.push_back(matchJumpCode);
		return true;
	}
	else if (relation == FIX_INSERTION_GAP)
	{
		return fixInsertionGap(firstFragNo, secondFragNo, segMapRangeStart, segMapRangeEnd, segmentLocInRead, segMapLoc, 
			read, chrom, segmentNum, cigarInfo, //chromString, 
			readString, readLength, indexInfo
			);
	}
	else if (relation == FIX_DELETION_NEIGHBOUR)
	{
		int frag2Length = 
		getFragmentLength(secondFragNo, segmentNum, segmentLocInRead, segMapRangeStart, segMapRangeEnd, readLength);
		int deletionLength = segMapLoc[secondFragNo-1] - segMapLoc[firstFragNo-1];
		Jump_Code deletionJumpCode(deletionLength, "D");
		Jump_Code matchJumpCode(frag2Length+1, "M");
		cigarInfo->jump_code.push_back(deletionJumpCode);
		cigarInfo->jump_code.push_back(matchJumpCode);
		return true;
	}
	else if (relation == FIX_DELETION_GAP)
	{
		////debugln("fix_deletion_gap ");
		return fixDeletionGap(firstFragNo, secondFragNo, segMapRangeStart, segMapRangeEnd, segmentLocInRead, segMapLoc, 
			read, chrom, segmentNum, cigarInfo, //chromString, 
			readString, readLength, indexInfo
			);			
	}
	else if (relation == FIX_SPLICE_NEIGHBOUR)
	{
		/*
		int frag2Length = 
		getFragmentLength(secondFragNo, segmentNum, segmentLocInRead, segMapRangeStart, segMapRangeEnd);
		int deletionLength = segMapLoc[secondFragNo-1] - segMapLoc[firstFragNo-1];
		Jump_Code deletionJumpCode(deletionLength, "N");
		Jump_Code matchJumpCode(frag2Length+1, "M");
		cigarInfo->jump_code.push_back(deletionJumpCode);
		cigarInfo->jump_code.push_back(matchJumpCode);*/
		return fixSpliceNeighbour(firstFragNo, secondFragNo, segMapRangeStart, segMapRangeEnd, segmentLocInRead, segMapLoc, 
			read, chrom, segmentNum, cigarInfo, readString, readLength, indexInfo);
	}
	else if (relation == FIX_SPLICE_GAP)
	{
		return fixSpliceGap(firstFragNo, secondFragNo, segMapRangeStart, segMapRangeEnd, segmentLocInRead, segMapLoc, 
			read, chrom, segmentNum, cigarInfo, readString, readLength, indexInfo//, strand
			);
	}
	else
	{
		////debugln("wrong relation");
		return FIX_NO_RELATIONSHIP;
	}
}