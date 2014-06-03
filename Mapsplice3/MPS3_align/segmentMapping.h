/* FIXME - THIS NEEDS TO BE SORTED OUT LATER KLM 5/29/14
bool mapMainSecondLevelForTargetMapping_compressedIndex(char *read,
	unsigned int* sa,
	BYTE* lcpCompress, unsigned int* child, BYTE* verifyChild,
	char* chrom, int readLength, unsigned int indexSize,
	unsigned int midPartMapPosForLongHead,
	unsigned int midPartMapPosForLongHeadInSecondLevelIndex,
	int* targetMappingAlignNum, unsigned int* targetMappingAlignLoc,
	Index_Info* indexInfo)
{
	//cout << "start target mapping " << endl;
	//cout << "read: " << endl;
	//cout << read << endl;
    unsigned int norSegmentNum = 0;
   	unsigned int *norSegmentLength = (unsigned int*)malloc(SEGMENTNUM * sizeof(unsigned int));
   	unsigned int *norSegmentLocInRead = (unsigned int*)malloc(SEGMENTNUM * sizeof(unsigned int));
   	unsigned int *norSegmentAlignNum = (unsigned int*)malloc(SEGMENTNUM * sizeof(unsigned int));
   	//Denote the start location of an alignment for a segment, start from 1;
   	unsigned int *norSegmentAlignLoc = (unsigned int*)malloc(SEGMENTNUM * CANDALILOC * sizeof(unsigned int));
   	unsigned int valLength;
   	//cout << "finish to set space " << endl;
	bool mapMain = false;
	unsigned int stop_loc = 0; // location in one segment for iterations
	unsigned int stop_loc_overall = 0; //location in the whole read for iterations
	unsigned int segment_num = 0;
	unsigned int segment_length = 0;
	unsigned int segment_length_max = 0;//used to compare with segment_length for eache segment to get the maximum length segment
	unsigned int segment_align_SArange[2] = {0,0};//legal align location in SA
	unsigned int segment_align_rangeNum = 0;
	unsigned int read_length = readLength; //READ_LENGTH;
	unsigned int interval_begin, interval_end;
	unsigned int n = indexInfo->indexSize;//size of SA
	unsigned int norAlignLoc;
	unsigned int align_chr_location;
	unsigned int align_chr;
	valLength = 0;
	char* read_local = read;
	//cout << "start to check every base......" << endl;
	while (stop_loc_overall < read_length) //- 15)
	{
		//cout << "stop_loc_overall = " << stop_loc_overall << endl;
		segment_num++;
		bool queryFound = true;

   	 	if((*read_local != 'A')&&(*read_local != 'C')&&(*read_local != 'G')&&(*read_local != 'T'))
   	 	{
   	 		queryFound = false;
   	 		//align_length[0] ++;
   	 		stop_loc = 0;
   	 		segment_align_SArange[0] = 1;
   	 		segment_align_SArange[1] = 0;
   	 		segment_align_rangeNum = 0;
   	 		queryFound = false;
   	 	}
   	 	unsigned int lcp_length = 0;
   	 	unsigned int start = 0, end = n-1;
   	 	unsigned int Min;
   	 	unsigned int c = 0;

   	 	//getFirstIntervalCompress(*read_local, &interval_begin, &interval_end, child, verifyChild);

   	 	//getFirstInterval(*read_local, &interval_begin, &interval_end, next);

   	 	getFirstInterval(*read_local, &interval_begin, &interval_end, child, verifyChild);

   	 	segment_align_SArange[0] = interval_begin;
   	 	segment_align_SArange[1] = interval_end;
   	 	segment_align_rangeNum = interval_end - interval_begin + 1;

   	 	//cout << endl << "Another Segment ... " << endl << "readheadChar = " << (*read_local) << " interval_begin = " << interval_begin <<
   	 	//" interval_end = " << interval_end << endl;

   	 	unsigned int iterateNum = 0;//debug;
   	 	while((c + stop_loc_overall< read_length) && (queryFound == true))
   	 	{
   	 		iterateNum++;
   	 		if(iterateNum>read_length)
   	 		{
   	 			return false;
   	 		}
   	 		unsigned int c_old = c;

			if(interval_begin != interval_end)
			{

 				//lcp_length = getlcp(interval_begin, interval_end, lcp, up, down);
	 			lcp_length = getlcp(interval_begin, interval_end, lcpCompress, child, verifyChild);

				Min = min(lcp_length, read_length - stop_loc_overall);

				unsigned int loc_pos = 0;
            	for(loc_pos = 0; loc_pos < Min - c_old; loc_pos++)
            	{
            		queryFound = (*(read_local+c_old+loc_pos) == *(chrom+sa[interval_begin]+c_old+loc_pos));
            		if (!queryFound)
            		{
            			break;
            		}
            	}
            	//cout << "queryFound = " << queryFound << endl;
            	if(!queryFound)
            	{
            		stop_loc = c_old + loc_pos;
            		break;
            	}

            	c = Min;
            	if(*(read_local+c) == 'N')
            	{
            		queryFound = false;
            		stop_loc = c;
            		break;
            	}
				start = interval_begin; end = interval_end;
				if (c + stop_loc_overall == read_length)
				{
					break;
				}
				//cout << "to get interval" << endl;
				unsigned int interval_begin_ori = interval_begin;
				unsigned int interval_end_ori = interval_end;
		    	//getIntervalCompress(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next,
		    	//	chrom, verifyChild);

		    	//getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, up, down, next,
		    	//	chrom);
			   	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next,
			   		chrom, verifyChild);
		    	//cout << "interval_begin != interval_end " << endl << "char = " << (*(read_local+c))
		    	//	<< " interval_begin = " << interval_begin << " interval_end = " << interval_end << endl;

		    	if(interval_begin > interval_end)
		    	{
		    		queryFound = false;
		    		stop_loc = c-1;
          			segment_align_SArange[0] = interval_begin_ori;
            		segment_align_SArange[1] = interval_end_ori;
            		segment_align_rangeNum = interval_end_ori - interval_begin_ori + 1;
		    		break;
		    	}
		    	else
		    	{
          			segment_align_SArange[0] = interval_begin;
            		segment_align_SArange[1] = interval_end;
            		segment_align_rangeNum = interval_end - interval_begin + 1;
		    	}
			}//end if
			else
			{
				//cout << "interval_begin == interval_end " << endl;
				unsigned int loc_pos = 0;
            	for(loc_pos = 0; loc_pos < read_length - c - stop_loc_overall; loc_pos++)
            	{
            		queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
            		if (!queryFound)
            			break;
            	}

	    		if(queryFound)
	    		{
	    			//align_length[101] ++;
	    		}
	    		else
	    		{
	    			stop_loc = c+loc_pos;
	    			//align_length[stop_loc] ++;
	    		}
          		segment_align_SArange[0] = interval_begin;
            	segment_align_SArange[1] = interval_end;
            	segment_align_rangeNum = interval_end - interval_begin + 1;

	    		break;
	    	}
		} //end while
		///////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////SEGMENT MAP RESULT////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////

    	if (queryFound && (interval_end >= interval_begin))
    	{
    		(norSegmentNum) ++;
    		//debug("segmentNum = "); debugln(*norSegmentNum);
    		unsigned int tmpSegLength = read_length - stop_loc_overall;

			if(tmpSegLength >= minValSegLength)
			{
				valLength = valLength + tmpSegLength;
			}

    		norSegmentLength[norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;

    		*(norSegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
    		*(norSegmentAlignNum + norSegmentNum - 1) = segment_align_rangeNum;


			for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,(unsigned int)100); alignment_num++)
			{
    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num)
    				= sa[segment_align_SArange[0] + alignment_num] + 1;
    				//- midPartMapPosForLongHeadInSecondLevelIndex + midPartMapPosForLongHead;
    		}

			segment_length = read_length-stop_loc_overall;

			break;
		}
		else
		{
			(norSegmentNum) ++;
			//debug("segmentNum = "); debugln(*norSegmentNum);
			if(norSegmentNum > (int)(read_length/5))
			{
				//debugln("map error, too many segments, there may exist too many Ns");
				return false;
			}
			norSegmentLength[norSegmentNum - 1] = stop_loc;

			if(stop_loc >= minValSegLength )
			{
				valLength = valLength + stop_loc;
			}

			norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
			norSegmentAlignNum[norSegmentNum - 1] = segment_align_rangeNum;

			for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,(unsigned int)100); alignment_num++)
		    {
    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num)
    				= sa[segment_align_SArange[0] + alignment_num] + 1;
    				//- midPartMapPosForLongHeadInSecondLevelIndex + midPartMapPosForLongHead;
    		}

			unsigned int stop_loc_overall_ori = stop_loc_overall;
			read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
			stop_loc_overall = stop_loc_overall + stop_loc + 1;

    		segment_length = stop_loc;

		}
   	}

	#ifdef DEBUG
    cout << "# of Segments = " << norSegmentNum << endl;//; debugln();

   	for(unsigned int k1 = 0; k1 < norSegmentNum; k1++)
   	{
  		cout << "segment "<< k1+1 << ": " << norSegmentLocInRead[k1] << "~"
  		<< (norSegmentLocInRead[k1] + norSegmentLength[k1] - 1)
  		<< "  Length: " << norSegmentLength[k1] << " Num: " << norSegmentAlignNum[k1] << endl;
      	if((norSegmentLength[k1]>=1)&&(norSegmentAlignNum[k1] < 40))
      	{
      		cout << "\tAlign Location: " << endl;
      		for(unsigned int k2 = 0; k2 < norSegmentAlignNum[k1]; k2++)
      		{
      			unsigned int locationInCurrentIndex = *(norSegmentAlignLoc + k1*CANDALILOC + k2);
      			getChrLocation(locationInCurrentIndex
      				- midPartMapPosForLongHeadInSecondLevelIndex + midPartMapPosForLongHead,
      				&align_chr, &align_chr_location);
      			cout << "\t" << "InSecondLevelIndex: " << locationInCurrentIndex
      			//<< *(norSegmentAlignLoc + k1*CANDALILOC + k2)
      			<<  " "
      			<< chr_name[align_chr] << " " << align_chr_location << endl;
      		}
      	}
   	}
	cout << ("segment_length_max = ") << segment_length_max << endl;
	#endif

	mapMain = (((norSegmentNum) == 1) && (norSegmentLength[0]==readLength) && ((norSegmentAlignNum[0])<40));
	//cout << "mapMain: " << mapMain << endl;
	if(mapMain)
	{
		*targetMappingAlignNum = (int)(norSegmentAlignNum[0]);
		for(int tmp = 0; tmp < (*targetMappingAlignNum); tmp++)
		{
			targetMappingAlignLoc[tmp] = norSegmentAlignLoc[tmp] + midPartMapPosForLongHead - midPartMapPosForLongHeadInSecondLevelIndex;
		}
	}
	//debugln("mapMain ended!!!");
   	free(norSegmentLength);
   	free(norSegmentLocInRead);
   	free(norSegmentAlignNum);
   	//Denote the start location of an alignment for a segment, start from 1;
   	free(norSegmentAlignLoc);

	return mapMain;
}*/
