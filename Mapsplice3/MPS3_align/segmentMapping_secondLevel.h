
void getFirstInterval(char ch, unsigned int *interval_begin, unsigned int *interval_end, unsigned int* child_next)
{
	if (ch == 'C')
	{
		*interval_begin = child_next[0];
		*interval_end = child_next[*interval_begin] - 1; 
	}
	else if(ch == 'A') //ch == 'A'
	{
		*interval_begin = 0;
		*interval_end = child_next[0] - 1;
	}
	else //ch == 'G' or ch == 'T'
	{
		if (ch == 'G')
		{
			*interval_begin = child_next[child_next[0]];
			*interval_end = child_next[*interval_begin] - 1;
		}
		else
		{
			*interval_begin = child_next[child_next[child_next[0]]];
			*interval_end = child_next[*interval_begin] - 1;
		}
	}
}

void getInterval_old(unsigned int start, unsigned int end, unsigned int position, char ch, unsigned int *interval_begin, unsigned int *interval_end, 
	unsigned int* sa, unsigned int* child_up, unsigned int* child_down, unsigned int* child_next, char* chrom)
{   
	unsigned int index_begin;
	unsigned int pos;
	*interval_end = 0;
	*interval_begin = 1;
	if(ch == 'C')
	{
		if ((start < child_up[end + 1]) && (end >= child_up[end +1]))
       		index_begin = child_up[end+1];
		else index_begin = child_down[start];

		if(chrom[sa[index_begin]+position] == 'C')
		{
			*interval_begin = index_begin;
			*interval_end = child_next[index_begin]-1;
			if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
		}	
		else if(chrom[sa[start]+position] == 'C')
		{
			*interval_begin = start;
			*interval_end = index_begin - 1;
			if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
		}		
	}
	else if(ch == 'A') // ch == 'A'
	{
		if ((start < child_up[end + 1]) && (end >= child_up[end +1]))
       		index_begin = child_up[end+1];
		else index_begin = child_down[start];

		if(chrom[sa[start]+position] == 'A')
		{
			*interval_begin = start;
			*interval_end = index_begin - 1;
			if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
		}	
	}
	else // ch == 'G' or ch == 'T'
	{	
		if(ch == 'G')
		{
			//cout << "ch = " << "G" << endl;
			if ((start < child_up[end + 1]) && (end >= child_up[end +1]))
       			index_begin = child_up[end+1];
			else index_begin = child_down[start];

			;
			pos = child_next[index_begin];
			/*cout << "index_begin " << index_begin << endl;
			cout << "pos " << pos << endl;
			cout << "chrom[sa[pos]+position] " << chrom[sa[pos]+position] << endl;
			cout << "chrom[sa[start]+position] " << chrom[sa[start]+position] << endl;
			cout << "chrom[sa[index_begin]+position] " << chrom[sa[index_begin]+position] << endl;

			cout << "start" << start << endl;*/


			if(chrom[sa[pos]+position] == 'G')
			{
				//cout << "here 1" << endl;
				*interval_begin = pos;
				//cout << "interval_begin " << pos << endl; 
				*interval_end = child_next[pos] - 1;
				//cout << "interval_end " << child_next[pos] - 1 << endl;;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}
			else if(chrom[sa[start]+position] == 'G')
			{
				//cout << "here 2" << endl;
				*interval_begin = start;
				*interval_end = index_begin - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}	
			/*if(chrom[sa[start]+position] == 'G')
			{
				//cout << "here 2" << endl;
				*interval_begin = start;
				//cout << "interval_begin" << start << endl;
				*interval_end = index_begin - 1;
				//cout << 
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}	
			else if(chrom[sa[pos]+position] == 'G')
			{
				cout << "here 1" << endl;
				*interval_begin = pos;
				cout << "interval_begin " << pos << endl; 
				*interval_end = child_next[pos] - 1;
				cout << "interval_end " << child_next[pos] - 1 << endl;;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}*/
			else if(chrom[sa[index_begin]+position] == 'G')
			{
				//cout << "here 3" << endl;
				*interval_begin = index_begin;
				*interval_end = child_next[index_begin]-1;

				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}	
			else
			{
				cout << "here 4" << endl;
			}			
		}
		else //(ch == 'T')
		{
			if ((start < child_up[end + 1]) && (end >= child_up[end +1]))
       			index_begin = child_up[end+1];
			else index_begin = child_down[start];

			if(chrom[sa[start]+position] == 'T')
			{
				*interval_begin = start;
				*interval_end = index_begin - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}	
			else if(chrom[sa[index_begin]+position] == 'T')
			{
				*interval_begin = index_begin;
				*interval_end = child_next[index_begin]-1;

				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}		
			else 
			{
				pos = child_next[index_begin];

				if(chrom[sa[pos]+position] == 'T')
				{
					*interval_begin = pos;
					*interval_end = child_next[pos] - 1;
					if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
				}
				else
				{
					pos = child_next[pos];
					if(chrom[sa[pos]+position] == 'T')
					{
						*interval_begin = pos;
						*interval_end = child_next[pos] - 1;
						if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
					}
				}						
			}
		}
	}			
}

void getInterval(unsigned int start, unsigned int end, unsigned int position, char ch, unsigned int *interval_begin, unsigned int *interval_end, 
	unsigned int* sa, unsigned int* child_up, unsigned int* child_down, unsigned int* child_next, char* chrom)
{   
	unsigned int index_begin;
	unsigned int pos;
	*interval_end = 0;
	*interval_begin = 1;
	if(ch == 'C')
	{
		if ((start < child_up[end + 1]) && (end >= child_up[end +1]))
       		index_begin = child_up[end+1];
		else index_begin = child_down[start];

		if(chrom[sa[start]+position] == 'C')
		{
			*interval_begin = start;
			*interval_end = index_begin - 1;
			if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
		}	
		else if(chrom[sa[index_begin]+position] == 'C')
		{
			*interval_begin = index_begin;
			*interval_end = child_next[index_begin]-1;
			if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
		}
		else
		{
			//cout << "error in getInterval C" << endl;
		}		
	}
	else if(ch == 'A') // ch == 'A'
	{
		if ((start < child_up[end + 1]) && (end >= child_up[end +1]))
       		index_begin = child_up[end+1];
		else index_begin = child_down[start];

		if(chrom[sa[start]+position] == 'A')
		{
			*interval_begin = start;
			*interval_end = index_begin - 1;
			if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
		}
		else
		{
			//cout << "error in getInterval A" << endl;
		}	
	}
	else // ch == 'G' or ch == 'T'
	{	
		if(ch == 'G')
		{
			if ((start < child_up[end + 1]) && (end >= child_up[end +1]))
       			index_begin = child_up[end+1];
			else index_begin = child_down[start];

			;
			pos = child_next[index_begin];

			if(chrom[sa[start]+position] == 'G')
			{
				*interval_begin = start;
				*interval_end = index_begin - 1;

				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}	
			else if(chrom[sa[index_begin]+position] == 'G')
			{
				*interval_begin = index_begin;
				*interval_end = child_next[index_begin]-1;

				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}
			else if(chrom[sa[pos]+position] == 'G')
			{
				*interval_begin = pos; 
				*interval_end = child_next[pos] - 1;

				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}	
			else
			{
				//cout << "error in getInterval G" << endl;
			}			
		}
		else //(ch == 'T')
		{
			if ((start < child_up[end + 1]) && (end >= child_up[end +1]))
       			index_begin = child_up[end+1];
			else index_begin = child_down[start];

			if(chrom[sa[start]+position] == 'T')
			{
				*interval_begin = start;
				*interval_end = index_begin - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}	
			else if(chrom[sa[index_begin]+position] == 'T')
			{
				*interval_begin = index_begin;
				*interval_end = child_next[index_begin]-1;

				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}		
			else 
			{
				pos = child_next[index_begin];

				if(chrom[sa[pos]+position] == 'T')
				{
					*interval_begin = pos;
					*interval_end = child_next[pos] - 1;
					if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
				}
				else
				{
					pos = child_next[pos];
					if(chrom[sa[pos]+position] == 'T')
					{
						*interval_begin = pos;
						*interval_end = child_next[pos] - 1;
						if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
					}
				}						
			}
		}
	}			
}

unsigned int getlcp(unsigned int start, unsigned int end, unsigned int* lcp, unsigned int* child_up, unsigned int* child_down)
{
	//getLcpBegin = clock();
	unsigned int lcp_length;
	
	if ((start < child_up[end + 1]) && (end >= child_up[end +1]))
        lcp_length = lcp[child_up[end+1]];
	else 
		lcp_length = lcp[child_down[start]];
		
	//getLcpEnd = clock();
	//getLcpCost = getLcpCost + getLcpEnd - getLcpBegin;
	return lcp_length;
}

bool mapMainSecondLevel(char *read, unsigned int* sa, unsigned int * lcp, unsigned int* up, unsigned int* down,
	unsigned int* next, char* chrom, unsigned int* norSegmentNum, 
	unsigned int* norSegmentLength, unsigned int* norSegmentLocInRead, unsigned int* norSegmentAlignNum, 
	unsigned int* norSegmentAlignLoc, unsigned int* valLength, int readLength, unsigned int indexSize, 
	unsigned int midPartMapPosForLongHead, unsigned int midPartMapPosForLongHeadInSecondLevelIndex,
	Index_Info* indexInfo)
{
	//input : read, sa, up, down, next, chrom; 
	//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc 
	//cout << "sa[0] = " << sa[0] << endl;

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
	//unsigned int align_length[102] = {0}; 
	unsigned int n = indexInfo->indexSize;//size of SA
	unsigned int norAlignLoc;
	unsigned int align_chr_location;
	unsigned int align_chr;
	*valLength = 0;
	//debugln("start mapMain Function!!!");
	char* read_local = read;

	//cout << "start to map " << endl;
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
   	 	getFirstInterval(*read_local, &interval_begin, &interval_end, next);
   	 	segment_align_SArange[0] = interval_begin;
   	 	segment_align_SArange[1] = interval_end;
   	 	segment_align_rangeNum = interval_end - interval_begin + 1;

   	 	//cout << endl << "Another Segment ... " << endl << "readheadChar = " << (*read_local) << " interval_begin = " << interval_begin << 
   	 	//" interval_end = " << interval_end << endl;

   	 	unsigned int iterateNum = 0;//debug;
   	 	while((c < read_length) && (queryFound == true))
   	 	{
   	 		iterateNum++;
   	 		if(iterateNum>read_length)
   	 		{
   	 			return false;
   	 		}
   	 		unsigned int c_old = c;
			
			if(interval_begin != interval_end)
			{ 
           		//Xinan: COMPRESS INDEX
           		//lcp_length = getlcp(interval_begin, interval_end, lcp, child_up, child_down);

 				//lcp_length = getlcpCompress(interval_begin, interval_end, lcpCompress, child, verifyChild//child_up, child_down
 				//	);

 				lcp_length = getlcp(interval_begin, interval_end, lcp, up, down);
				//lcp_length ++;
				//cout << "lcp_length = " << lcp_length << endl;
				//if(lcp_length < 1)
				//{
				//	lcp_length ++;
				//}
				//cout << "lcp_length = " << lcp_length << endl;
				Min = min(lcp_length, read_length);
				
				unsigned int loc_pos = 0;
            	for(loc_pos = 0; loc_pos < Min - c_old; loc_pos++)
            	{
            		queryFound = (*(read_local+c_old+loc_pos) == *(chrom+sa[interval_begin]+c_old+loc_pos));
            		if (!queryFound)
            		{	
            			//cout << "char in read = " << (*(read_local+c_old+loc_pos)) << endl
            			//	<< "char in chrom = " << (*(chrom+sa[interval_begin]+c_old+loc_pos)) << endl;
            			//cout << "interval_begin = " << interval_begin << endl;
            			//cout << "sa[interval_begin] = " << sa[interval_begin] << endl;
            			//cout << "queryFound failed at " << loc_pos << endl;
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
				if (c == read_length)
				{				
					break;			
				}	
				//cout << "to get interval" << endl;
				unsigned int interval_begin_ori = interval_begin;
				unsigned int interval_end_ori = interval_end;
		    	//getIntervalCompress(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next, 
		    	//	chrom, verifyChild);

		    	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, up, down, next, 
		    		chrom);		    	
		    	
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
            	for(loc_pos = 0; loc_pos < read_length - c; loc_pos++)
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
    		(*norSegmentNum) ++;
    		//debug("segmentNum = "); debugln(*norSegmentNum);
    		unsigned int tmpSegLength = read_length - stop_loc_overall;

			if(tmpSegLength >= minValSegLength)
			{
				*valLength = *valLength + tmpSegLength;
			}

    		norSegmentLength[*norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;

    		*(norSegmentLocInRead + *norSegmentNum - 1) = stop_loc_overall + 1;
    		*(norSegmentAlignNum + *norSegmentNum - 1) = segment_align_rangeNum;				
			

			for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,100); alignment_num++)
			{    			
    			*(norSegmentAlignLoc + (*norSegmentNum - 1) * CANDALILOC + alignment_num) 
    				= sa[segment_align_SArange[0] + alignment_num] + 1;
    				//- midPartMapPosForLongHeadInSecondLevelIndex + midPartMapPosForLongHead;
    		}

			segment_length = read_length-stop_loc_overall;

			break;
		}
		else 
		{    
			(*norSegmentNum) ++;
			//debug("segmentNum = "); debugln(*norSegmentNum);
			if(*norSegmentNum > (int)(read_length/5))
			{
				//debugln("map error, too many segments, there may exist too many Ns");
				return false;
			}
			norSegmentLength[*norSegmentNum - 1] = stop_loc;

			if(stop_loc >= minValSegLength )
			{
				*valLength = *valLength + stop_loc;
			}

			norSegmentLocInRead[*norSegmentNum - 1] = stop_loc_overall + 1;
			norSegmentAlignNum[*norSegmentNum - 1] = segment_align_rangeNum;

			for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,100); alignment_num++)
		    {    			
    			*(norSegmentAlignLoc + (*norSegmentNum - 1) * CANDALILOC + alignment_num) 
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
   	//cout << " " << endl;
    cout << "# of Segments = " << *norSegmentNum << endl;//; debugln();

   	for(unsigned int k1 = 0; k1 < *norSegmentNum; k1++)
   	{
  		cout << "segment "<< k1+1 << ": " << norSegmentLocInRead[k1] << "~"   
  		<< (norSegmentLocInRead[k1] + norSegmentLength[k1] - 1) 
  		<< "  Length: " << norSegmentLength[k1] << " Num: " << norSegmentAlignNum[k1] << endl;
      	if((norSegmentLength[k1]>=1)&&(norSegmentAlignNum[k1] < 40))
      	{
      		cout << "\tAlign Location: " << endl;
      		for(unsigned int k2 = 0; k2 < norSegmentAlignNum[k1]; k2++)
      		{
      			/*//getChrLocation(*(norSegmentAlignLoc + k1*CANDALILOC + k2), &align_chr, &align_chr_location);
      			unsigned int locationInCurrentIndex = *(norSegmentAlignLoc + k1*CANDALILOC + k2);
      			cout << "\t" 
      			<< locationInCurrentIndex << endl;*/
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
	
	mapMain = true;
	//debugln("mapMain ended!!!");
	return mapMain;
}

int LimitRange = 20;
 
bool mapMainSecondLevelForHead(char *read, unsigned int* sa, unsigned int * lcp, unsigned int* up, unsigned int* down,
	unsigned int* next, char* chrom, unsigned int* norSegmentNum, 
	unsigned int* norSegmentLength, unsigned int* norSegmentLocInRead, unsigned int* norSegmentAlignNum, 
	unsigned int* norSegmentAlignLoc, unsigned int* valLength, int readLength, unsigned int indexSize, 
	unsigned int midPartMapPosForLongHead, unsigned int midPartMapPosForLongHeadInSecondLevelIndex,
	int* mappedLengthWithLimitRange, unsigned int* rangeStart, unsigned int* rangeEnd, Index_Info* indexInfo
	)
{
	//input : read, sa, up, down, next, chrom; 
	//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc 
	//cout << "sa[0] = " << sa[0] << endl;

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
	//unsigned int align_length[102] = {0}; 
	unsigned int n = indexInfo->indexSize;//size of SA
	unsigned int norAlignLoc;
	unsigned int align_chr_location;
	unsigned int align_chr;
	*valLength = 0;
	bool getMappedLengthWithLimitRange = false;
	//debugln("start mapMain Function!!!");
	char* read_local = read;

	//cout << "start to map " << endl;
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
   	 	getFirstInterval(*read_local, &interval_begin, &interval_end, next);
   	 	segment_align_SArange[0] = interval_begin;
   	 	segment_align_SArange[1] = interval_end;
   	 	segment_align_rangeNum = interval_end - interval_begin + 1;

   	 	//cout << endl << "Another Segment ... " << endl << "readheadChar = " << (*read_local) << " interval_begin = " << interval_begin << 
   	 	//" interval_end = " << interval_end << endl;

   	 	unsigned int iterateNum = 0;//debug;
   	 	while((c < read_length) && (queryFound == true))
   	 	{
   	 		iterateNum++;
   	 		if(iterateNum>read_length)
   	 		{
   	 			return false;
   	 		}
   	 		unsigned int c_old = c;
			
			if(interval_begin != interval_end)
			{ 
           		//Xinan: COMPRESS INDEX
           		//lcp_length = getlcp(interval_begin, interval_end, lcp, child_up, child_down);

 				//lcp_length = getlcpCompress(interval_begin, interval_end, lcpCompress, child, verifyChild//child_up, child_down
 				//	);

 				lcp_length = getlcp(interval_begin, interval_end, lcp, up, down);
 				if((interval_end - interval_begin + 1 < LimitRange) && (!getMappedLengthWithLimitRange))
 				{
 					getMappedLengthWithLimitRange = true;
 					(*rangeStart) = interval_begin;
 					(*rangeEnd) = interval_end;
 					(*mappedLengthWithLimitRange) = lcp_length;
 					#ifdef DEBUG
 					cout << "..................   mappedLength = " << lcp_length << "  ................." << endl;
 					cout << "..................   range = " << interval_end - interval_begin + 1 << "  ................." << endl;
 					#endif
 				}
				//lcp_length ++;
				//cout << "lcp_length = " << lcp_length << endl;
				//if(lcp_length < 1)
				//{
				//	lcp_length ++;
				//}
				//cout << "lcp_length = " << lcp_length << endl;
				Min = min(lcp_length, read_length);
				
				unsigned int loc_pos = 0;
            	for(loc_pos = 0; loc_pos < Min - c_old; loc_pos++)
            	{
            		queryFound = (*(read_local+c_old+loc_pos) == *(chrom+sa[interval_begin]+c_old+loc_pos));
            		if (!queryFound)
            		{	
            			//cout << "char in read = " << (*(read_local+c_old+loc_pos)) << endl
            			//	<< "char in chrom = " << (*(chrom+sa[interval_begin]+c_old+loc_pos)) << endl;
            			//cout << "interval_begin = " << interval_begin << endl;
            			//cout << "sa[interval_begin] = " << sa[interval_begin] << endl;
            			//cout << "queryFound failed at " << loc_pos << endl;
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
				if (c == read_length)
				{				
					break;			
				}	
				//cout << "to get interval" << endl;
				unsigned int interval_begin_ori = interval_begin;
				unsigned int interval_end_ori = interval_end;
		    	//getIntervalCompress(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next, 
		    	//	chrom, verifyChild);

		    	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, up, down, next, 
		    		chrom);		    	
		    	
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
            	for(loc_pos = 0; loc_pos < read_length - c; loc_pos++)
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
    		(*norSegmentNum) ++;
    		//debug("segmentNum = "); debugln(*norSegmentNum);
    		unsigned int tmpSegLength = read_length - stop_loc_overall;

			if(tmpSegLength >= minValSegLength)
			{
				*valLength = *valLength + tmpSegLength;
			}

    		norSegmentLength[*norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;

    		*(norSegmentLocInRead + *norSegmentNum - 1) = stop_loc_overall + 1;
    		*(norSegmentAlignNum + *norSegmentNum - 1) = segment_align_rangeNum;				
			

			for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,100); alignment_num++)
			{    			
    			*(norSegmentAlignLoc + (*norSegmentNum - 1) * CANDALILOC + alignment_num) 
    				= sa[segment_align_SArange[0] + alignment_num] + 1;
    				//- midPartMapPosForLongHeadInSecondLevelIndex + midPartMapPosForLongHead;
    		}

			segment_length = read_length-stop_loc_overall;

			break;
		}
		else 
		{    
			(*norSegmentNum) ++;
			//debug("segmentNum = "); debugln(*norSegmentNum);
			if(*norSegmentNum > (int)(read_length/5))
			{
				//debugln("map error, too many segments, there may exist too many Ns");
				return false;
			}
			norSegmentLength[*norSegmentNum - 1] = stop_loc;

			if(stop_loc >= minValSegLength )
			{
				*valLength = *valLength + stop_loc;
			}

			norSegmentLocInRead[*norSegmentNum - 1] = stop_loc_overall + 1;
			norSegmentAlignNum[*norSegmentNum - 1] = segment_align_rangeNum;

			for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,100); alignment_num++)
		    {    			
    			*(norSegmentAlignLoc + (*norSegmentNum - 1) * CANDALILOC + alignment_num) 
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
   	//cout << " " << endl;
    cout << "# of Segments = " << *norSegmentNum << endl;//; debugln();

   	for(unsigned int k1 = 0; k1 < *norSegmentNum; k1++)
   	{
  		cout << "segment "<< k1+1 << ": " << norSegmentLocInRead[k1] << "~"   
  		<< (norSegmentLocInRead[k1] + norSegmentLength[k1] - 1) 
  		<< "  Length: " << norSegmentLength[k1] << " Num: " << norSegmentAlignNum[k1] << endl;
      	if((norSegmentLength[k1]>=1)&&(norSegmentAlignNum[k1] < 40))
      	{
      		cout << "\tAlign Location: " << endl;
      		for(unsigned int k2 = 0; k2 < norSegmentAlignNum[k1]; k2++)
      		{
      			/*//getChrLocation(*(norSegmentAlignLoc + k1*CANDALILOC + k2), &align_chr, &align_chr_location);
      			unsigned int locationInCurrentIndex = *(norSegmentAlignLoc + k1*CANDALILOC + k2);
      			cout << "\t" 
      			<< locationInCurrentIndex << endl;*/
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
	
	mapMain = true;
	//debugln("mapMain ended!!!");
	return mapMain;
}


bool mapMainSecondLevelForTargetMapping(char *read, unsigned int* sa, unsigned int * lcp, unsigned int* up, unsigned int* down,
	unsigned int* next, char* chrom, int readLength, unsigned int indexSize, unsigned int midPartMapPosForLongHead, unsigned int midPartMapPosForLongHeadInSecondLevelIndex,
	int* targetMappingAlignNum, unsigned int* targetMappingAlignLoc, Index_Info* indexInfo)
{
	//cout << "start target mapping " << endl;
	//cout << "read: " << endl;
	//cout << read << endl;
    unsigned int norSegmentNum = 0;  	
    //cout << "*norSegmentNum = " << endl; cout << (*norSegmentNum) << endl;
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
   	 	getFirstInterval(*read_local, &interval_begin, &interval_end, next);
   	 	segment_align_SArange[0] = interval_begin;
   	 	segment_align_SArange[1] = interval_end;
   	 	segment_align_rangeNum = interval_end - interval_begin + 1;

   	 	//cout << endl << "Another Segment ... " << endl << "readheadChar = " << (*read_local) << " interval_begin = " << interval_begin << 
   	 	//" interval_end = " << interval_end << endl;

   	 	unsigned int iterateNum = 0;//debug;
   	 	while((c < read_length) && (queryFound == true))
   	 	{
   	 		iterateNum++;
   	 		if(iterateNum>read_length)
   	 		{
   	 			return false;
   	 		}
   	 		unsigned int c_old = c;
			
			if(interval_begin != interval_end)
			{ 

 				lcp_length = getlcp(interval_begin, interval_end, lcp, up, down);

				Min = min(lcp_length, read_length);
				
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
				if (c == read_length)
				{				
					break;			
				}	
				//cout << "to get interval" << endl;
				unsigned int interval_begin_ori = interval_begin;
				unsigned int interval_end_ori = interval_end;
		    	//getIntervalCompress(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next, 
		    	//	chrom, verifyChild);

		    	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, up, down, next, 
		    		chrom);		    	
		    	
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
            	for(loc_pos = 0; loc_pos < read_length - c; loc_pos++)
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
			

			for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,100); alignment_num++)
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

			for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,100); alignment_num++)
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
      			/*//getChrLocation(*(norSegmentAlignLoc + k1*CANDALILOC + k2), &align_chr, &align_chr_location);
      			unsigned int locationInCurrentIndex = *(norSegmentAlignLoc + k1*CANDALILOC + k2);
      			cout << "\t" 
      			<< locationInCurrentIndex << endl;*/
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
}
