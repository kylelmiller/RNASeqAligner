#include <stdlib.h>
#include <stdio.h>

using namespace std;

class Seg2ndOri_Info
{
public:
	unsigned int segmentNum;
	unsigned int norSegmentLength[SEGMENTNUM];
	unsigned int norSegmentLocInRead[SEGMENTNUM];
	unsigned int norSegmentAlignNum[SEGMENTNUM];
	unsigned int norSegmentAlignLoc[SEGMENTNUM * CANDALILOC];

	int longSegMinLength;

	Seg2ndOri_Info()
	{}

	string int_to_str(int numerical)
	{
			char c[100];
			sprintf(c,"%d",numerical);
			string str(c);
			return str;
	}

	//input : read, sa, up, down, next, chrom;
	//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc
	bool mapMainSecondLevel(char *read, unsigned int* sa, unsigned int * lcp, 
		unsigned int* up, unsigned int* down,
		unsigned int* next, char* chrom, 
		int readLength, unsigned int indexSize, 
		Index_Info* indexInfo)
	{
		bool mapMain = false;
		int norSegmentNum = 0;
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
		char* read_local = read;

		// start to map
		while (stop_loc_overall < read_length) // 15
		{
			segment_num++;
			bool queryFound = true;

	   	 	if((*read_local != 'A')&&(*read_local != 'C')&&(*read_local != 'G')&&(*read_local != 'T')) 
	   	 	{
	   	 		queryFound = false;
	   	 		stop_loc = 1;
	   	 		segment_align_SArange[0] = 1;
	   	 		segment_align_SArange[1] = 0;
	   	 		segment_align_rangeNum = 0;
	   	 		queryFound = false;   	

	   	 		if(norSegmentNum >= SEGMENTNUM)
	   	 		{
	   	 			segmentNum = SEGMENTNUM;
	   	 			return false;
	   	 		}
	   	 		
	   	 		(norSegmentNum) ++;

	   	 		norSegmentLength[norSegmentNum - 1] = 1;
				norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[norSegmentNum - 1] = 0;
				
				stop_loc = 1;	
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;   	 		
	   	 		continue;
	   	 	}

	   	 	//cout << "not N" << endl;
	   	 	unsigned int lcp_length = 0;
	   	 	unsigned int start = 0, end = n-1;
	   	 	unsigned int Min;
	   	 	unsigned int c = 0;
	 	 	
	   	 	getFirstInterval(*read_local, &interval_begin, &interval_end, next);
	   	 	segment_align_SArange[0] = interval_begin;
	   	 	segment_align_SArange[1] = interval_end;
	   	 	segment_align_rangeNum = interval_end - interval_begin + 1;

	   	 	unsigned int iterateNum = 0;
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
	           		//Xinan: COMPRESS INDEX

	 				lcp_length = getlcp(interval_begin, interval_end, lcp, up, down);

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

					unsigned int interval_begin_ori = interval_begin;
					unsigned int interval_end_ori = interval_end;

			    	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, up, down, next, 
			    		chrom);		    	
			    	
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

			/////////////////////////////////////////SEGMENT MAP RESULT////////////////////////////////////
			
	    	if (queryFound && (interval_end >= interval_begin)) 
	    	{
	   	 		if(norSegmentNum >= SEGMENTNUM)
	   	 		{
	   	 			segmentNum = SEGMENTNUM;
	   	 			return false;
	   	 		}
	    		(norSegmentNum) ++;

	    		unsigned int tmpSegLength = read_length - stop_loc_overall;

	    		norSegmentLength[norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;

	    		*(norSegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
	    		*(norSegmentAlignNum + norSegmentNum - 1) = segment_align_rangeNum;				
				
				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum, (unsigned int)CANDALILOC); alignment_num++)
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
				//cout << "!queryFound or end <= begin" << endl << "segNum: " << (norSegmentNum) << endl;
	   	 		if(norSegmentNum >= SEGMENTNUM)
	   	 		{
	   	 			segmentNum = SEGMENTNUM;
	   	 			return false;
	   	 		}

				(norSegmentNum) ++;

				norSegmentLength[norSegmentNum - 1] = stop_loc;

				norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[norSegmentNum - 1] = segment_align_rangeNum;

				for (unsigned int alignment_num = 0; alignment_num < std::min(segment_align_rangeNum, (unsigned int)CANDALILOC); alignment_num++)
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
		
		segmentNum = (norSegmentNum);

		return true;
	}

	//input : read, sa, up, down, next, chrom;
	//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc
	bool mapMainSecondLevel_compressedIndex(
		char *read, unsigned int* sa, 
		BYTE* lcpCompress, 
		unsigned int* child, char* chrom, 
		BYTE* verifyChild, int readLength, 
		Index_Info* indexInfo)
	{
		
		unsigned int norSegmentNum = 0;
		bool mapMain = false;	
		unsigned int stop_loc = 0; // location in one segment for iterations
		unsigned int stop_loc_overall = 0; //location in the whole read for iterations
		unsigned int segment_length = 0;
		unsigned int segment_num = 0;
		unsigned int segment_align_SArange[2] = {0,0};//legal align location in SA
		unsigned int segment_align_rangeNum = 0;
		unsigned int read_length = readLength; //READ_LENGTH;
		unsigned int interval_begin, interval_end;
		unsigned int n = (indexInfo->indexSize);//size of SA
		//*valLength = 0;
		char* read_local = read;
		while (stop_loc_overall < read_length) //- 15)
		{
			segment_num++;
			bool queryFound = true;

	   	 	if((*read_local != 'A')&&(*read_local != 'C')&&(*read_local != 'G')&&(*read_local != 'T')) 
	   	 	{
	   	 		queryFound = false;
	   	 		//align_length[0] ++;
	   	 		stop_loc = 1;
	   	 		segment_align_SArange[0] = 1;
	   	 		segment_align_SArange[1] = 0;
	   	 		segment_align_rangeNum = 0;
	   	 		queryFound = false;   	

	   	 		if(norSegmentNum >= SEGMENTNUM)
	   	 		{
	   	 			segmentNum = SEGMENTNUM;
	   	 			return false;
	   	 		}
	   	 		
	   	 		(norSegmentNum) ++;

	   	 		norSegmentLength[norSegmentNum - 1] = 1;
				norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[norSegmentNum - 1] = 0;
				
				//cout << "norSegmentNum: " << norSegmentNum << endl;
				//cout << "norSegmentLength: " << norSegmentLength[norSegmentNum - 1] << endl;
				//cout << "norSegmentLocInRead: "<< norSegmentLocInRead[norSegmentNum - 1] << endl;
				//cout << "norSegmentAlignNum: " << norSegmentAlignNum[norSegmentNum - 1] << endl << endl;
				stop_loc = 1;	
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;   	 		
	   	 		continue;
	   	 	}
	   	 	unsigned int lcp_length = 0;
	   	 	unsigned int start = 0, end = n-1;
	   	 	unsigned int Min;
	   	 	unsigned int c = 0;
	   	 	//cout << "char: " << (*read_local) << endl;
	   	 	getFirstInterval(*read_local, &interval_begin, &interval_end, child, verifyChild);
	   	 	//cout << "new Segment: " << segment_num << endl;
	   	 	//cout << "char: " << (*read_local) << " firstInterval: " << interval_begin << " ~ " << interval_end << endl;
	   	 	segment_align_SArange[0] = interval_begin;
	   	 	segment_align_SArange[1] = interval_end;
	   	 	segment_align_rangeNum = interval_end - interval_begin + 1;

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
	           		//Xinan: COMPRESS INDEX
	 				lcp_length = getlcp(interval_begin, interval_end, lcpCompress, child, verifyChild);

	 				//cout << "lcp_length: " << lcp_length << endl;
					//Min = min(lcp_length, read_length);
					Min = min(lcp_length, read_length - stop_loc_overall);
					unsigned int loc_pos = 0;
	            	for(loc_pos = 0; loc_pos < Min - c_old; loc_pos++)
	            	{
	            		queryFound = (*(read_local+c_old+loc_pos) == *(chrom+sa[interval_begin]+c_old+loc_pos));
	            		if (!queryFound)
	            			break;
	            	}

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
					if (c + stop_loc_overall== read_length)
						break;

					unsigned int interval_begin_ori = interval_begin;
					unsigned int interval_end_ori = interval_end;
			    	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child, chrom, verifyChild);

			    	//cout << "char: " << *(read_local+c) << " interval: " << interval_begin << " ~ " << interval_end << endl;
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
					unsigned int loc_pos = 0;
	            	for(loc_pos = 0; loc_pos < read_length - c- stop_loc_overall; loc_pos++)
	            	{
	            		queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
	            		if (!queryFound)
	            			break;
	            	}

		    		if(!queryFound)
		    			stop_loc = c+loc_pos;

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
	   	 		if(norSegmentNum > SEGMENTNUM)
	   	 		{
	   	 			segmentNum = SEGMENTNUM;
	   	 			return false;
	   	 		}
	   	 		if(norSegmentNum > (int)(read_length/5))
	   	 		{
	   	 			segmentNum = (int)(read_length/5);
	   	 			return false;
	   	 		}

	    		
	    		unsigned int tmpSegLength = read_length - stop_loc_overall;

	    		norSegmentLength[norSegmentNum - 1] = tmpSegLength;

	    		*(norSegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
	    		*(norSegmentAlignNum + norSegmentNum - 1) = segment_align_rangeNum;				
				
				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum, (unsigned int)CANDALILOC); alignment_num++)
				{    			
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				segment_length = read_length-stop_loc_overall;
				
				break;
			}
			else 
			{    
				(norSegmentNum) ++;
	   	 		if(norSegmentNum > SEGMENTNUM)
	   	 		{
	   	 			segmentNum = SEGMENTNUM;
	   	 			return false;
	   	 		}
				
				if(norSegmentNum > (int)(read_length/5))
				{
					segmentNum = (int)(read_length/5);
					//debugln("map error, too many segments, there may exist too many Ns");
					return false;
				}
				norSegmentLength[norSegmentNum - 1] = stop_loc;

				#ifdef SEGLENGTH
					segmentLength1[stop_loc]++;
					segmentLength2[stop_loc]++;
				#endif

				norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[norSegmentNum - 1] = segment_align_rangeNum;

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum, (unsigned int)CANDALILOC); alignment_num++)
			    {    			
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				unsigned int stop_loc_overall_ori = stop_loc_overall;
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;

	    		segment_length = stop_loc;
 
			}		
	   	}

		
		segmentNum = (norSegmentNum);
		mapMain = true;
		//debugln("mapMain ended!!!");
		return mapMain;
	}
};

class Seg_Info
{
public: 
	unsigned int segmentNum;
	unsigned int norSegmentLength[SEGMENTNUM];
	unsigned int norSegmentLocInRead[SEGMENTNUM];
	unsigned int norSegmentAlignNum[SEGMENTNUM];
	unsigned int norSegmentAlignLoc[SEGMENTNUM * CANDALILOC];

	int longSegMinLength;

	Seg_Info()
	{
		longSegMinLength = 20;
		segmentNum = 0;
	}

	Seg_Info(Seg2ndOri_Info* seg2ndOriInfo, int mapPosIntervalStart, int mapPosIntervalEnd,
		int chrPosStartIn2ndLevelIndex, Index_Info* indexInfo, const string& chromNameStr)
	{
		longSegMinLength = 18;
		segmentNum = seg2ndOriInfo->segmentNum;
		//cout << "segmentNum: " << segmentNum << endl;
		for(int tmpSeg = 0; tmpSeg < segmentNum; tmpSeg++)
		{
			//cout << "tmpSeg: " << tmpSeg << endl;

			norSegmentLength[tmpSeg] = seg2ndOriInfo->norSegmentLength[tmpSeg];
			norSegmentLocInRead[tmpSeg] = seg2ndOriInfo->norSegmentLocInRead[tmpSeg];
			int tmpSegCandiNum = 0;

			if(seg2ndOriInfo->norSegmentAlignNum[tmpSeg] <= CANDALILOC)
			{
				for(int tmpSegCandi = 0; tmpSegCandi < seg2ndOriInfo->norSegmentAlignNum[tmpSeg];
					tmpSegCandi++)
				{

					int tmpLoc = *(seg2ndOriInfo->norSegmentAlignLoc + tmpSeg*CANDALILOC + tmpSegCandi)
						+ chrPosStartIn2ndLevelIndex - 1;

					if((tmpLoc <= mapPosIntervalEnd)||(tmpLoc >= mapPosIntervalStart))
					{
						*(norSegmentAlignLoc + tmpSeg*CANDALILOC + tmpSegCandi)
							= indexInfo->getWholeGenomeLocation(chromNameStr, tmpLoc);
						tmpSegCandiNum ++;
					}
				}
			}
			else
			{
				tmpSegCandiNum = 0;
			}
			norSegmentAlignNum[tmpSeg] = tmpSegCandiNum;
		}
	}		

	bool checkSegLongOrNot(int segGroupNO)
	{
		return (norSegmentLength[segGroupNO] >= longSegMinLength);
	}
	
	int checkSegRelation(int segNO_1, int segNO_1_candiNO, 
		int segNO_2, int segNO_2_candiNO)
	{
		unsigned int segEndNum1 = segNO_1;
		unsigned int segStartNum2 = segNO_2;

		unsigned int alignLoc1 = *(norSegmentAlignLoc + segNO_1 * CANDALILOC + segNO_1_candiNO) 
			- norSegmentLocInRead[segNO_1] + 1;
		unsigned int alignLoc2 = *(norSegmentAlignLoc + segNO_2 * CANDALILOC + segNO_2_candiNO) 
			- norSegmentLocInRead[segNO_2] + 1;

		if (segEndNum1 < segStartNum2)
		{
			if ((segStartNum2 - segEndNum1) == 1)
			{
				if(alignLoc1 == alignLoc2)
				{
					return FIX_MATCH;
				}
				else if (alignLoc1 > alignLoc2)
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
				if(alignLoc1 == alignLoc2)
				{
					return FIX_MATCH;
				}
				else if (alignLoc1 > alignLoc2)
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

	int distanceBetweenSegment(int segNO_1, int segNO_1_candiNO, 
		int segNO_2, int segNO_2_candiNO)
	{
		int segDist = 1000000;
		unsigned int alignLoc1 = *(norSegmentAlignLoc + segNO_1 * CANDALILOC + segNO_1_candiNO) 
			- norSegmentLocInRead[segNO_1] + 1;
		unsigned int alignLoc2 = *(norSegmentAlignLoc + segNO_2 * CANDALILOC + segNO_2_candiNO) 
			- norSegmentLocInRead[segNO_2] + 1;

		unsigned int segmentDistance;
		if(alignLoc2 >= alignLoc1)
		{
			segmentDistance = alignLoc2 - alignLoc1;
		}
		else
		{
			segmentDistance = alignLoc1 - alignLoc2;
		}
		//unsigned int segmentDistance = alignLoc2 - alignLoc1;
		if((segmentDistance < 300000))
		{
			segDist = segmentDistance;
		} 		
		return segDist;
	}

	int getFirstLongSegNO()
	{
		for(int tmpSegNO = 0; tmpSegNO < segmentNum; tmpSegNO++)
			if((norSegmentLength[tmpSegNO] >= longSegMinLength)&&(norSegmentAlignNum[tmpSegNO] <= SEGMENTNUM))
				return tmpSegNO;

		return -1;
	}

	bool mapMain_SegInfo_preIndex(char *read, unsigned int* sa, BYTE* lcpCompress, 
		unsigned int* child, char* chrom, 
		unsigned int* valLength, BYTE* verifyChild, int readLength, Index_Info* indexInfo,
		int* PreIndexMappedLengthArray, unsigned int* PreIndexIntervalStartArray,
		unsigned int* PreIndexIntervalEndArray)
	{
		//cout << "mapMain_SegInfo_preIndex function starts ...... "<< endl; 
		//input : read, sa, up, down, next, chrom; 
		//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc 
		//cout << "readLength: " << readLength << endl;
		//*(read+readLength) = 'Y';
		unsigned int norSegmentNum;

		string readStringStr = read;

		(norSegmentNum) = 0;
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
		unsigned int n = (indexInfo->indexSize);//size of SA
		unsigned int norAlignLoc;
		unsigned int align_chr_location;
		unsigned int align_chr;
		*valLength = 0;
		char* read_local = read;

		string KmerToSearchInIndexStr;

		while (stop_loc_overall < read_length) //- 15)
		{
			segment_num++;
			bool queryFound = true;

			/////////////////////////////////////////////////////////////////////////
			/////////////////////////   K-mer search  ///////////////////////////////

			bool KmerSearchFound = false;
			int KmerMappedLength;
			unsigned int KmerSearchIndexIntervalStart = 0;
			unsigned int KmerSearchIndexIntervalEnd = 0;			
			
			//cout << "$$-1\nnew segment ! segment_num: " << segment_num << endl;
			//cout << "stop_loc_overall: " << stop_loc_overall << endl << endl;

			if(read_length - stop_loc_overall < INDEX_KMER_LENGTH)
			{
				//cout << "read_length - stop_loc_overall: " << read_length-stop_loc_overall << endl;
				//cout << "< INDEX_KMER_LENGTH" << endl;
				KmerSearchFound = false;
			}
			else
			{
				//cout << "read_length - stop_loc_overall: " << read_length-stop_loc_overall << endl;
				//cout << "> INDEX_KMER_LENGTH" << endl;
				KmerToSearchInIndexStr = readStringStr.substr(stop_loc_overall, INDEX_KMER_LENGTH);
				//cout << "KmerToSearchInIndexStr: " << KmerToSearchInIndexStr << endl;
				KmerSearchFound = this->getIndexInterval_PreIndex(KmerToSearchInIndexStr, &KmerMappedLength,
					&KmerSearchIndexIntervalStart, &KmerSearchIndexIntervalEnd, PreIndexMappedLengthArray, 
					PreIndexIntervalStartArray,	PreIndexIntervalEndArray);
				//cout << "KmerMappedLength: " << KmerMappedLength << endl;

			}

			//cout << "$$-2\nKmerSearchFound: " << KmerSearchFound << endl << endl;

			///////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////

	   	 	unsigned int lcp_length = 0;
	   	 	unsigned int start = 0, end = n-1;
	   	 	unsigned int Min;
	   	 	unsigned int c = 0;
	   	 	unsigned int iterateNum = 0;//debug;

	   	 	if(!KmerSearchFound)
	   	 	{	
		   	 	if((*read_local != 'A')&&(*read_local != 'C')&&(*read_local != 'G')&&(*read_local != 'T')) 
		   	 	{
		   	 		//cout << "read_local: " << (*read_local) << endl;
		   	 		queryFound = false;
		   	 		stop_loc = 0;
		   	 		segment_align_SArange[0] = 1;
		   	 		segment_align_SArange[1] = 0;
		   	 		segment_align_rangeNum = 0;
		   	 		queryFound = false;   	 			
		   	 		
		   	 		if(norSegmentNum >= SEGMENTNUM)
		   	 		{
		   	 			segmentNum = SEGMENTNUM;
		   	 			return false;
		   	 		}
		   	 		
		   	 		(norSegmentNum) ++;

		   	 		norSegmentLength[norSegmentNum - 1] = 1;
					norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
					norSegmentAlignNum[norSegmentNum - 1] = 0;
					
					//cout << "norSegmentNum: " << norSegmentNum << endl;
					//cout << "norSegmentLength: " << norSegmentLength[norSegmentNum - 1] << endl;
					//cout << "norSegmentLocInRead: "<< norSegmentLocInRead[norSegmentNum - 1] << endl;
					//cout << "norSegmentAlignNum: " << norSegmentAlignNum[norSegmentNum - 1] << endl << endl;
					stop_loc = 1;	
					read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
					stop_loc_overall = stop_loc_overall + stop_loc + 1;   	 		
		   	 		continue;		   	 		
		   	 	}
	   	 		lcp_length = 0;
	   	 		start = 0, end = n-1;
	   	 		c = 0;
		   	 	getFirstInterval(*read_local, &interval_begin, &interval_end, child, verifyChild);
		   	 	segment_align_SArange[0] = interval_begin;
		   	 	segment_align_SArange[1] = interval_end;
		   	 	segment_align_rangeNum = interval_end - interval_begin + 1;
		   	 	iterateNum = 0;
	   	 	}
	   	 	else // K-mer found in preIndex base
	   	 	{
	   	 		lcp_length = 0;
	   	 		start = 0, end = n-1;
	   	 		c = 0;
	   	 		interval_begin = KmerSearchIndexIntervalStart;
	   	 		interval_end = KmerSearchIndexIntervalEnd;
	    	 	segment_align_SArange[0] = interval_begin;
	   	 		segment_align_SArange[1] = interval_end;
	   	 		segment_align_rangeNum = interval_end - interval_begin + 1;
	   	 		iterateNum = 0;//debug;  	 		
	   	 	}

	   	 	//cout << "$$-3\ninterval_begin: " << interval_begin << endl << 
	   	 	//	"interval_end: " << interval_end << endl << endl;

	   	 	while((c + stop_loc_overall < read_length) && (queryFound == true))
	   	 	{
	   	 		//cout << "$$-4\n c: " << c << endl << endl;
	   	 		iterateNum++;
	   	 		if(iterateNum + stop_loc_overall >read_length)
	   	 		{
	   	 			return false;
	   	 		}
	   	 		unsigned int c_old = c;

	   	 		//cout << "c_old: " << c_old << endl;
			    //cout << "interval_begin: " << interval_begin << endl;
			   	//cout << "interval_end: " << interval_end << endl;
				
				if(interval_begin != interval_end)
				{ 
					//cout << "$$-5\ninterval_begin != interval_end" << endl;
	           		//Xinan: COMPRESS INDEX
	           		//lcp_length = getlcp(interval_begin, interval_end, lcp, child_up, child_down);
	 				lcp_length = getlcp(interval_begin, interval_end, lcpCompress, child, verifyChild//child_up, child_down
	 					);
					//cout << "lcp_length: " << lcp_length << endl;
					//cout << "read_length - stop_loc_overall: " << read_length - stop_loc_overall << endl;
					//Min = min(lcp_length, read_length);
					Min = min(lcp_length, read_length - stop_loc_overall);
					//cout << "Min: " << Min << endl;
					unsigned int loc_pos = 0;

					//cout << "c_old: " << c_old << endl;
					//////////////////////////// changed for preIndex  ///////////////////////////////
	    	   	 	//if(!KmerSearchFound)
	   	 			//{
		            	for(loc_pos = 0; loc_pos < Min - c_old; loc_pos++)
		            	{
		            		//cout << "loc_pos: " << loc_pos << endl;
		            		//cout << "...*(read_local+c_old+loc_pos): " << *(read_local+c_old+loc_pos) << endl;
		            		//cout << "...*(chrom+sa[interval_begin]+c_old+loc_pos): " << *(chrom+sa[interval_begin]+c_old+loc_pos) << endl;
		            		queryFound = (*(read_local+c_old+loc_pos) == *(chrom+sa[interval_begin]+c_old+loc_pos));
		            		//cout << "...queryFound: " << queryFound << endl;
		            		if (!queryFound)
		            		{	
		            			//cout << "......queryFound false! " << endl;
		            			//cout << "......stop at loc_pos: " << loc_pos << endl;
		            			break;
		            		}
		            	}
		            
		            //cout << "$$-A\nstop at loc_pos: " << loc_pos << endl;
		            //cout << "---- queryFound: " << queryFound << endl;
	            	
	            	if(!queryFound)
	            	{
	            		//cout << "$$-B\n! if (!queryFound)" << endl;
	            		//cout << "---- c_old: " << c_old << 
	            		//	 "---- loc_pos: " << loc_pos << endl;

	            		stop_loc = c_old + loc_pos;
	            		
	            		//cout << "---- stop_loc: " << stop_loc << endl;
	            		
	            		break;
	            	}
	            	
	            	c = Min;
	            	//cout << "** c = Min: " << c << endl << endl;

	            	if(*(read_local+c) == 'N')
	            	{
	            		queryFound = false; 
	            		stop_loc = c;
	            		break;
	            	}
					start = interval_begin; end = interval_end;
					//cout << "start: " << start << " end: " << end << endl << endl;

					if (c + stop_loc_overall == read_length)
					{	
						//cout << "** c==read_length  break !!!" << endl << endl; 			
						break;			
					}	
					//cout << " +++ c: " << c << endl << endl;
					unsigned int interval_begin_ori = interval_begin;
					unsigned int interval_end_ori = interval_end;
					//cout << "*(read_local+c): " << *(read_local+c) << endl;
			    	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next, 
			    		chrom, verifyChild);
			    	//cout << "+++ interval_begin: " << interval_begin << endl;
			    	//cout << "+++ interval_end: " << interval_end << endl << endl;
			    	if(interval_begin > interval_end)
			    	{
			    		//cout << "$$-8\ninterval_begin > interval_end" << endl;
			    		queryFound = false;
			    		stop_loc = c-1;
	          			segment_align_SArange[0] = interval_begin_ori;
	            		segment_align_SArange[1] = interval_end_ori;
	            		
	            		//cout << endl << "segment_align_SArange[0]: " << segment_align_SArange[0] 
	            		//	<< " segment_align_SArange[1]: " << segment_align_SArange[1] << endl;
	            		//cout << "........ break happens !!!!" << endl;
	            		segment_align_rangeNum = interval_end_ori - interval_begin_ori + 1;		    			
			    		break;
			    	}
			    	else
			    	{
			    		//cout << "$$-9\ninterval_begin <= interval_end" << endl;
	          			segment_align_SArange[0] = interval_begin;
	            		segment_align_SArange[1] = interval_end;
	            		
	            		//cout << endl << "segment_align_SArange[0]: " << segment_align_SArange[0] 
	            		//	<< " segment_align_SArange[1]: " << segment_align_SArange[1] << endl;

	            		segment_align_rangeNum = interval_end - interval_begin + 1;
			    	}
				}//end if
				else 
				{
					unsigned int loc_pos;

					// Compares the read to the reference, character by character
					for(loc_pos = 0; loc_pos < read_length - c - stop_loc_overall; loc_pos++)
					{
						queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
						if (!queryFound)
							break;
					}

					// If we didn't match the entire read, note where we stopped
		    		if(!queryFound)
		    			stop_loc = c+loc_pos;

		    		// Store the alignment information
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
			//cout << "queryFound: " << queryFound << endl;
			//cout << "interval_begin: " << interval_begin << endl;
			//cout << "interval_end: " << interval_end << endl;

	    	if (queryFound && (interval_end >= interval_begin)) 
	    	{
	    		//cout << "\n$$-01\nqueryFound && (interval_end >= interval_begin)" << endl;
	    		(norSegmentNum) ++;
	    		if(norSegmentNum > SEGMENTNUM)
	    		{
	    			segmentNum = (SEGMENTNUM);
	    			return false;
	    		}
	   	 		if(norSegmentNum > (int)(read_length/5))
	   	 		{
	   	 			segmentNum = (int)(read_length/5);
	   	 			return false;
	   	 		}
	    		//cout << "segmentNum: " << norSegmentNum << endl;
	    		unsigned int tmpSegLength = read_length - stop_loc_overall;

				if(tmpSegLength >= minValSegLength)
				{
					*valLength = *valLength + tmpSegLength;
				}

	    		norSegmentLength[norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;

	    		*(norSegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
	    		*(norSegmentAlignNum + norSegmentNum - 1) = segment_align_rangeNum;				
				
	    		//cout << "...SegLength: " << tmpSegLength << endl;
	    		//cout << "...SegLocInRead: " << stop_loc_overall + 1 << endl;
	    		//cout << "...SegAlignment_num: " << segment_align_rangeNum << endl << endl;

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum, (unsigned int)CANDALILOC); alignment_num++)
				{    			
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				//segment_length = read_length-stop_loc_overall;

				break;
			}
			else 
			{    
				(norSegmentNum) ++;
				//cout << "\n$$-02\nelse !" << endl;
				//cout << "segmentNum: " << norSegmentNum << endl;
				if((norSegmentNum > (int)(read_length/5)) || (norSegmentNum > SEGMENTNUM))
				{
					//debugln("map error, too many segments, there may exist too many Ns");
					if((int)(read_length/5) > SEGMENTNUM)
						segmentNum = (SEGMENTNUM);
					else
						segmentNum = (int)(read_length/5);

					return false;
				}
				//cout << "stop_loc: " << stop_loc << endl;
				norSegmentLength[norSegmentNum - 1] = stop_loc;

				if(stop_loc >= minValSegLength )
				{
					*valLength = *valLength + stop_loc;
				}

				norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[norSegmentNum - 1] = segment_align_rangeNum;

	    		//cout << "...SegLength: " << stop_loc << endl;
	    		//cout << "...SegLocInRead: " << stop_loc_overall + 1 << endl;
	    		//cout << "...SegAlignment_num: " << segment_align_rangeNum << endl;

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum, (unsigned int)CANDALILOC); alignment_num++)
			    {    			
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				unsigned int stop_loc_overall_ori = stop_loc_overall;
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				//cout << "stop_loc + 1: " << stop_loc +1 << endl;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;
				//cout << "stop_loc_overall: " << stop_loc_overall << endl << endl;
	    		//segment_length = stop_loc;
	    		//cout << "segment_length: " << segment_length << endl;
			}		
	   	}

		segmentNum = (norSegmentNum);
		mapMain = true;
		//cout << "mapMain_SegInfo_preIndex function ends ...... "<< endl; 
		//debugln("mapMain ended!!!");
		return mapMain;
	}

	bool getIndexInterval_PreIndex(const string& readPreStr, int* mappedLength, 
		unsigned int* indexIntervalStart, unsigned int* indexIntervalEnd,
		int* PreIndexMappedLengthArray, unsigned int* PreIndexIntervalStartArray,
		unsigned int* PreIndexIntervalEndArray)
	{
		if (readPreStr.find("N") != readPreStr.npos)
			return false;

		int preIndexStrSize = readPreStr.length();

		unsigned int preIndexNO = 0;

		int baseForCount = 1;
		for(int tmp = preIndexStrSize - 1; tmp >= 0; tmp--)
		{
			preIndexNO = preIndexNO + baseChar2intArray[(readPreStr.at(tmp) - 'A')] * baseForCount;
			baseForCount = baseForCount * 4;
		}

		(*mappedLength) = PreIndexMappedLengthArray[preIndexNO];
		(*indexIntervalStart) = PreIndexIntervalStartArray[preIndexNO];
		(*indexIntervalEnd) = PreIndexIntervalEndArray[preIndexNO];

		return true;
	}

};

class Path_Info
{
public:

	vector< vector< pair<int, int> > > PathVec_seg; // vector <vector <seg_group, seg_candi> >
	vector< bool > PathValidBoolVec;

	vector < pair<int, pair<int, int> > > validPathVec_toPair; // vector < chromNameInt, <alignPosStart, alignPosEnd> >
	vector <int> validPathVec; // vector <validPathNO in PathVec_seg>

	vector< pair<int, Splice_Info*> > fixedPathVec;
	vector< int > fixedPathMismatchVec;

	vector< bool > PathFixedBoolVec;

	vector < pair< pair<int, int>, Splice_Info*> > finalPathVec;

	Path_Info()
	{}

	int pathValidNumInt()
	{
		int pathValidNum = 0;
		for(int tmp = 0; tmp < PathValidBoolVec.size(); tmp++)
		{
			if(PathValidBoolVec[tmp])
				pathValidNum++;
		}
		return pathValidNum;
	}

	void memoryFree()
	{
		for(int tmp = 0; tmp < fixedPathVec.size(); tmp++)
		{
			delete(fixedPathVec[tmp].second);
		}
		for(int tmp = 0; tmp < finalPathVec.size(); tmp++)
		{
			delete(finalPathVec[tmp].second);
		}
	}

	void getFinalPath_extend2HeadTail(Index_Info* indexInfo, Seg_Info* segInfo, int readLength, const string& readSeq_inProcess)
	{
		//cout << "start to get Final path" << endl;

		
		if(segInfo->segmentNum < 1)
		{
			return;
		}
		//int readLength = segInfo->norSegmentLocInRead[segInfo->segmentNum - 1] + segInfo->norSegmentLength[segInfo->segmentNum - 1] - 1;
		for(int tmpPath = 0; tmpPath < fixedPathVec.size(); tmpPath++)
		{
			int mismatchNumToAdd = 0;
			int fixedPathNO = fixedPathVec[tmpPath].first;

			int tmpPath1stSegGroupNO = (PathVec_seg[fixedPathNO])[0].first;
			int tmpPath1stSegCandiNO = (PathVec_seg[fixedPathNO])[0].second;

			int tmpPathElementSize = (PathVec_seg[fixedPathNO]).size();
			int tmpPathLastSegGroupNO = (PathVec_seg[fixedPathNO])[tmpPathElementSize - 1].first;

			unsigned int PathMapPos = *(segInfo->norSegmentAlignLoc + tmpPath1stSegGroupNO * CANDALILOC + tmpPath1stSegCandiNO);
			unsigned int tmpChrNameInt, tmpChrPosInt;
			indexInfo->getChrLocation(PathMapPos, &tmpChrNameInt, &tmpChrPosInt);

			
			int tmpPathFinalMapPos = tmpChrPosInt;
			int tmpPathFinalMapChr = tmpChrNameInt;
			//string tmpChrNameStr = indexInfo->chrNameStr[tmpChrNameInt];
			//int tmpSegmentMapPos = tmpChrPosInt;
			//cout << "tmpPathFinalMapChr: " << tmpPathFinalMapChr << endl;
			//cout << "tmpPathFinalMapPos: " << tmpPathFinalMapPos << endl;

			Splice_Info* tmpSpliceInfo = new Splice_Info();
			tmpSpliceInfo->jump_code.clear(); 

			//cout << "here" << endl;

			int tmpUnfixedHeadLength = segInfo->norSegmentLocInRead[tmpPath1stSegGroupNO] - 1;
			//cout << "tmpUnfixedHeadLength " << tmpUnfixedHeadLength << endl;
			string readSubSeqInProcess_head = readSeq_inProcess.substr(0, tmpUnfixedHeadLength);
			
			bool scoreStringBool_head; 
			string chromSubSeqInProcess_head; 

			size_t max_mismatch_head = (tmpUnfixedHeadLength)/8 + 1;
			size_t mismatch_bits_head = 0;
			size_t comb_bits_head = 0;

			//cout << "tmpUnfixedHeadLength " << tmpUnfixedHeadLength << endl;

			if(tmpPathFinalMapPos - tmpUnfixedHeadLength - 1 < 0)
			{
				scoreStringBool_head = false;
			}
			else
			{
				chromSubSeqInProcess_head = (indexInfo->chromStr[tmpChrNameInt]).substr(tmpPathFinalMapPos - tmpUnfixedHeadLength - 1, tmpUnfixedHeadLength);

				scoreStringBool_head = score_string(readSubSeqInProcess_head, chromSubSeqInProcess_head, max_mismatch_head, mismatch_bits_head, comb_bits_head);
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
			else
			{

			}

			tmpSpliceInfo->appendJumpCode(fixedPathVec[tmpPath].second);


			//////////////////////////// add last jump code /////////////////////////////////////////////////////
			int endMappedBaseMapPos = (fixedPathVec[tmpPath].second)->getEndBaseMapPos_jump_code(tmpPathFinalMapPos);
			//int readLength = segInfo->norSegmentLocInRead[segInfo->segmentNum - 1] + segInfo->norSegmentLength[segInfo->segmentNum - 1] - 1;
			int tmpUnfixedTailLength = readLength - (segInfo->norSegmentLocInRead[tmpPathLastSegGroupNO] 
										+ segInfo->norSegmentLength[tmpPathLastSegGroupNO] - 1);
		
			//cout << "tmpUnfixedTailLength: " << tmpUnfixedTailLength << endl;
			//cout << "readLength - tmpUnfixedTailLength: " << readLength - tmpUnfixedTailLength << endl;
			string readSubSeqInProcess_tail 
				= readSeq_inProcess.substr(readLength - tmpUnfixedTailLength, tmpUnfixedTailLength);

			bool scoreStringBool_tail;
			string chromSubSeqInProcess_tail;

			size_t max_mismatch_tail = (tmpUnfixedTailLength)/8 + 1;
			size_t mismatch_bits_tail = 0;
			size_t comb_bits_tail = 0;


			if(endMappedBaseMapPos + tmpUnfixedTailLength >= indexInfo->chromLength[tmpChrNameInt])
			{
				//cout << "here 1" << endl;
				scoreStringBool_tail = false;
			}
			else
			{
				//cout << "here 2" << endl;
				//cout << "endMappedBaseMapPos: " << endMappedBaseMapPos << endl;
				//cout << "tmpUnf"
				chromSubSeqInProcess_tail = (indexInfo->chromStr[tmpChrNameInt]).substr(endMappedBaseMapPos, tmpUnfixedTailLength);

				scoreStringBool_tail 
					= score_string(readSubSeqInProcess_tail, chromSubSeqInProcess_tail, 
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
			else
			{

			}
			//if(scoreStringBool_head || scoreStringBool_tail)
			//{
			//	cout << "extend to head/tail successfully" << endl;
			//}

			tmpSpliceInfo->getFinalJumpCode();			
			
			if((tmpUnfixedHeadLength > 0)&&(scoreStringBool_head))
			{
				tmpPathFinalMapPos = tmpPathFinalMapPos - tmpUnfixedHeadLength;
			}

			//cout << "tmpPath: " << tmpPath + 1 << endl;

			int oldMismatchNum = fixedPathMismatchVec[tmpPath];
			int newMismatchNum = oldMismatchNum + mismatchNumToAdd;

			//cout << "oldMismatchNum: " << oldMismatchNum << endl;
			//cout << "newMismatchNum: " << newMismatchNum << endl;


			fixedPathMismatchVec[tmpPath] = newMismatchNum;

			finalPathVec.push_back(pair< pair<int, int>, Splice_Info*> (pair<int,int> (tmpPathFinalMapChr, tmpPathFinalMapPos), tmpSpliceInfo) );
			//delete(tmpSpliceInfo);
		}		
	}

	void addNewSegGroupToCurrentPathInfoVec(Seg_Info* segInfo, int segGroupNO)
	{
		bool longSegBool = segInfo->checkSegLongOrNot(segGroupNO);

		int currentPathNum = PathVec_seg.size();
		for(int tmpSegCandiLoc = 0; tmpSegCandiLoc < segInfo->norSegmentAlignNum[segGroupNO]; tmpSegCandiLoc++)
		{
			this->addNewSegCandiToCurrentPathInfoVec_matchIndelUniqueSpliceMultiPath(segInfo, segGroupNO, tmpSegCandiLoc, longSegBool, currentPathNum);
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

	int minDistanceWithCurrentPath(Seg_Info* segInfo, int segGroupNO, int segCandiNO, int currentPathNum, int* minSegNumGap)
	{
		int minDistance_value = 900000;
		int minDistance_path = 1000;
		int tmpMinSegNumGap = 100;
		for(int tmpPathElementNO = 0; tmpPathElementNO < currentPathNum; tmpPathElementNO++)
		{
			int tmpPathElementVecSize = (PathVec_seg[tmpPathElementNO]).size();
			int tmpPathLastSegGroupNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).first;
			int tmpPathLastSegCandiNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).second;
			int tmpSegDistance = segInfo->distanceBetweenSegment(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, 
				segGroupNO, segCandiNO);
			int tmpSegNumGap = segGroupNO - tmpPathLastSegGroupNO;
			
			if(tmpSegNumGap <= 0)
			{
				continue;
			}

			if(tmpSegDistance < SPLICEDISTANCE)
			{
				if(tmpSegNumGap < tmpMinSegNumGap)
				{
					tmpMinSegNumGap = tmpSegNumGap;
					minDistance_value = tmpSegDistance;
					minDistance_path = tmpPathElementNO;
				}
				else if(tmpSegNumGap == tmpMinSegNumGap)
				{
					if(tmpSegDistance < minDistance_value)
					{
						minDistance_value = tmpSegDistance;
						minDistance_path = tmpPathElementNO;
					}
					else
					{

					}
				}
				else
				{

				}
			}

		}
		(*minSegNumGap) = tmpMinSegNumGap;
		return minDistance_value;		
	}

	void addNewSegCandiToCurrentPathInfoVec_matchIndelUniqueSpliceMultiPath(Seg_Info* segInfo, 
		int segGroupNO, int segCandiNO, bool longSegBool, int currentPathNum)
	{
		bool relatedToSomePath = false;
		int minSegNumGap = 0;
		int minDistance = this->minDistanceWithCurrentPath(segInfo, segGroupNO, segCandiNO, currentPathNum, &minSegNumGap);
		//cout << "minSegNumGap: " << minSegNumGap << " " << segGroupNO+1 << "," << segCandiNO+1 << endl;
		//cout << "min Distance = " << minDistance << endl;
		
		//int currentPathNum = PathVec_seg.size();
		if(minDistance <= MAX_DELETION_LENGTH)
		{
			//cout << "PathVec_seg.size(): " << PathVec_seg.size() << endl;

			for(int tmpPathElementNO = 0; tmpPathElementNO < currentPathNum; tmpPathElementNO++)
			{
				//cout << tmpPathElementNO << endl;
				//cout << " ... PathVec_seg.size(): " << PathVec_seg.size() << endl;
				int tmpPathElementVecSize = (PathVec_seg[tmpPathElementNO]).size();
				int tmpPathLastSegGroupNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).first;
				int tmpPathLastSegCandiNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).second;
				int tmpSegDistance = segInfo->distanceBetweenSegment(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, 
					segGroupNO, segCandiNO);
				int tmpSegNumGap = segGroupNO - tmpPathLastSegGroupNO;
				if((tmpSegDistance <= MAX_DELETION_LENGTH)&&(tmpSegNumGap <= minSegNumGap))
				{
					//cout << "tmpPathElementNO: " << tmpPathElementNO << endl;
					PathValidBoolVec[tmpPathElementNO] = false;
					this->copyOldPath_AddSegCandi_AddNewPath(tmpPathElementNO, segGroupNO, segCandiNO);
					relatedToSomePath = true;
				}
			}
		}
		else
		{
			for(int tmpPathElementNO = 0; tmpPathElementNO < currentPathNum; tmpPathElementNO ++)
			{

				int tmpPathElementVecSize = (PathVec_seg[tmpPathElementNO]).size();
				int tmpPathLastSegGroupNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).first;
				int tmpPathLastSegCandiNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).second;
				int tmpSegDistance = segInfo->distanceBetweenSegment(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, 
					segGroupNO, segCandiNO);
				int tmpSegNumGap = segGroupNO - tmpPathLastSegGroupNO;
				if((tmpSegDistance < SPLICEDISTANCE)&&(tmpSegNumGap <= minSegNumGap))
				{
					relatedToSomePath = true;
					
					PathValidBoolVec[tmpPathElementNO] = false;
					this->copyOldPath_AddSegCandi_AddNewPath(tmpPathElementNO, segGroupNO, segCandiNO);					
				}
			}
		}

		if((!relatedToSomePath)&&
			(longSegBool))
		{
			vector< pair<int,int> > tmpPathElementVec;
			tmpPathElementVec.push_back( pair<int,int> (segGroupNO, segCandiNO) );	
			PathVec_seg.push_back(tmpPathElementVec);
			PathValidBoolVec.push_back(true);		
			vector< pair<int,int> >().swap(tmpPathElementVec);
		}
	}

	bool getPossiPathFromSeg(Seg_Info* segInfo)
	{
		bool possiPathExists = false;
		int firstLongSegNO = segInfo->getFirstLongSegNO();
		
		if(firstLongSegNO < 0)
		{
			return false;
		}
		
			for(int tmpSegCandiLoc = 0; tmpSegCandiLoc < segInfo->norSegmentAlignNum[firstLongSegNO]; tmpSegCandiLoc++)
			{
				vector< pair<int,int> > tmpPathElementVec;
				tmpPathElementVec.push_back( pair<int,int> (firstLongSegNO, tmpSegCandiLoc) );
				PathVec_seg.push_back(tmpPathElementVec);
				PathValidBoolVec.push_back(true);
				vector< pair<int,int> > ().swap(tmpPathElementVec);
			}

		for(int tmpSegNO = firstLongSegNO + 1; tmpSegNO < segInfo->segmentNum; tmpSegNO ++)
		{
			if(//(!segInfo->checkSegLongOrNot(tmpSegNO)) && 
				(segInfo->norSegmentAlignNum[tmpSegNO]) > CANDALILOC)
				continue;
			this->addNewSegGroupToCurrentPathInfoVec(segInfo, tmpSegNO);
		}

		possiPathExists = true;
		return possiPathExists;
	}
};
