#include <stdlib.h>
#include <stdio.h>

using namespace std;

//unsigned int baseForCountArray[14] = {1, 4, 16, 64, 256, 1024, 4096, 16384, 0};

/*inline int baseChar2int(char baseChar)
{
	return baseChar2intArray[(baseChar - 'A')];
}*/

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

	bool mapMainSecondLevel(char *read, unsigned int* sa, unsigned int * lcp, 
		unsigned int* up, unsigned int* down,
		unsigned int* next, char* chrom, 
		int readLength, unsigned int indexSize, 
		//unsigned int midPartMapPosForLongHead, unsigned int midPartMapPosForLongHeadInSecondLevelIndex,
		Index_Info* indexInfo)
	{
		//input : read, sa, up, down, next, chrom; 
		//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc 
		//cout << "sa[0] = " << sa[0] << endl;
		//*(read + readLength) = 'Y';
		bool mapMain = false;
		int norSegmentNum = 0;	
		//(norSegmentNum) = 0;
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
		//*valLength = 0;
		//debugln("start mapMain Function!!!");
		char* read_local = read;

		//cout << "start to map " << endl;
		while (stop_loc_overall < read_length) //- 15)
		{
			//cout << "stop_loc_overall = " << stop_loc_overall << endl;
			segment_num++;
			bool queryFound = true;

   	 		//cout << "*read_local: " << (*read_local) << endl; 
   	 		//cout << "stop_loc_overall: " << stop_loc_overall << endl;

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

	   	 	//cout << "not N" << endl;
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
				
				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
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

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
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

		mapMain = true;
		//debugln("mapMain ended!!!");
		return mapMain;
	}

	bool mapMainSecondLevel_compressedIndex(
		char *read, unsigned int* sa, 
		BYTE* lcpCompress, 
		unsigned int* child, char* chrom, 
		BYTE* verifyChild, int readLength, 
		Index_Info* indexInfo)
	{
		//input : read, sa, up, down, next, chrom; 
		//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc 
		
		//cout << "start to mapMainSecondLevel_compressedIndex " << endl;
		//*(read + readLength) = 'Y';
		unsigned int norSegmentNum;

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
	           		//lcp_length = getlcp(interval_begin, interval_end, lcp, child_up, child_down);
	 				lcp_length = getlcp(interval_begin, interval_end, lcpCompress, child, verifyChild//child_up, child_down
	 					);
	 				//cout << "lcp_length: " << lcp_length << endl;
					//Min = min(lcp_length, read_length);
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
					if (c + stop_loc_overall== read_length)
					{				
						break;			
					}	
					unsigned int interval_begin_ori = interval_begin;
					unsigned int interval_end_ori = interval_end;
			    	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next, 
			    		chrom, verifyChild);
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

		    		if(queryFound) 
		    		{
		    		}
		    		else 
		    		{ 
		    			stop_loc = c+loc_pos;
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

				//if(tmpSegLength >= minValSegLength)
				//{
				//	*valLength = *valLength + tmpSegLength;
				//}

	    		norSegmentLength[norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;

	    		*(norSegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
	    		*(norSegmentAlignNum + norSegmentNum - 1) = segment_align_rangeNum;				
				

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
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

				//if(stop_loc >= minValSegLength )
				//{
				//	*valLength = *valLength + stop_loc;
				//}

				norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[norSegmentNum - 1] = segment_align_rangeNum;

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
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

	string segInfoStr(Index_Info* indexInfo, int chrPosStartIn2ndLevelIndex, string chrNameStr)
	{
		string segInfoStr;
		segInfoStr += "\noriginal 2nd level segment Info: \n";
		segInfoStr = segInfoStr + "SegmentNum: " + int_to_str(segmentNum) + "\n";
		unsigned int align_chr, align_chr_location;
	   	for(unsigned int k1 = 0; k1 < segmentNum; k1++)
	   	{
	  		segInfoStr = segInfoStr + "... segment " + int_to_str(k1+1) + ": " + int_to_str(norSegmentLocInRead[k1]) + "~"   
	  		+ int_to_str(norSegmentLocInRead[k1] + norSegmentLength[k1] - 1) + 
	  		"  Length: " + int_to_str(norSegmentLength[k1]) + " Num: " + int_to_str(norSegmentAlignNum[k1]) + "\n";

	      	if(//(norSegmentLength[k1]>=10)&&
	      		(norSegmentAlignNum[k1] <= CANDALILOC))
	      	{
	      		//segInfoStr += "...... Align Location: \n";
	      		for(unsigned int k2 = 0; k2 < norSegmentAlignNum[k1]; k2++)
	      		{
	      			segInfoStr = segInfoStr + "\t" + int_to_str(k2+1) +
	      			//<< *(norSegmentAlignLoc + k1*CANDALILOC + k2) 
	      			+ ". " 
	      			+ chrNameStr//(indexInfo->chrNameStr)[align_chr] 
	      			+ " " + int_to_str((int)((*(norSegmentAlignLoc + k1*CANDALILOC + k2)) + chrPosStartIn2ndLevelIndex - 1)) 
	      			+ "\n";
	      		}
	      	}
	   	}
		return segInfoStr;//+"\n";
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

		//cout << "segNO_1: " << segNO_1 << "segNO_1_candiNO: " << segNO_1_candiNO << " alignLoc1: " << alignLoc1 << endl;
		//cout << "segNO_2: " << segNO_2 << "segNO_2_candiNO: " << segNO_2_candiNO << " alignLoc2: " << alignLoc2 << endl;

		//cout << "alignLoc1: " << alignLoc1 << endl;
		//cout << "alignLoc2: " << alignLoc2 << endl;

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

	bool checkSegRelatedOrNot(int segNO_1, int segNO_1_candiNO, 
		int segNO_2, int segNO_2_candiNO)
	{
		int segRelation = checkSegRelation(segNO_1, segNO_1_candiNO, segNO_2, segNO_2_candiNO);
		//cout << endl << "SegRelation: " << segRelation << endl;
		if((segRelation == FIX_TOO_FAR)||(segRelation == FIX_TOO_CLOSE)||(segRelation == FIX_NO_RELATIONSHIP))
		{
			return false;
		} 
		else
		{
			return true;
		}
	}

	int getFirstLongSegNO()
	{
		bool longSegExists = false;
		for(int tmpSegNO = 0; tmpSegNO < segmentNum; tmpSegNO++)
		{
			if((norSegmentLength[tmpSegNO] >= longSegMinLength)&&(norSegmentAlignNum[tmpSegNO] <= SEGMENTNUM))
			{
				longSegExists = true;
				return tmpSegNO;
			}
		}
		if(!longSegExists)
			return -1;
	}

	bool mapMain_SegInfo_originalSize(char *read, unsigned int* sa, BYTE* lcpCompress, //unsigned int* child, 
		char* chrom, 
		unsigned int* valLength,  int readLength, Index_Info* indexInfo, 
		unsigned int* child_up, 
		unsigned int* child_down, 
		unsigned int* child_next)
	{
		//input : read, sa, up, down, next, chrom; 
		//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc 
		//*(read + readLength) = 'Y';
		unsigned int norSegmentNum;

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
	   	 	getFirstInterval_originalSize(*read_local, &interval_begin, &interval_end, //child, verifyChild
	   	 		//child_up, child_down, 
	   	 		child_next);
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
	           		//lcp_length = getlcp(interval_begin, interval_end, lcp, child_up, child_down);
	 				lcp_length = getlcp_originalSize(interval_begin, interval_end, lcpCompress, //child, verifyChild//child_up, child_down
	 					child_up, child_down, child_next);
	 				//cout << "lcp_length: " << lcp_length << endl;
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
			    	getInterval_originalSize(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, //child,//child_up, child_down, child_next, 
			    		chrom, //verifyChild
			    		child_up, child_down, child_next);
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
	            	for(loc_pos = 0; loc_pos < read_length - c - stop_loc_overall; loc_pos++)
	            	{
	            		queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
	            		if (!queryFound)
	            			break;
	            	}

		    		if(queryFound) 
		    		{
		    		}
		    		else 
		    		{ 
		    			stop_loc = c+loc_pos;
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

				if(tmpSegLength >= minValSegLength)
				{
					*valLength = *valLength + tmpSegLength;
				}

	    		norSegmentLength[norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;
	    		#ifdef SEGLENGTH
	    			segmentLength1[stop_loc]++;
	    		#endif
	    		*(norSegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
	    		*(norSegmentAlignNum + norSegmentNum - 1) = segment_align_rangeNum;				
				

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
				{    			
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				segment_length = read_length-stop_loc_overall;
				
				#ifdef DEBUG
				if(segment_length_max < segment_length)
				{
					norAlignLoc = align_chr_location - stop_loc_overall;					
					segment_length_max = segment_length;
				}
				#endif

				break;
			}
			else 
			{  
	   	 		if(norSegmentNum >= SEGMENTNUM)
	   	 		{
	   	 			segmentNum = SEGMENTNUM;
	   	 			return false;
	   	 		}  
				(norSegmentNum) ++;
				if(norSegmentNum > (int)(read_length/5))
				{
					//debugln("map error, too many segments, there may exist too many Ns");
					return false;
				}
				norSegmentLength[norSegmentNum - 1] = stop_loc;

				#ifdef SEGLENGTH
					segmentLength1[stop_loc]++;
					segmentLength2[stop_loc]++;
				#endif

				if(stop_loc >= minValSegLength )
				{
					*valLength = *valLength + stop_loc;
				}

				norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[norSegmentNum - 1] = segment_align_rangeNum;

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum, CANDALILOC); alignment_num++)
			    {    			
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				unsigned int stop_loc_overall_ori = stop_loc_overall;
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;

	    		segment_length = stop_loc;
				
				#ifdef DEBUG 		
				if (segment_length_max < segment_length)
				{
					norAlignLoc = align_chr_location - stop_loc_overall_ori;
					segment_length_max = segment_length;
				}
				#endif
			}		
	   	}
		
		segmentNum = (norSegmentNum);
		mapMain = true;
		return mapMain;
	}

	bool mapMain_SegInfo(char *read, unsigned int* sa, 
		BYTE* lcpCompress, unsigned int* child, 
		char* chrom, 
		unsigned int* valLength, BYTE* verifyChild, 
		int readLength, Index_Info* indexInfo)
	{
		//input : read, sa, up, down, next, chrom; 
		//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc 
		//*(read + readLength) = 'Y';
		unsigned int norSegmentNum;

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
	           		//lcp_length = getlcp(interval_begin, interval_end, lcp, child_up, child_down);
	 				lcp_length = getlcp(interval_begin, interval_end, lcpCompress, child, verifyChild//child_up, child_down
	 					);
	 				//cout << "lcp_length: " << lcp_length << endl;
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
					if (c + stop_loc_overall== read_length)
					{				
						break;			
					}	
					unsigned int interval_begin_ori = interval_begin;
					unsigned int interval_end_ori = interval_end;
			    	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next, 
			    		chrom, verifyChild);
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
	            	for(loc_pos = 0; loc_pos < read_length - c - stop_loc_overall; loc_pos++)
	            	{
	            		queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
	            		if (!queryFound)
	            			break;
	            	}

		    		if(queryFound) 
		    		{
		    		}
		    		else 
		    		{ 
		    			stop_loc = c+loc_pos;
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
	   	 		if(norSegmentNum >= SEGMENTNUM)
	   	 		{
	   	 			segmentNum = SEGMENTNUM;
	   	 			return false;
	   	 		}
	    		(norSegmentNum) ++;
	    		unsigned int tmpSegLength = read_length - stop_loc_overall;

				if(tmpSegLength >= minValSegLength)
				{
					*valLength = *valLength + tmpSegLength;
				}

	    		norSegmentLength[norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;
	    		#ifdef SEGLENGTH
	    			segmentLength1[stop_loc]++;
	    		#endif
	    		*(norSegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
	    		*(norSegmentAlignNum + norSegmentNum - 1) = segment_align_rangeNum;				
				

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
				{    			
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				segment_length = read_length-stop_loc_overall;
				
				#ifdef DEBUG
				if(segment_length_max < segment_length)
				{
					norAlignLoc = align_chr_location - stop_loc_overall;					
					segment_length_max = segment_length;
				}
				#endif

				break;
			}
			else 
			{    
	   	 		if(norSegmentNum >= SEGMENTNUM)
	   	 		{
	   	 			segmentNum = SEGMENTNUM;
	   	 			return false;
	   	 		}
				(norSegmentNum) ++;
				if(norSegmentNum > (int)(read_length/5))
				{
					//debugln("map error, too many segments, there may exist too many Ns");
					return false;
				}
				norSegmentLength[norSegmentNum - 1] = stop_loc;

				#ifdef SEGLENGTH
					segmentLength1[stop_loc]++;
					segmentLength2[stop_loc]++;
				#endif

				if(stop_loc >= minValSegLength )
				{
					*valLength = *valLength + stop_loc;
				}

				norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[norSegmentNum - 1] = segment_align_rangeNum;

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
			    {    			
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				unsigned int stop_loc_overall_ori = stop_loc_overall;
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;

	    		segment_length = stop_loc;
				
				#ifdef DEBUG 		
				if (segment_length_max < segment_length)
				{
					norAlignLoc = align_chr_location - stop_loc_overall_ori;
					segment_length_max = segment_length;
				}
				#endif
			}		
	   	}
		//cout << "mapMain_SegInfo function starts 3 ..." << endl;
		/*#ifdef DEBUG
	   	//cout << " " << endl;
	    cout << "# of Segments = " << norSegmentNum << endl;//; debugln();

	   	for(unsigned int k1 = 0; k1 < norSegmentNum; k1++)
	   	{
	  		cout << "segment "<< k1+1 << ": " << norSegmentLocInRead[k1] << "~"   
	  		<< (norSegmentLocInRead[k1] + norSegmentLength[k1] - 1) 
	  		<< "  Length: " << norSegmentLength[k1] << " Num: " << norSegmentAlignNum[k1] << endl;
	      	if((norSegmentLength[k1]>=10)&&(norSegmentAlignNum[k1] < CANDALILOC))
	      	{
	      		cout << "\tAlign Location: " << endl;
	      		for(unsigned int k2 = 0; k2 < norSegmentAlignNum[k1]; k2++)
	      		{
	      			getChrLocation(*(norSegmentAlignLoc + k1*CANDALILOC + k2), &align_chr, &align_chr_location);
	      			cout << "\t" 
	      			//<< *(norSegmentAlignLoc + k1*CANDALILOC + k2) 
	      			<<  " " 
	      			<< chr_name[align_chr] << " " << align_chr_location << endl;
	      		}
	      	}
	   	}
	   	cout << "mapMain_SegInfo function starts 4 ..." << endl;
		cout << ("segment_length_max = ") << segment_length_max << endl;
		#endif */
		
		segmentNum = (norSegmentNum);
		mapMain = true;
		//debugln("mapMain ended!!!");
		return mapMain;
	}
	/*
	bool mapMain_SegInfo_toCheck(char *read, unsigned int* sa, 
		BYTE* lcpCompress, unsigned int* child, char* chrom, 
		unsigned int* valLength, BYTE* verifyChild, 
		int readLength, Index_Info* indexInfo)
	{
		//input : read, sa, up, down, next, chrom; 
		//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc 
		unsigned int norSegmentNum;

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

		while (stop_loc_overall < read_length) //- 15)
		{
			segment_num++;
			bool queryFound = true;

	   	 	if((*read_local != 'A')&&(*read_local != 'C')&&(*read_local != 'G')&&(*read_local != 'T')) 
	   	 	{
	   	 		queryFound = false;
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
	   	 	getFirstInterval(*read_local, &interval_begin, &interval_end, child, verifyChild);
	   	 	segment_align_SArange[0] = interval_begin;
	   	 	segment_align_SArange[1] = interval_end;
	   	 	segment_align_rangeNum = interval_end - interval_begin + 1;

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
	 				lcp_length = getlcp(interval_begin, interval_end, lcpCompress, child, verifyChild//child_up, child_down
	 					);
					Min = min(lcp_length, read_length);
					
					unsigned int loc_pos = 0;
	            	
	            	for(loc_pos = 0; loc_pos < Min - c_old; loc_pos++)//******* Xinan changed 
	            	{
	            		queryFound = (*(read_local+c_old+loc_pos) == *(chrom+sa[interval_begin]+c_old+loc_pos));
	            		if (!queryFound)
	            		{	
	            			//loc_pos = loc_pos - 1;
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
					if (c == read_length)
					{				
						break;			
					}	
					unsigned int interval_begin_ori = interval_begin;
					unsigned int interval_end_ori = interval_end;
			    	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next, 
			    		chrom, verifyChild);

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
	            	for(loc_pos = 0; loc_pos < read_length - c; loc_pos++)
	            	{
	            		queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
	            		if (!queryFound)
	            			break;
	            	}

		    		if(queryFound) 
		    		{
		    		}
		    		else 
		    		{ 
		    			stop_loc = c+loc_pos;
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
	    		unsigned int tmpSegLength = read_length - stop_loc_overall;

				if(tmpSegLength >= minValSegLength)
				{
					*valLength = *valLength + tmpSegLength;
				}

	    		norSegmentLength[norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;
	    		#ifdef SEGLENGTH
	    			segmentLength1[stop_loc]++;
	    		#endif
	    		*(norSegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
	    		*(norSegmentAlignNum + norSegmentNum - 1) = segment_align_rangeNum;				
				

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
				{    			
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				segment_length = read_length-stop_loc_overall;
				
				#ifdef DEBUG
				if(segment_length_max < segment_length)
				{
					norAlignLoc = align_chr_location - stop_loc_overall;					
					segment_length_max = segment_length;
				}
				#endif

				break;
			}
			else 
			{    
				(norSegmentNum) ++;

				norSegmentLength[norSegmentNum - 1] = stop_loc;

				norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[norSegmentNum - 1] = segment_align_rangeNum;

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
			    {    			
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				//unsigned int stop_loc_overall_ori = stop_loc_overall;
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;

	    		segment_length = stop_loc;
				
			}		
	   	}
		//cout << "mapMain_SegInfo function starts 3 ..." << endl;
		#ifdef DEBUG
	   	//cout << " " << endl;
	    cout << "# of Segments = " << norSegmentNum << endl;//; debugln();

	   	for(unsigned int k1 = 0; k1 < norSegmentNum; k1++)
	   	{
	  		cout << "segment "<< k1+1 << ": " << norSegmentLocInRead[k1] << "~"   
	  		<< (norSegmentLocInRead[k1] + norSegmentLength[k1] - 1) 
	  		<< "  Length: " << norSegmentLength[k1] << " Num: " << norSegmentAlignNum[k1] << endl;
	      	if((norSegmentLength[k1]>=10)&&(norSegmentAlignNum[k1] < CANDALILOC))
	      	{
	      		cout << "\tAlign Location: " << endl;
	      		for(unsigned int k2 = 0; k2 < norSegmentAlignNum[k1]; k2++)
	      		{
	      			getChrLocation(*(norSegmentAlignLoc + k1*CANDALILOC + k2), &align_chr, &align_chr_location);
	      			cout << "\t" 
	      			//<< *(norSegmentAlignLoc + k1*CANDALILOC + k2) 
	      			<<  " " 
	      			<< chr_name[align_chr] << " " << align_chr_location << endl;
	      		}
	      	}
	   	}
	   	cout << "mapMain_SegInfo function starts 4 ..." << endl;
		cout << ("segment_length_max = ") << segment_length_max << endl;
		#endif 
		
		segmentNum = (norSegmentNum);
		mapMain = true;
		//debugln("mapMain ended!!!");
		return mapMain;
	}

	bool mapMain_SegInfo_preIndex_NotReady(char *read, unsigned int* sa, BYTE* lcpCompress, 
		unsigned int* child, char* chrom, 
		unsigned int* valLength, BYTE* verifyChild, int readLength, Index_Info* indexInfo,
		int* PreIndexMappedLengthArray, unsigned int* PreIndexIntervalStartArray,
		unsigned int* PreIndexIntervalEndArray)
	{
		//cout << "mapMain_SegInfo_preIndex function starts ...... "<< endl; 
		//input : read, sa, up, down, next, chrom; 
		//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc 
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
			
			//cout << "stop_loc_overall: " << stop_loc_overall << endl;

			if(read_length - stop_loc_overall < INDEX_KMER_LENGTH)
			{
				KmerSearchFound = false;
			}
			else
			{

				KmerToSearchInIndexStr = readStringStr.substr(stop_loc_overall, INDEX_KMER_LENGTH);
				//cout << "KmerToSearchInIndexStr: " << KmerToSearchInIndexStr << endl;
				KmerSearchFound = this->getIndexInterval_PreIndex(KmerToSearchInIndexStr, &KmerMappedLength,
					&KmerSearchIndexIntervalStart, &KmerSearchIndexIntervalEnd, PreIndexMappedLengthArray, 
					PreIndexIntervalStartArray,	PreIndexIntervalEndArray);
				//cout << "KmerMappedLength: " << KmerMappedLength << endl;

			}

			//cout << "KmerSearchFound: " << KmerSearchFound << endl;

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
		   	 		queryFound = false;
		   	 		stop_loc = 0;
		   	 		segment_align_SArange[0] = 1;
		   	 		segment_align_SArange[1] = 0;
		   	 		segment_align_rangeNum = 0;
		   	 		queryFound = false;   	 			
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

	   	 	//cout << "interval_begin: " << interval_begin << endl << 
	   	 	//	"interval_end: " << interval_end << endl;

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
	 				lcp_length = getlcp(interval_begin, interval_end, lcpCompress, child, verifyChild//child_up, child_down
	 					);
					Min = min(lcp_length, read_length);
					
					unsigned int loc_pos = 0;

					//////////////////////////// changed for preIndex  ///////////////////////////////
	    	   	 	if(!KmerSearchFound)
	   	 			{
		            	for(loc_pos = 0; loc_pos < Min - c_old; loc_pos++)
		            	{
		            		queryFound = (*(read_local+c_old+loc_pos) == *(chrom+sa[interval_begin]+c_old+loc_pos));
		            		if (!queryFound)
		            		{	
		            			break;
		            		}
		            	}
	            	}
	            	else
	            	{
	            		queryFound = true;
		            	for(loc_pos = KmerMappedLength-1; loc_pos < Min - c_old; loc_pos++)
		            	{
		            		queryFound = (*(read_local+c_old+loc_pos) == *(chrom+sa[interval_begin]+c_old+loc_pos));
		            		if (!queryFound)
		            		{	
		            			break;
		            		}
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
					if (c == read_length)
					{				
						break;			
					}	
					unsigned int interval_begin_ori = interval_begin;
					unsigned int interval_end_ori = interval_end;
			    	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next, 
			    		chrom, verifyChild);

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

	    	   	 	if(!KmerSearchFound)
	   	 			{
		            	for(loc_pos = 0; loc_pos < read_length - c; loc_pos++)
		            	{
		            		queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
		            		if (!queryFound)
		            			break;
		            	}
		            }
		            else
		            {
		            	queryFound = true;
		            	for(loc_pos = KmerMappedLength-1; loc_pos < read_length - c; loc_pos++)
		            	{
		            		queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
		            		if (!queryFound)
		            			break;
		            	}
		            }


		    		if(queryFound) 
		    		{
		    		}
		    		else 
		    		{ 
		    			stop_loc = c+loc_pos;
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
	    		unsigned int tmpSegLength = read_length - stop_loc_overall;

				if(tmpSegLength >= minValSegLength)
				{
					*valLength = *valLength + tmpSegLength;
				}

	    		norSegmentLength[norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;

	    		*(norSegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
	    		*(norSegmentAlignNum + norSegmentNum - 1) = segment_align_rangeNum;				
				

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
				{    			
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				segment_length = read_length-stop_loc_overall;

				break;
			}
			else 
			{    
				(norSegmentNum) ++;
				if(norSegmentNum > (int)(read_length/5))
				{
					//debugln("map error, too many segments, there may exist too many Ns");
					return false;
				}
				norSegmentLength[norSegmentNum - 1] = stop_loc;

				if(stop_loc >= minValSegLength )
				{
					*valLength = *valLength + stop_loc;
				}

				norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[norSegmentNum - 1] = segment_align_rangeNum;

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
			    {    			
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				unsigned int stop_loc_overall_ori = stop_loc_overall;
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;

	    		segment_length = stop_loc;

			}		
	   	}
		//cout << "mapMain_SegInfo function starts 3 ..." << endl;
		#ifdef DEBUG
	   	//cout << " " << endl;
	    cout << "# of Segments = " << norSegmentNum << endl;//; debugln();

	   	for(unsigned int k1 = 0; k1 < norSegmentNum; k1++)
	   	{
	  		cout << "segment "<< k1+1 << ": " << norSegmentLocInRead[k1] << "~"   
	  		<< (norSegmentLocInRead[k1] + norSegmentLength[k1] - 1) 
	  		<< "  Length: " << norSegmentLength[k1] << " Num: " << norSegmentAlignNum[k1] << endl;
	      	if((norSegmentLength[k1]>=10)&&(norSegmentAlignNum[k1] < CANDALILOC))
	      	{
	      		cout << "\tAlign Location: " << endl;
	      		for(unsigned int k2 = 0; k2 < norSegmentAlignNum[k1]; k2++)
	      		{
	      			getChrLocation(*(norSegmentAlignLoc + k1*CANDALILOC + k2), &align_chr, &align_chr_location);
	      			cout << "\t" 
	      			//<< *(norSegmentAlignLoc + k1*CANDALILOC + k2) 
	      			<<  " " 
	      			<< chr_name[align_chr] << " " << align_chr_location << endl;
	      		}
	      	}
	   	}
	   	cout << "mapMain_SegInfo function starts 4 ..." << endl;
		cout << ("segment_length_max = ") << segment_length_max << endl;
		#endif 
		
		segmentNum = (norSegmentNum);
		mapMain = true;
		//cout << "mapMain_SegInfo_preIndex function ends ...... "<< endl; 
		//debugln("mapMain ended!!!");
		return mapMain;
	}
	*/
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
		//unsigned int segment_length = 0; 
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
					//cout << "$$-C\nelse! interval_begin == interval_end" << endl;
					unsigned int loc_pos = 0;
					//cout << "read_length: " << read_length << endl;
					//cout << "read_length - c - stop_loc_overall " << read_length - c - stop_loc_overall << endl; 
		            	//for(loc_pos = 0; loc_pos < read_length - c; loc_pos++) // readLength - c -x = 21
		            	for(loc_pos = 0; loc_pos < read_length - c - stop_loc_overall; loc_pos++)
		            	{
		            		//cout << "loc_pos: " << loc_pos << endl;
		            		//cout << "*(read_local+c+loc_pos): " << *(read_local+c+loc_pos) << endl;
		            		//cout << "*(chrom+sa[interval_begin]+c+loc_pos)" << *(chrom+sa[interval_begin]+c+loc_pos) << endl; 
		            		queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
		            		if (!queryFound)
		            			break;
		            	}

		            //cout << "..... queryFound: " << queryFound << endl;
		            //cout << "..... stop at loc_pos: " << loc_pos << endl; 
		    		if(queryFound) 
		    		{
		    			//cout << "queryFound !!!!!!!!!!!!!!" << endl;
		    		}
		    		else 
		    		{ 
		    			//cout << " not query Found !!!!!!!!!!!" << endl;
		    			//cout << "c: " << c << " loc_pos: " << loc_pos << endl;
		    			stop_loc = c+loc_pos;
		    			//cout << "stop_loc: " << stop_loc << endl;
		    		}	
	          		segment_align_SArange[0] = interval_begin;
	            	segment_align_SArange[1] = interval_end;
	            	segment_align_rangeNum = interval_end - interval_begin + 1;   	
			    	//cout << "/////// interval_begin: " << interval_begin << endl;
			    	//cout << "/////// interval_end: " << interval_end << endl;
			    	//cout << "break happens !!!!!!!!!!!!" << endl;
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

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
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

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
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

	string segInfoStr(Index_Info* indexInfo)
	{
		string segInfoStr;
		segInfoStr += "\nsegment Info: \n";
		unsigned int align_chr, align_chr_location;
	   	for(unsigned int k1 = 0; k1 < segmentNum; k1++)
	   	{
	  		segInfoStr = segInfoStr + "... segment " + int_to_str(k1+1) + ": " + int_to_str(norSegmentLocInRead[k1]) + "~"   
	  		+ int_to_str(norSegmentLocInRead[k1] + norSegmentLength[k1] - 1) + 
	  		"  Length: " + int_to_str(norSegmentLength[k1]) + " Num: " + int_to_str(norSegmentAlignNum[k1]) + "\n";

	      	if(//(norSegmentLength[k1]>=10)&&
	      		(norSegmentAlignNum[k1] < CANDALILOC))
	      	{
	      		//segInfoStr += "...... Align Location: \n";
	      		for(unsigned int k2 = 0; k2 < norSegmentAlignNum[k1]; k2++)
	      		{
	      			indexInfo->getChrLocation(*(norSegmentAlignLoc + k1*CANDALILOC + k2), &align_chr, &align_chr_location);
	      			segInfoStr = segInfoStr + "\t" + int_to_str(k2+1) +
	      			//<< *(norSegmentAlignLoc + k1*CANDALILOC + k2) 
	      			+ ". " 
	      			+ (indexInfo->chrNameStr)[align_chr] + " " + int_to_str(align_chr_location) + "\n";
	      		}
	      	}
	   	}
		return segInfoStr;//+"\n";
	}

	/*bool mapMain_SegInfo_preIndex_toCheck(char *read, unsigned int* sa, BYTE* lcpCompress, 
		unsigned int* child, char* chrom, 
		unsigned int* valLength, BYTE* verifyChild, int readLength, Index_Info* indexInfo,
		int* PreIndexMappedLengthArray, unsigned int* PreIndexIntervalStartArray,
		unsigned int* PreIndexIntervalEndArray)
	{
		//cout << "mapMain_SegInfo_preIndex function starts ...... "<< endl; 
		//input : read, sa, up, down, next, chrom; 
		//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc 
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
			
			//cout << "stop_loc_overall: " << stop_loc_overall << endl;

			if(read_length - stop_loc_overall < INDEX_KMER_LENGTH)
			{
				KmerSearchFound = false;
			}
			else
			{

				KmerToSearchInIndexStr = readStringStr.substr(stop_loc_overall, INDEX_KMER_LENGTH);
				//cout << "KmerToSearchInIndexStr: " << KmerToSearchInIndexStr << endl;
				KmerSearchFound = this->getIndexInterval_PreIndex(KmerToSearchInIndexStr, &KmerMappedLength,
					&KmerSearchIndexIntervalStart, &KmerSearchIndexIntervalEnd, PreIndexMappedLengthArray, 
					PreIndexIntervalStartArray,	PreIndexIntervalEndArray);
				//cout << "KmerMappedLength: " << KmerMappedLength << endl;

			}

			//cout << "KmerSearchFound: " << KmerSearchFound << endl;

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
		   	 		queryFound = false;
		   	 		stop_loc = 0;
		   	 		segment_align_SArange[0] = 1;
		   	 		segment_align_SArange[1] = 0;
		   	 		segment_align_rangeNum = 0;
		   	 		queryFound = false;   	 			
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

	   	 	//cout << "interval_begin: " << interval_begin << endl << 
	   	 	//	"interval_end: " << interval_end << endl;

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
	 				lcp_length = getlcp(interval_begin, interval_end, lcpCompress, child, verifyChild//child_up, child_down
	 					);
					Min = min(lcp_length, read_length);
					
					unsigned int loc_pos = 0;

					//////////////////////////// changed for preIndex  ///////////////////////////////
	    	   	 	if(!KmerSearchFound)
	   	 			{
		            	for(loc_pos = 0; loc_pos < Min - c_old; loc_pos++)
		            	{
		            		queryFound = (*(read_local+c_old+loc_pos) == *(chrom+sa[interval_begin]+c_old+loc_pos));
		            		if (!queryFound)
		            		{	
		            			break;
		            		}
		            	}
	            	}
	            	else
	            	{
		            	for(loc_pos = KmerMappedLength-1; loc_pos < Min - c_old; loc_pos++)
		            	{
		            		queryFound = (*(read_local+c_old+loc_pos) == *(chrom+sa[interval_begin]+c_old+loc_pos));
		            		if (!queryFound)
		            		{	
		            			break;
		            		}
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
					if (c == read_length)
					{				
						break;			
					}	
					unsigned int interval_begin_ori = interval_begin;
					unsigned int interval_end_ori = interval_end;
			    	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next, 
			    		chrom, verifyChild);

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

	    	   	 	if(!KmerSearchFound)
	   	 			{
		            	for(loc_pos = 0; loc_pos < read_length - c; loc_pos++)
		            	{
		            		queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
		            		if (!queryFound)
		            			break;
		            	}
		            }
		            else
		            {
		            	for(loc_pos = KmerMappedLength-1; loc_pos < read_length - c; loc_pos++)
		            	{
		            		queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
		            		if (!queryFound)
		            			break;
		            	}
		            }


		    		if(queryFound) 
		    		{
		    		}
		    		else 
		    		{ 
		    			stop_loc = c+loc_pos;
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
	    		unsigned int tmpSegLength = read_length - stop_loc_overall;

				if(tmpSegLength >= minValSegLength)
				{
					*valLength = *valLength + tmpSegLength;
				}

	    		norSegmentLength[norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;

	    		*(norSegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
	    		*(norSegmentAlignNum + norSegmentNum - 1) = segment_align_rangeNum;				
				

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
				{    			
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				segment_length = read_length-stop_loc_overall;

				break;
			}
			else 
			{    
				(norSegmentNum) ++;
				if(norSegmentNum > (int)(read_length/5))
				{
					//debugln("map error, too many segments, there may exist too many Ns");
					return false;
				}
				norSegmentLength[norSegmentNum - 1] = stop_loc;

				if(stop_loc >= minValSegLength )
				{
					*valLength = *valLength + stop_loc;
				}

				norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[norSegmentNum - 1] = segment_align_rangeNum;

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
			    {    			
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				unsigned int stop_loc_overall_ori = stop_loc_overall;
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;

	    		segment_length = stop_loc;

			}		
	   	}
		//cout << "mapMain_SegInfo function starts 3 ..." << endl;
		#ifdef DEBUG
	   	//cout << " " << endl;
	    cout << "# of Segments = " << norSegmentNum << endl;//; debugln();

	   	for(unsigned int k1 = 0; k1 < norSegmentNum; k1++)
	   	{
	  		cout << "segment "<< k1+1 << ": " << norSegmentLocInRead[k1] << "~"   
	  		<< (norSegmentLocInRead[k1] + norSegmentLength[k1] - 1) 
	  		<< "  Length: " << norSegmentLength[k1] << " Num: " << norSegmentAlignNum[k1] << endl;
	      	if((norSegmentLength[k1]>=10)&&(norSegmentAlignNum[k1] < CANDALILOC))
	      	{
	      		cout << "\tAlign Location: " << endl;
	      		for(unsigned int k2 = 0; k2 < norSegmentAlignNum[k1]; k2++)
	      		{
	      			getChrLocation(*(norSegmentAlignLoc + k1*CANDALILOC + k2), &align_chr, &align_chr_location);
	      			cout << "\t" 
	      			//<< *(norSegmentAlignLoc + k1*CANDALILOC + k2) 
	      			<<  " " 
	      			<< chr_name[align_chr] << " " << align_chr_location << endl;
	      		}
	      	}
	   	}
	   	cout << "mapMain_SegInfo function starts 4 ..." << endl;
		cout << ("segment_length_max = ") << segment_length_max << endl;
		#endif 
		
		segmentNum = (norSegmentNum);
		mapMain = true;
		//cout << "mapMain_SegInfo_preIndex function ends ...... "<< endl; 
		//debugln("mapMain ended!!!");
		return mapMain;
	}*/

	unsigned int getPreIndexNO(const string& readPreStr)
	{
		int preIndexStrSize = readPreStr.length();

		unsigned int preIndexNO = 0;

		int baseForCount = 1;

		for(int tmp = preIndexStrSize - 1; tmp >= 0; tmp--)
		{
			//char tmpChar = readPreStr.at(tmp);
			//unsigned int tmpArrayIndex = tmpChar - 'A';
			//baseChar2int(tmpChar);
			//preIndexNO = preIndexNO + baseChar2int(tmpChar) * baseForCount;
			preIndexNO = preIndexNO + baseChar2intArray[(readPreStr.at(tmp) - 'A')] * baseForCount;
			baseForCount = baseForCount * 4;
		}		
		return preIndexNO;
	}

	bool getIndexInterval_PreIndex(const string& readPreStr, int* mappedLength, 
		unsigned int* indexIntervalStart, unsigned int* indexIntervalEnd,
		int* PreIndexMappedLengthArray, unsigned int* PreIndexIntervalStartArray,
		unsigned int* PreIndexIntervalEndArray)
	{
		////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////	
		bool getIndexInterval = false;

		if (readPreStr.find("N") != readPreStr.npos)
			return false;

		int preIndexStrSize = readPreStr.length();

		unsigned int preIndexNO = 0;

		int baseForCount = 1;
		for(int tmp = preIndexStrSize - 1; tmp >= 0; tmp--)
		{
			//if()
			//char tmpChar = readPreStr.at(tmp);
			//unsigned int tmpArrayIndex = tmpChar - 'A';
			//baseChar2int(tmpChar);
			//preIndexNO = preIndexNO + baseChar2int(tmpChar) * baseForCount;
			//preIndexNO = preIndexNO + baseChar2intArray[(readPreStr.at(tmp) - 'A')] * baseForCount;
			//preIndexNO = preIndexNO + baseCharCount2intArray[preIndexStrSize - 1 - tmp][(readPreStr.at(tmp) - 'A')];

			//preIndexNO = preIndexNO + baseCharCount2intArray[preIndexStrSize - 1 - tmp][('C' - 'A')];
			//preIndexNO = preIndexNO + baseChar2intArray[tmp];
			preIndexNO = preIndexNO + baseChar2intArray[(readPreStr.at(tmp) - 'A')] * baseForCount;
			baseForCount = baseForCount * 4;
		}

		//if(preIndexNO < 100)
		//{
			(*mappedLength) = PreIndexMappedLengthArray[preIndexNO];
			(*indexIntervalStart) = PreIndexIntervalStartArray[preIndexNO];
			(*indexIntervalEnd) = PreIndexIntervalEndArray[preIndexNO];

		//	getIndexInterval = true;
		//}
		//else
		//{
		//	getIndexInterval = false;
		//}*/	

		getIndexInterval = true;
		return getIndexInterval;
	}

};

/*class LocalSeg_Info
{
public:
	int secondLevelIndexNO;
	int interval_begin;
	int interval_end;
	bool crossLocalIndex;

	unsigned int finalSegmentNum;
	unsigned int finalSegmentLength[SEGMENTNUM];
	unsigned int finalSegmentLocInRead[SEGMENTNUM];
	unsigned int finalSegmentAlignNum[SEGMENTNUM];
	unsigned int egmentAlignLoc[SEGMENTNUM * CANDALILOC];

	unsigned int valSegmentNum;
	unsigned int valSegmentLength[SEGMENTNUM];
	unsigned int valSegmentLocInRead[SEGMENTNUM];
	unsigned int valSegmentAlignNum[SEGMENTNUM];
	unsigned int valSegmentAlignLoc[SEGMENTNUM * CANDALILOC];

	string chromName;
	int chromMapPos;
	int secondLevelIndexSize;
	int longSegMinLength;

	bool filterLocalSegWithInterval()
	{
		finalSegmentNum = valSegmentNum;
		for(int tmp = 0; tmp < valSegemntNum; tmp++)
		{
			finalSegmentLength[tmp] = valSegmentLength[tmp];
			finalSegmentLocInRead[tmp] = 
		}
	}

	LocalSeg_Info(Index_Info* indexInfo)
	{
		longSegMinLength = 15;
		secondLevelIndexSize = indexInfo->secondLevelIndexNormalSize;
	}
	LocalSeg_Info(int minLengthOfLongSeg, Index_Info* indexInfo)
	{
		longSegMinLength = minLengthOfLongSeg;
		secondLevelIndexSize = indexInfo->secondLevelIndexNormalSize;
	}

	LocalSeg_Info(string chromNameStr, int chrMapPos, Index_Info* indexInfo)
	{
		longSegMinLength = 15;
		chromName = chromNameStr;
		chromMapPos = chrMapPos;
		secondLevelIndexSize = indexInfo->secondLevelIndexNormalSize;

		secondLevelIndexNO = chromMapPos/secondLevelIndexSize; ///start from 0
	}

	bool crossLocalIndexOrNot()
	{}

	//bool mapMain_SegInfo(char *read, unsigned int* sa, BYTE* lcpCompress, unsigned int* child, char* chrom, 
	//	unsigned int* valLength, BYTE* verifyChild, int readLength, Index_Info* indexInfo)

	bool mapMain_SegInfo_localCompressedIndex(char *read, unsigned int* sa, 
		BYTE* lcpCompress, unsigned int* child, char* chrom, 
		unsigned int* valLength, BYTE* verifyChild, int readLength, 
		Index_Info* indexInfo, int indexSize)
	{
		//input : read, sa, up, down, next, chrom; 
		//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc 
		unsigned int norSegmentNum = 0;

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
		//unsigned int n = (indexInfo->indexSize);//size of SA
		unsigned int n = indexSize;
		unsigned int norAlignLoc;
		unsigned int align_chr_location;
		unsigned int align_chr;
		*valLength = 0;
		char* read_local = read;
		while (stop_loc_overall < read_length) //- 15)
		{
			segment_num++;
			bool queryFound = true;

	   	 	if((*read_local != 'A')&&(*read_local != 'C')&&(*read_local != 'G')&&(*read_local != 'T')) 
	   	 	{
	   	 		queryFound = false;
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
	   	 	getFirstInterval(*read_local, &interval_begin, &interval_end, child, verifyChild);
	   	 	//cout << "new Segment: " << segment_num << endl;
	   	 	//cout << "char: " << (*read_local) << " firstInterval: " << interval_begin << " ~ " << interval_end << endl;
	   	 	segment_align_SArange[0] = interval_begin;
	   	 	segment_align_SArange[1] = interval_end;
	   	 	segment_align_rangeNum = interval_end - interval_begin + 1;

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
	 				lcp_length = getlcp(interval_begin, interval_end, lcpCompress, child, verifyChild//child_up, child_down
	 					);
	 				//cout << "lcp_length: " << lcp_length << endl;
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
					unsigned int interval_begin_ori = interval_begin;
					unsigned int interval_end_ori = interval_end;
			    	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next, 
			    		chrom, verifyChild);
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
	            	for(loc_pos = 0; loc_pos < read_length - c; loc_pos++)
	            	{
	            		queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
	            		if (!queryFound)
	            			break;
	            	}

		    		if(queryFound) 
		    		{
		    		}
		    		else 
		    		{ 
		    			stop_loc = c+loc_pos;
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
	    		unsigned int tmpSegLength = read_length - stop_loc_overall;

				if(tmpSegLength >= minValSegLength)
				{
					*valLength = *valLength + tmpSegLength;
				}

	    		SegmentLength[norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;

	    		*(SegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
	    		*(SegmentAlignNum + norSegmentNum - 1) = segment_align_rangeNum;				
				

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
				{    			
	    			*(SegmentAlignLoc + (SegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				segment_length = read_length-stop_loc_overall;

				break;
			}
			else 
			{    
				(norSegmentNum) ++;
				if(norSegmentNum > (int)(read_length/5))
				{
					//debugln("map error, too many segments, there may exist too many Ns");
					return false;
				}
				SegmentLength[norSegmentNum - 1] = stop_loc;

				if(stop_loc >= minValSegLength )
				{
					*valLength = *valLength + stop_loc;
				}

				valSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				valSegmentAlignNum[norSegmentNum - 1] = segment_align_rangeNum;

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
			    {    			
	    			*(valSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				unsigned int stop_loc_overall_ori = stop_loc_overall;
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;

	    		segment_length = stop_loc;
				
			}		
	   	}
		//cout << "mapMain_SegInfo function starts 3 ..." << endl;
		#ifdef DEBUG
	   	//cout << " " << endl;
	    cout << "# of Segments = " << norSegmentNum << endl;//; debugln();

	   	for(unsigned int k1 = 0; k1 < norSegmentNum; k1++)
	   	{
	  		cout << "segment "<< k1+1 << ": " << norSegmentLocInRead[k1] << "~"   
	  		<< (norSegmentLocInRead[k1] + norSegmentLength[k1] - 1) 
	  		<< "  Length: " << norSegmentLength[k1] << " Num: " << norSegmentAlignNum[k1] << endl;
	      	if((norSegmentLength[k1]>=10)&&(norSegmentAlignNum[k1] < CANDALILOC))
	      	{
	      		cout << "\tAlign Location: " << endl;
	      		for(unsigned int k2 = 0; k2 < norSegmentAlignNum[k1]; k2++)
	      		{
	      			getChrLocation(*(norSegmentAlignLoc + k1*CANDALILOC + k2), &align_chr, &align_chr_location);
	      			cout << "\t" 
	      			//<< *(norSegmentAlignLoc + k1*CANDALILOC + k2) 
	      			<<  " " 
	      			<< chr_name[align_chr] << " " << align_chr_location << endl;
	      		}
	      	}
	   	}
	   	cout << "mapMain_SegInfo function starts 4 ..." << endl;
		cout << ("segment_length_max = ") << segment_length_max << endl;
		#endif 
		
		valSegmentNum = (norSegmentNum);
		mapMain = true;
		//debugln("mapMain ended!!!");
		return mapMain;
	}




	string segInfoStr(//Index_Info* indexInfo
		const string& chrNameStr, int secondLevelIndexNO)
	{
		string segInfoStr;
		segInfoStr += "\nsegment Info: \n";
		unsigned int align_chr, align_chr_location;
	   	for(unsigned int k1 = 0; k1 < segmentNum; k1++)
	   	{
	  		segInfoStr = segInfoStr + "... segment " + int_to_str(k1+1) + ": " + int_to_str(norSegmentLocInRead[k1]) + "~"   
	  		+ int_to_str(norSegmentLocInRead[k1] + norSegmentLength[k1] - 1) + 
	  		"  Length: " + int_to_str(norSegmentLength[k1]) + " Num: " + int_to_str(norSegmentAlignNum[k1]) + "\n";

	      	if(//(norSegmentLength[k1]>=10)&&
	      		(norSegmentAlignNum[k1] < CANDALILOC))
	      	{
	      		//segInfoStr += "...... Align Location: \n";
	      		for(unsigned int k2 = 0; k2 < norSegmentAlignNum[k1]; k2++)
	      		{
	      			indexInfo->getChrLocation(*(norSegmentAlignLoc + k1*CANDALILOC + k2), &align_chr, &align_chr_location);
	      			segInfoStr = segInfoStr + "\t" + int_to_str(k2+1) +
	      			//<< *(norSegmentAlignLoc + k1*CANDALILOC + k2) 
	      			+ ". " 
	      			+ (indexInfo->chrNameStr)[align_chr] + " " + int_to_str(align_chr_location) + "\n";
	      		}
	      	}
	   	}
		return segInfoStr;//+"\n";
	}
};*/

class Path_Info
{
public:
	//vector< vector<int> > PathVec_frag; // vector < vector <mapLebal > >
	//vector< vector<int> > PathVec_frag_relation; // size = PathVec_frag.size(); 
				//PathVec_frag[tmp].size() == PathVec_frag_relation[tmp].size() + 1
	vector< vector< pair<int, int> > > PathVec_seg; // vector <vector <seg_group, seg_candi> >
	vector< bool > PathValidBoolVec;

	vector < pair<int, pair<int, int> > > validPathVec_toPair; // vector < chromNameInt, <alignPosStart, alignPosEnd> >
	vector <int> validPathVec; // vector <validPathNO in PathVec_seg>

	vector< pair<int, Splice_Info*> > fixedPathVec;
	vector< int > fixedPathMismatchVec;

	vector< bool > PathFixedBoolVec;

	//vector< bool > PathValidBoolVec_afterPairing;

	vector < pair< pair<int, int>, Splice_Info*> > finalPathVec;

	//
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

		//vector< vector< pair<int, int> > >().swap(PathVec_seg);
		//vector< bool >().swap( PathValidBoolVec);
		//vector < pair<int, pair<int, int> > >().swap(validPathVec_toPair);
		//vector <int>().swap(validPathVec);
		//vector< pair<int, Splice_Info*> >().swap(fixedPathVec);
		//vector< int >().swap(fixedPathMismatchVec);
		//vector< bool >().swap(PathFixedBoolVec);
		//vector< bool >().swap(PathValidBoolVec_afterPairing);
		//vector < pair< pair<int, int>, Splice_Info*> >().swap(finalPathVec);

	}

	void getPathInfoForPair(Index_Info* indexInfo, Seg_Info* segInfo)
	{
		for(int tmp = 0; tmp < PathVec_seg.size(); tmp++)
		{
			if(PathValidBoolVec[tmp])
			{
				validPathVec.push_back(tmp);
			}
		}

		for(int tmp = 0; tmp < validPathVec.size(); tmp++)
		{
			int tmpPath = validPathVec[tmp];

			int firstSegGroupNO = (PathVec_seg[tmpPath])[0].first;
			int firstSegCandiNO = (PathVec_seg[tmpPath])[0].second;

			int lastSegGroupNO = (PathVec_seg[tmpPath])[PathVec_seg[tmpPath].size() - 1].first;
			int lastSegCandiNO = (PathVec_seg[tmpPath])[PathVec_seg[tmpPath].size() - 1].second;

			unsigned int mapPos_start_tmp = *(segInfo->norSegmentAlignLoc + firstSegGroupNO * CANDALILOC + firstSegCandiNO);
			unsigned int mapPos_end_tmp = *(segInfo->norSegmentAlignLoc + lastSegGroupNO * CANDALILOC + lastSegCandiNO) 
				+ segInfo->norSegmentLength[lastSegGroupNO] - 1;
			
			unsigned int tmpChrNameInt, tmpChrPosInt;
			indexInfo->getChrLocation(mapPos_start_tmp, &tmpChrNameInt, &tmpChrPosInt);

			int mapChrNameInt = tmpChrNameInt;
			int mapChrPos_start = tmpChrPosInt;
			int mapChrPos_end = (int)(mapPos_end_tmp - mapPos_start_tmp) + mapChrPos_start;

			validPathVec_toPair.push_back(pair < int, pair<int, int> > (mapChrNameInt, pair<int,int> (mapChrPos_start, mapChrPos_end)));
		}
	}

	void getFinalPath(Index_Info* indexInfo, Seg_Info* segInfo, int readLength)
	{
		//cout << "start to get Final path" << endl;
		if(segInfo->segmentNum < 1)
		{
			return;
		}
		//int readLength = segInfo->norSegmentLocInRead[segInfo->segmentNum - 1] + segInfo->norSegmentLength[segInfo->segmentNum - 1] - 1;
		for(int tmpPath = 0; tmpPath < fixedPathVec.size(); tmpPath++)
		{

			int fixedPathNO = fixedPathVec[tmpPath].first;

			int tmpPath1stSegGroupNO = (PathVec_seg[fixedPathNO])[0].first;
			int tmpPath1stSegCandiNO = (PathVec_seg[fixedPathNO])[0].second;

			int tmpPathElementSize = (PathVec_seg[fixedPathNO]).size();
			int tmpPathLastSegGroupNO = (PathVec_seg[fixedPathNO])[tmpPathElementSize - 1].first;

			unsigned int PathMapPos = *(segInfo->norSegmentAlignLoc + tmpPath1stSegGroupNO * CANDALILOC + tmpPath1stSegCandiNO);
			unsigned int tmpChrNameInt, tmpChrPosInt;
			indexInfo->getChrLocation(PathMapPos, &tmpChrNameInt, &tmpChrPosInt);
			//string tmpChrNameStr = indexInfo->chrNameStr[tmpChrNameInt];
			//int tmpSegmentMapPos = tmpChrPosInt;
			int tmpUnfixedHeadLength = segInfo->norSegmentLocInRead[tmpPath1stSegGroupNO] - 1;
			//int readLength = segInfo->norSegmentLocInRead[segInfo->segmentNum - 1] + segInfo->norSegmentLength[segInfo->segmentNum - 1] - 1;
			int tmpUnfixedTailLength = readLength - (segInfo->norSegmentLocInRead[tmpPathLastSegGroupNO] 
										+ segInfo->norSegmentLength[tmpPathLastSegGroupNO] - 1);
			//cout << "start to get JumpCode" << endl;
			//cout << tmpUnfixedHeadLength << endl << tmpUnfixedTailLength << endl;
			


			//cout << "finish getting jumpCode" << endl;
			Splice_Info* tmpSpliceInfo = new Splice_Info();
			tmpSpliceInfo->jump_code.clear(); 

			if(tmpUnfixedHeadLength > 0)
			{
				Jump_Code tmpHeadJumpCode(tmpUnfixedHeadLength, "S");
				tmpSpliceInfo->jump_code.push_back(tmpHeadJumpCode);
			}

			tmpSpliceInfo->appendJumpCode(fixedPathVec[tmpPath].second);

			if(tmpUnfixedTailLength > 0)
			{
				Jump_Code tmpTailJumpCode(tmpUnfixedTailLength, "S");
				tmpSpliceInfo->jump_code.push_back(tmpTailJumpCode);
			}

			tmpSpliceInfo->getFinalJumpCode();
			
			int tmpPathFinalMapPos = tmpChrPosInt;
			int tmpPathFinalMapChr = tmpChrNameInt;
			
			finalPathVec.push_back(pair< pair<int, int>, Splice_Info*> (pair<int,int> (tmpPathFinalMapChr, tmpPathFinalMapPos), tmpSpliceInfo) );
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

	string finalFixedPathStr(Index_Info* indexInfo)
	{
		string finalFixedfixedPathStr = "final Fixed Path Info: \n";
		for(int tmpPath = 0; tmpPath < finalPathVec.size(); tmpPath++)
		{

			int tmpChrNameInt = (finalPathVec[tmpPath].first).first;
			string tmpChrNameStr = indexInfo->chrNameStr[tmpChrNameInt];
			int tmpSegmentMapPos = (finalPathVec[tmpPath].first).second;		
			
			finalFixedfixedPathStr = finalFixedfixedPathStr + int_to_str(tmpPath+1) + " ... fixed Path  " + ": "		
				+ tmpChrNameStr + " " + int_to_str(tmpSegmentMapPos) + " " + (finalPathVec[tmpPath].second)->printFinalJumpCode();
			finalFixedfixedPathStr += "\n";
		}

		return finalFixedfixedPathStr;
	}

	string fixedPathVecStr(Index_Info* indexInfo, Seg_Info* segInfo)
	{
		string fixedPathStr = "fixed Path Info: \n";


		for(int tmpPath = 0; tmpPath < fixedPathVec.size(); tmpPath++)
		{

			int fixedPathNO = fixedPathVec[tmpPath].first;

			int tmpPath1stSegGroupNO = (PathVec_seg[fixedPathNO])[0].first;
			int tmpPath1stSegCandiNO = (PathVec_seg[fixedPathNO])[0].second;
			unsigned int PathMapPos = *(segInfo->norSegmentAlignLoc + tmpPath1stSegGroupNO * CANDALILOC + tmpPath1stSegCandiNO);
			unsigned int tmpChrNameInt, tmpChrPosInt;
			indexInfo->getChrLocation(PathMapPos, &tmpChrNameInt, &tmpChrPosInt);
			string tmpChrNameStr = indexInfo->chrNameStr[tmpChrNameInt];
			int tmpSegmentMapPos = tmpChrPosInt;			
			
			fixedPathStr = fixedPathStr + int_to_str(tmpPath+1) + " ... fixed Path " + int_to_str(fixedPathNO+1) + ": "			
				+ tmpChrNameStr + " " + int_to_str(tmpSegmentMapPos) + " " + (fixedPathVec[tmpPath].second)->printJumpCode();
			fixedPathStr += "\n";
		}

		return fixedPathStr;
	}

	/*void getFragInfoAndPossiblePath(Seg_Info* segInfo, Frag_Info* fragInfo)
	{
		vector<bool> mapLabelCheckedVec(fragInfo->mapLabelNum);
		for(int tmp = 0; tmp < fragInfo->mapLabelNum; tmp++)
		{
			mapLabelCheckedVec[tmp] = false; 
		}
		for(int tmpMapLabelNO = 0; 
			tmpMapLabelNO < fragInfo->mapLabelNum;
			tmpMapLabelNO++)
		{
			if(mapLabelCheckedVec[tmpMapLabelNO])
			{
				continue;
			}

		}
	}*/

	void addNewSegGroupToCurrentPathInfoVec(Seg_Info* segInfo, int segGroupNO)
	{
		bool longSegBool = segInfo->checkSegLongOrNot(segGroupNO);
		//cout << "start to add segGroupNO: " << segGroupNO << endl;
		int currentPathNum = PathVec_seg.size();
		for(int tmpSegCandiLoc = 0; tmpSegCandiLoc < segInfo->norSegmentAlignNum[segGroupNO]; tmpSegCandiLoc++)
		{
			//this->addNewSegCandiToCurrentPathInfoVec_uniquePath(segInfo, segGroupNO, tmpSegCandiLoc, longSegBool, currentPathNum);
			//this->addNewSegCandiToCurrentPathInfoVec_allPath(segInfo, segGroupNO, tmpSegCandiLoc, longSegBool, currentPathNum);
			//cout << "start to add segCandiNO: " << tmpSegCandiLoc << endl;
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

	int chooseBestCurrentPathForNewSegCandi(Seg_Info* segInfo, int segGroupNO, int segCandiNO, int currentPathNum)
	{
		//int bestCurrentPath = 1000000;
		int minDistance_value = 900000;
		int minDistance_path = 1000;
		for(int tmpPathElementNO = 0; tmpPathElementNO < currentPathNum; tmpPathElementNO++)
		{
			int tmpPathElementVecSize = (PathVec_seg[tmpPathElementNO]).size();
			int tmpPathLastSegGroupNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).first;
			int tmpPathLastSegCandiNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).second;
			int tmpSegDistance = segInfo->distanceBetweenSegment(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, 
				segGroupNO, segCandiNO);
			if(tmpSegDistance <= minDistance_value)
			{
				minDistance_value = tmpSegDistance;
				minDistance_path = tmpPathElementNO;
			}	
		}
		return minDistance_path;
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
				/*if(relatedOrNot)
				{
					relatedToSomePath = true;
					
					PathValidBoolVec[tmpPathElementNO] = false;
					this->copyOldPath_AddSegCandi_AddNewPath(tmpPathElementNO, segGroupNO, segCandiNO);
					relatedToSomePath = true;
					//(PathVec_seg[tmpPathElementNO]).push_back(pair<int,int> (segGroupNO, segCandiNO));
				}*/
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

	void addNewSegCandiToCurrentPathInfoVec_uniquePath(Seg_Info* segInfo, int segGroupNO, int segCandiNO, bool longSegBool, int currentPathNum)
	{
		int targetPath = this->chooseBestCurrentPathForNewSegCandi(segInfo, segGroupNO, segCandiNO, currentPathNum);
		if(targetPath <= 100)
		{
			PathValidBoolVec[targetPath] = false;
			this->copyOldPath_AddSegCandi_AddNewPath(targetPath, segGroupNO, segCandiNO);			
		}
		else if(longSegBool)
		{
			vector< pair<int,int> > tmpPathElementVec;
			tmpPathElementVec.push_back( pair<int,int> (segGroupNO, segCandiNO) );	
			PathVec_seg.push_back(tmpPathElementVec);
			PathValidBoolVec.push_back(true);		
			vector< pair<int,int> >().swap(tmpPathElementVec);
		}
		else
		{

		}
	}

	void addNewSegCandiToCurrentPathInfoVec_allPath(Seg_Info* segInfo, int segGroupNO, int segCandiNO, bool longSegBool, int currentPathNum)
	{
		bool relatedToSomePath = false;
		//int currentPathNum = PathVec_seg.size();
		for(int tmpPathElementNO = 0; tmpPathElementNO < currentPathNum; tmpPathElementNO ++)
		{

			int tmpPathElementVecSize = (PathVec_seg[tmpPathElementNO]).size();
			int tmpPathLastSegGroupNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).first;
			int tmpPathLastSegCandiNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).second;
			bool relatedOrNot = segInfo->checkSegRelatedOrNot(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, 
				segGroupNO, segCandiNO);
			//cout << endl << "relatedOrNot: " << relatedOrNot << endl << "group1: " << tmpPathLastSegGroupNO 
			//<< " candi1: " << tmpPathLastSegCandiNO << endl << "group2: " << segGroupNO << " candi2: " << segCandiNO << endl;
			if(relatedOrNot)
			{
				relatedToSomePath = true;
				
				PathValidBoolVec[tmpPathElementNO] = false;
				this->copyOldPath_AddSegCandi_AddNewPath(tmpPathElementNO, segGroupNO, segCandiNO);
				relatedToSomePath = true;
				//(PathVec_seg[tmpPathElementNO]).push_back(pair<int,int> (segGroupNO, segCandiNO));
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

	bool getPossiPathFromSeg_incompleteEnd(Seg_Info* segInfo)
	{
		bool possiPathExists = false;
		int firstLongSegNO = -1;//segInfo->getFirstLongSegNO();// no limit on starting from first long segment
			
		for(int tmp = 0; tmp < segInfo->segmentNum; tmp++)	
		{
			if((segInfo->norSegmentAlignNum[tmp] < 10) || (segInfo->norSegmentLength[tmp] > 20))
			{
				firstLongSegNO = tmp;
			}
		}

		if(firstLongSegNO < 0)//||()
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

	string possiPathStr()
	{
		string possiPathStr = "Path Info: \n";
		for(int tmpPath = 0; tmpPath < PathVec_seg.size(); tmpPath++)
		{
			possiPathStr = possiPathStr + "... Path " + int_to_str(tmpPath+1) + ": ";
			for(int tmpPathElement = 0; tmpPathElement < PathVec_seg[tmpPath].size(); tmpPathElement++)
			{
				int tmpSegGroupNO = ((PathVec_seg[tmpPath])[tmpPathElement]).first;
				int tmpSegCandiNO = ((PathVec_seg[tmpPath])[tmpPathElement]).second;
				possiPathStr = possiPathStr + int_to_str(tmpSegGroupNO+1) + "," + int_to_str(tmpSegCandiNO+1) + "--";
			}
			if(PathValidBoolVec[tmpPath])
			{
				possiPathStr += " -- true";
			}
			else
			{
				possiPathStr += " -- false";
			}
			possiPathStr += "\n";
		}
		return possiPathStr;
	}

};


class Path_Pair_Info
{
public:

	vector < pair<int, int> > pathPair; // vec <int, int> (pathNO in pathNorVec, pathNO in pathRcmVec); 

	vector < pair<int, int> > finalPathPair;

	Path_Pair_Info()
	{

	}

	void getFinalPathPair(Path_Info* pathInfoNor, Path_Info* pathInfoRcm)
	{
		for(int tmp = 0; tmp < pathPair.size(); tmp++)
		{
			int pathPair_nor = pathPair[tmp].first;
			int pathPair_rcm = pathPair[tmp].second;
			if(((pathInfoNor->PathFixedBoolVec)[pathPair_nor])
				&&((pathInfoRcm->PathFixedBoolVec)[pathPair_rcm]))
			{
				finalPathPair.push_back(pair<int,int> (pathPair_nor, pathPair_rcm));
			} 
		}
	}

	bool checkPairingOrNot(Path_Info* pathInfoNor, Path_Info* pathInfoRcm, 
		int pathInValidPathVec_nor, int pathInValidPathVec_rcm)
	{
		bool pairBool = false;
		int pathNor_chrNameInt = (pathInfoNor->validPathVec_toPair)[pathInValidPathVec_nor].first;
		int pathRcm_chrNameInt = (pathInfoRcm->validPathVec_toPair)[pathInValidPathVec_rcm].first;
		if(pathNor_chrNameInt == pathRcm_chrNameInt)
		{
			int pathNor_chrPosStart = ((pathInfoNor->validPathVec_toPair)[pathInValidPathVec_nor].second).first;
			int pathNor_chrPosEnd = ((pathInfoNor->validPathVec_toPair)[pathInValidPathVec_nor].second).second;
			int pathRcm_chrPosStart = ((pathInfoRcm->validPathVec_toPair)[pathInValidPathVec_rcm].second).first;
			int pathRcm_chrPosEnd = ((pathInfoRcm->validPathVec_toPair)[pathInValidPathVec_rcm].second).second;

			if(pathRcm_chrPosStart >= pathNor_chrPosEnd)
			{
				int distance  = pathRcm_chrPosStart - pathNor_chrPosEnd;
				if(distance < PAIR_READ_DISTANCE_MAX)
				{
					pairBool = true;
				}
				else
				{
					pairBool = false;
				}
			}
			else
			{
				if((pathNor_chrPosStart - PAIR_READ_DISTANCE_MAX < pathRcm_chrPosStart) 
					&& (pathNor_chrPosEnd < pathRcm_chrPosEnd + PAIR_READ_DISTANCE_MAX))
				{
					pairBool = true;
				}
				else
				{
					pairBool = false;
				}
			}
		}
		else
		{
			pairBool = false;
		}
		return pairBool;
	}

	void pairingPath(Index_Info* indexInfo, Seg_Info* segInfoNor, Seg_Info* segInfoRcm, Path_Info* pathInfoNor, Path_Info* pathInfoRcm)
	{
		//cout << "start to get path info " << endl;
		pathInfoNor->getPathInfoForPair(indexInfo, segInfoNor);
		pathInfoRcm->getPathInfoForPair(indexInfo, segInfoRcm);
		//cout << "start to pair path " << endl;

		//cout << "(pathInfoNor->validPathVec).size(): " << (pathInfoNor->validPathVec).size() << endl;

		for(int tmp = 0; tmp < (pathInfoNor->validPathVec).size(); tmp++)
		{
			//cout << "tmp: " << tmp << endl;
			int tmpPathNor = (pathInfoNor->validPathVec)[tmp];

			for(int tmp2 = 0; tmp2 < (pathInfoRcm->validPathVec).size(); tmp2++)
			{
				//cout << "tmp2: " << tmp2 << endl;
				int tmpPathRcm = (pathInfoRcm->validPathVec)[tmp2];
				bool pairOrNot = this->checkPairingOrNot(pathInfoNor, pathInfoRcm, tmp, tmp2);
				//cout << "pairOrNot: " << pairOrNot << endl;
				if(pairOrNot)
				{
					pathPair.push_back(pair<int, int>(tmpPathNor, tmpPathRcm));
				}
			}
		}
		//cout << "finish pairingPath function ... " << endl;
	}	

	string pairingPathStr()
	{
		//cout << "start pairingPathStr\npathPair.size(): " << pathPair.size() << endl;
		string pairingPathStr = "## Paired Path Info: \n";
		for(int tmp = 0; tmp < pathPair.size(); tmp++)
		{
			//cout << "tmp: " << tmp << endl;
			pairingPathStr = pairingPathStr + ".... Pair " + int_to_str(tmp+1) + ": " 
				+ "Path_Nor " + int_to_str(pathPair[tmp].first + 1) + " -- " 
				+ "Path_Rcm " + int_to_str(pathPair[tmp].second + 1) + "\n";
			//cout << "pairingPathStr" << pairingPathStr << endl;
		} 
		//cout << "finish pairingPathStr" << endl;
		return pairingPathStr;
	}

	string finalPairingPathStr()
	{
		string pairingPathStr = "## Final Paired Path Info: \n";
		for(int tmp = 0; tmp < finalPathPair.size(); tmp++)
		{
			pairingPathStr = pairingPathStr + ".... Final Pair " + int_to_str(tmp+1) + ": " 
				+ "Path_Nor " + int_to_str(finalPathPair[tmp].first + 1) + " -- " 
				+ "Path_Rcm " + int_to_str(finalPathPair[tmp].second + 1) + "\n";
		} 
		return pairingPathStr;
	}
};
