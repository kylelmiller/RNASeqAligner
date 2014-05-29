
#ifndef __SEG_INFO_H_INCLUDED__
#define __SEG_INFO_H_INCLUDED__

#include <stdlib.h>
#include <stdio.h>

#include "secondLevelChromosome.h"
#include "chromosome.h"
#include "segment.h"

using namespace std;

class Seg2ndOri_Info
{
private:

	int _longSegMinLength;
	vector<Segment*> _segments;

	/*
	 * Returns the current segment
	 */
	Segment* getCurrentSegment()
	{
		return _segments.size() != 0
			? _segments[_segments.size() - 1]
			: NULL;
	}

	void addNewSegment()
	{
		_segments.push_back(new Segment());
	}

	int getNumberOfSegments()
	{
		return _segments.size();
	}

public:
	Seg2ndOri_Info()
	{
		_longSegMinLength = 20;
	}

	~Seg2ndOri_Info()
	{
		for (vector<Segment*>::iterator it = _segments.begin(); it != _segments.end(); ++it)
			delete *it;
		_segments.clear();
	}

	/*
	 * FIX ME - 5/28/14 KLM
	 * This will take the place of mapMainSecondLevel_compressedIndex
	 */
	bool mapMainSecondLevel_compressedIndex(Read read, SecondLevelChromosome* secondLevelChromosome)
	{
		unsigned int stop_loc = 0; // location in one segment for iterations
		unsigned int stop_loc_overall = 0; //location in the whole read for iterations
		unsigned int segment_length = 0;
		unsigned int segment_align_SArange[2] = {0,0};//legal align location in SA
		unsigned int segment_align_rangeNum = 0;
		unsigned int interval_begin, interval_end;
		char* read_local = (char*)read.getSequence().c_str();

		while (stop_loc_overall < read.length())
		{
			bool queryFound = true;

	   	 	if((*read_local != 'A')&&(*read_local != 'C')&&(*read_local != 'G')&&(*read_local != 'T')) 
	   	 	{
	   	 		queryFound = false;
	   	 		stop_loc = 1;
	   	 		segment_align_SArange[0] = 1;
	   	 		segment_align_SArange[1] = 0;
	   	 		segment_align_rangeNum = 0;
	   	 		queryFound = false;   	

	   	 		if(getNumberOfSegments() >= SEGMENTNUM)
	   	 			return false;
	   	 		
	   	 		Segment* segment = getCurrentSegment();
	   	 		segment->setLength(1);
	   	 		segment->setLocationInRead(stop_loc_overall + 1);
	   	 		segment->setAlignmentNumber(0);
				addNewSegment();

				stop_loc = 1;	
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;   	 		
	   	 		continue;
	   	 	}
	   	 	unsigned int lcp_length = 0;
	   	 	unsigned int start = 0;
	   	 	unsigned int end = secondLevelChromosome->getIndexInfo()->indexSize - 1;
	   	 	unsigned int minimum;
	   	 	unsigned int walkSteps = 0;

	   	 	secondLevelChromosome->getFirstInterval(*read_local, &interval_begin, &interval_end);

	   	 	segment_align_SArange[0] = interval_begin;
	   	 	segment_align_SArange[1] = interval_end;
	   	 	segment_align_rangeNum = interval_end - interval_begin + 1;

	   	 	while(walkSteps + stop_loc_overall < read.length() && queryFound)
	   	 	{
				if(interval_begin != interval_end)
				{
					minimum = min(
						secondLevelChromosome->getLcpLength(interval_begin, interval_end),
						read.length() - stop_loc_overall);

					unsigned int loc_pos;
	            	for(loc_pos = 0; loc_pos < minimum - walkSteps; loc_pos++)
	            	{
	            		queryFound = (*(read_local+walkSteps+loc_pos) ==
							*(secondLevelChromosome->getReference() + secondLevelChromosome->getSuffixArray()[interval_begin]+walkSteps+loc_pos));

	            		if (!queryFound)
	            			break;
	            	}

	            	if(!queryFound)
	            	{
	            		stop_loc = walkSteps + loc_pos;
	            		break;
	            	}
	            	
	            	walkSteps = minimum;
	            	if(*(read_local+walkSteps) == 'N')
	            	{
	            		queryFound = false; 
	            		stop_loc = walkSteps;
	            		break;
	            	}

					start = interval_begin;
					end = interval_end;

					if (walkSteps + stop_loc_overall == read.length())
						break;

					unsigned int interval_begin_ori = interval_begin;
					unsigned int interval_end_ori = interval_end;
					secondLevelChromosome->getInterval(
						start,
						end,
						walkSteps,
						*(read_local+walkSteps),
						&interval_begin,
						&interval_end);

			    	if(interval_begin > interval_end)
			    	{
			    		queryFound = false;
			    		stop_loc = walkSteps - 1;
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
				} // if(interval_begin != interval_end)
				else 
				{
					unsigned int loc_pos = 0;
	            	for(loc_pos = 0; loc_pos < minimum - walkSteps; loc_pos++)
	            	{
	            		queryFound = (*(read_local+walkSteps+loc_pos) ==
							*(secondLevelChromosome->getReference() +
									secondLevelChromosome->getSuffixArray()[interval_begin]+walkSteps+loc_pos));

	            		if (!queryFound)
	            			break;
	            	}

		    		if(!queryFound)
		    			stop_loc =  walkSteps + loc_pos;

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

   	 		if(getNumberOfSegments() > SEGMENTNUM || getNumberOfSegments() > read.length() / 5)
				return false;

	    	if (queryFound && (interval_end >= interval_begin)) 
	    	{
	   	 		Segment* segment = getCurrentSegment();
	   	 		segment->setLength(read.length() - stop_loc_overall);
	   	 		segment->setLocationInRead(stop_loc_overall + 1);
	   	 		segment->setAlignmentNumber(segment_align_rangeNum);
				addNewSegment();
				
				for (unsigned int i=0; i < min(segment_align_rangeNum, (unsigned int)CANDALILOC); i++)
					segment->setAlignmentLocation(secondLevelChromosome->getSuffixArray()[segment_align_SArange[0] + i] + 1, i);
				break;
			}
			else 
			{    
				#ifdef SEGLENGTH
					segmentLength1[stop_loc]++;
					segmentLength2[stop_loc]++;
				#endif

				Segment* segment = getCurrentSegment();
				segment->setLength(stop_loc);
				segment->setLocationInRead(stop_loc_overall + 1);
				segment->setAlignmentNumber(segment_align_rangeNum);
				addNewSegment();

				for (unsigned int i=0; i<min(segment_align_rangeNum, (unsigned int)CANDALILOC); i++)
					segment->setAlignmentLocation(secondLevelChromosome->getSuffixArray()[segment_align_SArange[0] + i] + 1,i);

				// if read-align failed at some location, then restart from that location
				unsigned int stop_loc_overall_ori = stop_loc_overall;
				read_local = read_local + stop_loc + 1;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;

	    		segment_length = stop_loc;
			}		
	   	}
		return true;
	}
};

class Seg_Info
{
private:

	int _longSegMinLength;
	vector<Segment*> _segments;

	/*
	 * Returns the current segment
	 */
	Segment* getCurrentSegment()
	{
		return _segments.size() != 0
			? _segments[_segments.size() - 1]
			: NULL;
	}

	void addNewSegment()
	{
		_segments.push_back(new Segment());
	}

public:

	// FIX ME THESE METHODS NEED TO GO AWAY KLM 5/29/14
	int getNumberOfSegments()
	{
		return _segments.size();
	}

	Segment* getSegment(int index)
	{
		return _segments[index];
	}
	// END OF METHODS TO DELETE

	Seg_Info()
	{
		_longSegMinLength = 20;
	}

	Seg_Info(Seg2ndOri_Info* other, int mapPosIntervalStart, int mapPosIntervalEnd,
		int chrPosStartIn2ndLevelIndex, Index_Info* indexInfo, const string& chromNameStr)
	{
		_longSegMinLength = 18;

		/* FIX ME - KLM 5/29/14 This needs to get sorted out later
		for(int i=0; i < segmentNum; i++)
		{
			_segments[i] = other->_segments[i];
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
		*/
	}		

	bool checkSegLongOrNot(int index)
	{
		return _segments[index]->isLong();
	}
	
	int checkSegRelation(int firstSegmentIndex, int firstCandidateNumber,
		int secondSegmentIndex, int secondCandidateNumber)
	{
		unsigned int segEndNum1 = firstSegmentIndex;
		unsigned int segStartNum2 = secondSegmentIndex;

		unsigned int alignLoc1 = _segments[firstSegmentIndex]->getAlignmentLocation(firstCandidateNumber);
		unsigned int alignLoc2 = _segments[secondSegmentIndex]->getAlignmentLocation(secondCandidateNumber);

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

	int distanceBetweenSegment(int firstSegmentIndex, int firstCandidateNumber,
			int secondSegmentIndex, int secondCandidateNumber)
	{
		unsigned int alignLoc1 = _segments[firstSegmentIndex]->getAlignmentLocation(firstCandidateNumber);
		unsigned int alignLoc2 = _segments[secondSegmentIndex]->getAlignmentLocation(secondCandidateNumber);

		unsigned int segmentDistance = alignLoc2 >= alignLoc1
			? alignLoc2 - alignLoc1
			: alignLoc1 - alignLoc2;

		return segmentDistance < 300000 ? segmentDistance : 1000000;
	}

	int getFirstLongSegNO()
	{
		for(int i=0; i<getNumberOfSegments(); i++)
			if(_segments[i]->getLength() >= _longSegMinLength && _segments[i]->getAlignmentNumber() <= SEGMENTNUM)
				return i;

		return -1;
	}

	bool mapMain_SegInfo_preIndex(Read read, Chromosome* chrom)
	{
		unsigned int segmentNumber = 0; // current segment number
		unsigned int stop_loc = 0; // location in one segment for iterations
		unsigned int stop_loc_overall = 0; //location in the whole read for iterations
		unsigned int segment_align_SArange[2] = {0,0};//legal align location in SA
		unsigned int segment_align_rangeNum = 0;
		unsigned int interval_begin, interval_end;
		char* localRead = (char*)read.getSequence().c_str();

		while (stop_loc_overall < read.length())
		{
			/////   K-mer search  /////
			bool KmerSearchFound = false;
			int KmerMappedLength;
			unsigned int KmerSearchIndexIntervalStart = 0;
			unsigned int KmerSearchIndexIntervalEnd = 0;

			// If we hit our kmer limit then we fail
			if(read.length() - stop_loc_overall >= INDEX_KMER_LENGTH)
				KmerSearchFound = chrom->getMappedLocation(
					read,
					stop_loc_overall,
					&KmerMappedLength,
					&KmerSearchIndexIntervalStart,
					&KmerSearchIndexIntervalEnd);

	   	 	unsigned int start = 0;
			unsigned int end = chrom->getIndexInfo()->indexSize - 1;
	   	 	unsigned int minimum = 0;

	   	 	// If we don't fine a match then skip 1 nucleotide, create
	   	 	// a new segment and attempt to match another 14 nucleotides
	   	 	if(!KmerSearchFound)
	   	 	{
	   	 		// if the nucleotide is an 'N'
		   	 	if(	(*localRead != 'A') &&
					(*localRead != 'C') &&
					(*localRead != 'G') &&
					(*localRead != 'T'))
		   	 	{
		   	 		if(segmentNumber >= SEGMENTNUM)
		   	 			return false;

		   	 		segment_align_SArange[0] = 1;
		   	 		segment_align_SArange[1] = 0;

		   	 		Segment* segment = getCurrentSegment();
		   	 		segment->setLength(1);
		   	 		segment->setLocationInRead(stop_loc_overall + 1);
		   	 		segment->setAlignmentNumber(0);
					addNewSegment();

					stop_loc = 1;
					localRead = localRead + stop_loc + 1;
					stop_loc_overall = stop_loc_overall + stop_loc + 1;
		   	 		continue;
		   	 	}

		   	 	chrom->getFirstInterval(*localRead, &interval_begin, &interval_end);
	   	 	}
	   	 	else // K-mer found in preIndex base
	   	 	{
	   	 		interval_begin = KmerSearchIndexIntervalStart;
	   	 		interval_end = KmerSearchIndexIntervalEnd;
	   	 	}

	   	 	segment_align_SArange[0] = interval_begin;
	   	 	segment_align_SArange[1] = interval_end;
   	 		segment_align_rangeNum = interval_end - interval_begin + 1;

   	 		while(minimum + stop_loc_overall < read.length())
			{
   	   	 		unsigned int oldMinimum = minimum;

				if(interval_begin != interval_end)
				{
					minimum = min(
						chrom->getLcp(interval_begin, interval_end),
						read.length() - stop_loc_overall);

		   	 		unsigned int walkSteps;
					for(walkSteps=0;
						walkSteps < minimum - oldMinimum &&
						read.getSequence().at(oldMinimum + walkSteps) == chrom->getReference()[chrom->getSuffixArray()[interval_begin] + oldMinimum + walkSteps];
						walkSteps++)
					{
					}

					// If this is true then we know we stopped walking because we hit
					// an unmatched nucleotide
					if(walkSteps != minimum - oldMinimum)
					{
						stop_loc = oldMinimum + walkSteps;
						break;
					}
				}
				else // (interval_begin == interval_end)
				{

					unsigned int walkSteps;
					// Compares the read to the reference, character by character
					for(walkSteps=0;
						walkSteps<read.length() - minimum - stop_loc_overall &&
						read.getSequence().at(*localRead + minimum + walkSteps) == chrom->getReference()[chrom->getSuffixArray()[interval_begin] + minimum + walkSteps];
						walkSteps++)
					{
					}

					// If we didn't match the entire read, note where we stopped
		    		if(walkSteps < read.length() - minimum - stop_loc_overall)
		    			stop_loc = minimum + walkSteps;

		    		// Store the alignment information
	          		segment_align_SArange[0] = interval_begin;
	            	segment_align_SArange[1] = interval_end;
	            	segment_align_rangeNum = interval_end - interval_begin + 1;
		    		break;
				}
			} // while(minimum + stop_loc_overall < read.length())

    		if(getNumberOfSegments() > SEGMENTNUM || getNumberOfSegments() > read.length() / 5)
   	 			return false;

   	 		// Segment Mapping
	    	if (interval_end >= interval_begin)
	    	{
	   	 		Segment* segment = getCurrentSegment();
	   	 		segment->setLength(read.length() - stop_loc_overall);
	   	 		segment->setLocationInRead(stop_loc_overall + 1);
	   	 		segment->setAlignmentNumber(segment_align_rangeNum);
				for (unsigned int i=0; i < min(segment_align_rangeNum, (unsigned int)CANDALILOC); i++)
					segment->setAlignmentLocation(chrom->getSuffixArray()[segment_align_SArange[0] + i] + 1, i);

				addNewSegment();
				break;
			}
			else
			{
				// if read-align failed at some location, then restart from that location
				Segment* segment = getCurrentSegment();
				segment->setLength(stop_loc);
				segment->setLocationInRead(stop_loc_overall + 1);
				segment->setAlignmentNumber(segment_align_rangeNum);

				for (unsigned int i=0; i<min(segment_align_rangeNum, (unsigned int)CANDALILOC); i++)
					segment->setAlignmentLocation(chrom->getSuffixArray()[segment_align_SArange[0] + i] + 1,i);

				addNewSegment();

				localRead = localRead + stop_loc + 1;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;
			}
		} // while(stop_loc_overall < read.length())

		return true;
	}
/*
	bool mapMain_SegInfo_preIndex(char *read, unsigned int* sa, BYTE* lcpCompress, 
		unsigned int* child, char* chrom, 
		unsigned int* valLength, BYTE* verifyChild, int readLength, Index_Info* indexInfo,
		int* PreIndexMappedLengthArray, unsigned int* PreIndexIntervalStartArray,
		unsigned int* PreIndexIntervalEndArray)
	{
		string readStringStr = read;

		unsigned int norSegmentNum = 0;
		unsigned int stop_loc = 0; // location in one segment for iterations
		unsigned int stop_loc_overall = 0; //location in the whole read for iterations
		unsigned int segment_length = 0;
		unsigned int segment_length_max = 0;//used to compare with segment_length for each segment to get the maximum length segment
		unsigned int segment_align_SArange[2] = {0,0};//legal align location in SA
		unsigned int segment_align_rangeNum = 0;
		unsigned int interval_begin, interval_end;
		*valLength = 0;
		char* read_local = read;

		while (stop_loc_overall < readLength) //- 15)
		{
			bool queryFound = true;

			/////   K-mer search  /////
			bool KmerSearchFound = false;
			int KmerMappedLength;
			unsigned int KmerSearchIndexIntervalStart = 0;
			unsigned int KmerSearchIndexIntervalEnd = 0;

			// If we hit our kmer limit then we fail
			if(readLength - stop_loc_overall < INDEX_KMER_LENGTH)
				KmerSearchFound = false;
			else
				KmerSearchFound = this->getIndexInterval_PreIndex(
					readStringStr.substr(stop_loc_overall, INDEX_KMER_LENGTH),
					&KmerMappedLength,
					&KmerSearchIndexIntervalStart,
					&KmerSearchIndexIntervalEnd,
					PreIndexMappedLengthArray,
					PreIndexIntervalStartArray,
					PreIndexIntervalEndArray);

	   	 	unsigned int start = 0,
				end = indexInfo->indexSize - 1,
				c = 0,
				iterateNum = 0;
	   	 	unsigned int Min;

	   	 	if(!KmerSearchFound)
	   	 	{	
		   	 	if(	(*read_local != 'A') &&
					(*read_local != 'C') &&
					(*read_local != 'G') &&
					(*read_local != 'T'))
		   	 	{
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
		   	 		
		   	 		norSegmentNum++;

		   	 		norSegmentLength[norSegmentNum - 1] = 1;
					norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
					norSegmentAlignNum[norSegmentNum - 1] = 0;

					stop_loc = 1;	
					read_local = read_local + stop_loc + 1;
					stop_loc_overall = stop_loc_overall + stop_loc + 1;   	 		
		   	 		continue;		   	 		
		   	 	}

		   	 	getFirstInterval(*read_local, &interval_begin, &interval_end, child, verifyChild);
		   	 	segment_align_SArange[0] = interval_begin;
		   	 	segment_align_SArange[1] = interval_end;
		   	 	segment_align_rangeNum = interval_end - interval_begin + 1;
	   	 	}
	   	 	else // K-mer found in preIndex base
	   	 	{
	   	 		interval_begin = KmerSearchIndexIntervalStart;
	   	 		interval_end = KmerSearchIndexIntervalEnd;
	    	 	segment_align_SArange[0] = interval_begin;
	   	 		segment_align_SArange[1] = interval_end;
	   	 		segment_align_rangeNum = interval_end - interval_begin + 1;
	   	 	}

	   	 	while((c + stop_loc_overall < readLength) && queryFound == true)
	   	 	{
	   	 		iterateNum++;
	   	 		if(iterateNum + stop_loc_overall > readLength)
	   	 			return false;

	   	 		unsigned int c_old = c;

				if(interval_begin != interval_end)
				{ 
					Min = min(
						getlcp(interval_begin, interval_end, lcpCompress, child, verifyChild),
						readLength - stop_loc_overall);
					c = Min;

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

	            	if(*(read_local+c) == 'N')
	            	{
	            		queryFound = false; 
	            		stop_loc = c;
	            		break;
	            	}
					start = interval_begin; end = interval_end;

					if (c + stop_loc_overall == readLength)
						break;

					unsigned int interval_begin_ori = interval_begin;
					unsigned int interval_end_ori = interval_end;

			    	getInterval(start, end, c, *(read_local+c), &interval_begin,
						&interval_end, sa, child, chrom, verifyChild);

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
				}
				else // (interval_begin == interval_end)
				{
					unsigned int loc_pos;

					// Compares the read to the reference, character by character
					for(loc_pos = 0; loc_pos < readLength - c - stop_loc_overall; loc_pos++)
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

			///  SEGMENT MAP RESULT ///

	    	if (queryFound && (interval_end >= interval_begin)) 
	    	{

	    		norSegmentNum++;
	    		if(norSegmentNum > SEGMENTNUM)
	    		{
	    			segmentNum = (SEGMENTNUM);
	    			return false;
	    		}
	   	 		if(norSegmentNum > (int)(readLength / 5))
	   	 		{
	   	 			segmentNum = (int)(readLength / 5);
	   	 			return false;
	   	 		}

	    		unsigned int tmpSegLength = readLength - stop_loc_overall;

				if(tmpSegLength >= minValSegLength)
					*valLength = *valLength + tmpSegLength;

	    		norSegmentLength[norSegmentNum - 1] = tmpSegLength;

	    		*(norSegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
	    		*(norSegmentAlignNum + norSegmentNum - 1) = segment_align_rangeNum;

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum, (unsigned int)CANDALILOC); alignment_num++)
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;

				break;
			}
			else 
			{    
				norSegmentNum++;

				if(norSegmentNum > (int)(readLength / 5) || norSegmentNum > SEGMENTNUM)
				{
					segmentNum = (int)(readLength / 5) > SEGMENTNUM
						? SEGMENTNUM
						: (int)(readLength / 5);

					return false;
				}

				norSegmentLength[norSegmentNum - 1] = stop_loc;

				if(stop_loc >= minValSegLength )
					*valLength = *valLength + stop_loc;

				norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[norSegmentNum - 1] = segment_align_rangeNum;

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum, (unsigned int)CANDALILOC); alignment_num++)
			    {    			
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				unsigned int stop_loc_overall_ori = stop_loc_overall;
				read_local = read_local + stop_loc + 1;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;
			}		
	   	}

		segmentNum = norSegmentNum;
		return true;
	}
*/
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

	~Path_Info()
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

	// FIX ME - THIS METHOD NEEDS TO BE CLEANED UP KLM 5/29/14
	void getFinalPath_extend2HeadTail(Index_Info* indexInfo, Seg_Info* segInfo, Read readSeq_inProcess)
	{
		//cout << "start to get Final path" << endl;
		if(segInfo->getNumberOfSegments() < 1)
			return;

		for(int tmpPath = 0; tmpPath < fixedPathVec.size(); tmpPath++)
		{
			int mismatchNumToAdd = 0;
			int fixedPathNO = fixedPathVec[tmpPath].first;

			int tmpPath1stSegGroupNO = (PathVec_seg[fixedPathNO])[0].first;
			int tmpPath1stSegCandiNO = (PathVec_seg[fixedPathNO])[0].second;

			Segment* currentSegment = segInfo->getSegment(tmpPath1stSegGroupNO);

			int tmpPathElementSize = (PathVec_seg[fixedPathNO]).size();
			int tmpPathLastSegGroupNO = (PathVec_seg[fixedPathNO])[tmpPathElementSize - 1].first;

			unsigned int PathMapPos = currentSegment->getAlignmentLocation(tmpPath1stSegCandiNO);
			unsigned int tmpChrNameInt, tmpChrPosInt;
			indexInfo->getChrLocation(PathMapPos, &tmpChrNameInt, &tmpChrPosInt);
			
			int tmpPathFinalMapPos = tmpChrPosInt;
			int tmpPathFinalMapChr = tmpChrNameInt;

			Splice_Info* tmpSpliceInfo = new Splice_Info();
			tmpSpliceInfo->jump_code.clear(); 

			int tmpUnfixedHeadLength = currentSegment->getLocationInRead();
			string readSubSeqInProcess_head = readSeq_inProcess.getSequence().substr(0, tmpUnfixedHeadLength);
			
			bool scoreStringBool_head; 
			string chromSubSeqInProcess_head; 

			size_t max_mismatch_head = (tmpUnfixedHeadLength)/8 + 1;
			size_t mismatch_bits_head = 0;
			size_t comb_bits_head = 0;

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

			tmpSpliceInfo->appendJumpCode(fixedPathVec[tmpPath].second);

			//////////////////////////// add last jump code /////////////////////////////////////////////////////
			int endMappedBaseMapPos = (fixedPathVec[tmpPath].second)->getEndBaseMapPos_jump_code(tmpPathFinalMapPos);
			int tmpUnfixedTailLength =
				readSeq_inProcess.length() -
				(currentSegment->getLocationInRead() + currentSegment->getLength() - 1);
		
			string readSubSeqInProcess_tail 
				= readSeq_inProcess.getSequence().substr(readSeq_inProcess.length() - tmpUnfixedTailLength, tmpUnfixedTailLength);

			bool scoreStringBool_tail;
			string chromSubSeqInProcess_tail;

			size_t max_mismatch_tail = (tmpUnfixedTailLength)/8 + 1;
			size_t mismatch_bits_tail = 0;
			size_t comb_bits_tail = 0;

			if(endMappedBaseMapPos + tmpUnfixedTailLength >= indexInfo->chromLength[tmpChrNameInt])
			{
				scoreStringBool_tail = false;
			}
			else
			{
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

			tmpSpliceInfo->getFinalJumpCode();			
			
			if((tmpUnfixedHeadLength > 0)&&(scoreStringBool_head))
			{
				tmpPathFinalMapPos = tmpPathFinalMapPos - tmpUnfixedHeadLength;
			}

			int oldMismatchNum = fixedPathMismatchVec[tmpPath];
			int newMismatchNum = oldMismatchNum + mismatchNumToAdd;


			fixedPathMismatchVec[tmpPath] = newMismatchNum;

			finalPathVec.push_back(pair< pair<int, int>, Splice_Info*> (pair<int,int> (tmpPathFinalMapChr, tmpPathFinalMapPos), tmpSpliceInfo) );
		}		
	}

	void addNewSegGroupToCurrentPathInfoVec(Segment* segment)
	{
		for(int i=0; i<segment->getAlignmentNumber(); i++)
		{
			this->addNewSegCandiToCurrentPathInfoVec_matchIndelUniqueSpliceMultiPath(
				segment,
				i,
				segment->isLong(),
				PathVec_seg.size());
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

	// FIX ME - LEFT OFF HERE!!!!
	int minDistanceWithCurrentPath(Segment* segment, int segCandiNO, int currentPathNum, int* minSegNumGap)
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

	void addNewSegCandiToCurrentPathInfoVec_matchIndelUniqueSpliceMultiPath(Segment* segment,
			int segCandiNO, bool longSegBool, int currentPathNum)
	{
		bool relatedToSomePath = false;
		int minSegNumGap = 0;
		int minDistance = minDistanceWithCurrentPath(segment, segCandiNO, currentPathNum, &minSegNumGap);

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
			return false;

		Segment* currentSegment = segInfo->getSegment(firstLongSegNO);

		for(int i=0; i<currentSegment->getAlignmentNumber(); i++)
		{
			vector< pair<int,int> > tmpPathElementVec;
			tmpPathElementVec.push_back( pair<int,int> (firstLongSegNO, i));
			PathVec_seg.push_back(tmpPathElementVec);
			PathValidBoolVec.push_back(true);
			vector< pair<int,int> > ().swap(tmpPathElementVec);
		}

		for(int i=firstLongSegNO + 1; i < segInfo->getNumberOfSegments(); i++)
		{
			if(segInfo->getSegment(i)->getAlignmentNumber() > CANDALILOC)
				continue;

			this->addNewSegGroupToCurrentPathInfoVec(segInfo->getSegment(i));
		}

		return true;
	}
};

#endif
