
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

	void addNewSegment(unsigned int length, unsigned int locationInRead, unsigned int alignmentNumber)
	{
		_segments.push_back(new Segment(length, locationInRead, alignmentNumber));
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
	   	 		
	   	 		addNewSegment(1, stop_loc_overall + 1, 0);

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
				addNewSegment(read.length() - stop_loc_overall, stop_loc_overall, segment_align_rangeNum);
				Segment* segment = getCurrentSegment();
				
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

				addNewSegment(stop_loc, stop_loc_overall, segment_align_rangeNum);
				Segment* segment = getCurrentSegment();

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

	void addNewSegment(unsigned int length, unsigned int locationInRead, unsigned int alignmentNumber)
	{
		_segments.push_back(new Segment(length, locationInRead, alignmentNumber));
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

	int getGroupNumber(Segment* value)
	{
		for(int i=0; i<getNumberOfSegments();i++)
			if(getSegment(i)==value)
				return i;
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

	bool isLong(int index)
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
			Segment* segment, int secondCandidateNumber)
	{
		unsigned int alignLoc1 = _segments[firstSegmentIndex]->getAlignmentLocation(firstCandidateNumber);
		unsigned int alignLoc2 = segment->getAlignmentLocation(secondCandidateNumber);

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

	   	 	// If we don't find a match then skip 1 nucleotide, create
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

					addNewSegment(1, stop_loc_overall, 0);

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
					if(walkSteps >= minimum - oldMinimum)
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
						read.getSequence().at(minimum + walkSteps) == chrom->getReference()[chrom->getSuffixArray()[interval_begin] + minimum + walkSteps];
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
				addNewSegment(read.length() - stop_loc_overall, stop_loc_overall, segment_align_rangeNum);
	   	 		Segment* segment = getCurrentSegment();
				for (unsigned int i=0; i < min(segment_align_rangeNum, (unsigned int)CANDALILOC); i++)
					segment->setAlignmentLocation(chrom->getSuffixArray()[segment_align_SArange[0] + i] + 1, i);
				break;
			}
			else
			{
				// if read-align failed at some location, then restart from that location
				addNewSegment(stop_loc, stop_loc_overall, segment_align_rangeNum);
				Segment* segment = getCurrentSegment();

				for (unsigned int i=0; i<min(segment_align_rangeNum, (unsigned int)CANDALILOC); i++)
					segment->setAlignmentLocation(chrom->getSuffixArray()[segment_align_SArange[0] + i] + 1,i);

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
#endif
