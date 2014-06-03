
#ifndef __SEG_INFO_H_INCLUDED__
#define __SEG_INFO_H_INCLUDED__

#include <stdlib.h>
#include <stdio.h>

#include "secondLevelChromosome.h"
#include "chromosome.h"
#include "segment.h"
#include "read.h"

using namespace std;

class MappedRead
{
private:

	/*
	 * Constants
	 */
	const unsigned int CHARACTER_SKIP_DISTANCE = 1;

	/*
	 * Member Variables
	 */
	int _longSegMinLength;
	Read* _read;
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

	void mapRead(Chromosome* chrom)
	{
		unsigned int currentStartLocation = 0;
		unsigned int suffixArrayRange = 0;
		unsigned int interval_begin, interval_end;
		char* localRead = (char*)_read->getSequence().c_str();

		// We will continue to attempt to map as long as there is a segment left
		// that meets our minimum segment length
		while (_read->getLength() > currentStartLocation && chrom->couldContainSegment(_read->getLength() - currentStartLocation))
		{
			/////   K-mer search  /////
			int KmerMappedLength;
			unsigned int KmerSearchIndexIntervalStart = 0;
			unsigned int KmerSearchIndexIntervalEnd = 0;

	   	 	// If we don't find a match then skip 1 nucleotide, create
	   	 	// a new segment and attempt to match another 14 nucleotides
	   	 	if(!chrom->getMappedLocation(
				_read,
				currentStartLocation,
				&KmerMappedLength,
				&KmerSearchIndexIntervalStart,
				&KmerSearchIndexIntervalEnd))
	   	 	{
	   	 		// Didn't find anything, skip and attempt to map again
				localRead += CHARACTER_SKIP_DISTANCE;
				currentStartLocation += CHARACTER_SKIP_DISTANCE;
				continue;
	   	 	}

	   	 	// The range of indices in the suffix array which matches our substring
			interval_begin = KmerSearchIndexIntervalStart;
			interval_end = KmerSearchIndexIntervalEnd;
   	 		suffixArrayRange = interval_end - interval_begin + 1;

   	 		unsigned int walkLimit = interval_begin == interval_end
				? _read->getLength() - currentStartLocation
				: min(chrom->getLcp(interval_begin, interval_end), _read->getLength() - currentStartLocation);

			unsigned int walkSteps;
			for(walkSteps=0;
				walkSteps < walkLimit &&
				_read->getSequence()[currentStartLocation + walkSteps] == chrom->getReference()[chrom->getSuffixArray()[interval_begin] + walkSteps];
				walkSteps++)
			{
			}

			// This next line is wrong it should be this
			//addNewSegment(walkSteps, currentStartLocation, suffixArrayRange);
			// Needs to be changed in concert with the rest of the code
			addNewSegment(_read->getLength() - currentStartLocation, currentStartLocation, suffixArrayRange);

			Segment* segment = getCurrentSegment();
			for (unsigned int i=0; i < min(suffixArrayRange, (unsigned int)CANDALILOC); i++)
				segment->setAlignmentLocation(chrom->getSuffixArray()[interval_begin + i] + 1, i);

			// if read-align failed at some location, then restart from that location
			localRead += walkSteps;
			currentStartLocation += walkSteps;

			// FIXME - KLM 6/2/14: So, this is wrong but it's how it originally was.
			// We only map a 14 length segment then call it a day, rather than trying to map
			// as many segments as we can. This needs to be changed in conjunction with
			// the rest of the code
			break;

		} // while(stop_loc_overall < read->length())
	}

	/*
	 * FIXME - 5/28/14 KLM
	 * This will take the place of mapMainSecondLevel_compressedIndex
	 */
	bool mapRead(SecondLevelChromosome* secondLevelChromosome)
	{
		unsigned int stop_loc = 0; // location in one segment for iterations
		unsigned int stop_loc_overall = 0; //location in the whole read for iterations
		unsigned int segment_length = 0;
		unsigned int segment_align_SArange[2] = {0,0};//legal align location in SA
		unsigned int segment_align_rangeNum = 0;
		unsigned int interval_begin, interval_end;
		char* read_local = (char*)_read->getSequence().c_str();

		while (stop_loc_overall < _read->getLength())
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

	   	 	while(walkSteps + stop_loc_overall < _read->getLength() && queryFound)
	   	 	{
				if(interval_begin != interval_end)
				{
					minimum = min(
						secondLevelChromosome->getLcpLength(interval_begin, interval_end),
						_read->getLength() - stop_loc_overall);

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

					if (walkSteps + stop_loc_overall == _read->getLength())
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

   	 		if(getNumberOfSegments() > SEGMENTNUM || getNumberOfSegments() > _read->getLength() / 5)
				return false;

	    	if (queryFound && (interval_end >= interval_begin))
	    	{
				addNewSegment(_read->getLength() - stop_loc_overall, stop_loc_overall, segment_align_rangeNum);
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
public:

	// FIXME THESE METHODS NEED TO GO AWAY KLM 5/29/14
	int getNumberOfSegments()
	{
		return _segments.size();
	}

	Read* getRead()
	{
		return _read;
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

		return -1;
	}
	// END OF METHODS TO DELETE

	MappedRead(Read* read, Chromosome* chrom)
	{
		_longSegMinLength = 20;
		_read = new Read(read);
		mapRead(chrom);
	}

	MappedRead(Read* read, SecondLevelChromosome* secondLevelChromosome)
	{
		_longSegMinLength = 20;
		_read = new Read(read);
		mapRead(secondLevelChromosome);
	}

	~MappedRead()
	{
		for (vector<Segment*>::iterator it = _segments.begin(); it != _segments.end(); ++it)
			delete *it;
		_segments.clear();
		delete _read;
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

	Segment* getFirstLongSegment()
	{
		for(int i=0; i<getNumberOfSegments(); i++)
			if(_segments[i]->getLength() >= _longSegMinLength && _segments[i]->getAlignmentNumber() <= SEGMENTNUM)
				return _segments[i];

		return NULL;
	}

	int getFirstLongSegNO()
	{
		for(int i=0; i<getNumberOfSegments(); i++)
			if(_segments[i]->getLength() >= _longSegMinLength && _segments[i]->getAlignmentNumber() <= SEGMENTNUM)
				return i;

		return -1;
	}
};
#endif
