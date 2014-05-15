#include <string>
#include <string.h>
#include "splice_info.h"

using namespace std;

class Unfixed_Tail
{
public:
	string readName;
	string alignDirection; 
	unsigned int midPartMapPosInWholeGenome;
	string otherCigarString;
	int unfixedTailLength;
	int midPartLength;
	string readSeqOriginal;
	//int readLength;

	int midPartMapPosInChr; // map pos of 1st base in midPart
	string midPartMapChrName;
	int midPartMapChrInt;

	vector<int> possiSJposInRead; // start from 1, end at unfixedHeadLength + bufferLength

	vector<int> possiGTAGpos; // only AG end checked
	vector<int> possiCTACpos; // only AC end checked

	vector< pair<int,int> > GTAGsjPos; // sequence around the SJ has been checked including short anchor and midpart
	vector< pair<int,int> > CTACsjPos; // sequence around the SJ has been checked 

	vector< pair<int,int> > SJposFromRemapping; // // <posInRead, SJsize>

	Unfixed_Tail()
	{//readLength = 100;
	}

	Unfixed_Tail(int tail_length, int mapPos, string mapChrName)
	{
		unfixedTailLength = tail_length;
		midPartMapPosInChr = midPartMapPosInChr;
		midPartMapChrName = mapChrName;
	}

	~Unfixed_Tail()
	{}

	void getUnfixedTailInfoFromRecord(PE_Read_Info* readInfo, bool end1, Alignment_Info* alignInfo, Index_Info* indexInfo)
	{
		if(end1)
		{
			readName = (readInfo->readInfo_pe1).readName;
			readSeqOriginal = (readInfo->readInfo_pe1).readSeq;
		}
		else
		{
			readName = (readInfo->readInfo_pe2).readName;
			readSeqOriginal = (readInfo->readInfo_pe2).readSeq;
		}
		alignDirection = alignInfo->alignDirection;

		midPartMapChrName = alignInfo->alignChromName;
		//cout << "midpartMapChrPos: " << 
		midPartMapPosInChr = alignInfo->getEndMatchedPosInChr();
		midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);

		otherCigarString = alignInfo->otherJumpCodeVec2StrForTail();

		int jumpCodeVecSize = (alignInfo->cigarStringJumpCode).size();
		unfixedTailLength = (alignInfo->cigarStringJumpCode)[jumpCodeVecSize-1].len;
		midPartLength = (alignInfo->cigarStringJumpCode)[jumpCodeVecSize-2].len;
		midPartMapPosInWholeGenome = indexInfo->getWholeGenomeLocation(
			indexInfo->convertStringToInt(midPartMapChrName), midPartMapPosInChr);
	}

	void getUnfixedTailInfoFromRecordWithAlignInfoType(PE_Read_Info* readInfo, int alignInfoType, Alignment_Info* alignInfo, Index_Info* indexInfo)
	{
		if((alignInfoType == 1) || (alignInfoType == 2))
		{
			readName = (readInfo->readInfo_pe1).readName;
			readSeqOriginal = (readInfo->readInfo_pe1).readSeq;
		}
		else
		{
			readName = (readInfo->readInfo_pe2).readName;
			readSeqOriginal = (readInfo->readInfo_pe2).readSeq;
		}
		alignDirection = alignInfo->alignDirection;

		midPartMapChrName = alignInfo->alignChromName;
		//cout << "midpartMapChrPos: " << 
		midPartMapPosInChr = alignInfo->getEndMatchedPosInChr();
		midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);

		otherCigarString = alignInfo->otherJumpCodeVec2StrForTail();

		int jumpCodeVecSize = (alignInfo->cigarStringJumpCode).size();
		unfixedTailLength = (alignInfo->cigarStringJumpCode)[jumpCodeVecSize-1].len;
		midPartLength = (alignInfo->cigarStringJumpCode)[jumpCodeVecSize-2].len;
		midPartMapPosInWholeGenome = indexInfo->getWholeGenomeLocation(
			indexInfo->convertStringToInt(midPartMapChrName), midPartMapPosInChr);
	}


	void getPossibleSJpos(const string& readSeqWithDirection, const string& chromSeq, Index_Info* indexInfo)
	{
		int bufferLength = 4;
		if(bufferLength > midPartLength)
		{
			bufferLength = midPartLength;
		}
		//cout << "readLength: " << readLength << endl;
		//cout << "bufferlength: " << bufferLength << endl;
		int readLength = readSeqWithDirection.length();
		string pendingReadSeq = readSeqWithDirection.substr(
			readLength - unfixedTailLength - bufferLength, unfixedTailLength + bufferLength);
		//string pendingChroSeq = chromSeq.substr(
		//	midPartMapPosInWholeGenome - bufferLength, unfixedTailLength + bufferLength);

		//cout << "midPartChrMapPos: " << midPartMapPosInChr << endl;
		int midPartMapChrNameInt = indexInfo->convertStringToInt(midPartMapChrName);
		//cout << "midPartMapChrNameInt: " << midPartMapChrNameInt << endl;

		//cout << "midPartMapPosInChr-bufferLength" << midPartMapPosInChr-bufferLength << endl;

		//cout << "midPartMapPosInChr + unfixedTailLength " <<midPartMapPosInChr + unfixedTailLength << endl;
		//cout << "(indexInfo->chromLength)[midPartMapChrNameInt]: " << (indexInfo->chromLength)[midPartMapChrNameInt] << endl; 
 		if((midPartMapPosInChr-bufferLength < 0)||
			(midPartMapPosInChr-bufferLength + unfixedTailLength + bufferLength > (indexInfo->chromLength)[midPartMapChrNameInt] ))
		{
			return;
		}	

		string pendingChroSeq = indexInfo->chromStr[midPartMapChrNameInt].substr(
									midPartMapPosInChr-bufferLength, unfixedTailLength + bufferLength);

		//cout << "pendingReadSeq: " << pendingReadSeq << endl;
		//cout << "pendingChroSeq: " << pendingChroSeq << endl;
		const string SJendStrAG = "GT";
		const string SJendStrAC = "CT";	
		//search for "AG"
		int startSearchPos = 0;
		int foundPos = 0;

		int otherPartInReadLength = readLength - unfixedTailLength;
		for(int tmp = 0; ; tmp++)
		{
			foundPos = pendingChroSeq.find(SJendStrAG, startSearchPos);
			if(foundPos == pendingChroSeq.npos)
				break;
			if(foundPos > bufferLength)
			{
				size_t max_append_mismatch = (foundPos - bufferLength)/10 + 1;
				size_t mismatch_bits = 0;
				bool matchBool = score_string(pendingReadSeq.substr(bufferLength, foundPos - bufferLength),
												pendingChroSeq.substr(bufferLength, foundPos - bufferLength),
												max_append_mismatch, mismatch_bits);//append first
				if(matchBool)
					possiGTAGpos.push_back(foundPos + otherPartInReadLength - bufferLength + 1);
			}
			else
			{
				possiGTAGpos.push_back(foundPos + otherPartInReadLength - bufferLength + 1);
			}
			startSearchPos = foundPos + 1;
		}

		//search for "AC"
		startSearchPos = 0;
		foundPos = 0;
		for(int tmp = 0; ; tmp++)
		{
			foundPos = pendingChroSeq.find(SJendStrAC, startSearchPos);
			if(foundPos == pendingChroSeq.npos)
				break;
			if(foundPos > bufferLength)
			{
				size_t max_append_mismatch = (foundPos - bufferLength)/10 + 1;
				size_t mismatch_bits = 0;
				bool matchBool = score_string(pendingReadSeq.substr(bufferLength, foundPos - bufferLength),
												pendingChroSeq.substr(bufferLength, foundPos - bufferLength),
												max_append_mismatch, mismatch_bits);//append first
				if(matchBool)
					possiCTACpos.push_back(foundPos + otherPartInReadLength - bufferLength + 1);
			}
			else
			{
				possiCTACpos.push_back(foundPos + otherPartInReadLength - bufferLength + 1);
			}
			startSearchPos = foundPos + 1;
		}
	}

	bool SJsearchInSJhash(SJhash_Info* SJinfo, const string& readSeqWithDirection, Index_Info* indexInfo)
	//without use of areaHash
	{
		bool SJfoundInSJhash = false;

		int readLength = readSeqWithDirection.length();
		/*
		int mapPosArea = midPartMapPosInChr/(SJinfo.areaSize);

		vector<int> possibleSJposInArea;

		if((SJinfo.SJendPosAreaHash)[midPartMapChrInt].find(mapPosArea) 
			== (SJinfo.SJendPosAreaHash)[midPartMapChrInt].end())
		{
			SJfoundInSJhash = false;
			return SJfoundInSJhash;
		}
		else // area found 
		{
			for(int tmp = 0; tmp < (((SJinfo.SJendPosAreaHash)[midPartMapChrInt].find(mapPosArea))->second).size(); tmp++)
			{
				int tmpPossibleSJposInChr = (((SJinfo.SJendPosAreaHash)[midPartMapChrInt].find(mapPosArea))->second)[tmp];
			}
		}*/
		//vector < pair<int, int> > SJfoundInHash; // <posInRead, SJsize>
		int buffer = 4;
		string readPendingStr = readSeqWithDirection.substr(readLength - unfixedTailLength - buffer - 1, unfixedTailLength + buffer + 1);
		
		for(int tmp = -buffer; tmp < unfixedTailLength; tmp++)
		{
			int tmpSJdonerEndPosInRead = readLength - unfixedTailLength + tmp;
			int tmpSJacceptorStartPosInRead = tmpSJdonerEndPosInRead + 1;
			int tmpSJdonerEndPosInChr = midPartMapPosInChr + tmp;

			if(tmpSJdonerEndPosInChr > (indexInfo->chromLength)[midPartMapChrInt] - 2)
			{
				continue;
			}

			if( ((SJinfo->SJintHashNormal)[midPartMapChrInt]).find(tmpSJdonerEndPosInChr)
				== ((SJinfo->SJintHashNormal)[midPartMapChrInt]).end() )
			{}
			else
			{
				for(set<int>::iterator tmp2 = ((((SJinfo->SJintHashNormal)[midPartMapChrInt]).find(tmpSJdonerEndPosInChr))->second).begin(); 
					tmp2 != ((((SJinfo->SJintHashNormal)[midPartMapChrInt]).find(tmpSJdonerEndPosInChr))->second).end(); 
					tmp2++)
				{
					int tmpSJacceptorStartPosInChr 
						= *tmp2;
					//string readPendingStr = readSeqWithDirection.substr(tmpSJdonerEndPosInRead + buffer - 1, readLength - (tmpSJdonerEndPosInRead + buffer) + 1)
					string chromDonerEndStr;
					string chromAcceptorStartStr;
					string chromPendingStr;
					size_t max_append_mismatch;
					size_t mismatch_bits;
					bool matchBool;

					if((midPartMapPosInChr + tmpSJdonerEndPosInRead - readLength + unfixedTailLength > ((indexInfo->chromLength)[midPartMapChrInt] - 1))
						||
						(tmpSJacceptorStartPosInChr  + readLength - tmpSJacceptorStartPosInRead > ((indexInfo->chromLength)[midPartMapChrInt] - 1)))
					{
						matchBool = false;
					}	
					else
					{
						chromDonerEndStr = indexInfo->chromStr[midPartMapChrInt].substr(
										midPartMapPosInChr - buffer - 1, tmpSJdonerEndPosInRead - readLength + unfixedTailLength + buffer + 1);
				
						
						chromAcceptorStartStr = indexInfo->chromStr[midPartMapChrInt].substr(
										tmpSJacceptorStartPosInChr - 1, readLength - tmpSJacceptorStartPosInRead + 1);
						
						chromPendingStr = chromDonerEndStr + chromAcceptorStartStr;

						
						max_append_mismatch = (unfixedTailLength + buffer + 1)/10 + 1;
						
						mismatch_bits = 0;
						
						matchBool = score_string(readPendingStr, chromPendingStr,
													max_append_mismatch, mismatch_bits);//append first
						//cout << "readPendingStr: "  << endl << readPendingStr << endl; 
						//cout << "chroPendingStr: "  << endl << chromPendingStr << endl;
					}
					
					if(matchBool)
					{
						SJfoundInSJhash = true;
						//cout << "SJ found: at " << tmpSJdonerEndPosInRead << " SJsize: " << tmpSJacceptorStartPosInChr - tmpSJdonerEndPosInChr - 1 << endl;
						SJposFromRemapping.push_back(pair <int, int > (tmpSJdonerEndPosInRead, tmpSJacceptorStartPosInChr - tmpSJdonerEndPosInChr - 1));
					}
				} 
			}

		}


		return SJfoundInSJhash;
	}
};