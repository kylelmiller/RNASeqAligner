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

	vector<int> possiGTAGpos_mismatch; // only AG end checked
	vector<int> possiCTACpos_mismatch; // only AC end checked

	vector< pair<int,int> > GTAGsjPos; // sequence around the SJ has been checked including short anchor and midpart
	vector< pair<int,int> > CTACsjPos; // sequence around the SJ has been checked 

	vector< int > GTAGsjPos_mismatch; // sequence around the SJ has been checked including short anchor and midpart
	vector< int > CTACsjPos_mismatch; // sequence around the SJ has been checked 

	vector< pair<int,int> > SJposFromRemappingVec; // // <posInRead, SJsize>
	vector< int > SJposFromRemappingVec_mismatch;

	vector< pair<int,int> > SJposFromRemappingVec_candi;
	vector< int > SJposFromRemappingVec_candi_mismatch;   


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

	int countMismatchNumInBufferSeq(int buffer, const string& readSeqWithDirection, Index_Info* indexInfo)
	{
		int countNum = 0;
		int chromNameInt = indexInfo->convertStringToInt(midPartMapChrName);
		int readLength = readSeqWithDirection.length();
		for(int tmp = 0; tmp < buffer + 1; tmp ++)
		{
			if(readSeqWithDirection.at(readLength - unfixedTailLength - buffer - 1 + tmp ) 
				!= (indexInfo->chromStr[chromNameInt]).at(midPartMapPosInChr - buffer - 1 + tmp) )
			{
				countNum ++;
			}
			else
			{}
		}
		return countNum;
	}

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
				size_t comb_bits = 0;
				bool matchBool = score_string(pendingReadSeq.substr(bufferLength, foundPos - bufferLength),
												pendingChroSeq.substr(bufferLength, foundPos - bufferLength),
												max_append_mismatch, mismatch_bits, comb_bits);//append first
				if(matchBool)
				{
					possiGTAGpos.push_back(foundPos + otherPartInReadLength - bufferLength + 1);
					possiGTAGpos_mismatch.push_back(mismatch_bits);
				}
			}
			else
			{
				possiGTAGpos.push_back(foundPos + otherPartInReadLength - bufferLength + 1);
				possiGTAGpos_mismatch.push_back(0);
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
				size_t comb_bits = 0;
				bool matchBool = score_string(pendingReadSeq.substr(bufferLength, foundPos - bufferLength),
												pendingChroSeq.substr(bufferLength, foundPos - bufferLength),
												max_append_mismatch, mismatch_bits, comb_bits);//append first
				if(matchBool)
				{
					possiCTACpos.push_back(foundPos + otherPartInReadLength - bufferLength + 1);
					possiCTACpos_mismatch.push_back(mismatch_bits);
				}
			}
			else
			{
				possiCTACpos.push_back(foundPos + otherPartInReadLength - bufferLength + 1);
				possiCTACpos_mismatch.push_back(0);
			}
			startSearchPos = foundPos + 1;
		}
	}

	bool SJsearchInSJhash_areaStringHash(SJhash_Info* SJinfo, 
		const string& readSeqWithDirection, 
		Index_Info* indexInfo, int areaSize)
	{
		bool SJfoundInSJhash = false;

		int readLength = readSeqWithDirection.length();

		int buffer = 4;

		if(buffer > midPartLength - 1)
			buffer = midPartLength - 1;

		/*if(readLength - unfixedTailLength - buffer - 1 < 0)
		{
			buffer = readLength - unfixedTailLength - 1;
		}*/

		string readPendingStr = readSeqWithDirection.substr(readLength - unfixedTailLength - buffer - 1, unfixedTailLength + buffer + 1);
		int chromNameInt = indexInfo->convertStringToInt(midPartMapChrName);

		int areaNOmin = (int)((midPartMapPosInChr-buffer)/areaSize);
		int areaNOmax = (int)((midPartMapPosInChr+unfixedTailLength-1)/areaSize);

		vector<int> SJdonerSiteVec;

		for(int tmpArea = areaNOmin; tmpArea <= areaNOmax; tmpArea ++)
		{
			SJareaHashIter tmpSJareaHashIter
				= ((SJinfo->SJstartPosAreaHash)[chromNameInt]).find(tmpArea);
			if(tmpSJareaHashIter != ((SJinfo->SJstartPosAreaHash)[chromNameInt]).end())
			{
				for(set<int>::iterator intSetIter = (tmpSJareaHashIter->second).begin();
					intSetIter != (tmpSJareaHashIter->second).end(); intSetIter ++)
				{
					int tmpSJdonerPos = (*intSetIter);
					if( (tmpSJdonerPos >= (midPartMapPosInChr-buffer)) && (tmpSJdonerPos <= (midPartMapPosInChr+unfixedTailLength-1)) )
					{
						SJdonerSiteVec.push_back(tmpSJdonerPos);
					}
				}
			}
			else
			{}
		}

		if(SJdonerSiteVec.size() > 0)
		{
			string readPendingStr = readSeqWithDirection.substr(readLength - unfixedTailLength - buffer - 1, unfixedTailLength + buffer + 1);
			for(int tmp = 0; tmp < SJdonerSiteVec.size(); tmp ++)
			{
				vector<int> tmpAcceptorSiteVec;

				int tmpSJdonerSite = SJdonerSiteVec[tmp];
				int tmpTailLength = unfixedTailLength + midPartMapPosInChr - tmpSJdonerSite;

				SplicePosHashIter tmpPosHashIter
					= (SJinfo->spliceJunctionNormal)[chromNameInt].find(tmpSJdonerSite);

				if(tmpPosHashIter != (SJinfo->spliceJunctionNormal)[chromNameInt].end())
				{
					if(tmpTailLength >= (SJinfo->anchorStringLength))
					{
						string tmpAnchorString = readSeqWithDirection.substr(readLength - tmpTailLength, (SJinfo->anchorStringLength));
					
						SpliceEndStrHashIter tmpEndStrHashIter = (tmpPosHashIter->second).find(tmpAnchorString);
						if(tmpEndStrHashIter != (tmpPosHashIter->second).end())
						{
							for(set<int>::iterator tmpIntSetIter = (tmpEndStrHashIter->second).begin();
								tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter++)
							{
								tmpAcceptorSiteVec.push_back(*tmpIntSetIter);
							}
						}
						else
						{
							tmpEndStrHashIter = (tmpPosHashIter->second).find("*");
							if(tmpEndStrHashIter != (tmpPosHashIter->second).end())
							{
								for(set<int>::iterator tmpIntSetIter = (tmpEndStrHashIter->second).begin();
									tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter++)
								{
									tmpAcceptorSiteVec.push_back(*tmpIntSetIter);
								}
							}
							else
							{
								cout << "error! * should be found in endStrHash" << endl;
							}

						}
					}
					else
					{
						SpliceEndStrHashIter tmpEndStrHashIter = (tmpPosHashIter->second).find("*");
						if(tmpEndStrHashIter != (tmpPosHashIter->second).end())
						{
							for(set<int>::iterator tmpIntSetIter = (tmpEndStrHashIter->second).begin();
								tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter ++ )
							{
								tmpAcceptorSiteVec.push_back(*tmpIntSetIter);
							}
						}
						else
						{
							cout << "error! * should be found in endStrHash" << endl;
						}
					}
				}
				else
				{
					cout << "error in SJsearchInAreaAndStringHash ! tmpSJdonerSite should be found in hash ! " << endl;

				}

				for(int tmpVecNO = 0; tmpVecNO < tmpAcceptorSiteVec.size(); tmpVecNO ++)
				{
					int tmpSJdonerEndPosInRead = readLength - tmpTailLength;
					int tmpSJacceptorStartPosInRead = tmpSJdonerEndPosInRead + 1;

					int tmpSJdonerEndPosInChr = tmpSJdonerSite;
					int tmpSJacceptorStartPosInChr = tmpAcceptorSiteVec[tmpVecNO];

					string chromDonerEndStr;
					string chromAcceptorStartStr;
					string chromPendingStr;
					size_t max_append_mismatch;
					size_t mismatch_bits;
					size_t comb_bits;
					bool matchBool;

					if(
						(midPartMapPosInChr + tmpSJdonerEndPosInRead - readLength + unfixedTailLength 
							<= ((indexInfo->chromLength)[midPartMapChrInt] - 1))
						&&
						(tmpSJacceptorStartPosInChr  + readLength - tmpSJacceptorStartPosInRead 
							<= ((indexInfo->chromLength)[midPartMapChrInt] - 1))
						)
					{
						chromDonerEndStr = (indexInfo->chromStr[chromNameInt]).substr(midPartMapPosInChr - buffer - 1,
							tmpSJdonerEndPosInRead - readLength + unfixedTailLength + buffer + 1);
						chromAcceptorStartStr = (indexInfo->chromStr[chromNameInt]).substr(
							tmpSJacceptorStartPosInChr - 1, readLength - tmpSJacceptorStartPosInRead + 1);
						chromPendingStr = chromDonerEndStr + chromAcceptorStartStr;
			
						max_append_mismatch = (unfixedTailLength)/10 + 1;			
						mismatch_bits = 0;
						comb_bits = 0;
						matchBool = score_string(readPendingStr, chromPendingStr,
														max_append_mismatch, mismatch_bits, comb_bits);//append first
							//cout << "readPendingStr: "  << endl << readPendingStr << endl; 
							//cout << "chroPendingStr: "  << endl << chromPendingStr << endl;
					}
					else
					{
						matchBool = false;
					}

					if(matchBool)
					{
						SJfoundInSJhash = true;
						//cout << "SJ found: at " << tmpSJdonerEndPosInRead << " SJsize: " << tmpSJacceptorStartPosInChr - tmpSJdonerEndPosInChr - 1 << endl;
						SJposFromRemappingVec_candi.push_back(pair <int, int > (tmpSJdonerEndPosInRead, tmpSJacceptorStartPosInChr - tmpSJdonerEndPosInChr - 1));
						SJposFromRemappingVec_candi_mismatch.push_back(mismatch_bits);
					}					

				}
			}
		}
		else
		{}

		if(SJfoundInSJhash)
		{
			//int mismatchNumInBufferSeq = this->countMismatchNumInBufferSeq(buffer, readSeqWithDirection, indexInfo);
			int mismatchNumInBufferSeq = 0;
			this->filterSJposFromRemapping_candi(mismatchNumInBufferSeq);
		}
		
		return SJfoundInSJhash;
	}

	void filterSJposFromRemapping_candi(int mismatchNumInBufferSeq)//, Index_Info* indexInfo)
	{
		int currentBestSJposNO = 0;
		int currentBestSJposNO_SJdistance = (SJposFromRemappingVec_candi[0]).second;
		int currentBestSJposNO_mismatch = SJposFromRemappingVec_candi_mismatch[0];
		for(int tmpSJposNO = 1; tmpSJposNO < SJposFromRemappingVec_candi.size(); tmpSJposNO ++)
		{
			int tmpSJpos_mismatch = SJposFromRemappingVec_candi_mismatch[tmpSJposNO]; 
			int tmpSJpos_SJdistance = (SJposFromRemappingVec_candi[tmpSJposNO]).second;

			if(tmpSJpos_mismatch < currentBestSJposNO_mismatch)
			{
				currentBestSJposNO = tmpSJposNO;
				currentBestSJposNO_mismatch = tmpSJpos_mismatch;
				currentBestSJposNO_SJdistance = tmpSJpos_SJdistance;
			}
			else if(tmpSJpos_mismatch == currentBestSJposNO_mismatch)
			{
				if(tmpSJpos_SJdistance < currentBestSJposNO_SJdistance)
				{
					currentBestSJposNO = tmpSJposNO;
					currentBestSJposNO_mismatch = tmpSJpos_mismatch;
					currentBestSJposNO_SJdistance = tmpSJpos_SJdistance;
				}
				else
				{}
			}
			else
			{}
		}

		int SJposFromRemapping_first = (SJposFromRemappingVec_candi[currentBestSJposNO]).first;
		int SJposFromRemapping_second = (SJposFromRemappingVec_candi[currentBestSJposNO]).second;
		int SJposFromRemapping_mismatch = SJposFromRemappingVec_candi_mismatch[currentBestSJposNO];

		SJposFromRemappingVec.push_back(pair<int, int> (SJposFromRemapping_first, SJposFromRemapping_second) );
		SJposFromRemappingVec_mismatch.push_back(SJposFromRemapping_mismatch - mismatchNumInBufferSeq);
	}

	bool SJsearchInSJhash(SJhash_Info* SJinfo, const string& readSeqWithDirection, Index_Info* indexInfo)
	//without use of areaHash
	{
		bool SJfoundInSJhash = false;

		int readLength = readSeqWithDirection.length();

		int buffer = 4;

		if(readLength - unfixedTailLength - buffer - 1 < 0)
		{
			buffer = readLength - unfixedTailLength - 1;
		}

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
						SJposFromRemappingVec.push_back(pair <int, int > (tmpSJdonerEndPosInRead, tmpSJacceptorStartPosInChr - tmpSJdonerEndPosInChr - 1));
					}
				} 
			}

		}


		return SJfoundInSJhash;
	}
};