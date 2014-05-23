#include <string>
#include <string.h>
#include "splice_info.h"

using namespace std;

class Unfixed_Head
{
public:
	string readName;
	string alignDirection; 
	unsigned int midPartMapPosInWholeGenome;
	string otherCigarString;
	int unfixedHeadLength;
	int midPartLength;
	string readSeqOriginal;
	//int readLength;

	int midPartMapPosInChr; // map pos of 1st base in midPart
	string midPartMapChrName;
	int midPartMapChrInt;

	vector<int> possiSJposInRead; // start from 1, end at unfixedHeadLength + bufferLength

	vector<int> possiGTAGpos; // only AG end checked
	vector<int> possiCTACpos; // only AC end checked

	vector<int> possiGTAGpos_mismatch;
	vector<int> possiCTACpos_mismatch;

	vector< pair<int,int> > GTAGsjPos; // sequence around the SJ has been checked including short anchor and midpart
	vector< pair<int,int> > CTACsjPos; // sequence around the SJ has been checked 

	vector< int > GTAGsjPos_mismatch; // sequence around the SJ has been checked including short anchor and midpart
	vector< int > CTACsjPos_mismatch; // sequence around the SJ has been checked 

	vector< pair<int, int> > SJposFromRemappingVec;
	vector< int > SJposFromRemappingVec_mismatch;
	
	vector< pair<int, int> > SJposFromRemappingVec_candi;
	vector< int > SJposFromRemappingVec_candi_mismatch;	


	Unfixed_Head()
	{//readLength = 100;
	}

	Unfixed_Head(int head_length, int mapPos, string mapChrName)
	{
		unfixedHeadLength = head_length;
		midPartMapPosInChr = midPartMapPosInChr;
		midPartMapChrName = mapChrName;
	}

	~Unfixed_Head()
	{}

	int countMismatchNumInBufferSeq(int buffer, const string& readSeqWithDirection, Index_Info* indexInfo)
	{
		int countNum = 0;
		int chromNameInt = indexInfo->convertStringToInt(midPartMapChrName);
		for(int tmp = 0; tmp < buffer + 1; tmp ++)
		{
			if(readSeqWithDirection.at(unfixedHeadLength + 1 + tmp - 1) != indexInfo->chromStr[chromNameInt].at(midPartMapPosInChr + tmp - 1) )
			{
				countNum ++;
			}
			else
			{}
		}
		return countNum;
	}

	void getUnfixedHeadInfoFromRecord(PE_Read_Info* readInfo, bool end1, Alignment_Info* alignInfo, Index_Info* indexInfo)
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
		midPartMapPosInChr = alignInfo->alignChromPos;

		midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);
		//cout << "midPartMapChrInt: " << midPartMapChrInt << endl;

		otherCigarString = alignInfo->otherJumpCodeVec2Str();
		unfixedHeadLength = (alignInfo->cigarStringJumpCode)[0].len;
		midPartLength = (alignInfo->cigarStringJumpCode)[1].len;

		//cout << "start to get whole genome loc" << endl;
		//cout << "midPartMapPosInChr: " << midPartMapPosInChr << endl;
		//cout << "midPartMapChrInt: " << midPartMapChrInt << endl;
		midPartMapPosInWholeGenome = indexInfo->getWholeGenomeLocation(
			(unsigned int)midPartMapChrInt, (unsigned int)midPartMapPosInChr);
		//cout << "midPartMapPosInWholeGenome: " << midPartMapPosInWholeGenome << endl;
	}

	void getUnfixedHeadInfoFromRecordWithAlignInfoType(PE_Read_Info* readInfo, int alignInfoType, Alignment_Info* alignInfo, Index_Info* indexInfo)
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
		midPartMapPosInChr = alignInfo->alignChromPos;

		midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);
		//cout << "midPartMapChrInt: " << midPartMapChrInt << endl;

		otherCigarString = alignInfo->otherJumpCodeVec2Str();
		unfixedHeadLength = (alignInfo->cigarStringJumpCode)[0].len;
		midPartLength = (alignInfo->cigarStringJumpCode)[1].len;

		//cout << "start to get whole genome loc" << endl;
		//cout << "midPartMapPosInChr: " << midPartMapPosInChr << endl;
		//cout << "midPartMapChrInt: " << midPartMapChrInt << endl;
		midPartMapPosInWholeGenome = indexInfo->getWholeGenomeLocation(
			(unsigned int)midPartMapChrInt, (unsigned int)midPartMapPosInChr);
		//cout << "midPartMapPosInWholeGenome: " << midPartMapPosInWholeGenome << endl;
	}
	void fromRecordStringToClass(char* unfixedHeadRecordString, int readLength)
	{
		//cout << "unfixedHeadRecordString: " << endl;
		//cout << unfixedHeadRecordString << endl;

		char readNameCharForUnfixedHead[100]; 
		char alignDirectionCharForUnfixedHead;
		unsigned int midPartMapPosForUnfixedHead;
		int firstJumpCodeLengthForUnfixedHead; // type = S
		int secondJumpCodeLengthForUnfixedHead;
		char otherCigarStringCharForUnfixedHead[50]; 
		char readSeqCharForUnfixedHead[600];
		
		sscanf(unfixedHeadRecordString, "%s\t%s\t%d\t%s\t%d\t%d\t%s",
			readNameCharForUnfixedHead, &alignDirectionCharForUnfixedHead,
			&midPartMapPosForUnfixedHead, otherCigarStringCharForUnfixedHead,
			&firstJumpCodeLengthForUnfixedHead, &secondJumpCodeLengthForUnfixedHead,
			readSeqCharForUnfixedHead);

		readName = readNameCharForUnfixedHead;
		alignDirection = alignDirectionCharForUnfixedHead;
		midPartMapPosInWholeGenome = midPartMapPosForUnfixedHead;// + firstJumpCodeLengthForUnfixedHead;
		otherCigarString = otherCigarStringCharForUnfixedHead;
		unfixedHeadLength = firstJumpCodeLengthForUnfixedHead;
		midPartLength = secondJumpCodeLengthForUnfixedHead;
		readSeqOriginal = readSeqCharForUnfixedHead;
		//int readLength = 100;
		readSeqOriginal = readSeqOriginal.substr(0, readLength);
		//cout << "unfixedHeadInfo: " << endl;
		//cout << readName << endl << alignDirection << endl << midPartMapPosInWholeGenome << endl << unfixedHeadLength << endl << readSeqOriginal << endl;
	}

	void getPossibleSJpos(const string& readSeqWithDirection,/* const string& chromSeq,*/ Index_Info* indexInfo)
	{
		int bufferLength = 4;
		if(bufferLength > midPartLength)
		{
			bufferLength = midPartLength;
		}

		string pendingReadSeq = readSeqWithDirection.substr(0, unfixedHeadLength+bufferLength);
		//string pendingChroSeq = chromSeq.substr(midPartMapPosInWholeGenome-unfixedHeadLength-1, unfixedHeadLength+bufferLength);		
		int pendingChroSeqStartPos = midPartMapPosInChr - unfixedHeadLength-1;
		if((pendingChroSeqStartPos < 0)||(pendingChroSeqStartPos + unfixedHeadLength+bufferLength > (indexInfo->chromLength[midPartMapChrInt])))
		{
			/*cout << endl << "readName: " << readName << endl;
			cout << "pendingChroSeqStartPos: " << pendingChroSeqStartPos << endl;
			cout << "midPartMapChrInt: " << midPartMapChrInt << endl;
			cout << "midPartMapChrName: " << midPartMapChrName << endl;
			cout << "chromSize: " << indexInfo->chromStr[midPartMapChrInt].length() << endl;
			cout << "realChromSize: " << (indexInfo->chrEndPosInGenome)[midPartMapChrInt]-(indexInfo->chrEndPosInGenome)[midPartMapChrInt-1]-1 << endl;
			cout << endl;*/
			return;
		}	
		string pendingChroSeq = indexInfo->chromStr[indexInfo->convertStringToInt(midPartMapChrName)].substr(
									midPartMapPosInChr - unfixedHeadLength-1,
									unfixedHeadLength+bufferLength);
		//cout << "pendingReadSeq: " << pendingReadSeq << endl;
		//cout << "pendingChroSeq: " << pendingChroSeq << endl;

		string SJendStrAG = "AG";
		string SJendStrAC = "AC";

		//cout << "searching for AGs" << endl;
		//search for "AG" 
		int startSearchPos = 0;
		int foundPos = 0;
		for(int tmp = 0; ; tmp++)
		{
			foundPos = pendingChroSeq.find(SJendStrAG, startSearchPos);
			if(foundPos == pendingChroSeq.npos)
				break;
			if(unfixedHeadLength-foundPos-2 > 0)
			{
				size_t max_append_mismatch = (unfixedHeadLength-foundPos-2)/10 + 1;
				size_t mismatch_bits = 0;
				size_t comb_bits = 0;
				bool matchBool = score_string(pendingReadSeq.substr(foundPos+2, unfixedHeadLength-foundPos-2), 
												pendingChroSeq.substr(foundPos+2, unfixedHeadLength-foundPos-2), 
												max_append_mismatch, mismatch_bits, comb_bits);//append first
				if(matchBool)
				{
					possiGTAGpos.push_back(foundPos+3);
					possiGTAGpos_mismatch.push_back(mismatch_bits);
				}
			}
			else
			{
				possiGTAGpos.push_back(foundPos+3);	
				possiGTAGpos_mismatch.push_back(0);
			}

			startSearchPos = foundPos+1;
			//cout << "found at foundPos " << endl;
		}

		//cout << "searching for ACs" << endl;
		//search for "AC" 
		startSearchPos = 0;
		foundPos = 0;
		for(int tmp = 0; ; tmp++)
		{
			foundPos = pendingChroSeq.find(SJendStrAC, startSearchPos);
			if(foundPos == pendingChroSeq.npos)
				break;			
			if(unfixedHeadLength-foundPos-2 > 0)
			{
				size_t max_append_mismatch = (unfixedHeadLength-foundPos-2)/10;
				size_t mismatch_bits = 0;
				size_t comb_bits = 0;
				bool matchBool = score_string(pendingReadSeq.substr(foundPos+2, unfixedHeadLength-foundPos-2), 
												pendingChroSeq.substr(foundPos+2, unfixedHeadLength-foundPos-2), 
												max_append_mismatch, mismatch_bits, comb_bits);//append first
				if(matchBool)
				{
					possiCTACpos.push_back(foundPos+3);
					possiCTACpos_mismatch.push_back(mismatch_bits);
				}
			}
			else
			{
				possiCTACpos.push_back(foundPos+3);	
				possiCTACpos_mismatch.push_back(0);
			}

			startSearchPos = foundPos+1;
			//cout << "found at foundPos " << endl;
		}
		return;
	}

	/*void getPossibleSJpos_test(const string& readSeqWithDirection, Index_Info* indexInfo)
	{
		int bufferLength = 4;
		if(bufferLength > midPartLength)
		{
			bufferLength = midPartLength;
		}
		//cout << "size of read: " << readSeqWithDirection.length()<< endl;
		//cout << "size of genome: " << chromSeq.length()<< endl;
		string pendingReadSeq = readSeqWithDirection.substr(0, unfixedHeadLength+bufferLength);
		//string pendingChroSeq = chromSeq.substr(midPartMapPosInWholeGenome-unfixedHeadLength-1, unfixedHeadLength+bufferLength);
		
		int pendingChroSeqStartPos = midPartMapPosInChr - unfixedHeadLength-1;
		if((pendingChroSeqStartPos < 0)||(pendingChroSeqStartPos + unfixedHeadLength+bufferLength > indexInfo->chromStr[midPartMapChrInt].length()))
		{
			//cout << "readName: " << readName << endl;
			//cout << "pendingChroSeqStartPos: " << pendingChroSeqStartPos << endl;
			//cout << "midPartMapChrInt: " << midPartMapChrInt << endl;
			//cout << "midPartMapChrName: " << midPartMapChrName << " int: " << indexInfo->convertStringToInt(midPartMapChrName) << endl;
			//cout << "chromSize: " << indexInfo->chromStr[midPartMapChrInt].length() << endl;
			//cout << "realChromSize: " << (indexInfo->chrEndPosInGenome)[midPartMapChrInt]-(indexInfo->chrEndPosInGenome)[midPartMapChrInt-1]-1 << endl << endl;
			return;
		}	
		string pendingChroSeq = indexInfo->chromStr[midPartMapChrInt].substr(
									pendingChroSeqStartPos,
									unfixedHeadLength+bufferLength);
		//cout << "pendingReadSeq: " << pendingReadSeq << endl;
		//cout << "pendingChroSeq: " << pendingChroSeq << endl;

		string SJendStrAG = "AG";
		string SJendStrAC = "AC";

		//cout << "searching for AGs" << endl;
		//search for "AG" 
		int startSearchPos = 0;
		int foundPos = 0;
		for(int tmp = 0; ; tmp++)
		{
			foundPos = pendingChroSeq.find(SJendStrAG, startSearchPos);
			if(foundPos == pendingChroSeq.npos)
				break;
			if(unfixedHeadLength-foundPos-2 > 0)
			{
				size_t max_append_mismatch = (unfixedHeadLength-foundPos-2)/10 + 1;
				size_t mismatch_bits = 0;
				bool matchBool = score_string(pendingReadSeq.substr(foundPos+2, unfixedHeadLength-foundPos-2), 
												pendingChroSeq.substr(foundPos+2, unfixedHeadLength-foundPos-2), 
												max_append_mismatch, mismatch_bits);//append first
				if(matchBool)
					possiGTAGpos.push_back(foundPos+3);
			}
			else
				possiGTAGpos.push_back(foundPos+3);	

			startSearchPos = foundPos+1;
			//cout << "found at foundPos " << endl;
		}

		//cout << "searching for ACs" << endl;
		//search for "AC" 
		startSearchPos = 0;
		foundPos = 0;
		for(int tmp = 0; ; tmp++)
		{
			foundPos = pendingChroSeq.find(SJendStrAC, startSearchPos);
			if(foundPos == pendingChroSeq.npos)
				break;			
			if(unfixedHeadLength-foundPos-2 > 0)
			{
				size_t max_append_mismatch = (unfixedHeadLength-foundPos-2)/10;
				size_t mismatch_bits = 0;
				bool matchBool = score_string(pendingReadSeq.substr(foundPos+2, unfixedHeadLength-foundPos-2), 
												pendingChroSeq.substr(foundPos+2, unfixedHeadLength-foundPos-2), 
												max_append_mismatch, mismatch_bits);//append first
				if(matchBool)
					possiCTACpos.push_back(foundPos+3);
			}
			else
				possiCTACpos.push_back(foundPos+3);	

			startSearchPos = foundPos+1;
			//cout << "found at foundPos " << endl;
		}
		return;
	}*/


	void printSJ()
	{
		cout << "...... GTAG sj ......" << endl;
		for(int tmp = 0; tmp < GTAGsjPos.size(); tmp++)
		{
			cout << "SJ:" << tmp << " headLength: " << GTAGsjPos[tmp].first - 1 
			<< " SJdistance: " <<  GTAGsjPos[tmp].second << endl;
		}
		cout << "...... CTAC sj ......" << endl;
		for(int tmp = 0; tmp < CTACsjPos.size(); tmp++)
		{
			cout << "SJ:" << tmp << " headLength: " << CTACsjPos[tmp].first - 1 
			<< " SJdistance: " <<  CTACsjPos[tmp].second << endl;
		}
	}

	void printUnfixedHeadInfo()
	{
		cout << readName << endl << alignDirection << endl << midPartMapPosInWholeGenome << endl << unfixedHeadLength << readSeqOriginal << endl; 
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

		//if(unfixedHeadLength + buffer >= readLength)
		//	buffer = readLength - unfixedHeadLength - 1;

		//string readPendingStr
		//set<int> areaNOset;
		int chromNameInt = indexInfo->convertStringToInt(midPartMapChrName);

		int areaNOmin = (int)((midPartMapPosInChr - unfixedHeadLength + 1)/areaSize);
		int areaNOmax = (int)((midPartMapPosInChr + buffer)/areaSize);

		//int areaCandidateNum = areaNOmax - areaNOmin + 1;

		vector<int> SJacceptorSiteVec;

		for(int tmpArea = areaNOmin; tmpArea <= areaNOmax; tmpArea ++)
		{
			SJareaHashIter tmpSJareaHashIter 
				= ((SJinfo->SJendPosAreaHash)[chromNameInt]).find(tmpArea);
			if(tmpSJareaHashIter != ((SJinfo->SJendPosAreaHash)[chromNameInt]).end())
			{
				for(set<int>::iterator intSetIter = (tmpSJareaHashIter->second).begin();
					intSetIter != (tmpSJareaHashIter->second).end(); intSetIter ++)
				{
					int tmpSJacceptorPos = (*intSetIter);
					if( (tmpSJacceptorPos >= midPartMapPosInChr - unfixedHeadLength + 1) 
						&& (tmpSJacceptorPos <= midPartMapPosInChr + buffer) )
					{
						SJacceptorSiteVec.push_back(tmpSJacceptorPos);
					}
				}
			}
			else
			{}
		}
		
		if(SJacceptorSiteVec.size() > 0)
		{
			string readPendingStr = readSeqWithDirection.substr(
				0, unfixedHeadLength + buffer + 1);
			for(int tmp = 0; tmp < SJacceptorSiteVec.size(); tmp++)
			{
				vector<int> tmpDonerSiteVec;

				int tmpSJacceptorSite = SJacceptorSiteVec[tmp];
				int tmpHeadLength = unfixedHeadLength + tmpSJacceptorSite - midPartMapPosInChr;
				
				SplicePosHashIter tmpPosHashIter 
					= (SJinfo->spliceJunctionReverse)[chromNameInt].find(tmpSJacceptorSite);
				
				if(tmpPosHashIter != ((SJinfo->spliceJunctionReverse)[chromNameInt]).end())
				{
					if( tmpHeadLength >= (SJinfo->anchorStringLength) )
					{
						string tmpAnchorString = readSeqWithDirection.substr(tmpHeadLength - (SJinfo->anchorStringLength) , (SJinfo->anchorStringLength));

						SpliceEndStrHashIter tmpEndStrHashIter = (tmpPosHashIter->second).find(tmpAnchorString);
						if(tmpEndStrHashIter != (tmpPosHashIter->second).end())
						{
							for(set<int>::iterator tmpIntSetIter = (tmpEndStrHashIter->second).begin(); 
								tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter ++)
							{
								tmpDonerSiteVec.push_back(*tmpIntSetIter);
							}
						}
						else
						{
							//cout << "error! * should be found in endStrHash" << endl;
							tmpEndStrHashIter = (tmpPosHashIter->second).find("*");
							if(tmpEndStrHashIter != (tmpPosHashIter->second).end())
							{
								for(set<int>::iterator tmpIntSetIter = (tmpEndStrHashIter->second).begin(); 
									tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter ++)
								{
									tmpDonerSiteVec.push_back(*tmpIntSetIter);
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
								tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter ++)
							{
								tmpDonerSiteVec.push_back(*tmpIntSetIter);
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
					cout << "error in SJsearchInAreaAndStringHash ! tmpSJacceptorSite should be found in hash ! " << endl;
				}

				for(int tmpVecNO = 0; tmpVecNO < tmpDonerSiteVec.size(); tmpVecNO ++)
				{
					//int tmpSJdonerSite = tmpDonerSiteVec[tmpVec];
					//int tmpDonerLocInRead = tmpHeadLength;
					//int tmpAcceptorLocInRead = tmpDonerLocInRead + 1;
					int tmpSJdonerEndPosInRead = tmpHeadLength;
					int tmpSJacceptorStartPosInRead = tmpSJdonerEndPosInRead + 1;

					int tmpSJdonerEndPosInChr = tmpDonerSiteVec[tmpVecNO];
					int tmpSJacceptorStartPosInChr = tmpSJacceptorSite;

					string chromDonerEndStr;
					string chromAcceptorStartStr;
					string chromPendingStr;
					size_t max_append_mismatch;
					size_t mismatch_bits;
					size_t comb_bits;
					bool matchBool;

					if(tmpSJdonerEndPosInChr - tmpSJdonerEndPosInRead >= 0)
					{
						chromDonerEndStr = indexInfo->chromStr[chromNameInt].substr(
										tmpSJdonerEndPosInChr - tmpSJdonerEndPosInRead, tmpSJdonerEndPosInRead);

						chromAcceptorStartStr = indexInfo->chromStr[chromNameInt].substr(
										tmpSJacceptorStartPosInChr - 1, unfixedHeadLength + 1 - tmpSJacceptorStartPosInRead + buffer + 1);

						chromPendingStr = chromDonerEndStr + chromAcceptorStartStr;
						
						max_append_mismatch = (unfixedHeadLength)/10 + 1;
						
						mismatch_bits = 0;
						
						comb_bits = 0;

						matchBool = score_string(readPendingStr, chromPendingStr,
													max_append_mismatch, mismatch_bits, comb_bits);//append first
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

	void filterSJposFromRemapping_candi(int mismatchNumInBufferSeq)
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

		//int countNum = this->countMismatchNumInBufferSeq()

		SJposFromRemappingVec_mismatch.push_back(SJposFromRemapping_mismatch - mismatchNumInBufferSeq);
	}

	bool SJsearchInSJhash(SJhash_Info* SJinfo, const string& readSeqWithDirection, 
		Index_Info* indexInfo)
	//without use of areaHash
	{
		bool SJfoundInSJhash = false;

		int buffer = 4;
		string readPendingStr = readSeqWithDirection.substr(0, unfixedHeadLength + buffer + 1);
		
		for(int tmp = 1; tmp < unfixedHeadLength + buffer + 1; tmp++)
		{
			int tmpSJdonerEndPosInRead = tmp;
			int tmpSJacceptorStartPosInRead = tmpSJdonerEndPosInRead + 1;
			int tmpSJacceptorStartPosInChr = tmpSJacceptorStartPosInRead + midPartMapPosInChr - unfixedHeadLength - 1;

			//cout << "chr: " << midPartMapChrInt << endl
			//	<< "donerEndInRead: " << tmpSJdonerEndPosInRead << endl 
			//	<< "acceptorStartPosInChr: " << tmpSJacceptorStartPosInChr << endl; 
			if(tmpSJacceptorStartPosInChr < 2)
			{
				continue;
			}

			if( ((SJinfo->SJintHashReverse)[midPartMapChrInt]).find(tmpSJacceptorStartPosInChr)
				== ((SJinfo->SJintHashReverse)[midPartMapChrInt]).end() )
			{}
			else
			{
			//	cout << "found in SJintHashReverse" << endl;
				for(set<int>::iterator tmp2 = ((((SJinfo->SJintHashReverse)[midPartMapChrInt]).find(tmpSJacceptorStartPosInChr))->second).begin(); 
					tmp2 != ((((SJinfo->SJintHashReverse)[midPartMapChrInt]).find(tmpSJacceptorStartPosInChr))->second).end(); 
					tmp2++)
				{
					int tmpSJdonerEndPosInChr 
						= *tmp2;
					//string readPendingStr = readSeqWithDirection.substr(tmpSJdonerEndPosInRead + buffer - 1, readLength - (tmpSJdonerEndPosInRead + buffer) + 1)
					string chromDonerEndStr;
					string chromAcceptorStartStr;
					string chromPendingStr;
					size_t max_append_mismatch;
					size_t mismatch_bits;
					bool matchBool;

					if(tmpSJdonerEndPosInChr - tmpSJdonerEndPosInRead >= 0)
					{
						chromDonerEndStr = indexInfo->chromStr[midPartMapChrInt].substr(
										tmpSJdonerEndPosInChr - tmpSJdonerEndPosInRead, tmpSJdonerEndPosInRead);

						chromAcceptorStartStr = indexInfo->chromStr[midPartMapChrInt].substr(
										tmpSJacceptorStartPosInChr - 1, unfixedHeadLength + 1 - tmpSJacceptorStartPosInRead + buffer + 1);

						chromPendingStr = chromDonerEndStr + chromAcceptorStartStr;

						
						max_append_mismatch = (unfixedHeadLength + buffer + 1)/10 + 1;
						
						mismatch_bits = 0;
						
						matchBool = score_string(readPendingStr, chromPendingStr,
													max_append_mismatch, mismatch_bits);//append first
					}
					else
					{
						matchBool = false;
					}
					//cout << "readPendingStr: "  << endl << readPendingStr << endl; 
					//cout << "chroPendingStr: "  << endl << chromPendingStr << endl;

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