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

	vector< pair<int,int> > GTAGsjPos; // sequence around the SJ has been checked including short anchor and midpart
	vector< pair<int,int> > CTACsjPos; // sequence around the SJ has been checked 

	vector< pair<int,int> > SJposFromRemapping;
	
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
	}

	void getPossibleSJpos_test(const string& readSeqWithDirection,/* const string& chromSeq,*/ Index_Info* indexInfo)
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
	}


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

	bool SJsearchInAreaAndStringHash(SJhash_Info* SJinfo, 
		const string& readSeqWithDirection, 
		Index_Info* indexInfo, int areaSize)
	{
		bool SJfoundInSJhash = false;

		int readLength = readSeqWithDirection.length();
		int buffer = 4;
		if(unfixedHeadLength + buffer >= readLength)
			buffer = readLength - unfixedHeadLength - 1;

		//string readPendingStr
		//set<int> areaNOset;
		int chromNameInt = indexInfo->convertStringToInt(midPartMapChrName);

		int areaNOmin = (int)((midPartMapPosInChr - unfixedHeadLength - 1)/areaSize);
		int areaNOmax = (int)((midPartMapPosInChr + buffer)/areaSize);

		int areaCandidateNum = areaNOmax - areaNOmin + 1;

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
					if( (tmpSJacceptorPos >= midPartMapPosInChr - unfixedHeadLength - 1) 
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
				int tmpSJacceptorSite = SJacceptorSiteVec[tmp];

			}
		}
		else
		{}	
		
		return SJfoundInSJhash;
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
						SJposFromRemapping.push_back(pair <int, int > (tmpSJdonerEndPosInRead, tmpSJacceptorStartPosInChr - tmpSJdonerEndPosInChr - 1));
					}
				} 
			}

		}
		return SJfoundInSJhash;
	}

};