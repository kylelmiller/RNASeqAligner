#ifndef __ALIGN_INFO_H_INCLUDED__
#define __ALIGN_INFO_H_INCLUDED__

#include "path.h"
#include "pairedEndRead.h"
#include "splice_info.h"
#include "mappedRead.h"

// FIXME - KLM 5/30/14 THIS NEEDS TO GET REWORKED
class SpliceJunction_Alignment
{
public:

	string SJchrName;
	int SJdonerEnd;
	int SJacceptorStart;
	string flankString;

	SpliceJunction_Alignment(string chrName, int donerPos, int acceptorPos, Index_Info* indexInfo)
	{
		//cout << "chrName: " << chrName.length() << " " << chrName << endl;
		SJchrName = chrName;
		SJdonerEnd = donerPos;
		SJacceptorStart = acceptorPos;
		//cout << SJdonerEnd << " " << SJacceptorStart << endl;
		this->getFlankString(indexInfo);
	}

	void getFlankString(Index_Info* indexInfo)
	{
		int chromNO = indexInfo->convertStringToInt(SJchrName);
		//cout << "chromNO: " << chromNO << endl;
		flankString = (indexInfo->chromStr)[chromNO].substr(SJdonerEnd, 2) + (indexInfo->chromStr)[chromNO].substr(SJacceptorStart-3, 2);
		//cout << "flankString: " << flankString << endl;
	}
};


class SamFormat
{
public:

	int FLAG;
	int MAPQ;
	string RNEXT;
	int PNEXT;
	int TLEN;
	int IH;
	int HI;

};


class Alignment_Info
{
public:
	string alignDirection;

	unsigned int mapPosInWholeGenome;

	string alignChromName;
	int alignChromPos;
	vector<Jump_Code> cigarStringJumpCode;
	vector< pair< int, SpliceJunction_Alignment > > spliceJunctionVec; // <posInRead, SpliceJunction>
	int mismatchNum;
	string SJstrand;
	int mappedBaseNum;
	int endMatchedPosInChr;
	
	SamFormat currentSam;

	//vector< pair< int, int > > SJforCompareVec; 

	Alignment_Info()
	{}

	int unfixedHeadLength()
	{
		return cigarStringJumpCode[0].len;
	}

	int unfxiedTailLength()
	{
		int jumpCodeNum = cigarStringJumpCode.size();
		return cigarStringJumpCode[jumpCodeNum - 1].len;
	}
	bool unfixedHeadExistsBool()
	{
		int jumpCodeNum = cigarStringJumpCode.size();
		if(jumpCodeNum < 0)
		{
			return true;
		}

		bool unfixedHeadBool = (cigarStringJumpCode[0].type == "S");
		return unfixedHeadBool;
	}

	bool unfixedTailExistsBool()
	{
		int jumpCodeNum = cigarStringJumpCode.size();
		if(jumpCodeNum < 0)
		{
			return true;
		}

		bool unfixedHeadBool = (cigarStringJumpCode[jumpCodeNum-1].type == "S");
		return unfixedHeadBool;		
	}

	bool noUnfixedHeadTailBool()
	{
		int jumpCodeNum = cigarStringJumpCode.size();
		if(jumpCodeNum < 0)
		{
			return false;
		}

		bool unfixedHeadBool = (cigarStringJumpCode[0].type == "S");
		bool unfixedTailBool = (cigarStringJumpCode[jumpCodeNum - 1].type == "S");

		return  !unfixedHeadBool && !unfixedTailBool;
	}

	Alignment_Info(string alignDir, const string& mapChromName, 
		int mapChromPos, vector<Jump_Code>& cigarString, int mismatch, Index_Info* indexInfo)
	{
		//cout << "stop 9" << endl;
		alignDirection = alignDir;
		alignChromName = mapChromName;		
		alignChromPos = mapChromPos;
		//cout << "stop 10" << endl;
		cigarStringJumpCode = cigarString;
		//cout << "stop 11" << endl;
		this->jumpCodeVec2spliceJunctionVec(indexInfo);
		//cout << "stop 12" << endl;
		mismatchNum = mismatch;
		//cout << "stop 13" << endl;
		SJstrand = this->getStrandFromSJ();
		//cout << "stop 14" << endl;
		//currentSam.setSamFormat(this, this);
		//SJstrand
	}	

	Alignment_Info(const string& alignInfoStr, const string& alignDir)
	{
		alignDirection = alignDir;

		int tmpFieldStartPosInStr = 0;
		int tmpFieldEndPosInStr;
		string tmpFieldStr;

		tmpFieldEndPosInStr = alignInfoStr.find(",",tmpFieldStartPosInStr) - 1;
		tmpFieldStr = alignInfoStr.substr(tmpFieldStartPosInStr, 
				tmpFieldEndPosInStr - tmpFieldStartPosInStr + 1);
		alignChromName = tmpFieldStr;

		tmpFieldStartPosInStr = tmpFieldEndPosInStr + 2;
		tmpFieldEndPosInStr = alignInfoStr.find(",",tmpFieldStartPosInStr) - 1;
		tmpFieldStr = alignInfoStr.substr(tmpFieldStartPosInStr, 
				tmpFieldEndPosInStr - tmpFieldStartPosInStr + 1);
		alignChromPos = atoi(tmpFieldStr.c_str());

		tmpFieldStartPosInStr = tmpFieldEndPosInStr + 2;
		tmpFieldEndPosInStr = alignInfoStr.find(",",tmpFieldStartPosInStr) - 1;
		tmpFieldStr = alignInfoStr.substr(tmpFieldStartPosInStr, 
				tmpFieldEndPosInStr - tmpFieldStartPosInStr + 1);		
		this->jumpCodeStr2jumpCodeVec(tmpFieldStr);

		tmpFieldStartPosInStr = tmpFieldEndPosInStr + 2;
		tmpFieldEndPosInStr = alignInfoStr.find(",",tmpFieldStartPosInStr) - 1;
		tmpFieldStr = alignInfoStr.substr(tmpFieldStartPosInStr, 
				tmpFieldEndPosInStr - tmpFieldStartPosInStr + 1);
		mismatchNum = atoi(tmpFieldStr.c_str());
	}

	Alignment_Info* newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(int firstMatchLength, int spliceJunctionDistance, Index_Info* indexInfo, int newMismatchToAddInHead)
	{
		//cout << "start to get new alignInfo" << endl;
		string newAlignDir = alignDirection;
		string newAlignChromName = alignChromName;
		int newAlignChromPos = alignChromPos - spliceJunctionDistance - cigarStringJumpCode[0].len;
		
		//cout << "start to creat new jumpCode " << endl;
		vector<Jump_Code> newCigarStringJumpCode;
		//cout << "stop 0 " << endl;
		Jump_Code* tmpJumpCode = new Jump_Code(firstMatchLength, "M");
		//tmpJumpCode.len = firstMatchLength; tmpJumpCode.type = "M";
		//cout << "stop 1 " << endl;
		newCigarStringJumpCode.push_back(*tmpJumpCode);
		tmpJumpCode->len = spliceJunctionDistance; tmpJumpCode->type = "N";
		//cout << "stop 2 " << endl;
		newCigarStringJumpCode.push_back(*tmpJumpCode);
		tmpJumpCode->len = cigarStringJumpCode[0].len + 
			cigarStringJumpCode[1].len - firstMatchLength; tmpJumpCode->type = "M";
		//cout << "stop 3 " << endl;
		newCigarStringJumpCode.push_back(*tmpJumpCode);
		//cout << "stop 4 " << endl;
		for(int tmp = 2; tmp < cigarStringJumpCode.size(); tmp++)
		{
			newCigarStringJumpCode.push_back(cigarStringJumpCode[tmp]);
		}
		
		int newMismatch = mismatchNum + newMismatchToAddInHead;
		return new Alignment_Info(newAlignDir, newAlignChromName,
			newAlignChromPos, newCigarStringJumpCode, newMismatch, indexInfo);
	}

	Alignment_Info* newAlignInfoAfterAddTailInfo2SoftClippingAlignment(int lastMatchLength, int spliceJunctionDistance, Index_Info* indexInfo, int newMismatchToAddInTail)
	{
		//cout << "start to get new alignInfo" << endl;
		string newAlignDir = alignDirection;
		string newAlignChromName = alignChromName;
		int newAlignChromPos = alignChromPos;// - spliceJunctionDistance - cigarStringJumpCode[0].len;
		
		//cout << "start to creat new jumpCode " << endl;
		vector<Jump_Code> newCigarStringJumpCode;
		for(int tmp = 0; tmp < cigarStringJumpCode.size()-2; tmp++)
		{
			newCigarStringJumpCode.push_back(cigarStringJumpCode[tmp]);
		}

		int jumpCodeSize = cigarStringJumpCode.size();
		int penultimateMatchLength 
			= cigarStringJumpCode[jumpCodeSize-2].len + cigarStringJumpCode[jumpCodeSize-1].len - lastMatchLength; 
		Jump_Code* tmpJumpCode = new Jump_Code(penultimateMatchLength, "M");
		newCigarStringJumpCode.push_back(*tmpJumpCode);
		tmpJumpCode->len = spliceJunctionDistance; tmpJumpCode->type = "N";
		newCigarStringJumpCode.push_back(*tmpJumpCode);
		tmpJumpCode->len = lastMatchLength; tmpJumpCode->type = "M";
		newCigarStringJumpCode.push_back(*tmpJumpCode);

		int newMismatch = mismatchNum + newMismatchToAddInTail;

		Alignment_Info* newAlignInfo = new Alignment_Info(newAlignDir, newAlignChromName,
			newAlignChromPos, newCigarStringJumpCode, newMismatch, indexInfo);

		return newAlignInfo;	
	}

	void jumpCodeStr2jumpCodeVec(const string& jumpCodeStr)
	{
		int tmpJumpCodeLength;
		string tmpJumpCodeType;

		int jumpCodeStartPosInCigarStr = 0;
		int jumpCodeEndPosInCigarStr;
		
		string candidateJumpCodeType = "SMNID";
		while(1)
		{
			jumpCodeEndPosInCigarStr = 
				jumpCodeStr.find_first_of(candidateJumpCodeType, jumpCodeStartPosInCigarStr);
			if(jumpCodeEndPosInCigarStr == jumpCodeStr.npos)
				{break;}
			else
			{
				tmpJumpCodeLength = 
					atoi((jumpCodeStr.substr(jumpCodeStartPosInCigarStr, jumpCodeEndPosInCigarStr - jumpCodeStartPosInCigarStr)).c_str());
				tmpJumpCodeType = jumpCodeStr.substr(jumpCodeEndPosInCigarStr, 1);
				cigarStringJumpCode.push_back(Jump_Code(tmpJumpCodeLength, tmpJumpCodeType));
				jumpCodeStartPosInCigarStr = jumpCodeEndPosInCigarStr + 1;
			}
		}
	}

	void jumpCodeVec2spliceJunctionVec(Index_Info* indexInfo) 
	{
		int tmpPosInRead = 0; 
		int tmpDonerPosInChr = alignChromPos - 1;
		int tmpAcceptorPosInChr = 0;
		//cout << "start to get sjVec ..." << endl;
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			//cout << "JumpCode[" << tmp << "]: " << cigarStringJumpCode[tmp].len << " " << cigarStringJumpCode[tmp].type << endl;
			if(cigarStringJumpCode[tmp].type == "S")
			{
				tmpPosInRead += cigarStringJumpCode[tmp].len;
			}
			else if (cigarStringJumpCode[tmp].type == "M")
			{
				tmpPosInRead += cigarStringJumpCode[tmp].len;
				//cout << "tmpPosInRead: " << tmpPosInRead << endl;
				tmpDonerPosInChr += cigarStringJumpCode[tmp].len;
				//cout << "tmpDonerPosInChr: " << tmpDonerPosInChr << endl;
			}
			else if (cigarStringJumpCode[tmp].type == "N")
			{
				//cout << "tmpDonerPosInChr: " << tmpDonerPosInChr << endl;

				tmpAcceptorPosInChr = tmpDonerPosInChr + cigarStringJumpCode[tmp].len + 1;
				//cout << "tmpAcceptorPosInChr: " << tmpAcceptorPosInChr << endl;
				SpliceJunction_Alignment* tmpSJ = 
					new SpliceJunction_Alignment(alignChromName, tmpDonerPosInChr, tmpAcceptorPosInChr, indexInfo);
				//cout << "tmpPosInRead: " << tmpPosInRead << endl;
				//tmpSJ
				spliceJunctionVec.push_back(pair <int, SpliceJunction_Alignment> (tmpPosInRead+1, (*tmpSJ)));
				tmpDonerPosInChr += cigarStringJumpCode[tmp].len;
			}
			else if (cigarStringJumpCode[tmp].type == "I")
			{
				tmpPosInRead += cigarStringJumpCode[tmp].len;
				//tmpDonerPosInChr += cigarStringJumpCode[tmp].len;
			}
			else if (cigarStringJumpCode[tmp].type == "D")
			{
				tmpDonerPosInChr += cigarStringJumpCode[tmp].len;
			}
			else
			{
				cout << "other jumpCodeType" << endl;
			}
		}
	}

	string getStrandFromSJ()
	{	
		string tmpFlankString;
		string currentStrand = "N";
		for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
		{
			tmpFlankString = (spliceJunctionVec[tmp].second).flankString;
			if((tmpFlankString == "ATAC") || (tmpFlankString == "CTGC") || (tmpFlankString == "CTAC"))
			{
				if(currentStrand == "+")
					return "X";
				else
					currentStrand = "-";
			}
			else
			{
				if(currentStrand == "-")
					return "X";
				else
					return "+";
			}
		}
		return currentStrand;
	}

	string jumpCodeVec2Str()
	{
		string str;
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{	
			str += cigarStringJumpCode[tmp].toString();
		}
		return str;
	}

	string otherJumpCodeVec2Str() // only for Head
	{
		string str;
		for(int tmp = 1; tmp < cigarStringJumpCode.size(); tmp++)
		{	
			str += cigarStringJumpCode[tmp].toString();
		}
		return str;
	}

	string otherJumpCodeVec2StrForTail() // only for Tail
	{
		string str;
		for(int tmp = 0; tmp < cigarStringJumpCode.size()-1; tmp++)
		{	
			str += cigarStringJumpCode[tmp].toString();
		}
		return str;
	}

	int mappedLength()
	{
		int pos = 0;
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			if ((cigarStringJumpCode[tmp].type == "M")||(cigarStringJumpCode[tmp].type == "I"))
			{
				pos += cigarStringJumpCode[tmp].len;
			}
			else
			{
				continue;
				//cout << "error, unexpected jumpCodeType: " << cigarStringJumpCode[tmp].type << endl;
			}
		}
		return pos;
	}

	int getEndMatchedPosInChr()
	{
		int pos = 0;//alignChromPos;
		//int tmpJumpCodeLength = 0;
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			if ((cigarStringJumpCode[tmp].type == "S")||(cigarStringJumpCode[tmp].type == "I"))
			{
				continue;
			}
			else if ((cigarStringJumpCode[tmp].type == "M")||(cigarStringJumpCode[tmp].type == "N")||(cigarStringJumpCode[tmp].type == "D"))
			{
				pos += cigarStringJumpCode[tmp].len;
			}
			else
			{
				cout << "error, unexpected jumpCodeType: " << cigarStringJumpCode[tmp].type << endl;
			}
		}

		endMatchedPosInChr = (alignChromPos + pos - 1);
		return endMatchedPosInChr;
		//return (alignChromPos + pos - 1);
	}

	string getSamFormatStr (const string& readName, const string& readSeq,
		const string& qualitySeq)
	{
		string samString;

		//string QNAME;
		int FLAG = (this->currentSam).FLAG;
		string RNAME = alignChromName;
		int POS = alignChromPos;
		int MAPQ = (this->currentSam).MAPQ;
		string CIGAR = this->jumpCodeVec2Str();
		string RNEXT = (this->currentSam).RNEXT;
		int PNEXT = (this->currentSam).PNEXT;
		int TLEN = (this->currentSam).TLEN;
		int NM = mismatchNum;
		int IH = (this->currentSam).IH;
		int HI = (this->currentSam).HI;
		string XS = SJstrand;
		string XF;
		for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
		{
			XF = XF + (spliceJunctionVec[tmp].second).flankString + ","; 
		}		

		samString = readName + "\t" + Utilities::int_to_str(FLAG) + "\t" + RNAME + "\t"
			+ Utilities::int_to_str(POS) + "\t" + Utilities::int_to_str(MAPQ) + "\t" + CIGAR
			+ "\t" + RNEXT + "\t" + Utilities::int_to_str(PNEXT) + "\t"
			+ Utilities::int_to_str(TLEN) + "\t" + readSeq + "\t" + qualitySeq + "\tNM:i:"+
			Utilities::int_to_str(NM) + "\tIH:i:" + Utilities::int_to_str(IH) + "\tHI:i:" + Utilities::int_to_str(HI)
			+ "\tXS:A:" + XS + "\tXF:Z:" + XF;


		return samString;
	}

	int getFlag(bool pairEndReadsBool, bool bothEndMappedBool, bool readUnmappedBool, bool anotherEndReadUnmappedBool,
		bool mappedAsRcmBool, bool anotherEndReadMappedAsRcmBool, bool pairEnd_1_bool, bool pairEnd_2_bool, 
		bool notPrimaryAlignmentBool, bool notPassQualControlBool, bool PCRorOpticalDuplicateBool, bool suppleAlignmentBool)
	{
		int flagInt = 0;
		if(pairEndReadsBool)
			flagInt += 1;
		if(bothEndMappedBool)
			flagInt += 2;
		if(readUnmappedBool)
			flagInt += 4;
		if(anotherEndReadUnmappedBool)
			flagInt += 8;
		if(mappedAsRcmBool)
			flagInt += 16;
		if(anotherEndReadMappedAsRcmBool)
			flagInt += 32;
		if(pairEnd_1_bool)
			flagInt += 64;
		if(pairEnd_2_bool)
			flagInt += 128;
		if(notPrimaryAlignmentBool)
			flagInt += 256;
		if(notPassQualControlBool)
			flagInt += 512;
		if(PCRorOpticalDuplicateBool)
			flagInt += 1024;
		if(suppleAlignmentBool)
			flagInt += 2048;
	
		return flagInt;
	}

	int getFlag_unpaired(bool mapDirection,
		bool End1OrEnd2Bool )
	{
		bool pairEndsReadBool = true;
		bool bothEndMappedBool = false;
		bool readUnmappedBool = false;
		bool anotherEndReadUnmappedBool = true;

		bool mappedAsRcmBool, anotherEndReadMappedAsRcmBool;
		if(mapDirection)
		{
			mappedAsRcmBool = false;
			anotherEndReadMappedAsRcmBool = false;
		}
		else
		{
			mappedAsRcmBool = true;
			anotherEndReadMappedAsRcmBool = false;
		}

		bool pairEnd_1_bool, pairEnd_2_bool;
		if(End1OrEnd2Bool)
		{
			pairEnd_1_bool = true;
			pairEnd_2_bool = false;
		}
		else
		{
			pairEnd_1_bool = false;
			pairEnd_2_bool = true;
		}	

		bool notPrimaryAlignmentBool = false;
		bool notPassQualControlBool = false;
		bool PCRorOpticalDuplicateBool = false;
		bool suppleAlignmentBool = false;

		int flagInt = this->getFlag(pairEndsReadBool, bothEndMappedBool, readUnmappedBool, anotherEndReadUnmappedBool,
			mappedAsRcmBool, anotherEndReadMappedAsRcmBool, pairEnd_1_bool, pairEnd_2_bool, 
			notPrimaryAlignmentBool, notPassQualControlBool, PCRorOpticalDuplicateBool, suppleAlignmentBool); 

		return flagInt;		

	}

	int getFlag_unpaired_secondaryOrNot(bool mapDirection,
		bool End1OrEnd2Bool, bool SecondaryOrNot )
	{
		bool pairEndsReadBool = true;
		bool bothEndMappedBool = false;
		bool readUnmappedBool = false;
		bool anotherEndReadUnmappedBool = true;

		bool mappedAsRcmBool, anotherEndReadMappedAsRcmBool;
		if(mapDirection)
		{
			mappedAsRcmBool = false;
			anotherEndReadMappedAsRcmBool = false;
		}
		else
		{
			mappedAsRcmBool = true;
			anotherEndReadMappedAsRcmBool = false;
		}

		bool pairEnd_1_bool, pairEnd_2_bool;
		if(End1OrEnd2Bool)
		{
			pairEnd_1_bool = true;
			pairEnd_2_bool = false;
		}
		else
		{
			pairEnd_1_bool = false;
			pairEnd_2_bool = true;
		}	

		bool notPrimaryAlignmentBool = SecondaryOrNot;
		bool notPassQualControlBool = false;
		bool PCRorOpticalDuplicateBool = false;
		bool suppleAlignmentBool = false;

		int flagInt = this->getFlag(pairEndsReadBool, bothEndMappedBool, readUnmappedBool, anotherEndReadUnmappedBool,
			mappedAsRcmBool, anotherEndReadMappedAsRcmBool, pairEnd_1_bool, pairEnd_2_bool, 
			notPrimaryAlignmentBool, notPassQualControlBool, PCRorOpticalDuplicateBool, suppleAlignmentBool); 

		return flagInt;		

	}

	int getFlag_paired(bool mapDirection, 
		bool End1OrEnd2Bool)
	{
		//int flagInt = 4;

		bool pairEndsReadBool = true;
		bool bothEndMappedBool = true;
		bool readUnmappedBool = false;
		bool anotherEndReadUnmappedBool = false;
		
		bool mappedAsRcmBool, anotherEndReadMappedAsRcmBool;
		if(mapDirection)
		{
			mappedAsRcmBool = false;
			anotherEndReadMappedAsRcmBool = true;
		}
		else
		{
			mappedAsRcmBool = true;
			anotherEndReadMappedAsRcmBool = false;
		}

		bool pairEnd_1_bool, pairEnd_2_bool;
		if(End1OrEnd2Bool)
		{
			pairEnd_1_bool = true;
			pairEnd_2_bool = false;
		}
		else
		{
			pairEnd_1_bool = false;
			pairEnd_2_bool = true;
		}

		bool notPrimaryAlignmentBool = false;
		bool notPassQualControlBool = false;
		bool PCRorOpticalDuplicateBool = false;
		bool suppleAlignmentBool = false;

		int flagInt = this->getFlag(pairEndsReadBool, bothEndMappedBool, readUnmappedBool, anotherEndReadUnmappedBool,
			mappedAsRcmBool, anotherEndReadMappedAsRcmBool, pairEnd_1_bool, pairEnd_2_bool, 
			notPrimaryAlignmentBool, notPassQualControlBool, PCRorOpticalDuplicateBool, suppleAlignmentBool); 

		return flagInt;
	}

	int getFlag_paired_secondaryOrNot(bool mapDirection, 
		bool End1OrEnd2Bool, bool SecondaryOrNot)
	{
		//int flagInt = 4;

		bool pairEndsReadBool = true;
		bool bothEndMappedBool = true;
		bool readUnmappedBool = false;
		bool anotherEndReadUnmappedBool = false;
		
		bool mappedAsRcmBool, anotherEndReadMappedAsRcmBool;
		if(mapDirection)
		{
			mappedAsRcmBool = false;
			anotherEndReadMappedAsRcmBool = true;
		}
		else
		{
			mappedAsRcmBool = true;
			anotherEndReadMappedAsRcmBool = false;
		}

		bool pairEnd_1_bool, pairEnd_2_bool;
		if(End1OrEnd2Bool)
		{
			pairEnd_1_bool = true;
			pairEnd_2_bool = false;
		}
		else
		{
			pairEnd_1_bool = false;
			pairEnd_2_bool = true;
		}

		bool notPrimaryAlignmentBool = SecondaryOrNot;
		bool notPassQualControlBool = false;
		bool PCRorOpticalDuplicateBool = false;
		bool suppleAlignmentBool = false;

		int flagInt = this->getFlag(pairEndsReadBool, bothEndMappedBool, readUnmappedBool, anotherEndReadUnmappedBool,
			mappedAsRcmBool, anotherEndReadMappedAsRcmBool, pairEnd_1_bool, pairEnd_2_bool, 
			notPrimaryAlignmentBool, notPassQualControlBool, PCRorOpticalDuplicateBool, suppleAlignmentBool); 

		return flagInt;
	}

	string getSamFormatString(const string& readName, const string& readSeq)
	{
		string samString;

		//string QNAME;
		int FLAG;
		string RNAME;
		int POS;

		///////////////////// get FLAG, RNAME and POS ///////////////////////////
		if(alignDirection == "+")
		{
			RNAME = alignChromName;
			POS = alignChromPos;
			FLAG = 0;

		}
		else if(alignDirection == "-")
		{
			RNAME = alignChromName;
			POS = alignChromPos;			
			FLAG = 16;
		}
		else
		{
			RNAME = "*";
			POS = 0;
			FLAG = 4;
		}
		////////////////////////////////////////////////////////////////////////

		int MAPQ = 255;
		string CIGAR = this->jumpCodeVec2Str();
		string RNEXT = "*";
		int PNEXT = 0;
		int TLEN = 0;
		//string SEQ;
		string QUAL = "*";
		string strandStr = SJstrand;


		samString = readName + "\t" + 
			Utilities::int_to_str(FLAG) + "\t" + RNAME + "\t"
			+ Utilities::int_to_str(POS) + "\t" + Utilities::int_to_str(MAPQ) + "\t" + CIGAR
			+ "\t" + RNEXT + "\t" + Utilities::int_to_str(PNEXT) + "\t"
			+ Utilities::int_to_str(TLEN) + "\t" + readSeq + "\t" +
			QUAL + "\tNM:i:" + Utilities::int_to_str(mismatchNum) + "\tXS:A:" + strandStr + "\tXF:Z:";
		for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
		{
			samString = samString + (spliceJunctionVec[tmp].second).flankString + ","; 
		}

		return samString;
	}

	string getSamFormatString_paired(const string& readName, const string& readSeq, 
		Alignment_Info* mateAlignInfo, bool End1OrEnd2)
	{
		string samString;

		//return samString;

		int FLAG;
		string RNAME;
		int POS;

		///////////////////// get FLAG, RNAME and POS ///////////////////////////
		if(alignDirection == "+")
		{
			RNAME = alignChromName;
			POS = alignChromPos;
			FLAG = getFlag_paired(true, End1OrEnd2);//0;
		}
		else if(alignDirection == "-")
		{
			RNAME = alignChromName;
			POS = alignChromPos;			
			FLAG = getFlag_paired(false, End1OrEnd2); //16;
		}
		else
		{
			RNAME = "*";
			POS = 0;
			FLAG = 4;
		}
		////////////////////////////////////////////////////////////////////////

		int MAPQ = 255;
		string CIGAR = this->jumpCodeVec2Str();
		string RNEXT = "="; //mateAlignInfo->alignChromName; //"*";
		int PNEXT = mateAlignInfo->alignChromPos; // 0;
		int TLEN = 0;
		//string SEQ;
		string QUAL = "*";
		string strandStr = SJstrand;

		string FLAGstr; 
		string POSstr; 
		string MAPQstr; 
		string PNEXTstr; 
		string TLENstr; 
		string mismatchNumStr;
		
		FLAGstr = Utilities::int_to_str(FLAG);
		POSstr = Utilities::int_to_str(POS);
		MAPQstr = Utilities::int_to_str(MAPQ);
		PNEXTstr = Utilities::int_to_str(PNEXT);
		TLENstr = Utilities::int_to_str(TLEN);
		mismatchNumStr = Utilities::int_to_str(mismatchNum);

		samString = readName + "\t" 
			+ FLAGstr + "\t" 
			+ RNAME + "\t"
			+ POSstr + "\t" 
			+ MAPQstr + "\t" 
			+ CIGAR + "\t" 
			+ RNEXT + "\t" 
			+ PNEXTstr + "\t" 
			+ TLENstr + "\t" 
			+ readSeq + "\t" 
			+ QUAL 
			+ "\tNM:i:" 
			+ mismatchNumStr + "\tXS:A:" + strandStr + "\tXF:Z:"
			;


		for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
		{
			samString = samString + (spliceJunctionVec[tmp].second).flankString + ","; 
		}

		return samString;
	}

	string getSamFormatString_unpaired(const string& readName, const string& readSeq, 
		//Alignment_Info* mateAlignInfo, 
		bool End1OrEnd2, int IH_num, int HI_num)
	{
		string samString;

		int FLAG;
		string RNAME;
		int POS;
		///////////////////// get FLAG, RNAME and POS ///////////////////////////
		if(alignDirection == "+")
		{
			RNAME = alignChromName;
			POS = alignChromPos;
			FLAG = getFlag_unpaired(true, End1OrEnd2);//0;
		}
		else if(alignDirection == "-")
		{
			RNAME = alignChromName;
			POS = alignChromPos;			
			FLAG = getFlag_unpaired(false, End1OrEnd2); //16;
		}
		else
		{
			RNAME = "*";
			POS = 0;
			FLAG = 4;
		}
		////////////////////////////////////////////////////////////////////////
		int MAPQ = 255;
		string CIGAR = this->jumpCodeVec2Str();
		string RNEXT = "*"; //mateAlignInfo->alignChromName; //"*";
		int PNEXT = 0;//mateAlignInfo->alignChromPos; // 0;
		int TLEN = 0;
		//string SEQ;
		string QUAL = "*";
		string strandStr = SJstrand;

		string FLAGstr; 
		string POSstr; 
		string MAPQstr; 
		string PNEXTstr; 
		string TLENstr; 
		string mismatchNumStr;
		string IHstr;
		string HIstr;

		FLAGstr = Utilities::int_to_str(FLAG);
		POSstr = Utilities::int_to_str(POS);
		MAPQstr = Utilities::int_to_str(MAPQ);
		PNEXTstr = Utilities::int_to_str(PNEXT);
		TLENstr = Utilities::int_to_str(TLEN);
		mismatchNumStr = Utilities::int_to_str(mismatchNum);
		IHstr = Utilities::int_to_str(IH_num);
		HIstr = Utilities::int_to_str(HI_num);

		samString = readName + "\t" 
			+ FLAGstr + "\t" 
			+ RNAME + "\t"
			+ POSstr + "\t" 
			+ MAPQstr + "\t" 
			+ CIGAR + "\t" 
			+ RNEXT + "\t" 
			+ PNEXTstr + "\t" 
			+ TLENstr + "\t" 
			+ readSeq + "\t" 
			+ QUAL 
			+ "\tNM:i:" + mismatchNumStr 
			+ "\tIH:i:" + IHstr
			+ "\tHI:i:" + HIstr
			+ "\tXS:A:" + strandStr 
			+ "\tXF:Z:"
			;


		for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
		{
			samString = samString + (spliceJunctionVec[tmp].second).flankString + ","; 
		}

		return samString;

	}

	string getSamFormatString_paired(const string& readName, const string& readSeq, 
		Alignment_Info* mateAlignInfo, bool End1OrEnd2, int IH_num, int HI_num)
	{
		string samString;

		int FLAG;
		string RNAME;
		int POS;

		///////////////////// get FLAG, RNAME and POS ///////////////////////////
		if(alignDirection == "+")
		{
			RNAME = alignChromName;
			POS = alignChromPos;
			FLAG = getFlag_paired(true, End1OrEnd2);//0;
		}
		else if(alignDirection == "-")
		{
			RNAME = alignChromName;
			POS = alignChromPos;			
			FLAG = getFlag_paired(false, End1OrEnd2); //16;
		}
		else
		{
			RNAME = "*";
			POS = 0;
			FLAG = 4;
		}
		////////////////////////////////////////////////////////////////////////

		int MAPQ = 255;
		string CIGAR = this->jumpCodeVec2Str();
		string RNEXT = "="; //mateAlignInfo->alignChromName; //"*";
		int PNEXT = mateAlignInfo->alignChromPos; // 0;
		int TLEN = 0;
		//string SEQ;
		string QUAL = "*";
		string strandStr = SJstrand;

		string FLAGstr; 
		string POSstr; 
		string MAPQstr; 
		string PNEXTstr; 
		string TLENstr; 
		string mismatchNumStr;
		string IHstr;
		string HIstr;

		FLAGstr = Utilities::int_to_str(FLAG);
		POSstr = Utilities::int_to_str(POS);
		MAPQstr = Utilities::int_to_str(MAPQ);
		PNEXTstr = Utilities::int_to_str(PNEXT);
		TLENstr = Utilities::int_to_str(TLEN);
		mismatchNumStr = Utilities::int_to_str(mismatchNum);
		IHstr = Utilities::int_to_str(IH_num);
		HIstr = Utilities::int_to_str(HI_num);

		samString = readName + "\t" 
			+ FLAGstr + "\t" 
			+ RNAME + "\t"
			+ POSstr + "\t" 
			+ MAPQstr + "\t" 
			+ CIGAR + "\t" 
			+ RNEXT + "\t" 
			+ PNEXTstr + "\t" 
			+ TLENstr + "\t" 
			+ readSeq + "\t" 
			+ QUAL 
			+ "\tNM:i:" + mismatchNumStr 
			+ "\tIH:i:" + IHstr
			+ "\tHI:i:" + HIstr
			+ "\tXS:A:" + strandStr 
			+ "\tXF:Z:"
			;


		for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
		{
			samString = samString + (spliceJunctionVec[tmp].second).flankString + ","; 
		}

		return samString;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	string getSamFormatString_unpaired_secondaryOrNot(const string& readName, const string& readSeq, 
		//Alignment_Info* mateAlignInfo, 
		bool End1OrEnd2, int IH_num, int HI_num, bool secondaryOrNot)
	{
		string samString;

		int FLAG;
		string RNAME;
		int POS;
		///////////////////// get FLAG, RNAME and POS ///////////////////////////
		if(alignDirection == "+")
		{
			RNAME = alignChromName;
			POS = alignChromPos;
			FLAG = getFlag_unpaired_secondaryOrNot(true, End1OrEnd2, secondaryOrNot);//0;
		}
		else if(alignDirection == "-")
		{
			RNAME = alignChromName;
			POS = alignChromPos;			
			FLAG = getFlag_unpaired_secondaryOrNot(false, End1OrEnd2, secondaryOrNot); //16;
		}
		else
		{
			RNAME = "*";
			POS = 0;
			FLAG = 4;
		}
		////////////////////////////////////////////////////////////////////////
		int MAPQ = 255;
		string CIGAR = this->jumpCodeVec2Str();
		string RNEXT = "*"; //mateAlignInfo->alignChromName; //"*";
		int PNEXT = 0;//mateAlignInfo->alignChromPos; // 0;
		int TLEN = 0;
		//string SEQ;
		string QUAL = "*";
		string strandStr = SJstrand;

		string FLAGstr; 
		string POSstr; 
		string MAPQstr; 
		string PNEXTstr; 
		string TLENstr; 
		string mismatchNumStr;
		string IHstr;
		string HIstr;

		FLAGstr = Utilities::int_to_str(FLAG);
		POSstr = Utilities::int_to_str(POS);
		MAPQstr = Utilities::int_to_str(MAPQ);
		PNEXTstr = Utilities::int_to_str(PNEXT);
		TLENstr = Utilities::int_to_str(TLEN);
		mismatchNumStr = Utilities::int_to_str(mismatchNum);
		IHstr = Utilities::int_to_str(IH_num);
		HIstr = Utilities::int_to_str(HI_num);

		samString = readName + "\t" 
			+ FLAGstr + "\t" 
			+ RNAME + "\t"
			+ POSstr + "\t" 
			+ MAPQstr + "\t" 
			+ CIGAR + "\t" 
			+ RNEXT + "\t" 
			+ PNEXTstr + "\t" 
			+ TLENstr + "\t" 
			+ readSeq + "\t" 
			+ QUAL 
			+ "\tNM:i:" + mismatchNumStr 
			+ "\tIH:i:" + IHstr
			+ "\tHI:i:" + HIstr
			+ "\tXS:A:" + strandStr 
			+ "\tXF:Z:"
			;


		for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
		{
			samString = samString + (spliceJunctionVec[tmp].second).flankString + ","; 
		}

		return samString;

	}

	string getSamFormatString_paired_secondaryOrNot(const string& readName, const string& readSeq, 
		Alignment_Info* mateAlignInfo, bool End1OrEnd2, int IH_num, int HI_num, bool SecondaryOrNot, int templateLength)
	{
		string samString;

		int FLAG;
		string RNAME;
		int POS;

		int TLEN = 0;
		///////////////////// get FLAG, RNAME and POS ///////////////////////////
		if(alignDirection == "+")
		{
			RNAME = alignChromName;
			POS = alignChromPos;
			FLAG = getFlag_paired_secondaryOrNot(true, End1OrEnd2, SecondaryOrNot);//0;
			TLEN = templateLength;
		}
		else if(alignDirection == "-")
		{
			RNAME = alignChromName;
			POS = alignChromPos;			
			FLAG = getFlag_paired_secondaryOrNot(false, End1OrEnd2, SecondaryOrNot); //16;
			TLEN = 0 - templateLength;
		}
		else
		{
			RNAME = "*";
			POS = 0;
			FLAG = 4;
		}
		////////////////////////////////////////////////////////////////////////

		int MAPQ = 255;
		string CIGAR = this->jumpCodeVec2Str();
		string RNEXT = "="; //mateAlignInfo->alignChromName; //"*";
		int PNEXT = mateAlignInfo->alignChromPos; // 0;
		//int TLEN = 0;
		//string SEQ;
		string QUAL = "*";
		string strandStr = SJstrand;

		string FLAGstr; 
		string POSstr; 
		string MAPQstr; 
		string PNEXTstr; 
		string TLENstr; 
		string mismatchNumStr;
		string IHstr;
		string HIstr;

		FLAGstr = Utilities::int_to_str(FLAG);
		POSstr = Utilities::int_to_str(POS);
		MAPQstr = Utilities::int_to_str(MAPQ);
		PNEXTstr = Utilities::int_to_str(PNEXT);
		TLENstr = Utilities::int_to_str(TLEN);
		mismatchNumStr = Utilities::int_to_str(mismatchNum);
		IHstr = Utilities::int_to_str(IH_num);
		HIstr = Utilities::int_to_str(HI_num);

		samString = readName + "\t" 
			+ FLAGstr + "\t" 
			+ RNAME + "\t"
			+ POSstr + "\t" 
			+ MAPQstr + "\t" 
			+ CIGAR + "\t" 
			+ RNEXT + "\t" 
			+ PNEXTstr + "\t" 
			+ TLENstr + "\t" 
			+ readSeq + "\t" 
			+ QUAL 
			+ "\tNM:i:" + mismatchNumStr 
			+ "\tIH:i:" + IHstr
			+ "\tHI:i:" + HIstr
			+ "\tXS:A:" + strandStr 
			+ "\tXF:Z:"
			;


		for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
		{
			samString = samString + (spliceJunctionVec[tmp].second).flankString + ","; 
		}

		return samString;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	bool checkOverlapPairAlignment(Alignment_Info* secondAlignmentInfo) // parameter is another alignment 
	{
		bool SJinOverlapAreaNoConflict = true;

		int overlapStartPos = secondAlignmentInfo->alignChromPos;
		int overlapEndPos = this->endMatchedPosInChr;

		int tmpDonerPosInChr, tmpAcceptorPosInChr;
		vector< pair<int, int> > norAlignSJvecInOverlapArea;
		vector< pair<int, int> > rcmAlignSJvecInOverlapArea;
		for(int tmp = 0; tmp < (this->spliceJunctionVec).size(); tmp++)
		{
			if(((this->spliceJunctionVec)[tmp].second).SJdonerEnd <= overlapStartPos)
			{
				tmpDonerPosInChr = ((this->spliceJunctionVec)[tmp].second).SJdonerEnd;
				tmpAcceptorPosInChr = ((this->spliceJunctionVec)[tmp].second).SJacceptorStart;
				norAlignSJvecInOverlapArea.push_back(pair<int,int> (tmpDonerPosInChr, tmpAcceptorPosInChr));
			}
			else
			{}
		} 

		for(int tmp = 0; tmp < (secondAlignmentInfo->spliceJunctionVec).size(); tmp++)
		{
			if(((secondAlignmentInfo->spliceJunctionVec)[tmp].second).SJacceptorStart >= overlapEndPos)
			{
				tmpDonerPosInChr = ((secondAlignmentInfo->spliceJunctionVec)[tmp].second).SJdonerEnd;
				tmpAcceptorPosInChr = ((secondAlignmentInfo->spliceJunctionVec)[tmp].second).SJacceptorStart;
				rcmAlignSJvecInOverlapArea.push_back(pair<int,int> (tmpDonerPosInChr, tmpAcceptorPosInChr));
			}
			else
			{}
		} 

		if(norAlignSJvecInOverlapArea.size() == rcmAlignSJvecInOverlapArea.size())
		{
			for(int tmp = 0; tmp < norAlignSJvecInOverlapArea.size(); tmp++)
			{
				if((norAlignSJvecInOverlapArea[tmp].first == rcmAlignSJvecInOverlapArea[tmp].first)
					&&(norAlignSJvecInOverlapArea[tmp].second == rcmAlignSJvecInOverlapArea[tmp].second))
				{
					SJinOverlapAreaNoConflict = true;
				}
				else
				{
					return false;
				}
			}
		}
		else
		{
			return false;
		}

		return SJinOverlapAreaNoConflict;
	}

	bool checkOverlapPairAlignment_new(Alignment_Info* secondAlignmentInfo) // parameter is another alignment 
	{
		bool SJinOverlapAreaNoConflict = true;

		int overlapStartPos_1 = this->alignChromPos;
		int overlapStartPos_2 = secondAlignmentInfo->alignChromPos;
		int overlapEndPos_1 = this->endMatchedPosInChr;
		int overlapEndPos_2 = secondAlignmentInfo->endMatchedPosInChr;

		int overlapStartPos;
		int overlapEndPos;
		
		if(overlapStartPos_1 > overlapStartPos)
			overlapStartPos = overlapStartPos_1;
		else
			overlapStartPos = overlapStartPos_2;

		if(overlapEndPos_1 < overlapEndPos_2)
			overlapEndPos = overlapEndPos_1;
		else
			overlapEndPos = overlapEndPos_2;


		int tmpDonerPosInChr, tmpAcceptorPosInChr;
		vector< pair<int, int> > norAlignSJvecInOverlapArea;
		vector< pair<int, int> > rcmAlignSJvecInOverlapArea;
		for(int tmp = 0; tmp < (this->spliceJunctionVec).size(); tmp++)
		{
			if(((this->spliceJunctionVec)[tmp].second).SJdonerEnd <= overlapStartPos)
			{
				tmpDonerPosInChr = ((this->spliceJunctionVec)[tmp].second).SJdonerEnd;
				tmpAcceptorPosInChr = ((this->spliceJunctionVec)[tmp].second).SJacceptorStart;
				norAlignSJvecInOverlapArea.push_back(pair<int,int> (tmpDonerPosInChr, tmpAcceptorPosInChr));
			}
			else
			{}
		} 

		for(int tmp = 0; tmp < (secondAlignmentInfo->spliceJunctionVec).size(); tmp++)
		{
			if(((secondAlignmentInfo->spliceJunctionVec)[tmp].second).SJacceptorStart >= overlapEndPos)
			{
				tmpDonerPosInChr = ((secondAlignmentInfo->spliceJunctionVec)[tmp].second).SJdonerEnd;
				tmpAcceptorPosInChr = ((secondAlignmentInfo->spliceJunctionVec)[tmp].second).SJacceptorStart;
				rcmAlignSJvecInOverlapArea.push_back(pair<int,int> (tmpDonerPosInChr, tmpAcceptorPosInChr));
			}
			else
			{}
		} 

		if(norAlignSJvecInOverlapArea.size() == rcmAlignSJvecInOverlapArea.size())
		{
			for(int tmp = 0; tmp < norAlignSJvecInOverlapArea.size(); tmp++)
			{
				if((norAlignSJvecInOverlapArea[tmp].first == rcmAlignSJvecInOverlapArea[tmp].first)
					&&(norAlignSJvecInOverlapArea[tmp].second == rcmAlignSJvecInOverlapArea[tmp].second))
				{
					SJinOverlapAreaNoConflict = true;
				}
				else
				{
					return false;
				}
			}
		}
		else
		{
			return false;
		}

		return SJinOverlapAreaNoConflict;
	}
};

class PE_Read_Alignment_Info
{
public:

	vector<Alignment_Info*> norAlignmentInfo_PE_1;
	vector<Alignment_Info*> rcmAlignmentInfo_PE_1;
	vector<Alignment_Info*> norAlignmentInfo_PE_2;
	vector<Alignment_Info*> rcmAlignmentInfo_PE_2;

	vector< pair< int, vector<int> > > oriAlignPair_Nor1Rcm2;
	vector< pair< int, vector<int> > > oriAlignPair_Nor2Rcm1;

	vector< pair< int, int > > oriAlignPair_Nor1Rcm2_filtered;
	vector< pair< int, int > > oriAlignPair_Nor2Rcm1_filtered;	

	vector< pair< int, int > > finalAlignPair_Nor1Rcm2;
	vector< pair< int, int > > finalAlignPair_Nor2Rcm1;

	vector<bool> otherEndUnmappedBoolVec;

	bool betterNewOtherEndAlignInfoBool(
		//int oriAlignInfo_mappedLength, int newAlignInfo_mappedLength,
		//int oriAlignInfo_mismatchNum, int newAlignInfo_mismatchNum
		//Alignment_Info* fixedEndAlignInfo,
		Alignment_Info* otherEndAlignInfo_ori, Alignment_Info* otherEndAlignInfo_new)
	{
		int mappedLength_ori = otherEndAlignInfo_ori->mappedLength();
		int mappedLength_new = otherEndAlignInfo_new->mappedLength();

		int mismatchNum_ori = otherEndAlignInfo_ori->mismatchNum;
		int mismatchNum_new = otherEndAlignInfo_new->mismatchNum;

		int endMappedPos_ori = otherEndAlignInfo_ori->endMatchedPosInChr;
		int endmappedPos_new = otherEndAlignInfo_new->endMatchedPosInChr;

		if((mappedLength_ori - mismatchNum_ori) - (mappedLength_new - mismatchNum_new) < 0)
		{
			return true;
		}
		else if((mappedLength_ori - mismatchNum_ori) - (mappedLength_new - mismatchNum_new) > 0)
		{
			return false;
		}
		else // ((mappedLength_ori - mismatchNum_ori) - (mappedLength_new - mismatchNum_new) == 0)
		{
			if(endMappedPos_ori < endmappedPos_new)
			{
				return false;
			}
			else
			{
				return true;
			}
		}
	}

	int filterOriPair(bool oriPairNor1Rcm2OrNot, int oriPairVecNO)
	{
		int selectedAlignInfoNO = 0;

		if(oriPairNor1Rcm2OrNot)
		{
			int alignInfo_Nor1_NO = oriAlignPair_Nor1Rcm2[oriPairVecNO].first;
			selectedAlignInfoNO = (oriAlignPair_Nor1Rcm2[oriPairVecNO].second)[0];

			for(int tmp = 1; tmp < (oriAlignPair_Nor1Rcm2[oriPairVecNO].second).size(); tmp++)
			{
				int alignInfo_Rcm2_NO = (oriAlignPair_Nor1Rcm2[oriPairVecNO].second)[tmp];
				if( this->betterNewOtherEndAlignInfoBool(
					//norAlignmentInfo_PE_1[alignInfo_Nor1_NO], 
					rcmAlignmentInfo_PE_2[selectedAlignInfoNO], rcmAlignmentInfo_PE_2[alignInfo_Rcm2_NO] ) )
				{
					selectedAlignInfoNO = alignInfo_Rcm2_NO;
				}
			}
		}
		else
		{
			int alignInfo_Nor2_NO = oriAlignPair_Nor2Rcm1[oriPairVecNO].first;
			selectedAlignInfoNO = (oriAlignPair_Nor2Rcm1[oriPairVecNO].second)[0];

			for(int tmp = 1; tmp < (oriAlignPair_Nor2Rcm1[oriPairVecNO].second).size(); tmp++)
			{
				int alignInfo_Rcm1_NO = (oriAlignPair_Nor2Rcm1[oriPairVecNO].second)[tmp];
				if( this->betterNewOtherEndAlignInfoBool(
					//norAlignmentInfo_PE_2[alignInfo_Nor2_NO],
					rcmAlignmentInfo_PE_1[selectedAlignInfoNO], rcmAlignmentInfo_PE_1[alignInfo_Rcm1_NO]) )
				{
					selectedAlignInfoNO = alignInfo_Rcm1_NO;
				}
			}
		}
		
		return selectedAlignInfoNO;
	}

	void generatePeReadInfoAndPeAlignInfo_Fasta_toFixOneEndUnmapped_getline(const string& line1, const string& line2, 
		//const string& recordLine3, 
		const string& line4, const string& line5,
		//const string& recordLine6, 
		const string& line7, const string& line8,
		const string& line9, const string& line10, PairedEndRead* peReadInfo)
	{
		int Nor1Num = 0, Rcm1Num = 0, Nor2Num = 0, Rcm2Num = 0;
		string readNameStr_1, readNameStr_2;

		int startSearchPos = 0, foundSearchPos;	foundSearchPos = line1.find("\t", startSearchPos);
		readNameStr_1 = line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1);
			
		startSearchPos = foundSearchPos + 1; foundSearchPos = line1.find("\t", startSearchPos);
		Nor1Num = atoi((line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());		
			
		startSearchPos = foundSearchPos + 1; foundSearchPos = line1.find("\t", startSearchPos);		
		Rcm1Num = atoi((line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());			

		startSearchPos = 0; foundSearchPos = line4.find("\t", startSearchPos);
		readNameStr_2 = line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1);
				
		startSearchPos = foundSearchPos + 1; foundSearchPos = line4.find("\t", startSearchPos);
		Nor2Num = atoi((line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());		
				
		startSearchPos = foundSearchPos + 1; foundSearchPos = line4.find("\t", startSearchPos);		
		Rcm2Num = atoi((line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());	

		int readLength_1 = line2.length() - 1;
		int readLength_2 = line5.length() - 1;

		peReadInfo->setReadData(readNameStr_1, readNameStr_2, line2, line5);
		this->getPeReadAlignmentInfo(line7, line8, line9, line10,
			Nor1Num, Rcm1Num, Nor2Num, Rcm2Num);
	}		

	void generatePeReadInfoAndPeAlignInfo_Fasta_toFixIncompleteAlignment_getline(const string& line1, const string& line2,
		//const string& recordLine3, 
		const string& line4, const string& line5,
		//const string& recordLine6, 
		const string& line7, const string& line8,
		const string& line9, const string& line10, PairedEndRead* peReadInfo)
	{
		int Nor1Num = 0, Rcm1Num = 0, Nor2Num = 0, Rcm2Num = 0;
		string readNameStr_1, readNameStr_2;

		int startSearchPos = 0, foundSearchPos;	foundSearchPos = line1.find("\t", startSearchPos);
		readNameStr_1 = line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1);
	
		startSearchPos = foundSearchPos + 1; foundSearchPos = line1.find("\t", startSearchPos);
		Nor1Num = atoi((line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());		
		
		startSearchPos = foundSearchPos + 1; foundSearchPos = line1.find("\t", startSearchPos);		
		Rcm1Num = atoi((line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());			

		startSearchPos = 0; foundSearchPos = line4.find("\t", startSearchPos);
		readNameStr_2 = line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1);
				
		startSearchPos = foundSearchPos + 1; foundSearchPos = line4.find("\t", startSearchPos);
		Nor2Num = atoi((line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());		
				
		startSearchPos = foundSearchPos + 1; foundSearchPos = line4.find("\t", startSearchPos);		
		Rcm2Num = atoi((line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());	

		peReadInfo->setReadData(readNameStr_1, readNameStr_2, line2, line5);
		this->getPeReadAlignmentInfo(line7, line8, line9, line10,
			Nor1Num, Rcm1Num, Nor2Num, Rcm2Num);
	}

	~PE_Read_Alignment_Info()
	{
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
			delete norAlignmentInfo_PE_1[tmp];

		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
			delete rcmAlignmentInfo_PE_1[tmp];

		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
			delete norAlignmentInfo_PE_2[tmp];

		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
			delete rcmAlignmentInfo_PE_2[tmp];
	}

	bool finalPairExistsBool()
	{
		return finalAlignPair_Nor1Rcm2.size() + finalAlignPair_Nor2Rcm1.size() != 0;
	}

	bool oriPairExistsBool()
	{
		return oriAlignPair_Nor1Rcm2.size() + oriAlignPair_Nor2Rcm1.size() != 0;
	}

	bool alignInfoExistsBool()
	{
		int alignInfoNum = norAlignmentInfo_PE_1.size() + rcmAlignmentInfo_PE_1.size()
			+ norAlignmentInfo_PE_2.size() + rcmAlignmentInfo_PE_2.size();
		return alignInfoNum != 0;
	}

	bool allAlignmentInFinalPairCompleted()
	{
		for(int i=0; i<finalAlignPair_Nor1Rcm2.size(); i++)
		{
			if(!norAlignmentInfo_PE_1[(finalAlignPair_Nor1Rcm2[i].first)]->noUnfixedHeadTailBool())
				return false;

			if(!rcmAlignmentInfo_PE_2[(finalAlignPair_Nor1Rcm2[i].second)]->noUnfixedHeadTailBool())
				return false;
		}

		for(int i=0; i<finalAlignPair_Nor2Rcm1.size(); i++)
		{
			if(!norAlignmentInfo_PE_2[(finalAlignPair_Nor2Rcm1[i].first)]->noUnfixedHeadTailBool())
				return false;

			if(!rcmAlignmentInfo_PE_1[(finalAlignPair_Nor2Rcm1[i].second)]->noUnfixedHeadTailBool())
				return false;
		}

		return true;
	}

	bool allUnpairedAlignmentCompleted()
	{
		for(int i = 0; i < norAlignmentInfo_PE_1.size(); i++)
		{
			bool tmpBool = 
				norAlignmentInfo_PE_1[i]->noUnfixedHeadTailBool();
			if(!tmpBool)
				return false;
		}

		for(int i = 0; i < norAlignmentInfo_PE_2.size(); i++)
		{
			bool tmpBool = 
				norAlignmentInfo_PE_2[i]->noUnfixedHeadTailBool();
			if(!tmpBool)
				return false;			
		}

		for(int i = 0; i < rcmAlignmentInfo_PE_1.size(); i++)
		{
			bool tmpBool = 
				rcmAlignmentInfo_PE_1[i]->noUnfixedHeadTailBool();
			if(!tmpBool)
				return false;
		}

		for(int i = 0; i < rcmAlignmentInfo_PE_2.size(); i++)
		{
			bool tmpBool = 
				rcmAlignmentInfo_PE_2[i]->noUnfixedHeadTailBool();
			if(!tmpBool)
				return false;			
		}

		return true;
	}
	
	void getPeReadAlignmentInfo(const string& line5, const string& line6, const string& line7,
		const string& line8, int nor1Num, int rcm1Num, int nor2Num, int rcm2Num)
	{
		//Alignment_Info* tmpAlignmentInfo;
		//int tmp 
		int tmpAlignInfoStrStartPos;
		int tmpAlignInfoStrEndPos;
		string tmpAlignInfoStr;

		norAlignmentInfo_PE_1.clear();
		rcmAlignmentInfo_PE_1.clear();
		norAlignmentInfo_PE_2.clear();
		rcmAlignmentInfo_PE_2.clear();		
		//cout << 
		if(nor1Num < 40)
		{
			tmpAlignInfoStrStartPos = line5.find("\t", 0) + 1;
			for(int tmp = 0; tmp < nor1Num; tmp++)
			{
				tmpAlignInfoStrEndPos = line5.find("\t", tmpAlignInfoStrStartPos) - 1;
				tmpAlignInfoStr = line5.substr(tmpAlignInfoStrStartPos, 
					tmpAlignInfoStrEndPos - tmpAlignInfoStrStartPos + 1);
				//cout << "tmpAlignInfoStr: " << tmpAlignInfoStr << endl;
				Alignment_Info* tmpAlignInfo = new Alignment_Info(tmpAlignInfoStr, "+");
				norAlignmentInfo_PE_1.push_back(tmpAlignInfo);

				tmpAlignInfoStrStartPos = tmpAlignInfoStrEndPos + 2;
			}
		}
		
		if(rcm1Num < 40)
		{	
			tmpAlignInfoStrStartPos = line6.find("\t", 0) + 1;
			for(int tmp = 0; tmp < rcm1Num; tmp++)
			{
				tmpAlignInfoStrEndPos = line6.find("\t", tmpAlignInfoStrStartPos) - 1;
				tmpAlignInfoStr = line6.substr(tmpAlignInfoStrStartPos, 
					tmpAlignInfoStrEndPos - tmpAlignInfoStrStartPos + 1);
				//cout << "tmpAlignInfoStr: " << tmpAlignInfoStr << endl;
				Alignment_Info* tmpAlignInfo = new Alignment_Info(tmpAlignInfoStr, "-");
				rcmAlignmentInfo_PE_1.push_back(tmpAlignInfo);

				tmpAlignInfoStrStartPos = tmpAlignInfoStrEndPos + 2;
			}
		}	
		//cout << "norAlignmentInfo_PE_2.size():  " << norAlignmentInfo_PE_2.size() << endl;
		//cout << "norNum_2: " << nor2Num << endl;
		if(nor2Num < 40)
		{
			//cout << "line7: " << line7 << endl;
			tmpAlignInfoStrStartPos = line7.find("\t", 0) + 1;
			for(int tmp = 0; tmp < nor2Num; tmp++)
			{
				//cout << "start at:" << tmpAlignInfoStrStartPos << " end at:" 
				//	<< tmpAlignInfoStrEndPos << endl;
				tmpAlignInfoStrEndPos = line7.find("\t", tmpAlignInfoStrStartPos) - 1;
				tmpAlignInfoStr = line7.substr(tmpAlignInfoStrStartPos, 
					tmpAlignInfoStrEndPos - tmpAlignInfoStrStartPos + 1);
				//cout << "tmpAlignInfoStr: " << tmpAlignInfoStr << endl;
				Alignment_Info* tmpAlignInfo = new Alignment_Info(tmpAlignInfoStr, "+");
				norAlignmentInfo_PE_2.push_back(tmpAlignInfo);

				tmpAlignInfoStrStartPos = tmpAlignInfoStrEndPos + 2;
			}
		}
		//cout << "norAlignmentInfo_PE_2.size():  " << norAlignmentInfo_PE_2.size() << endl;
		//cout << "rcmAlignmentInfo_PE_2.size():  " << rcmAlignmentInfo_PE_2.size() << endl;
		//cout << "rcm2Num: " << rcm2Num << endl;
		if(rcm2Num < 40)
		{
			tmpAlignInfoStrStartPos = line8.find("\t", 0) + 1;
			for(int tmp = 0; tmp < rcm2Num; tmp++)
			{
				tmpAlignInfoStrEndPos = line8.find("\t", tmpAlignInfoStrStartPos) - 1;
				tmpAlignInfoStr = line8.substr(tmpAlignInfoStrStartPos, 
					tmpAlignInfoStrEndPos - tmpAlignInfoStrStartPos + 1);
				//cout << "tmpAlignInfoStr: " << tmpAlignInfoStr << endl;
				Alignment_Info* tmpAlignInfo = new Alignment_Info(tmpAlignInfoStr, "-");
				rcmAlignmentInfo_PE_2.push_back(tmpAlignInfo);

				tmpAlignInfoStrStartPos = tmpAlignInfoStrEndPos + 2;
			}
		}
		//cout << "rcmAlignmentInfo_PE_2.size():  " << rcmAlignmentInfo_PE_2.size() << endl;
	}

	void pushBackPathInfo2PeAlignInfo(Path* newPathInfo, bool End1OrEnd2,
		bool NorOrRcm, Index_Info* indexInfo) // Note: End1OrEnd2, NorOrRcm are both describing the other end alignment
	{
		//cout << "End1OrEnd2: " << End1OrEnd2 << endl;
		//cout << "NorOrRcm: " << NorOrRcm << endl;
		//cout << "PathSize: " << (newPathInfo->finalPathVec).size() << endl;
		if(End1OrEnd2 && NorOrRcm)
		{
			for(int tmpPath = 0; tmpPath < newPathInfo->finalPathVec.size(); tmpPath++)
			{
				int mapChromNameInt = (((newPathInfo->finalPathVec)[tmpPath]).first).first;
				int mapChromPosInt = (((newPathInfo->finalPathVec)[tmpPath]).first).second;
				int tmpMismatch = 0;

				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("-", indexInfo->chrNameStr[mapChromNameInt],
					mapChromPosInt, (((newPathInfo->finalPathVec)[tmpPath]).second)->final_jump_code,
					tmpMismatch, indexInfo);
				rcmAlignmentInfo_PE_2.push_back(tmpAlignmentInfo);
			}
		}
		else if(End1OrEnd2 && (!NorOrRcm))
		{
			for(int tmpPath = 0; tmpPath < newPathInfo->finalPathVec.size(); tmpPath++)
			{
				int mapChromNameInt = (((newPathInfo->finalPathVec)[tmpPath]).first).first;
				int mapChromPosInt = (((newPathInfo->finalPathVec)[tmpPath]).first).second;
				int tmpMismatch = 0;

				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->chrNameStr[mapChromNameInt],
					mapChromPosInt, (((newPathInfo->finalPathVec)[tmpPath]).second)->final_jump_code,
					tmpMismatch, indexInfo);
				norAlignmentInfo_PE_2.push_back(tmpAlignmentInfo);
			}
		}
		else if((!End1OrEnd2) && NorOrRcm)
		{
			for(int tmpPath = 0; tmpPath < newPathInfo->finalPathVec.size(); tmpPath++)
			{
				int mapChromNameInt = (((newPathInfo->finalPathVec)[tmpPath]).first).first;
				int mapChromPosInt = (((newPathInfo->finalPathVec)[tmpPath]).first).second;
				int tmpMismatch = 0;

				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("-", indexInfo->chrNameStr[mapChromNameInt],
					mapChromPosInt, (((newPathInfo->finalPathVec)[tmpPath]).second)->final_jump_code,
					tmpMismatch, indexInfo);
				rcmAlignmentInfo_PE_1.push_back(tmpAlignmentInfo);
			}
		}
		else
		{
			for(int tmpPath = 0; tmpPath < newPathInfo->finalPathVec.size(); tmpPath++)
			{
				int mapChromNameInt = (((newPathInfo->finalPathVec)[tmpPath]).first).first;
				int mapChromPosInt = (((newPathInfo->finalPathVec)[tmpPath]).first).second;
				int tmpMismatch = 0;

				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->chrNameStr[mapChromNameInt],
					mapChromPosInt, (((newPathInfo->finalPathVec)[tmpPath]).second)->final_jump_code, 
					tmpMismatch, indexInfo);
				norAlignmentInfo_PE_1.push_back(tmpAlignmentInfo);
			}
		}

	}

	PE_Read_Alignment_Info(Path* pathInfo_nor1, Path* pathInfo_rcm1, 
			Path* pathInfo_nor2, Path* pathInfo_rcm2, Index_Info* indexInfo)
	{
		//cout << "start PE_Read_Alignment_Info " << endl;
		for(int tmpPath = 0; tmpPath < pathInfo_nor1->finalPathVec.size(); tmpPath++)
		{
			int mapChromNameInt = (((pathInfo_nor1->finalPathVec)[tmpPath]).first).first;
			int mapChromPosInt = (((pathInfo_nor1->finalPathVec)[tmpPath]).first).second;
			int tmpMismatch = (pathInfo_nor1->fixedPathMismatchVec)[tmpPath];

			Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->chrNameStr[mapChromNameInt], 
				mapChromPosInt, (((pathInfo_nor1->finalPathVec)[tmpPath]).second)->final_jump_code, 
				tmpMismatch, indexInfo);
			norAlignmentInfo_PE_1.push_back(tmpAlignmentInfo);
		}
		//cout << "pe alignInfo Nor1 ends ..." << endl;
		for(int tmpPath = 0; tmpPath < pathInfo_rcm1->finalPathVec.size(); tmpPath++)
		{
			int mapChromNameInt = (((pathInfo_rcm1->finalPathVec)[tmpPath]).first).first;
			int mapChromPosInt = (((pathInfo_rcm1->finalPathVec)[tmpPath]).first).second;
			int tmpMismatch = (pathInfo_rcm1->fixedPathMismatchVec)[tmpPath];

			Alignment_Info* tmpAlignmentInfo = new Alignment_Info("-", indexInfo->chrNameStr[mapChromNameInt], 
				mapChromPosInt, (((pathInfo_rcm1->finalPathVec)[tmpPath]).second)->final_jump_code, 
				tmpMismatch, indexInfo);
			rcmAlignmentInfo_PE_1.push_back(tmpAlignmentInfo);
		}
		//cout << "pe alignInfo Rcm1 ends ..." << endl;
		for(int tmpPath = 0; tmpPath < pathInfo_nor2->finalPathVec.size(); tmpPath++)
		{
			int mapChromNameInt = (((pathInfo_nor2->finalPathVec)[tmpPath]).first).first;
			int mapChromPosInt = (((pathInfo_nor2->finalPathVec)[tmpPath]).first).second;
			int tmpMismatch = (pathInfo_nor2->fixedPathMismatchVec)[tmpPath];

			Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->chrNameStr[mapChromNameInt], 
				mapChromPosInt, (((pathInfo_nor2->finalPathVec)[tmpPath]).second)->final_jump_code, 
				tmpMismatch, indexInfo);
			norAlignmentInfo_PE_2.push_back(tmpAlignmentInfo);
		}
		//cout << "pe alignInfo Nor2 ends ..." << endl;
		for(int tmpPath = 0; tmpPath < pathInfo_rcm2->finalPathVec.size(); tmpPath++)
		{
			int mapChromNameInt = (((pathInfo_rcm2->finalPathVec)[tmpPath]).first).first;
			int mapChromPosInt = (((pathInfo_rcm2->finalPathVec)[tmpPath]).first).second;
			int tmpMismatch = (pathInfo_rcm2->fixedPathMismatchVec)[tmpPath];

			Alignment_Info* tmpAlignmentInfo = new Alignment_Info("-", indexInfo->chrNameStr[mapChromNameInt], 
				mapChromPosInt, (((pathInfo_rcm2->finalPathVec)[tmpPath]).second)->final_jump_code, 
				tmpMismatch, indexInfo);
			rcmAlignmentInfo_PE_2.push_back(tmpAlignmentInfo);
		}
		//cout << "pe alignInfo Rcm2 ends ..." << endl;
	}

	PE_Read_Alignment_Info()
	{}

	Alignment_Info* fixHeadTail_getAlignInfo(
		int replacedAlignInfoKind, int replacedAlignInfoNO)// 1. Nor_1; 2. Rcm_1; 3. Nor_2; 4. Rcm_2
	{
		if(replacedAlignInfoKind == 1)
		{
			return norAlignmentInfo_PE_1[replacedAlignInfoNO];
		}
		else if(replacedAlignInfoKind == 2)
		{
			return rcmAlignmentInfo_PE_1[replacedAlignInfoNO];
		}
		else if(replacedAlignInfoKind == 3)
		{
			return norAlignmentInfo_PE_2[replacedAlignInfoNO];
		}
		else if(replacedAlignInfoKind == 4)
		{
			return rcmAlignmentInfo_PE_2[replacedAlignInfoNO];
		}
		else
		{
			cout << "error in PE_Align_Info.fixHeadTail_getAlignInfo()" << endl;
		}
	}

	Alignment_Info* getAlignInfo(
		int replacedAlignInfoKind, int replacedAlignInfoNO)// 1. Nor_1; 2. Rcm_1; 3. Nor_2; 4. Rcm_2
	{
		if(replacedAlignInfoKind == 1)
		{
			return norAlignmentInfo_PE_1[replacedAlignInfoNO];
		}
		else if(replacedAlignInfoKind == 2)
		{
			return rcmAlignmentInfo_PE_1[replacedAlignInfoNO];
		}
		else if(replacedAlignInfoKind == 3)
		{
			return norAlignmentInfo_PE_2[replacedAlignInfoNO];
		}
		else if(replacedAlignInfoKind == 4)
		{
			return rcmAlignmentInfo_PE_2[replacedAlignInfoNO];
		}
		else
		{
			cout << "error in PE_Align_Info.fixHeadTail_getAlignInfo()" << endl;
		}
	}

	void fixHeadTail_replaceWithNewAlignInfo(Alignment_Info* newAlignmentInfo, 
		int replacedAlignInfoKind, int replacedAlignInfoNO)// 1. Nor_1; 2. Rcm_1; 3. Nor_2; 4. Rcm_2
	{
		if(replacedAlignInfoKind == 1)
		{
			norAlignmentInfo_PE_1[replacedAlignInfoNO] = newAlignmentInfo;
		}
		else if(replacedAlignInfoKind == 2)
		{
			rcmAlignmentInfo_PE_1[replacedAlignInfoNO] = newAlignmentInfo;
		}
		else if(replacedAlignInfoKind == 3)
		{
			norAlignmentInfo_PE_2[replacedAlignInfoNO] = newAlignmentInfo;
		}
		else if(replacedAlignInfoKind == 4)
		{
			rcmAlignmentInfo_PE_2[replacedAlignInfoNO] = newAlignmentInfo;
		}
		else
		{
			cout << "error in PE_Align_Info.fixHeadTail_replaceWithNewAlignInfo()" << endl;
		}
	}

	void fixHeadTail_addNewAlignInfo(Alignment_Info* newAlignmentInfo, 
		int addedAlignInfoKind)// 1. Nor_1; 2. Rcm_1; 3. Nor_2; 4. Rcm_2
	{
		if(addedAlignInfoKind == 1)
		{
			norAlignmentInfo_PE_1.push_back(newAlignmentInfo);
		}
		else if(addedAlignInfoKind == 2)
		{
			rcmAlignmentInfo_PE_1.push_back(newAlignmentInfo);
		}
		else if(addedAlignInfoKind == 3)
		{
			norAlignmentInfo_PE_2.push_back(newAlignmentInfo);
		}
		else if(addedAlignInfoKind == 4)
		{
			rcmAlignmentInfo_PE_2.push_back(newAlignmentInfo);
		}
		else
		{
			cout << "error in PE_Align_Info.fixHeadTail_addNewAlignInfo()" << endl;
		}
	}

	int fixHeadTail_getAlignInfoVecSize(int alignInfoKind)// 1. Nor_1; 2. Rcm_1; 3. Nor_2; 4. Rcm_2
	{
		int size = 100;
		if(alignInfoKind == 1)
		{
			size = norAlignmentInfo_PE_1.size();
		}
		else if(alignInfoKind == 2)
		{
			size = rcmAlignmentInfo_PE_1.size();
		}
		else if(alignInfoKind == 3)
		{
			size = norAlignmentInfo_PE_2.size();
		}
		else if(alignInfoKind == 4)
		{
			size = rcmAlignmentInfo_PE_2.size();
		}
		else
		{
			cout << "error in PE_Align_Info.fixHeadTail_getAlignInfoVecSize()" << endl;
		}		
		return size;
	}

	int getAlignInfoVecSize(int alignInfoKind)// 1. Nor_1; 2. Rcm_1; 3. Nor_2; 4. Rcm_2
	{
		//cout << "getAlignInfoVecSize starts ... " << endl;
		//cout << "alignInfoKind: " << alignInfoKind << endl;
		int size = 100;
		if(alignInfoKind == 1)
		{
			size = norAlignmentInfo_PE_1.size();
		}
		else if(alignInfoKind == 2)
		{
			size = rcmAlignmentInfo_PE_1.size();
		}
		else if(alignInfoKind == 3)
		{
			size = norAlignmentInfo_PE_2.size();
		}
		else if(alignInfoKind == 4)
		{
			size = rcmAlignmentInfo_PE_2.size();
		}
		else
		{
			cout << "error in PE_Align_Info.fixHeadTail_getAlignInfoVecSize()" << endl;
		}		
		//cout << "size: " << size << endl;
		return size;
	}

	void getEndMatchPosForEveryAlignment()
	{
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			norAlignmentInfo_PE_1[tmp]->getEndMatchedPosInChr();
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
		{
			rcmAlignmentInfo_PE_1[tmp]->getEndMatchedPosInChr();
		}
		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			norAlignmentInfo_PE_2[tmp]->getEndMatchedPosInChr();
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
		{
			rcmAlignmentInfo_PE_2[tmp]->getEndMatchedPosInChr();
		}
	}

	void pairingAlignment() 
	// get final PE alignment finalXXXAlignmentInfo_PE_X  
	// from original alignment info XXXAlignmentInfo_PE_X;
	{
		//cout << "start to pair ... " << endl;
		this->getEndMatchPosForEveryAlignment();
		Alignment_Info* tmpAlignInfo_1;
		Alignment_Info* tmpAlignInfo_2;
		bool newEntity = false;
		//cout << "start to pair Nor1Rcm2... " << endl;
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++) 
		//pair norAlignmentInfo_PE_1 & rcmAlignmentInfo_PE_2
		{
			newEntity = true;
			tmpAlignInfo_1 = norAlignmentInfo_PE_1[tmp];

			if(tmpAlignInfo_1 -> SJstrand == "X")
			{
				continue;
			}

			for(int tmp2 = 0; tmp2 < rcmAlignmentInfo_PE_2.size(); tmp2++)
			{

				tmpAlignInfo_2 = rcmAlignmentInfo_PE_2[tmp2];
				
				if((tmpAlignInfo_2 -> SJstrand == "X")
					||((tmpAlignInfo_1 -> SJstrand == "+")&&(tmpAlignInfo_2 -> SJstrand == "-"))
					||((tmpAlignInfo_1 -> SJstrand == "-")&&(tmpAlignInfo_2 -> SJstrand == "+")))
				{
					continue;
				}		

				if((tmpAlignInfo_1->alignChromName) == (tmpAlignInfo_2->alignChromName))
				{
					if(tmpAlignInfo_1->endMatchedPosInChr <= tmpAlignInfo_2->alignChromPos)
					{
						if( ((tmpAlignInfo_2->alignChromPos) - (tmpAlignInfo_1->endMatchedPosInChr))< PAIR_READ_DISTANCE_MAX)
						{
							if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (tmp, newTmpVec));
								newEntity = false;
							}
							else //tmp has already been in oriAlignPair_Nor1Rcm2
							{
								(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
							}
						}
						else
						{
							//two far away 
						}
					}
					else
					{
						if((tmpAlignInfo_1->alignChromPos <= tmpAlignInfo_2->alignChromPos)
							&&(tmpAlignInfo_1->endMatchedPosInChr <= tmpAlignInfo_2->endMatchedPosInChr))
						{
							//Note: In addition should check whether they cross the same SJs or not
							if(tmpAlignInfo_1->checkOverlapPairAlignment(tmpAlignInfo_2))
							{
								if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
								{
									vector<int> newTmpVec;
									newTmpVec.push_back(tmp2);
									oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (tmp, newTmpVec));
									newEntity = false;
								}
								else //tmp has already been in oriAlignPair_Nor1Rcm2
								{
									(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
								}
							}
							else
							{}
						}
						else 
						{
							// not in the correct order
						}
					}
				}
				else
				{

				}
			}
		}
		//cout << "start to pair Nor2Rcm1... " << endl;
		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			newEntity = true;
			tmpAlignInfo_1 = norAlignmentInfo_PE_2[tmp];
			
			if(tmpAlignInfo_1 -> SJstrand == "X")
			{
				continue;
			}

			for(int tmp2 = 0; tmp2 < rcmAlignmentInfo_PE_1.size(); tmp2++)
			{
				tmpAlignInfo_2 = rcmAlignmentInfo_PE_1[tmp2];

				if((tmpAlignInfo_2 -> SJstrand == "X")
					||((tmpAlignInfo_1 -> SJstrand == "+")&&(tmpAlignInfo_2 -> SJstrand == "-"))
					||((tmpAlignInfo_1 -> SJstrand == "-")&&(tmpAlignInfo_2 -> SJstrand == "+")))
				{
					continue;
				}					
				
				if((tmpAlignInfo_1->alignChromName) == (tmpAlignInfo_2->alignChromName))
				{
					if(tmpAlignInfo_1->endMatchedPosInChr <= tmpAlignInfo_2->alignChromPos)
					{
						if(((tmpAlignInfo_2->alignChromPos) - (tmpAlignInfo_1->endMatchedPosInChr)) < 500000)
						{
							if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor2Rcm1.push_back(pair<int, vector<int> > (tmp, newTmpVec));
								newEntity = false;
							}
							else //tmp has already been in oriAlignPair_Nor1Rcm2
							{
								(oriAlignPair_Nor2Rcm1[oriAlignPair_Nor2Rcm1.size()-1].second).push_back(tmp2);
								newEntity = false;
							}
						}
						else
						{
							//two far away 
						}
					}
					else
					{
						if((tmpAlignInfo_1->alignChromPos <= tmpAlignInfo_2->alignChromPos)
							&&(tmpAlignInfo_1->endMatchedPosInChr <= tmpAlignInfo_2->endMatchedPosInChr))
						{
							//Note: In addition should check whether they cross the same SJs or not
							if(tmpAlignInfo_1->checkOverlapPairAlignment(tmpAlignInfo_2))
							{
								if(newEntity) //tmp is not in oriAlignPair_Nor2Rcm1
								{
									vector<int> newTmpVec;
									newTmpVec.push_back(tmp2);
									oriAlignPair_Nor2Rcm1.push_back(pair<int, vector<int> > (tmp, newTmpVec));
									newEntity = false;
								}
								else //tmp has already been in oriAlignPair_Nor2Rcm1
								{
									(oriAlignPair_Nor2Rcm1[oriAlignPair_Nor2Rcm1.size()-1].second).push_back(tmp2);
									newEntity = false;
								}
							}
							else
							{}
						}
						else 
						{
							// not in the correct order
						}
					}
				}
				else
				{

				}
			}
		}
		///////////////// 1. only one end read mapped, another one unmapped ///////////////////////

		///////////////// 2. reads can be mapped under both directions //////////////////////////

		///////////////// 3. reads can be mapped to different places in one direction //////////////////
	}

	void chooseBestAlignment()
	{
		map < int, int > mapPos_Nor1Rcm2; // <alignChormPos of Nor1[#], # in oriAlignPair_Nor1Rcm2  >
		map < int, int > mapPos_Nor2Rcm1; // <alignChormPos of Nor2[#], # in oriAlignPair_Nor2Rcm1  >
		map < int, int >::iterator it;

		int tmpNor1NO, tmpNor2NO, tmpRcm1NO, tmpRcm2NO;
		int tmpNor1NO_2, tmpNor2NO_2;

		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2.size(); tmp++)
		{
			tmpNor1NO = oriAlignPair_Nor1Rcm2[tmp].first;
			it = mapPos_Nor1Rcm2.find(norAlignmentInfo_PE_1[tmpNor1NO]->alignChromPos);
			if(it == mapPos_Nor1Rcm2.end())
			{
				mapPos_Nor1Rcm2.insert(pair<int, int> (norAlignmentInfo_PE_1[tmpNor1NO]->alignChromPos, tmp));
			}
			else // the same start mapping pos
			{
				tmpNor1NO_2 = oriAlignPair_Nor1Rcm2[it->second].first;
				if((norAlignmentInfo_PE_1[tmpNor1NO]->endMatchedPosInChr)
					< (norAlignmentInfo_PE_1[tmpNor1NO_2]->endMatchedPosInChr))
				{
					mapPos_Nor1Rcm2.erase(it);
					mapPos_Nor1Rcm2.insert(pair<int, int> (norAlignmentInfo_PE_1[tmpNor1NO]->alignChromPos, tmp));
				}
				else
				{}
			}
		}

		for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1.size(); tmp++)
		{
			tmpNor2NO = oriAlignPair_Nor2Rcm1[tmp].first;
			it = mapPos_Nor2Rcm1.find(norAlignmentInfo_PE_2[tmpNor2NO]->alignChromPos);
			if(it == mapPos_Nor2Rcm1.end())
			{
				mapPos_Nor2Rcm1.insert(pair<int, int> (norAlignmentInfo_PE_2[tmpNor2NO]->alignChromPos, tmp));
			}
			else // the same start mapping pos
			{
				tmpNor2NO_2 = oriAlignPair_Nor2Rcm1[it->second].first;
				if((norAlignmentInfo_PE_2[tmpNor2NO]->endMatchedPosInChr)
					< (norAlignmentInfo_PE_2[tmpNor2NO_2]->endMatchedPosInChr))
				{
					mapPos_Nor2Rcm1.erase(it);
					mapPos_Nor2Rcm1.insert(pair<int, int> (norAlignmentInfo_PE_2[tmpNor2NO]->alignChromPos, tmp));
				}
			}
		}

		int currentBestRcm2NO = 0;
		int currentShortestPairDifference = 2147483647;
		int tmpNor1Rcm2NO;
		for(it = mapPos_Nor1Rcm2.begin(); it != mapPos_Nor1Rcm2.end(); it++)
		{
			currentShortestPairDifference = 2147483647;

			tmpNor1Rcm2NO = it->second;
			//cout << "tmpNor1Rcm2NO = " << tmpNor1Rcm2NO << endl;
			for(int tmp = 0; tmp < (oriAlignPair_Nor1Rcm2[tmpNor1Rcm2NO].second).size(); tmp++)
			{
				tmpRcm2NO = (oriAlignPair_Nor1Rcm2[tmpNor1Rcm2NO].second)[tmp];
				if(rcmAlignmentInfo_PE_2[tmpRcm2NO]->alignChromPos < currentShortestPairDifference ||
						(rcmAlignmentInfo_PE_2[tmpRcm2NO]->alignChromPos == currentShortestPairDifference &&
							rcmAlignmentInfo_PE_2[tmpRcm2NO]->endMatchedPosInChr
							< rcmAlignmentInfo_PE_2[currentBestRcm2NO]->endMatchedPosInChr))
					{
						currentBestRcm2NO = tmpRcm2NO;
						currentShortestPairDifference = rcmAlignmentInfo_PE_2[tmpRcm2NO]->alignChromPos;
					}
			}
			finalAlignPair_Nor1Rcm2.push_back( pair<int, int>
				(oriAlignPair_Nor1Rcm2[tmpNor1Rcm2NO].first, currentBestRcm2NO) );
		}

		int currentBestRcm1NO = 0;
		//currentShortestPairDifference = 2147483647;
		int tmpNor2Rcm1NO;
		for(it = mapPos_Nor2Rcm1.begin(); it != mapPos_Nor2Rcm1.end(); it++)
		{
			currentShortestPairDifference = 2147483647;
			tmpNor2Rcm1NO = it->second;
			for(int tmp = 0; tmp < (oriAlignPair_Nor2Rcm1[tmpNor2Rcm1NO].second).size(); tmp++)
			{
				tmpRcm1NO = (oriAlignPair_Nor2Rcm1[tmpNor2Rcm1NO].second)[tmp];
				if(rcmAlignmentInfo_PE_1[tmpRcm1NO]->alignChromPos < currentShortestPairDifference)
				{
					currentBestRcm1NO = tmpRcm1NO;
					//cout << "currentBestRcm1 = " << currentBestRcm1NO << endl;
					currentShortestPairDifference = rcmAlignmentInfo_PE_1[tmpRcm1NO]->alignChromPos;
				}
				else if(rcmAlignmentInfo_PE_1[tmpRcm1NO]->alignChromPos == currentShortestPairDifference)
				{
					if(rcmAlignmentInfo_PE_1[tmpRcm1NO]->endMatchedPosInChr
						< rcmAlignmentInfo_PE_1[currentBestRcm1NO]->endMatchedPosInChr)
					{
						currentBestRcm1NO = tmpRcm1NO;
						//cout << "currentBestRcm1 = " << currentBestRcm1NO << endl;
						currentShortestPairDifference = rcmAlignmentInfo_PE_1[tmpRcm1NO]->alignChromPos;
					}
					else
					{}
				}
				else
				{}
			}
			finalAlignPair_Nor2Rcm1.push_back( pair<int, int>
				(oriAlignPair_Nor2Rcm1[tmpNor2Rcm1NO].first, currentBestRcm1NO) );
		}

	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////// Output AlignInfo functions For Fasta Reads ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////// Output AlignInfo functions For Fasta Reads ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////


	string getTmpAlignInfo(PairedEndRead* pairedEndRead)
	{
		string readName_1 = pairedEndRead->getFirstRead()->getName();
		string readOriSeq_1 = pairedEndRead->getFirstRead()->getSequence();
		string readOriQualSeq_1 = pairedEndRead->getFirstRead()->getQuality();

		string readName_2 = pairedEndRead->getSecondRead()->getName();
		string readOriSeq_2 = pairedEndRead->getSecondRead()->getSequence();
		string readOriQualSeq_2 = pairedEndRead->getSecondRead()->getQuality();

		string tmpAlignInfoStr = "\n" + readName_1 + "\t"
			+ Utilities::int_to_str(norAlignmentInfo_PE_1.size()) + "\t"
			+ Utilities::int_to_str(rcmAlignmentInfo_PE_1.size()) + "\n"
			+ readOriSeq_1 + "\n" + readOriQualSeq_1 + "\n"
			+ readName_2 + "\t"
			+ Utilities::int_to_str(norAlignmentInfo_PE_2.size()) + "\t"
			+ Utilities::int_to_str(rcmAlignmentInfo_PE_2.size()) + "\n"
			+ readOriSeq_2 + "\n" + readOriQualSeq_2;

		Alignment_Info* tmpAlignInfo;
		string tmpNorAlignInfo_1;
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			if(norAlignmentInfo_PE_1.size() >= 40)
			{
				break;
			}
			tmpAlignInfo = norAlignmentInfo_PE_1[tmp];
			tmpNorAlignInfo_1 = tmpNorAlignInfo_1
				+ tmpAlignInfo->alignChromName + ","
				+ Utilities::int_to_str(tmpAlignInfo->alignChromPos) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ Utilities::int_to_str(tmpAlignInfo->mismatchNum) + ",\t";
		}
		string tmpRcmAlignInfo_1;
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
		{
			if(rcmAlignmentInfo_PE_1.size() >= 40)
			{
				break;
			}
			tmpAlignInfo = rcmAlignmentInfo_PE_1[tmp];
			tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1
				+ tmpAlignInfo->alignChromName + ","
				+ Utilities::int_to_str(tmpAlignInfo->alignChromPos) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ Utilities::int_to_str(tmpAlignInfo->mismatchNum) + ",\t";
		}
		string tmpNorAlignInfo_2;
		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			if(norAlignmentInfo_PE_2.size() >= 40)
			{
				break;
			}
			tmpAlignInfo = norAlignmentInfo_PE_2[tmp];
			tmpNorAlignInfo_2 = tmpNorAlignInfo_2
				+ tmpAlignInfo->alignChromName + ","
				+ Utilities::int_to_str(tmpAlignInfo->alignChromPos) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ Utilities::int_to_str(tmpAlignInfo->mismatchNum) + ",\t";
		}
		string tmpRcmAlignInfo_2;
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
		{
			if(rcmAlignmentInfo_PE_2.size() >= 40)
			{
				break;
			}
			tmpAlignInfo = rcmAlignmentInfo_PE_2[tmp];
			tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2
				+ tmpAlignInfo->alignChromName + ","
				+ Utilities::int_to_str(tmpAlignInfo->alignChromPos) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ Utilities::int_to_str(tmpAlignInfo->mismatchNum) + ",\t";
		}

		tmpAlignInfoStr = tmpAlignInfoStr + "\n"
			+ "Nor_1:\t" + tmpNorAlignInfo_1 + "\n"
			+ "Rcm_1:\t" + tmpRcmAlignInfo_1 + "\n"
			+ "Nor_2:\t" + tmpNorAlignInfo_2 + "\n"
			+ "Rcm_2:\t" + tmpRcmAlignInfo_2;
		return tmpAlignInfoStr;
	}

	string getTmpAlignInfoForFinalPair(PairedEndRead* pairedEndRead)
	{
		string readName_1 = pairedEndRead->getFirstRead()->getName();
		string readOriSeq_1 = pairedEndRead->getFirstRead()->getSequence();
		string readOriQualSeq_1 = pairedEndRead->getFirstRead()->getQuality();

		string readName_2 = pairedEndRead->getSecondRead()->getName();
		string readOriSeq_2 = pairedEndRead->getSecondRead()->getSequence();
		string readOriQualSeq_2 = pairedEndRead->getSecondRead()->getQuality();

		string tmpAlignInfoStr = "\n" + readName_1 + "\t"
			+ Utilities::int_to_str(finalAlignPair_Nor1Rcm2.size()) + "\t"
			+ Utilities::int_to_str(finalAlignPair_Nor2Rcm1.size()) + "\n"
			+ readOriSeq_1 + "\n" + readOriQualSeq_1 + "\n"
			+ readName_2 + "\t"
			+ Utilities::int_to_str(finalAlignPair_Nor2Rcm1.size()) + "\t"
			+ Utilities::int_to_str(finalAlignPair_Nor1Rcm2.size()) + "\n"
			+ readOriSeq_2 + "\n" + readOriQualSeq_2;

		Alignment_Info* tmpAlignInfo;

		string tmpNorAlignInfo_1;
		for(int tmpNor1Rcm2 = 0; tmpNor1Rcm2 < finalAlignPair_Nor1Rcm2.size(); tmpNor1Rcm2++)
		{
			if(finalAlignPair_Nor1Rcm2.size() >= 40)
				{break;}

			int tmp = finalAlignPair_Nor1Rcm2[tmpNor1Rcm2].first;
			tmpAlignInfo = norAlignmentInfo_PE_1[tmp];
			tmpNorAlignInfo_1 = tmpNorAlignInfo_1
				+ tmpAlignInfo->alignChromName + ","
				+ Utilities::int_to_str(tmpAlignInfo->alignChromPos) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ Utilities::int_to_str(tmpAlignInfo->mismatchNum) + ",\t";
		}
		string tmpRcmAlignInfo_1;
		for(int tmpNor2Rcm1 = 0; tmpNor2Rcm1 < finalAlignPair_Nor2Rcm1.size(); tmpNor2Rcm1++)
		{
			if(finalAlignPair_Nor2Rcm1.size() >= 40)
				{break;}
			int tmp = finalAlignPair_Nor2Rcm1[tmpNor2Rcm1].second;
			tmpAlignInfo = rcmAlignmentInfo_PE_1[tmp];
			tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1
				+ tmpAlignInfo->alignChromName + ","
				+ Utilities::int_to_str(tmpAlignInfo->alignChromPos) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ Utilities::int_to_str(tmpAlignInfo->mismatchNum) + ",\t";
		}
		string tmpNorAlignInfo_2;
		for(int tmpNor2Rcm1 = 0; tmpNor2Rcm1 < finalAlignPair_Nor2Rcm1.size(); tmpNor2Rcm1++)
		{
			if(finalAlignPair_Nor2Rcm1.size() >= 40)
				{break;}
			int tmp = finalAlignPair_Nor2Rcm1[tmpNor2Rcm1].first;
			tmpAlignInfo = norAlignmentInfo_PE_2[tmp];
			tmpNorAlignInfo_2 = tmpNorAlignInfo_2
				+ tmpAlignInfo->alignChromName + ","
				+ Utilities::int_to_str(tmpAlignInfo->alignChromPos) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ Utilities::int_to_str(tmpAlignInfo->mismatchNum) + ",\t";
		}
		string tmpRcmAlignInfo_2;
		for(int tmpNor1Rcm2 = 0; tmpNor1Rcm2 < finalAlignPair_Nor1Rcm2.size(); tmpNor1Rcm2++)
		{
			if(finalAlignPair_Nor1Rcm2.size() >= 40)
				{break;}
			int tmp = finalAlignPair_Nor1Rcm2[tmpNor1Rcm2].second;
			tmpAlignInfo = rcmAlignmentInfo_PE_2[tmp];
			tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2
				+ tmpAlignInfo->alignChromName + ","
				+ Utilities::int_to_str(tmpAlignInfo->alignChromPos) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ Utilities::int_to_str(tmpAlignInfo->mismatchNum) + ",\t";
		}
		tmpAlignInfoStr = tmpAlignInfoStr + "\n"
			+ "Nor_1:\t" + tmpNorAlignInfo_1 + "\n"
			+ "Rcm_1:\t" + tmpRcmAlignInfo_1 + "\n"
			+ "Nor_2:\t" + tmpNorAlignInfo_2 + "\n"
			+ "Rcm_2:\t" + tmpRcmAlignInfo_2;
		return tmpAlignInfoStr;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////// Output SAM functions For Fasta Reads ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////// Output SAM functions For Fasta Reads ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	string getSAMformatForFinalPair_secondaryOrNot(PairedEndRead* pairedEndRead)
	{
		string readName_1 = pairedEndRead->getFirstRead()->getName();
		string readSeq_1 = pairedEndRead->getFirstRead()->getSequence();

		string readName_2 = pairedEndRead->getSecondRead()->getName();
		string readSeq_2 = pairedEndRead->getSecondRead()->getSequence();

		string tmpSamStr;

		int IH_Nor1Rcm2 = finalAlignPair_Nor1Rcm2.size();
		int IH_Nor2Rcm1 = finalAlignPair_Nor2Rcm1.size();

		int IH_allPair = IH_Nor1Rcm2 + IH_Nor2Rcm1;

		for(int tmp = 0; tmp < finalAlignPair_Nor1Rcm2.size(); tmp++)
		{
			int HI_Nor1Rcm2_tmp = tmp + 1;
			int tmpNor1NO = finalAlignPair_Nor1Rcm2[tmp].first;
			int tmpRcm2NO = finalAlignPair_Nor1Rcm2[tmp].second;
			Alignment_Info* tmpAlignInfo_1 = norAlignmentInfo_PE_1[tmpNor1NO];
			Alignment_Info* tmpAlignInfo_2 = rcmAlignmentInfo_PE_2[tmpRcm2NO];

			int template_start = tmpAlignInfo_1->alignChromPos;
			int template_end = tmpAlignInfo_2->endMatchedPosInChr;

			int template_length = template_end - template_start + 1;

			tmpSamStr 
				= tmpSamStr + tmpAlignInfo_1->getSamFormatString_paired_secondaryOrNot(
					readName_1, readSeq_1, tmpAlignInfo_2, true, IH_allPair, HI_Nor1Rcm2_tmp, (tmp != 0), template_length);

			tmpSamStr += "\n";
			tmpSamStr 
				= tmpSamStr + tmpAlignInfo_2->getSamFormatString_paired_secondaryOrNot(
					readName_2, readSeq_2, tmpAlignInfo_1, false, IH_allPair, HI_Nor1Rcm2_tmp, (tmp != 0), template_length);
			tmpSamStr += "\n";
		}
		
		for(int tmp = 0; tmp < finalAlignPair_Nor2Rcm1.size(); tmp++)
		{
			int HI_Nor2Rcm1_tmp = tmp + 1
				+ IH_Nor1Rcm2;

			int tmpNor2NO = finalAlignPair_Nor2Rcm1[tmp].first;
			int tmpRcm1NO = finalAlignPair_Nor2Rcm1[tmp].second;
			Alignment_Info* tmpAlignInfo_1 = norAlignmentInfo_PE_2[tmpNor2NO];		
			Alignment_Info* tmpAlignInfo_2 = rcmAlignmentInfo_PE_1[tmpRcm1NO];
			
			int template_start = tmpAlignInfo_1->alignChromPos;
			int template_end = tmpAlignInfo_2->endMatchedPosInChr;

			int template_length = template_end - template_start + 1;		

			tmpSamStr 
				//= tmpAlignInfo_1 -> getSamFormatString(readName_2, readSeq_2);
				= tmpSamStr + tmpAlignInfo_1->getSamFormatString_paired_secondaryOrNot(
					readName_2, readSeq_2, tmpAlignInfo_2, false, IH_allPair, HI_Nor2Rcm1_tmp, ((tmp != 0)||(IH_Nor1Rcm2 > 0)), template_length);
			tmpSamStr += "\n";
			tmpSamStr 
				= tmpSamStr + tmpAlignInfo_2->getSamFormatString_paired_secondaryOrNot(
					readName_1, readSeq_1, tmpAlignInfo_1, true, IH_allPair, HI_Nor2Rcm1_tmp, ((tmp != 0)||(IH_Nor1Rcm2 > 0)), template_length);
			tmpSamStr += "\n";
		}
		string returnStr = tmpSamStr.substr(0, tmpSamStr.length()-1);
		return returnStr;
	}	

	string getSAMformatForUnpairedAlignments_secondaryOrNot(PairedEndRead* pairedEndRead)
	{
		string readName_1 = pairedEndRead->getFirstRead()->getName();
		string readSeq_1 = pairedEndRead->getFirstRead()->getSequence();
		string readOriQualSeq_1 = pairedEndRead->getFirstRead()->getQuality();

		string readName_2 = pairedEndRead->getSecondRead()->getName();
		string readSeq_2 = pairedEndRead->getSecondRead()->getSequence();
		string readOriQualSeq_2 = pairedEndRead->getSecondRead()->getQuality();

		string peAlignSamStr;

		int IH_Nor1 = norAlignmentInfo_PE_1.size();
		int IH_Rcm1 = rcmAlignmentInfo_PE_1.size();	

		int IH_1 = norAlignmentInfo_PE_1.size() + rcmAlignmentInfo_PE_1.size();	

		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); 
			tmp++)
		{
			int HI_Nor1_tmp = tmp + 1;
			string tmpSamStr = norAlignmentInfo_PE_1[tmp]->getSamFormatString_unpaired_secondaryOrNot(
				readName_1, readSeq_1, true, IH_1, HI_Nor1_tmp, (tmp!=0));
			//outputFile << tmpSamStr << endl;
			peAlignSamStr = peAlignSamStr + tmpSamStr + "\n";
		}

		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); 
			tmp++)
		{
			int HI_Rcm1_tmp = tmp + 1 + IH_Nor1;
			string tmpSamStr = rcmAlignmentInfo_PE_1[tmp]->getSamFormatString_unpaired_secondaryOrNot(
				readName_1, readSeq_1, true, IH_1, HI_Rcm1_tmp, ((tmp!=0)||(IH_Nor1 > 0)) );
			//outputFile << tmpSamStr << endl;
			peAlignSamStr = peAlignSamStr + tmpSamStr + "\n";
		}
		
		if((norAlignmentInfo_PE_1.size() + rcmAlignmentInfo_PE_1.size()) == 0)
		{
			//outputFile << readName_1 << "\t4\t*\t0\t255\t*\t*\t0\t0\t" << readSeq_1 << endl;
			peAlignSamStr = readName_1 + "\t69\t*\t0\t0\t*\t*\t0\t0\t" + readSeq_1 + "\t*\tIH:i:0\tHI:i:0\n";
		}

		int IH_Nor2 = norAlignmentInfo_PE_2.size();	
		int IH_Rcm2 = rcmAlignmentInfo_PE_2.size();

		int IH_2 = IH_Nor2 + IH_Rcm2;

		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); 
			tmp++)
		{
			int HI_Nor2_tmp = tmp + 1;
			string tmpSamStr = norAlignmentInfo_PE_2[tmp]->getSamFormatString_unpaired_secondaryOrNot(
				readName_2, readSeq_2, false, IH_2, HI_Nor2_tmp, (tmp!=0));
			//outputFile << tmpSamStr << endl;
			peAlignSamStr = peAlignSamStr + tmpSamStr + "\n";
		}


		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); 
			tmp++)
		{
			int HI_Rcm2_tmp = tmp + 1 + IH_Nor2;
			string tmpSamStr = rcmAlignmentInfo_PE_2[tmp]->getSamFormatString_unpaired_secondaryOrNot(
				readName_2, readSeq_2, false, IH_2, HI_Rcm2_tmp, ((tmp!=0)||(IH_Nor2>0)) );
			//outputFile << tmpSamStr << endl;
			peAlignSamStr = peAlignSamStr + tmpSamStr + "\n";
		}

		if((norAlignmentInfo_PE_2.size() + rcmAlignmentInfo_PE_2.size()) == 0)
		{
			//outputFile << readName_2 << "\t4\t*\t0\t255\t*\t*\t0\t0\t" << readSeq_2 << endl;
			peAlignSamStr = peAlignSamStr + readName_2 + "\t133\t*\t0\t0\t*\t*\t0\t0\t" + readSeq_2 + "\t*\tIH:i:0\tHI:i:0\n";
		}
		return peAlignSamStr.substr(0,peAlignSamStr.length()-1);
	}

	string getSAMformatForBothEndsUnmapped(PairedEndRead* pairedEndRead)
	{
		string readName_1 = pairedEndRead->getFirstRead()->getName();
		string readSeq_1 = pairedEndRead->getFirstRead()->getSequence();

		string readName_2 = pairedEndRead->getSecondRead()->getName();
		string readSeq_2 = pairedEndRead->getSecondRead()->getSequence();

		string peAlignSamStr;
		peAlignSamStr = readName_1 + "\t77\t*\t0\t0\t*\t*\t0\t0\t" + readSeq_1 + "\t*\tIH:i:0\tHI:i:0\n"
			+ readName_2 + "\t141\t*\t0\t0\t*\t*\t0\t0\t" + readSeq_2 + "\t*\tIH:i:0\tHI:i:0";
		return peAlignSamStr;
	}

};

#endif
