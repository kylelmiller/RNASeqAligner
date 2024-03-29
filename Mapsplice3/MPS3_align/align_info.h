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
		SJchrName = chrName;
		SJdonerEnd = donerPos;
		SJacceptorStart = acceptorPos;
		this->getFlankString(indexInfo);
	}

	void getFlankString(Index_Info* indexInfo)
	{
		int chromNO = indexInfo->convertStringToInt(SJchrName);
		flankString = indexInfo->getChromSequence(chromNO).substr(SJdonerEnd, 2) + indexInfo->getChromSequence(chromNO).substr(SJacceptorStart-3, 2);
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
	vector< pair< int, SpliceJunction_Alignment > > spliceJunctionVec;
	int mismatchNum;
	string SJstrand;
	int mappedBaseNum;
	int endMatchedPosInChr;
	
	SamFormat currentSam;

	Alignment_Info()
	{}

	bool unfixedHeadExists()
	{
		return cigarStringJumpCode.size() != 0 && cigarStringJumpCode[0]._type == "S";
	}

	bool unfixedTailExists()
	{
		return cigarStringJumpCode.size() != 0 && cigarStringJumpCode.back()._type == "S";
	}

	bool unfixedHeadAndTailExists()
	{
		return unfixedHeadExists() && unfixedTailExists();
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
		int newAlignChromPos = alignChromPos - spliceJunctionDistance - cigarStringJumpCode[0]._length;
		
		//cout << "start to creat new jumpCode " << endl;
		vector<Jump_Code> newCigarStringJumpCode;
		//cout << "stop 0 " << endl;
		Jump_Code* tmpJumpCode = new Jump_Code(firstMatchLength, "M");
		//tmpJumpCode.len = firstMatchLength; tmpJumpCode.type = "M";
		//cout << "stop 1 " << endl;
		newCigarStringJumpCode.push_back(*tmpJumpCode);
		tmpJumpCode->_length = spliceJunctionDistance; tmpJumpCode->_type = "N";
		//cout << "stop 2 " << endl;
		newCigarStringJumpCode.push_back(*tmpJumpCode);
		tmpJumpCode->_length = cigarStringJumpCode[0]._length + 
			cigarStringJumpCode[1]._length - firstMatchLength; tmpJumpCode->_type = "M";
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
			= cigarStringJumpCode[jumpCodeSize-2]._length + cigarStringJumpCode[jumpCodeSize-1]._length - lastMatchLength; 
		Jump_Code* tmpJumpCode = new Jump_Code(penultimateMatchLength, "M");
		newCigarStringJumpCode.push_back(*tmpJumpCode);
		tmpJumpCode->_length = spliceJunctionDistance; tmpJumpCode->_type = "N";
		newCigarStringJumpCode.push_back(*tmpJumpCode);
		tmpJumpCode->_length = lastMatchLength; tmpJumpCode->_type = "M";
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

	// FIXME - KLM 6/4/14 THIS DOES NOT BELONG HERE
	// THIS CAN'T BE IN ALIGN_INFO
	void jumpCodeVec2spliceJunctionVec(Index_Info* indexInfo) 
	{
		int tmpPosInRead = 0; 
		int tmpDonerPosInChr = alignChromPos - 1;
		int tmpAcceptorPosInChr = 0;
		//cout << "start to get sjVec ..." << endl;
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			//cout << "JumpCode[" << tmp << "]: " << cigarStringJumpCode[tmp].len << " " << cigarStringJumpCode[tmp].type << endl;
			if(cigarStringJumpCode[tmp]._type == "S")
			{
				tmpPosInRead += cigarStringJumpCode[tmp]._length;
			}
			else if (cigarStringJumpCode[tmp]._type == "M")
			{
				tmpPosInRead += cigarStringJumpCode[tmp]._length;
				//cout << "tmpPosInRead: " << tmpPosInRead << endl;
				tmpDonerPosInChr += cigarStringJumpCode[tmp]._length;
				//cout << "tmpDonerPosInChr: " << tmpDonerPosInChr << endl;
			}
			else if (cigarStringJumpCode[tmp]._type == "N")
			{
				//cout << "tmpDonerPosInChr: " << tmpDonerPosInChr << endl;

				tmpAcceptorPosInChr = tmpDonerPosInChr + cigarStringJumpCode[tmp]._length + 1;
				//cout << "tmpAcceptorPosInChr: " << tmpAcceptorPosInChr << endl;
				SpliceJunction_Alignment* tmpSJ = 
					new SpliceJunction_Alignment(alignChromName, tmpDonerPosInChr, tmpAcceptorPosInChr, indexInfo);
				//cout << "tmpPosInRead: " << tmpPosInRead << endl;
				//tmpSJ
				spliceJunctionVec.push_back(pair <int, SpliceJunction_Alignment> (tmpPosInRead+1, (*tmpSJ)));
				tmpDonerPosInChr += cigarStringJumpCode[tmp]._length;
			}
			else if (cigarStringJumpCode[tmp]._type == "I")
			{
				tmpPosInRead += cigarStringJumpCode[tmp]._length;
				//tmpDonerPosInChr += cigarStringJumpCode[tmp].len;
			}
			else if (cigarStringJumpCode[tmp]._type == "D")
			{
				tmpDonerPosInChr += cigarStringJumpCode[tmp]._length;
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
			if ((cigarStringJumpCode[tmp]._type == "M")||(cigarStringJumpCode[tmp]._type == "I"))
			{
				pos += cigarStringJumpCode[tmp]._length;
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
			if ((cigarStringJumpCode[tmp]._type == "S")||(cigarStringJumpCode[tmp]._type == "I"))
			{
				continue;
			}
			else if ((cigarStringJumpCode[tmp]._type == "M")||(cigarStringJumpCode[tmp]._type == "N")||(cigarStringJumpCode[tmp]._type == "D"))
			{
				pos += cigarStringJumpCode[tmp]._length;
			}
			else
			{
				cout << "error, unexpected jumpCodeType: " << cigarStringJumpCode[tmp]._type << endl;
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
		//cout << "start to check overlapArea SJ ..." << endl;
		bool SJinOverlapAreaNoConflict = true;

		int overlapStartPos = secondAlignmentInfo->alignChromPos;
		int overlapEndPos = this->endMatchedPosInChr;

		int tmpDonerPosInChr, tmpAcceptorPosInChr;
		vector< pair<int, int> > norAlignSJvecInOverlapArea;
		vector< pair<int, int> > rcmAlignSJvecInOverlapArea;
		for(int tmp = 0; tmp < (this->spliceJunctionVec).size(); tmp++)
		{
			if((((this->spliceJunctionVec)[tmp].second).SJdonerEnd >= overlapStartPos)
				&&(((this->spliceJunctionVec)[tmp].second).SJacceptorStart <= overlapEndPos))
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
			if((((secondAlignmentInfo->spliceJunctionVec)[tmp].second).SJacceptorStart <= overlapEndPos)
				&&(((secondAlignmentInfo->spliceJunctionVec)[tmp].second).SJdonerEnd >= overlapStartPos))
			{
				tmpDonerPosInChr = ((secondAlignmentInfo->spliceJunctionVec)[tmp].second).SJdonerEnd;
				tmpAcceptorPosInChr = ((secondAlignmentInfo->spliceJunctionVec)[tmp].second).SJacceptorStart;
				rcmAlignSJvecInOverlapArea.push_back(pair<int,int> (tmpDonerPosInChr, tmpAcceptorPosInChr));
			}
			else
			{}
		} 
		//cout << "norAlignSJvecInOverlapArea.size(): " << norAlignSJvecInOverlapArea.size() << endl;
		//cout << "rcmAlignSJvecInOverlapArea.size(): " << rcmAlignSJvecInOverlapArea.size() << endl;

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
		//cout << "start to check overlap SJ" << endl;
		bool SJinOverlapAreaNoConflict = true;

		int overlapStartPos_1 = this->alignChromPos;
		int overlapStartPos_2 = secondAlignmentInfo->alignChromPos;
		int overlapEndPos_1 = this->endMatchedPosInChr;
		int overlapEndPos_2 = secondAlignmentInfo->endMatchedPosInChr;

		int overlapStartPos;
		int overlapEndPos;
		
		if(overlapStartPos_1 > overlapStartPos_2)
			overlapStartPos = overlapStartPos_1;
		else
			overlapStartPos = overlapStartPos_2;

		if(overlapEndPos_1 < overlapEndPos_2)
			overlapEndPos = overlapEndPos_1;
		else
			overlapEndPos = overlapEndPos_2;

		//cout << "overlapStartPos: " << overlapStartPos << endl;
		//cout << "overlapEndPos: " << overlapEndPos << endl;

		int tmpDonerPosInChr, tmpAcceptorPosInChr;
		vector< pair<int, int> > norAlignSJvecInOverlapArea;
		vector< pair<int, int> > rcmAlignSJvecInOverlapArea;
		for(int tmp = 0; tmp < (this->spliceJunctionVec).size(); tmp++)
		{
			if((((this->spliceJunctionVec)[tmp].second).SJdonerEnd >= overlapStartPos)
				&&(((this->spliceJunctionVec)[tmp].second).SJacceptorStart <= overlapEndPos))
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
			if((((secondAlignmentInfo->spliceJunctionVec)[tmp].second).SJacceptorStart <= overlapEndPos)
				&&(((secondAlignmentInfo->spliceJunctionVec)[tmp].second).SJdonerEnd >= overlapStartPos))
			{
				tmpDonerPosInChr = ((secondAlignmentInfo->spliceJunctionVec)[tmp].second).SJdonerEnd;
				tmpAcceptorPosInChr = ((secondAlignmentInfo->spliceJunctionVec)[tmp].second).SJacceptorStart;
				rcmAlignSJvecInOverlapArea.push_back(pair<int,int> (tmpDonerPosInChr, tmpAcceptorPosInChr));
			}
			else
			{}
		} 
		//cout << "norAlignSJvecInOverlapArea.size(): " << norAlignSJvecInOverlapArea.size() << endl;
		//cout << "rcmAlignSJvecInOverlapArea.size(): " << rcmAlignSJvecInOverlapArea.size() << endl;
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

	vector< pair<int, int> > oriAlignPair_Nor1Rcm2_new;
	vector< string > oriAlignPairNew_chrNameVec_Nor1Rcm2;
	vector< int > oriAlignPairNew_startMapPosVec_Nor1Rcm2;
	vector< int > oriAlignPairNew_endMapPosVec_Nor1Rcm2;
	vector< int > oriAlignPairNew_mappedLength_Nor1Rcm2;
	vector< int > oriAlignPairNew_mismatchNum_Nor1Rcm2;
	vector< int > oriAlignPairNew_pairDistance_Nor1Rcm2;
	vector< double > oriAlignPairNew_score_Nor1Rcm2;
	//vector< int > oriAlignPairNew_score_Nor1Rcm2;

	vector< pair<int, int> > oriAlignPair_Nor2Rcm1_new; 
	vector< string > oriAlignPairNew_chrNameVec_Nor2Rcm1;
	vector< int > oriAlignPairNew_startMapPosVec_Nor2Rcm1;
	vector< int > oriAlignPairNew_endMapPosVec_Nor2Rcm1;
	vector< int > oriAlignPairNew_mappedLength_Nor2Rcm1;
	vector< int > oriAlignPairNew_mismatchNum_Nor2Rcm1;
	vector< int > oriAlignPairNew_pairDistance_Nor2Rcm1;
	vector< double > oriAlignPairNew_score_Nor2Rcm1;
	//vector< int > oriAlignPairNew_score_Nor1Rcm2;

	//vector< pair< int, int > > oriAlignPair_Nor1Rcm2_filtered;
	//vector< pair< int, int > > oriAlignPair_Nor2Rcm1_filtered;	

	vector< vector< int > > oriAlignPairGroupedByRegion_Nor1Rcm2;
	vector< vector< int > > oriAlignPairGroupedByRegion_Nor2Rcm1;

	vector< pair< int, int > > finalAlignPair_Nor1Rcm2;
	vector< pair< int, int > > finalAlignPair_Nor2Rcm1;


	vector<bool> otherEndUnmappedBoolVec;

	void oriAlignPair2oriAlignPairNew()
	{
		for(int tmp1 = 0; tmp1 < oriAlignPair_Nor1Rcm2.size(); tmp1++)
		{
			int tmpPair_NO_1 =  (oriAlignPair_Nor1Rcm2[tmp1].first);
			int tmpSize = (oriAlignPair_Nor1Rcm2[tmp1].second).size();
			for(int tmp2 = 0; tmp2 < tmpSize; tmp2++)
			{
				int tmpPair_NO_2 = (oriAlignPair_Nor1Rcm2[tmp1].second)[tmp2];

				string tmpChrName = (norAlignmentInfo_PE_1[tmpPair_NO_1])->alignChromName;
				int tmpChrMapPosStart_nor = (norAlignmentInfo_PE_1[tmpPair_NO_1])->alignChromPos;
				int tmpChrMapPosEnd_nor = (norAlignmentInfo_PE_1[tmpPair_NO_1])->endMatchedPosInChr;
				int tmpMappedLength_nor = (norAlignmentInfo_PE_1[tmpPair_NO_1])->mappedLength();
				int tmpMismatchNum_nor = (norAlignmentInfo_PE_1[tmpPair_NO_1])->mismatchNum;

				int tmpChrMapPosStart_rcm = (rcmAlignmentInfo_PE_2[tmpPair_NO_2])->alignChromPos;
				int tmpChrMapPosEnd_rcm = (rcmAlignmentInfo_PE_2[tmpPair_NO_2])->endMatchedPosInChr;
				int tmpMappedLength_rcm = (rcmAlignmentInfo_PE_2[tmpPair_NO_2])->mappedLength();
				int tmpMismatchNum_rcm = (rcmAlignmentInfo_PE_2[tmpPair_NO_2])->mismatchNum;

				int tmpChrMapPosStart, tmpChrMapPosEnd;
				if(tmpChrMapPosStart_nor <= tmpChrMapPosStart_rcm)
				{
					tmpChrMapPosStart = tmpChrMapPosStart_nor;
				}
				else
				{
					tmpChrMapPosStart = tmpChrMapPosStart_rcm;
				}

				if(tmpChrMapPosEnd_nor <= tmpChrMapPosEnd_rcm)
				{
					tmpChrMapPosEnd = tmpChrMapPosEnd_rcm;
				}
				else
				{
					tmpChrMapPosEnd = tmpChrMapPosEnd_nor;
				}
				int tmpMappedLength = tmpMappedLength_nor + tmpMappedLength_rcm;
				int tmpMismatchNum = tmpMismatchNum_nor + tmpMismatchNum_rcm;
				int tmpPairDistance = tmpChrMapPosStart_rcm - tmpChrMapPosEnd_nor;

				oriAlignPair_Nor1Rcm2_new.push_back(pair<int, int> (tmpPair_NO_1, tmpPair_NO_2) );
				oriAlignPairNew_chrNameVec_Nor1Rcm2.push_back(tmpChrName);
				oriAlignPairNew_startMapPosVec_Nor1Rcm2.push_back(tmpChrMapPosStart);
				oriAlignPairNew_endMapPosVec_Nor1Rcm2.push_back(tmpChrMapPosEnd);
				oriAlignPairNew_mappedLength_Nor1Rcm2.push_back(tmpMappedLength);
				oriAlignPairNew_mismatchNum_Nor1Rcm2.push_back(tmpMismatchNum);
				oriAlignPairNew_pairDistance_Nor1Rcm2.push_back(tmpPairDistance);
			}
		}

		for(int tmp1 = 0; tmp1 < oriAlignPair_Nor2Rcm1.size(); tmp1++)
		{
			int tmpPair_NO_1 = (oriAlignPair_Nor2Rcm1[tmp1]).first;
			int tmpSize = (oriAlignPair_Nor2Rcm1[tmp1].second).size();
			for(int tmp2 = 0; tmp2 < tmpSize; tmp2++)
			{
				int tmpPair_NO_2 = (oriAlignPair_Nor2Rcm1[tmp1].second)[tmp2];

				string tmpChrName = (norAlignmentInfo_PE_2[tmpPair_NO_1])->alignChromName;
				int tmpChrMapPosStart_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->alignChromPos;
				int tmpChrMapPosEnd_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->endMatchedPosInChr;
				int tmpMappedLength_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->mappedLength();
				int tmpMismatchNum_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->mismatchNum;

				int tmpChrMapPosStart_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->alignChromPos;
				int tmpChrMapPosEnd_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->endMatchedPosInChr;
				int tmpMappedLength_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->mappedLength();
				int tmpMismatchNum_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->mismatchNum;

				int tmpChrMapPosStart, tmpChrMapPosEnd;
				if(tmpChrMapPosStart_nor <= tmpChrMapPosStart_rcm)
				{
					tmpChrMapPosStart = tmpChrMapPosStart_nor;
				}
				else
				{
					tmpChrMapPosStart = tmpChrMapPosStart_rcm;
				}

				if(tmpChrMapPosEnd_nor <= tmpChrMapPosEnd_rcm)
				{
					tmpChrMapPosEnd = tmpChrMapPosEnd_rcm;
				}
				else
				{
					tmpChrMapPosEnd = tmpChrMapPosEnd_nor;
				}
				int tmpMappedLength = tmpMappedLength_nor + tmpMappedLength_rcm;
				int tmpMismatchNum = tmpMismatchNum_nor + tmpMismatchNum_rcm;
				int tmpPairDistance = tmpChrMapPosStart_rcm - tmpChrMapPosEnd_nor;

				oriAlignPair_Nor2Rcm1_new.push_back(pair<int, int> (tmpPair_NO_1, tmpPair_NO_2) );
				oriAlignPairNew_chrNameVec_Nor2Rcm1.push_back(tmpChrName);
				oriAlignPairNew_startMapPosVec_Nor2Rcm1.push_back(tmpChrMapPosStart);
				oriAlignPairNew_endMapPosVec_Nor2Rcm1.push_back(tmpChrMapPosEnd);
				oriAlignPairNew_mappedLength_Nor2Rcm1.push_back(tmpMappedLength);
				oriAlignPairNew_mismatchNum_Nor2Rcm1.push_back(tmpMismatchNum);
				oriAlignPairNew_pairDistance_Nor2Rcm1.push_back(tmpPairDistance);
			}
		}
	}

	void oriAlignPairGroupedByRegion()
	{
		vector < string > oriPairRegionChrVec_Nor1Rcm2;
		vector< pair <int, int> > oriPairRegionPosVec_Nor1Rcm2;

		vector < string > oriPairRegionChrVec_Nor2Rcm1;
		vector< pair <int, int> > oriPairRegionPosVec_Nor2Rcm1;

		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp ++)
		{
			int tmpPairStartMapPos = oriAlignPairNew_startMapPosVec_Nor1Rcm2[tmp];
			int tmpPairEndMapPos = oriAlignPairNew_endMapPosVec_Nor1Rcm2[tmp];
			string tmpPairChrName = oriAlignPairNew_chrNameVec_Nor1Rcm2[tmp];

			int regionVecSize = oriPairRegionPosVec_Nor1Rcm2.size();
			bool overlapRegionFound = false;
			for(int tmp2 = 0; tmp2 < regionVecSize; tmp2++)
			{
				string tmpRegionChr = oriPairRegionChrVec_Nor1Rcm2[tmp2];
				int tmpRegionStartPos = oriPairRegionPosVec_Nor1Rcm2[tmp2].first;
				int tmpRegionEndPos = oriPairRegionPosVec_Nor1Rcm2[tmp2].second;
				if ( ( tmpPairChrName == tmpRegionChr ) &&
					( !( (tmpPairStartMapPos >= tmpRegionEndPos) || ( tmpPairEndMapPos <= tmpRegionStartPos ) ) ) )
				{
					overlapRegionFound = true;					
					if(tmpPairStartMapPos < tmpRegionStartPos)
						oriPairRegionPosVec_Nor1Rcm2[tmp2].first = tmpPairStartMapPos;
					if(tmpPairEndMapPos > tmpRegionEndPos)
						oriPairRegionPosVec_Nor1Rcm2[tmp2].second = tmpPairEndMapPos;
					oriAlignPairGroupedByRegion_Nor1Rcm2[tmp2].push_back(tmp);
					break;
				}
			}

			if(!overlapRegionFound)
			{
				vector<int> newIntVec;
				newIntVec.push_back(tmp);
				oriAlignPairGroupedByRegion_Nor1Rcm2.push_back(newIntVec);
				oriPairRegionChrVec_Nor1Rcm2.push_back(tmpPairChrName);
				oriPairRegionPosVec_Nor1Rcm2.push_back(pair<int,int> (tmpPairStartMapPos, tmpPairEndMapPos));
			}
		}

		for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp ++)
		{
			int tmpPairStartMapPos = oriAlignPairNew_startMapPosVec_Nor2Rcm1[tmp];
			int tmpPairEndMapPos = oriAlignPairNew_endMapPosVec_Nor2Rcm1[tmp];
			string tmpPairChrName = oriAlignPairNew_chrNameVec_Nor2Rcm1[tmp];

			int regionVecSize = oriPairRegionPosVec_Nor2Rcm1.size();
			bool overlapRegionFound = false;
			for(int tmp2 = 0; tmp2 < regionVecSize; tmp2++)
			{
				string tmpRegionChr = oriPairRegionChrVec_Nor2Rcm1[tmp2];
				int tmpRegionStartPos = oriPairRegionPosVec_Nor2Rcm1[tmp2].first;
				int tmpRegionEndPos = oriPairRegionPosVec_Nor2Rcm1[tmp2].second;
				if ( ( tmpPairChrName == tmpRegionChr ) &&
					( !( (tmpPairStartMapPos >= tmpRegionEndPos) || ( tmpPairEndMapPos <= tmpRegionStartPos ) ) ) )
				{
					overlapRegionFound = true;					
					if(tmpPairStartMapPos < tmpRegionStartPos)
						oriPairRegionPosVec_Nor2Rcm1[tmp2].first = tmpPairStartMapPos;
					if(tmpPairEndMapPos > tmpRegionEndPos)
						oriPairRegionPosVec_Nor2Rcm1[tmp2].second = tmpPairEndMapPos;
					oriAlignPairGroupedByRegion_Nor2Rcm1[tmp2].push_back(tmp);
					break;
				}
			}

			if(!overlapRegionFound)
			{
				vector<int> newIntVec;
				newIntVec.push_back(tmp);
				oriAlignPairGroupedByRegion_Nor2Rcm1.push_back(newIntVec);
				oriPairRegionChrVec_Nor2Rcm1.push_back(tmpPairChrName);
				oriPairRegionPosVec_Nor2Rcm1.push_back(pair<int,int> (tmpPairStartMapPos, tmpPairEndMapPos));
			}
		}
	}

	void getScoreForEachPair_mapLen_mis_peDis()
	{
		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp++)
		{
			int tmpMappedLength = oriAlignPairNew_mappedLength_Nor1Rcm2[tmp];
			int tmpMismatchNum = oriAlignPairNew_mismatchNum_Nor1Rcm2[tmp];
			int tmpPairDistance = oriAlignPairNew_pairDistance_Nor1Rcm2[tmp];
			double tmpPairDistance_score;
			if(tmpPairDistance >= 500)
				tmpPairDistance_score = 0.1;
			else if(tmpPairDistance <= 0)
				tmpPairDistance_score = 0;
			else
				tmpPairDistance_score = (double)tmpPairDistance/5000;
			
			double tmpScore = tmpMappedLength - tmpMismatchNum - tmpPairDistance_score;

			oriAlignPairNew_score_Nor1Rcm2.push_back(tmpScore);
		}

		for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp++)
		{
			int tmpMappedLength = oriAlignPairNew_mappedLength_Nor2Rcm1[tmp];
			int tmpMismatchNum = oriAlignPairNew_mismatchNum_Nor2Rcm1[tmp];
			int tmpPairDistance = oriAlignPairNew_pairDistance_Nor2Rcm1[tmp];
			double tmpPairDistance_score;
			if(tmpPairDistance > 500)
				tmpPairDistance_score = 0.1;
			else if(tmpPairDistance <= 0)
				tmpPairDistance_score = 0;
			else
				tmpPairDistance_score = (double)tmpPairDistance/5000;
			double tmpScore = tmpMappedLength - tmpMismatchNum - tmpPairDistance_score;

			oriAlignPairNew_score_Nor2Rcm1.push_back(tmpScore);
		}
	}

	void getScoreForEachPair_mapLen_mis()
	{
		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp++)
		{
			int tmpMappedLength = oriAlignPairNew_mappedLength_Nor1Rcm2[tmp];
			int tmpMismatchNum = oriAlignPairNew_mismatchNum_Nor1Rcm2[tmp];
			double tmpScore = tmpMappedLength - tmpMismatchNum;

			oriAlignPairNew_score_Nor1Rcm2.push_back(tmpScore);
		}

		for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp++)
		{
			int tmpMappedLength = oriAlignPairNew_mappedLength_Nor2Rcm1[tmp];
			int tmpMismatchNum = oriAlignPairNew_mismatchNum_Nor2Rcm1[tmp];
			double tmpScore = tmpMappedLength - tmpMismatchNum;

			oriAlignPairNew_score_Nor2Rcm1.push_back(tmpScore);
		}
	}

	void chooseBestAlignment_selectRandomOneIfMulti()
	{
		this->pairingAlignment2OriPair();
		this->oriAlignPair2oriAlignPairNew();
		this->getScoreForEachPair_mapLen_mis_peDis();

		int selectedBestAlignmentNO = 0;
		bool selectedBestAlignment_Nor1Rcm2 = true;
		double tmpBestScore = 0;

		if(oriAlignPair_Nor1Rcm2_new.size() + oriAlignPair_Nor2Rcm1_new.size() == 0)
			return;

		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp ++)
		{
			double tmpScore = oriAlignPairNew_score_Nor1Rcm2[tmp];
			if(tmpScore > tmpBestScore)
			{
				tmpBestScore = tmpScore;
				selectedBestAlignmentNO = tmp;
				selectedBestAlignment_Nor1Rcm2 = true;
			}
		}

		for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp ++)
		{
			double tmpScore = oriAlignPairNew_score_Nor2Rcm1[tmp];
			if(tmpScore > tmpBestScore)
			{
				tmpBestScore = tmpScore;
				selectedBestAlignmentNO = tmp;
				selectedBestAlignment_Nor1Rcm2 = false;
			}
		}

		if(selectedBestAlignment_Nor1Rcm2)
		{
			finalAlignPair_Nor1Rcm2.push_back(oriAlignPair_Nor1Rcm2_new[selectedBestAlignmentNO]);
		}
		else
		{
			finalAlignPair_Nor2Rcm1.push_back(oriAlignPair_Nor2Rcm1_new[selectedBestAlignmentNO]);
		}
	}

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
			if(norAlignmentInfo_PE_1[(finalAlignPair_Nor1Rcm2[i].first)]->unfixedHeadAndTailExists())
				return false;

			if(rcmAlignmentInfo_PE_2[(finalAlignPair_Nor1Rcm2[i].second)]->unfixedHeadAndTailExists())
				return false;
		}

		for(int i=0; i<finalAlignPair_Nor2Rcm1.size(); i++)
		{
			if(norAlignmentInfo_PE_2[(finalAlignPair_Nor2Rcm1[i].first)]->unfixedHeadAndTailExists())
				return false;

			if(rcmAlignmentInfo_PE_1[(finalAlignPair_Nor2Rcm1[i].second)]->unfixedHeadAndTailExists())
				return false;
		}

		return true;
	}

	bool allUnpairedAlignmentCompleted()
	{
		for(int i = 0; i < norAlignmentInfo_PE_1.size(); i++)
			if(norAlignmentInfo_PE_1[i]->unfixedHeadAndTailExists())
				return false;

		for(int i = 0; i < norAlignmentInfo_PE_2.size(); i++)
			if(norAlignmentInfo_PE_2[i]->unfixedHeadAndTailExists())
				return false;

		for(int i = 0; i < rcmAlignmentInfo_PE_1.size(); i++)
			if(rcmAlignmentInfo_PE_1[i]->unfixedHeadAndTailExists())
				return false;

		for(int i = 0; i < rcmAlignmentInfo_PE_2.size(); i++)
			if(rcmAlignmentInfo_PE_2[i]->unfixedHeadAndTailExists())
				return false;

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

				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("-", indexInfo->getChromName(mapChromNameInt),
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

				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->getChromName(mapChromNameInt),
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

				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("-", indexInfo->getChromName(mapChromNameInt),
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

				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->getChromName(mapChromNameInt),
					mapChromPosInt, (((newPathInfo->finalPathVec)[tmpPath]).second)->final_jump_code, 
					tmpMismatch, indexInfo);
				norAlignmentInfo_PE_1.push_back(tmpAlignmentInfo);
			}
		}

	}

	PE_Read_Alignment_Info(Path* pathInfo_nor1, Path* pathInfo_rcm1, 
			Path* pathInfo_nor2, Path* pathInfo_rcm2, Index_Info* indexInfo)
	{
		// The final paths were found in the gap calls
		// now we loop of the final path of each of the types of mappings
		// This shouldn't know about different types of mappings though
		// it should just know about segments and where they are mapped
		for(int tmpPath = 0; tmpPath < pathInfo_nor1->finalPathVec.size(); tmpPath++)
		{
			int mapChromNameInt = (((pathInfo_nor1->finalPathVec)[tmpPath]).first).first;
			int mapChromPosInt = (((pathInfo_nor1->finalPathVec)[tmpPath]).first).second;
			int tmpMismatch = (pathInfo_nor1->fixedPathMismatchVec)[tmpPath];

			Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->getChromName(mapChromNameInt),
				mapChromPosInt, (((pathInfo_nor1->finalPathVec)[tmpPath]).second)->final_jump_code, 
				tmpMismatch, indexInfo);
			norAlignmentInfo_PE_1.push_back(tmpAlignmentInfo);
		}

		for(int tmpPath = 0; tmpPath < pathInfo_rcm1->finalPathVec.size(); tmpPath++)
		{
			int mapChromNameInt = (((pathInfo_rcm1->finalPathVec)[tmpPath]).first).first;
			int mapChromPosInt = (((pathInfo_rcm1->finalPathVec)[tmpPath]).first).second;
			int tmpMismatch = (pathInfo_rcm1->fixedPathMismatchVec)[tmpPath];

			Alignment_Info* tmpAlignmentInfo = new Alignment_Info("-", indexInfo->getChromName(mapChromNameInt),
				mapChromPosInt, (((pathInfo_rcm1->finalPathVec)[tmpPath]).second)->final_jump_code, 
				tmpMismatch, indexInfo);
			rcmAlignmentInfo_PE_1.push_back(tmpAlignmentInfo);
		}

		for(int tmpPath = 0; tmpPath < pathInfo_nor2->finalPathVec.size(); tmpPath++)
		{
			int mapChromNameInt = (((pathInfo_nor2->finalPathVec)[tmpPath]).first).first;
			int mapChromPosInt = (((pathInfo_nor2->finalPathVec)[tmpPath]).first).second;
			int tmpMismatch = (pathInfo_nor2->fixedPathMismatchVec)[tmpPath];

			Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->getChromName(mapChromNameInt),
				mapChromPosInt, (((pathInfo_nor2->finalPathVec)[tmpPath]).second)->final_jump_code, 
				tmpMismatch, indexInfo);
			norAlignmentInfo_PE_2.push_back(tmpAlignmentInfo);
		}

		for(int tmpPath = 0; tmpPath < pathInfo_rcm2->finalPathVec.size(); tmpPath++)
		{
			int mapChromNameInt = (((pathInfo_rcm2->finalPathVec)[tmpPath]).first).first;
			int mapChromPosInt = (((pathInfo_rcm2->finalPathVec)[tmpPath]).first).second;
			int tmpMismatch = (pathInfo_rcm2->fixedPathMismatchVec)[tmpPath];

			Alignment_Info* tmpAlignmentInfo = new Alignment_Info("-", indexInfo->getChromName(mapChromNameInt),
				mapChromPosInt, (((pathInfo_rcm2->finalPathVec)[tmpPath]).second)->final_jump_code, 
				tmpMismatch, indexInfo);
			rcmAlignmentInfo_PE_2.push_back(tmpAlignmentInfo);
		}
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
		for(int i=0; i<norAlignmentInfo_PE_1.size(); i++)
			norAlignmentInfo_PE_1[i]->getEndMatchedPosInChr();

		for(int i=0; i<rcmAlignmentInfo_PE_1.size(); i++)
			rcmAlignmentInfo_PE_1[i]->getEndMatchedPosInChr();

		for(int i=0; i<norAlignmentInfo_PE_2.size(); i++)
			norAlignmentInfo_PE_2[i]->getEndMatchedPosInChr();

		for(int i=0; i< rcmAlignmentInfo_PE_2.size(); i++)
			rcmAlignmentInfo_PE_2[i]->getEndMatchedPosInChr();
	}

	void pairingAlignment2OriPair()
	{
		this->getEndMatchPosForEveryAlignment();

		Alignment_Info* tmpAlignInfo_1;
		Alignment_Info* tmpAlignInfo_2;
		bool newEntity = false;

		for(int i= 0; i < norAlignmentInfo_PE_1.size(); i++)
		{
			newEntity = true;
			tmpAlignInfo_1 = norAlignmentInfo_PE_1[i];

			if(tmpAlignInfo_1->SJstrand == "X")
				continue;

			for(int tmp2 = 0; tmp2 < rcmAlignmentInfo_PE_2.size(); tmp2++)
			{
				tmpAlignInfo_2 = rcmAlignmentInfo_PE_2[tmp2];
				
				if((tmpAlignInfo_2 -> SJstrand == "X")
					||((tmpAlignInfo_1 -> SJstrand == "+")&&(tmpAlignInfo_2 -> SJstrand == "-"))
					||((tmpAlignInfo_1 -> SJstrand == "-")&&(tmpAlignInfo_2 -> SJstrand == "+")))
					continue;

				if((tmpAlignInfo_1->alignChromName) == (tmpAlignInfo_2->alignChromName))
				{
					if(tmpAlignInfo_1->endMatchedPosInChr <= tmpAlignInfo_2->alignChromPos)
					{
						if( ((tmpAlignInfo_2->alignChromPos) - (tmpAlignInfo_1->endMatchedPosInChr))< PAIR_READ_DISTANCE_MAX)
						{
							if(newEntity) //i is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (i, newTmpVec));
								newEntity = false;
							}
							else //i has already been in oriAlignPair_Nor1Rcm2
							{
								(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
							}
						}
					}
					else if((tmpAlignInfo_1->alignChromPos <= tmpAlignInfo_2->alignChromPos)
							&&(tmpAlignInfo_1->endMatchedPosInChr <= tmpAlignInfo_2->endMatchedPosInChr))
					{
						//cout << "parts of read overlap ..." << endl;
						//cout << "i: " << i << " i 2: " << tmp2 << endl;
						//Note: In addition should check whether they cross the same SJs or not
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
						{
							//cout << "overlap correct!" << endl;
							if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (i, newTmpVec));
								newEntity = false;
							}
							else //i has already been in oriAlignPair_Nor1Rcm2
							{
								(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
							}
						}
					}
					else if(
						(tmpAlignInfo_1->alignChromPos <= tmpAlignInfo_2->alignChromPos)
							&&(tmpAlignInfo_1->endMatchedPosInChr > tmpAlignInfo_2->endMatchedPosInChr)
							&&(tmpAlignInfo_2->unfixedTailExists()) // pe_2 read has unfixed tail
							 )
					{
						//cout << "type 3" << endl;
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
						{
							if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (i, newTmpVec));
								newEntity = false;
							}
							else //i has already been in oriAlignPair_Nor1Rcm2
							{
								(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
							}
						}
					}
					else if (
						(tmpAlignInfo_1->alignChromPos > tmpAlignInfo_2->alignChromPos)
							&&(tmpAlignInfo_1->endMatchedPosInChr <= tmpAlignInfo_2->endMatchedPosInChr)
							&&(tmpAlignInfo_1->unfixedHeadExists()) // pe_1 read has unfixed head
							)
					{
						//cout << "type 4" << endl;
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
						{
							if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (i, newTmpVec));
								newEntity = false;
							}
							else //i has already been in oriAlignPair_Nor1Rcm2
							{
								(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
							}
						}
					}
					else if (
							(tmpAlignInfo_1->unfixedHeadExists()) // pe_1 read has unfixed head
							&&(tmpAlignInfo_2->unfixedTailExists()) // pe_2 read has unfixed tail
							&&(tmpAlignInfo_1->alignChromPos > tmpAlignInfo_2->alignChromPos)
							&&(tmpAlignInfo_1->endMatchedPosInChr > tmpAlignInfo_2->endMatchedPosInChr)
							//&&(tmpAlignInfo_1->endMatchedPosInChr <= tmpAlignInfo_2->endMatchedPosInChr)
							&&( (tmpAlignInfo_1->endMatchedPosInChr)-(tmpAlignInfo_2->alignChromPos) < PAIR_READ_DISTANCE_MAX)							
							)
					{
						//cout << "type 5" << endl;
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
						{
							//cout << "overlap correct" << endl;
							if(newEntity) //i is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (i, newTmpVec));
								newEntity = false;
							}
							else //i has already been in oriAlignPair_Nor1Rcm2
							{
								(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
							}
						}
					}
				}
			}
		}
		//cout << "start to pair Nor2Rcm1... " << endl;
		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			//cout << "tmp: " << tmp << endl;
			newEntity = true;
			tmpAlignInfo_1 = norAlignmentInfo_PE_2[tmp];
			
			if(tmpAlignInfo_1 -> SJstrand == "X")
			{
				//cout << "SJstrand_1 = X" << endl;
				continue;
			}

			for(int tmp2 = 0; tmp2 < rcmAlignmentInfo_PE_1.size(); tmp2++)
			{
				//cout << "tmp2: " << tmp2 << endl;
				tmpAlignInfo_2 = rcmAlignmentInfo_PE_1[tmp2];

				if((tmpAlignInfo_2 -> SJstrand == "X")
					||((tmpAlignInfo_1 -> SJstrand == "+")&&(tmpAlignInfo_2 -> SJstrand == "-"))
					||((tmpAlignInfo_1 -> SJstrand == "-")&&(tmpAlignInfo_2 -> SJstrand == "+")))
				{
					//cout << "SJstrand_2 = X or conflict strand: " << tmpAlignInfo_1->SJstrand << ", " << tmpAlignInfo_2->SJstrand << endl;
					continue;
				}					
				
				if((tmpAlignInfo_1->alignChromName) == (tmpAlignInfo_2->alignChromName))
				{
					if(tmpAlignInfo_1->endMatchedPosInChr <= tmpAlignInfo_2->alignChromPos)
					{
						//cout << "type 1" << endl;
						if(((tmpAlignInfo_2->alignChromPos) - (tmpAlignInfo_1->endMatchedPosInChr)) < PAIR_READ_DISTANCE_MAX)
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
					else if((tmpAlignInfo_1->alignChromPos <= tmpAlignInfo_2->alignChromPos)
							&&(tmpAlignInfo_1->endMatchedPosInChr <= tmpAlignInfo_2->endMatchedPosInChr))
					{
						//cout << "parts of read overlap ..." << endl;
						//cout << "tmp: " << tmp << " tmp 2: " << tmp2 << endl;
						//Note: In addition should check whether they cross the same SJs or not
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
						{
							//cout << "overlap correct " << endl;
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
						{
							//cout << "overlap error !" << endl;
						}
					}
					else if(
						(tmpAlignInfo_1->alignChromPos <= tmpAlignInfo_2->alignChromPos)
							&&(tmpAlignInfo_1->endMatchedPosInChr > tmpAlignInfo_2->endMatchedPosInChr)
							&&(tmpAlignInfo_2->unfixedTailExists()) // pe_2 read has unfixed tail
							)
					{
						//cout << "type 3" << endl;
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
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
					else if (
						(tmpAlignInfo_1->alignChromPos > tmpAlignInfo_2->alignChromPos)
							&&(tmpAlignInfo_1->endMatchedPosInChr <= tmpAlignInfo_2->endMatchedPosInChr)
							&&(tmpAlignInfo_1->unfixedHeadExists()) // pe_1 read has unfixed head
							)
					{
						//cout << "type 4" << endl;
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
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
					else if (
							(tmpAlignInfo_1->unfixedHeadExists()) // pe_1 read has unfixed head
							&&(tmpAlignInfo_2->unfixedTailExists()) // pe_2 read has unfixed tail
							&&(tmpAlignInfo_1->alignChromPos > tmpAlignInfo_2->alignChromPos)
							&&(tmpAlignInfo_1->endMatchedPosInChr > tmpAlignInfo_2->endMatchedPosInChr)
							&&( (tmpAlignInfo_1->endMatchedPosInChr)-(tmpAlignInfo_2->alignChromPos) < PAIR_READ_DISTANCE_MAX)							
							)
					{
						//cout << "type 5" << endl;
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
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
					}
				}
			}
		}
		///////////////// 1. only one end read mapped, another one unmapped ///////////////////////

		///////////////// 2. reads can be mapped under both directions //////////////////////////

		///////////////// 3. reads can be mapped to different places in one direction //////////////////
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
