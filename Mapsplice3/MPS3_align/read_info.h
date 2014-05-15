#include <string>
#include <string.h>

using namespace std;

class Ori_Read_Info
{
public:
	string readName;
	int readSeqLength;
	string readSeq;
	string rcmReadSeq;
	string readQual;
	string rcmReadQual;
	//string readPeInfo;

	bool getSeqLength()
	{
		//return readSeq.length();
		readSeqLength = readSeq.length();
		if(readSeqLength <= 0)
		{
			return false;
		}	
		else
		{
			return true;
		}
	}
	/*
	bool getRcmRead()
	{
		for(int tmp = 0; )
	}
	*/

};

class PE_Read_Info
{
public:
	Ori_Read_Info readInfo_pe1;
	Ori_Read_Info readInfo_pe2; 

	PE_Read_Info()
	{}

	bool checkEnd1OrEnd2WithAlignInfoTypeNo(int alignInfoType)
	{
		if(alignInfoType <= 2)
			return true;
		else
			return false;
	}

	bool checkNorOrRcmWithAlignInfoTypeNo(int alignInfoType)
	{
		if((alignInfoType == 1)||(alignInfoType == 3))
			return true;
		else
			return false;
	}

	int checkReadLengthWithAlignInfoTypeNo(int alignInfoType)
	{
		if(alignInfoType <= 2)
			return (readInfo_pe1.readSeqLength);
		else
			return (readInfo_pe2.readSeqLength);
	}

	void get_PE_Read_Info(const string& readName1, const string& readName2,
		const string& readSeq1, const string& readSeq2)
	{
		(readInfo_pe1.readName) = readName1;
		//readInfo_pe1.readSeqLength = readSeq1.length();
		(readInfo_pe1.readSeq) = readSeq1;
		//readInfo_pe1.rcmReadSeq = rcmReadSeq1;
		//readInfo_pe1.readPeInfo = "1";

		(readInfo_pe2.readName) = readName2;
		//readInfo_pe2.readSeqLength = readSeq2.length();
		(readInfo_pe2.readSeq) = readSeq2;
		//readInfo_pe2.rcmReadSeq = rcmReadSeq2;
		//readInfo_pe2.readPeInfo = "2";
	}

	PE_Read_Info(const string& readName1, const string& readName2,
		const string& readSeq1, const string& readSeq2)
	{
		(readInfo_pe1.readName) = readName1;
		//readInfo_pe1.readSeqLength = readSeq1.length();
		(readInfo_pe1.readSeq) = readSeq1;
		//readInfo_pe1.rcmReadSeq = rcmReadSeq1;
		//readInfo_pe1.readPeInfo = "1";

		(readInfo_pe2.readName) = readName2;
		//readInfo_pe2.readSeqLength = readSeq2.length();
		(readInfo_pe2.readSeq) = readSeq2;
		//readInfo_pe2.rcmReadSeq = rcmReadSeq2;
		//readInfo_pe2.readPeInfo = "2";
	}

	/*void getPeReadInfo_Fasta(const string& readName_1, const string& readName_2, 
		const string& readSeq_1, const string& readSeq_2)
	{
		readInfo_pe1.readName = readName_1;
		readInfo_pe2.readName = readName_2;
		readInfo_pe1.readSeq = readSeq_1;
		readInfo_pe2.readSeq = readSeq_2;
		//readInfo_pe1.rcmReadSeq = convertStringToReverseComplement(readSeq_1);
		//readInfo_pe2.rcmReadSeq = convertStringToReverseComplement(readSeq_2);
		readInfo_pe1.readQual = "*";
		readInfo_pe2.readQual = "*";
		readInfo_pe1.rcmReadQual = "*";
		readInfo_pe2.rcmReadQual = "*";

		readInfo_pe1.readSeqLength = (readInfo_pe1.readSeq).length();
		readInfo_pe2.readSeqLength = (readInfo_pe2.readSeq).length();		
	}*/

	void getBothEndRcmReadSeq()
	{
		readInfo_pe1.rcmReadSeq = covertStringToReverseComplement(readInfo_pe1.readSeq);//, (readInfo_pe1.readSeq).length());
		readInfo_pe2.rcmReadSeq = covertStringToReverseComplement(readInfo_pe2.readSeq);//, (readInfo_pe2.readSeq).length());
	}

	void getFastaFormatReadInfo(const string& readName_1, const string& readName_2, 
		const string& readSeq_1, const string& readSeq_1_rcm, 
		const string& readSeq_2, const string& readSeq_2_rcm)
	{
		readInfo_pe1.readName = readName_1;
		readInfo_pe2.readName = readName_2;
		readInfo_pe1.readSeq = readSeq_1;
		readInfo_pe2.readSeq = readSeq_2;
		readInfo_pe1.rcmReadSeq = readSeq_1_rcm;
		readInfo_pe2.rcmReadSeq = readSeq_2_rcm;
		readInfo_pe1.readQual = "*";
		readInfo_pe2.readQual = "*";
		readInfo_pe1.rcmReadQual = "*";
		readInfo_pe2.rcmReadQual = "*";

		readInfo_pe1.readSeqLength = (readInfo_pe1.readSeq).length();
		readInfo_pe2.readSeqLength = (readInfo_pe2.readSeq).length();
		//int readSeqLength;
		//string readSeq;
		//string rcmReadSeq;
		//string readQual;
		//string rcmReadQual;
	}

	void getFastaFormatReadInfo(const string& readName_1, const string& readName_2, 
		const string& readSeq_1, const string& readSeq_2)
	{
		readInfo_pe1.readName = readName_1;
		readInfo_pe2.readName = readName_2;
		readInfo_pe1.readSeq = readSeq_1;
		readInfo_pe2.readSeq = readSeq_2;
		//readInfo_pe1.rcmReadSeq = convertStringToReverseComplement(readSeq_1);
		//readInfo_pe2.rcmReadSeq = convertStringToReverseComplement(readSeq_2);
		readInfo_pe1.readQual = "*";
		readInfo_pe2.readQual = "*";
		readInfo_pe1.rcmReadQual = "*";
		readInfo_pe2.rcmReadQual = "*";

		readInfo_pe1.readSeqLength = (readInfo_pe1.readSeq).length();
		readInfo_pe2.readSeqLength = (readInfo_pe2.readSeq).length();		
	}

	void getReverseComplementReadSeq(char* readChar, char* readChar_PE)
	{
		//readInfo_pe1.rcmReadSeq = convertStringToReverseComplement(readInfo_pe1.readSeq);
		//readInfo_pe2.rcmReadSeq = convertStringToReverseComplement(readInfo_pe2.readSeq);		
		readInfo_pe1.rcmReadSeq = 
			convertCharArrayToReverseCompletmentStr(readChar, readInfo_pe1.readSeqLength);
		readInfo_pe2.rcmReadSeq = 
			convertCharArrayToReverseCompletmentStr(readChar_PE, readInfo_pe2.readSeqLength);
	}

	void getRcmReadSeq()
	{
		readInfo_pe1.rcmReadSeq = revcomp(readInfo_pe1.readSeq);
		readInfo_pe2.rcmReadSeq = revcomp(readInfo_pe2.readSeq);
	}

	void getFastqFormatReadInfo(const string& readName_1, const string& readName_2, 
		const string& readSeq_1, const string& readSeq_2, 
		const string& readQualSeq_1, const string& readQualSeq_2)
	{
		readInfo_pe1.readName = readName_1;
		readInfo_pe2.readName = readName_2;
		readInfo_pe1.readSeq = readSeq_1;
		readInfo_pe2.readSeq = readSeq_2;
		readInfo_pe1.rcmReadSeq = convertStringToReverseComplement(readQualSeq_1);
		readInfo_pe2.rcmReadSeq = convertStringToReverseComplement(readQualSeq_2);
		readInfo_pe1.readQual = readQualSeq_1; //convertQualityScoreString2Reverse(readQualSeq_1);
		readInfo_pe2.readQual = readQualSeq_2; //convertQualityScoreString2Reverse(readQualSeq_2);
		readInfo_pe1.rcmReadQual = convertQualityScoreString2Reverse(readQualSeq_1);
		readInfo_pe2.rcmReadQual = convertQualityScoreString2Reverse(readQualSeq_2);

		readInfo_pe1.readSeqLength = (readInfo_pe1.readSeq).length();
		readInfo_pe2.readSeqLength = (readInfo_pe2.readSeq).length();		
	}

	PE_Read_Info(const string& readName1, const string& readName2,
		const string& readSeq1, const string& rcmReadSeq1, 
		const string& readSeq2, const string& rcmReadSeq2)
	{
		readInfo_pe1.readName = readName1;
		readInfo_pe1.readSeqLength = readSeq1.length();
		readInfo_pe1.readSeq = readSeq1;
		readInfo_pe1.rcmReadSeq = rcmReadSeq1;
		//readInfo_pe1.readPeInfo = "1";

		readInfo_pe2.readName = readName2;
		readInfo_pe2.readSeqLength = readSeq2.length();
		readInfo_pe2.readSeq = readSeq2;
		readInfo_pe2.rcmReadSeq = rcmReadSeq2;
		//readInfo_pe2.readPeInfo = "2";

	}

	string getIncompleteEndReadSeq(bool End1OrEnd2, bool NorOrRcm)
	{
		if(End1OrEnd2 && NorOrRcm)
		{
			return readInfo_pe2.rcmReadSeq;
		}
		else if(End1OrEnd2 && (!NorOrRcm))
		{
			return readInfo_pe2.readSeq;
		}
		else if((!End1OrEnd2) && NorOrRcm)
		{
			return readInfo_pe1.rcmReadSeq;
		}
		else
		{
			return readInfo_pe1.readSeq;
		}
	}

	string getReadSeq(bool End1OrEnd2, bool NorOrRcm)
	{
		if(End1OrEnd2 && NorOrRcm)
		{
			return readInfo_pe1.readSeq;
		}
		else if(End1OrEnd2 && (!NorOrRcm))
		{
			return readInfo_pe1.rcmReadSeq;
		}
		else if((!End1OrEnd2) && NorOrRcm)
		{
			return readInfo_pe2.readSeq;
		}
		else
		{
			return readInfo_pe2.rcmReadSeq;
		}
	}
	void printPEreadInfo()
	{
		cout << "end 1 read: " <<  readInfo_pe1.readName << endl;
		cout << "length: " << readInfo_pe1.readSeqLength << endl;
		cout << readInfo_pe1.readSeq << endl;

		cout << "end 2 read: " <<  readInfo_pe2.readName << endl;
		cout << "length: " << readInfo_pe2.readSeqLength << endl;
		cout << readInfo_pe2.readSeq << endl;		
	}

};



