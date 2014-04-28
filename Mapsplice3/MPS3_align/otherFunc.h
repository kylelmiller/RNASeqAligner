#include <string>
#include <string.h>

using namespace std;

char reverseComplement(char ch)
{
	if(ch == 'A')
		return 'T';
	else if(ch == 'C')
		return 'G';
	else if(ch == 'G')
		return 'C';
	else if(ch == 'T')
		return 'A';
	else 
		return 'N';
}

string reverseComplementStr(char ch)
{
	if(ch == 'A')
		return "T";
	else if(ch == 'C')
		return "G";
	else if(ch == 'G')
		return "C";
	else if(ch == 'T')
		return "A";
	else 
		return "N";	
}

string convertCharArrayToReverseCompletmentStr(char* readChar, int readLength)
{
	string rcmReadStr = "";
	for(int tmp = 0; tmp < readLength; tmp++)
	{
		rcmReadStr += reverseComplementStr(*(readChar + readLength - 1 - tmp));
	}
	return rcmReadStr;
}

string covertCharToReverseComplement(const string& Ori_Char)
{
	if(Ori_Char == "A")
	{
		return "T";
	}
	else if(Ori_Char == "T")
	{
		return "A";
	}
	else if(Ori_Char == "G")
	{
		return "C";
	}
	else if(Ori_Char == "C")
	{
		return "G";
	}
	else if(Ori_Char == "N")
	{
		return "N";
	}	
}

string covertStringToReverseComplement(const string& originalString)
{
	int stringLength = originalString.size();
	string resultString = covertCharToReverseComplement(originalString.substr(stringLength-1, 1));
	for (int tmp = 1; tmp < stringLength; tmp++)
	{
		resultString = resultString + covertCharToReverseComplement(
			originalString.substr(stringLength-1-tmp, 1));
	}
	return resultString;
}

string convertStringToReverseComplement(const string& originalString)
{
	int stringLength = originalString.size();
	string resultString = covertCharToReverseComplement(originalString.substr(stringLength-1, 1));
	for (int tmp = 1; tmp < stringLength; tmp++)
	{
		resultString = resultString + covertCharToReverseComplement(
			originalString.substr(stringLength-1-tmp, 1));
	}
	return resultString;
}

string convertQualityScoreString2Reverse(const string& originalQualityScoreString)
{
	int stringLength = originalQualityScoreString.size();
	string resultString = originalQualityScoreString.substr(stringLength-1, 1);//covertCharToReverseComplement(originalString.substr(stringLength-1, 1));
	for (int tmp = 1; tmp < stringLength; tmp++)
	{
		resultString = resultString + originalQualityScoreString.substr(stringLength-1-tmp, 1);
			//covertCharToReverseComplement(originalString.substr(stringLength-1-tmp, 1));
	}
	return resultString;
}

int extendBackInChromSeq(int readLoc, const string& readSeq, int chromLoc, const string& chromSeq, int extendBackLengthMax)
{
	int tmp = 1;
	for (tmp = 1; tmp <= extendBackLengthMax; tmp++)
	{
		if(readSeq.at(readLoc - tmp - 1) != chromSeq.at(chromLoc - tmp - 1))
		{
			return tmp - 1;
		}
	}
	return extendBackLengthMax;
}

class OtherFunc
{
public:
	vector<char> char2CharRcmVec;
	//vector<string> str2StrRcmVec;


	OtherFunc()
	{
		for(int tmp = 'A'-'A'; tmp <= 'Z' - 'A'; tmp++)
		{
			char2CharRcmVec.push_back('N');
		}
		char2CharRcmVec['A'-'A'] = 'T';
		char2CharRcmVec['C'-'A'] = 'G';
		char2CharRcmVec['G'-'A'] = 'C';
		char2CharRcmVec['T'-'A'] = 'A';

	}

	/*string stringRcm(const string& s)
	{
		string t;
		for(string::reverse_iterator iter = s.rbegin();
			iter != s.rend(); iter++)
		{
			t = t + char2CharRcmVec(*iter);
		}
		return t;		
	}*/

};