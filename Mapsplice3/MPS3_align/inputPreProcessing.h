#include <stdlib.h>
#include <string>
#include <string.h>
#include <map>

using namespace std;

class InputPreProcess
{
public:
	string validBaseCharStr;

	InputPreProcess()
	{
		validBaseCharStr = "ATGC";
	}

	string trimReadSeq(const string& inputReadSeq)
	{
		int validReadSeqStart;
		int validReadSeqEnd;
		
		int inputReadSeqLength = inputReadSeq.length();
		
		validReadSeqStart = inputReadSeq.find_first_of(validBaseCharStr);
		validReadSeqEnd = inputReadSeq.find_last_of(validBaseCharStr);
		return inputReadSeq.substr(validReadSeqStart, validReadSeqEnd - validReadSeqStart + 1);
	}
};