#ifndef __READ_INFO_H_INCLUDED__
#define __READ_INFO_H_INCLUDED__

#include <string>
#include <string.h>
#include <stdexcept>

using namespace std;

class Read
{
private:

	// Member variables
	string _name;
	string _sequence;
	string _quality;

	string reverseComplement(const char& character)
	{
		switch(character)
		{
			case 'A':
				return "T";
				break;
			case 'T':
				return "A";
				break;
			case 'G':
				return "C";
				break;
			case 'C':
				return "G";
				break;
			case 'N':
				return "N";
				break;
			default: // invalid character
				throw std::invalid_argument(string("Invalid character ") + character + string(" in read"));
				break;
		}
	}

	string getReverseComplementString()
	{
		string rc = "";
		for(int i = 0; i < length(); i++)
		{
			rc += reverseComplement(_sequence.at(length() - 1 - i));
		}
		return rc;
	}

	string getReverseComplementQuaility()
	{
		int length = _quality.size();
		string resultString = _quality.substr(length-1, 1);
		for (int i = 1; i < length; i++)
			resultString = resultString + _quality.substr(length - 1 - i, 1);

		return resultString;
	}

public:
	int length()
	{
		return _sequence.length();
	}

	/*
	 * Sets the name of this read
	 */
	void setName(string value)
	{
		_name = value;
	}

	/*
	 * Sets the value of this sequence
	 */
	void setSequence(string value)
	{
		_sequence = value;
	}

	/*
	 * Set quality score
	 */
	void setQualityScore(string value)
	{
		_quality = value;
	}

	string getName()
	{
		return _name;
	}

	/*
	 * Gets the nucleotide sequence of this read
	 */
	string getSequence()
	{
		return _sequence;
	}

	string getQuality()
	{
		return _quality;
	}

	/*
	 * returns the reverse complement of this read
	 */
	Read getReverseComplement()
	{
		Read reverseComplementRead;
		reverseComplementRead._name = _name;
		reverseComplementRead._sequence = getReverseComplementString();
		reverseComplementRead._quality = getReverseComplementQuaility();

		return reverseComplementRead;
	}

};

class PairedEndRead
{

public:

	Read firstPairedEndRead;
	Read secondPairedEndRead;

	PairedEndRead()
	{}
	/* FIX ME - DELETE LATER 5/28/14
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
			return firstPairedEndRead.length();
		else
			return secondPairedEndRead.length();
	}
	 */

	void setReadData(const string& readName1, const string& readName2,
		const string& readSeq1, const string& readSeq2)
	{
		firstPairedEndRead.setName(readName1);
		firstPairedEndRead.setSequence(readSeq1);
		firstPairedEndRead.setQualityScore("*");

		secondPairedEndRead.setName(readName2);
		secondPairedEndRead.setSequence(readSeq2);
		secondPairedEndRead.setQualityScore("*");
	}

	Read getFirstRead()
	{
		return firstPairedEndRead;
	}

	Read getSecondRead()
	{
		return secondPairedEndRead;
	}

	Read getFirstReadReverseComplement()
	{
		return firstPairedEndRead.getReverseComplement();
	}

	Read getSecondReadReverseComplement()
	{
		return secondPairedEndRead.getReverseComplement();
	}
};

#endif
