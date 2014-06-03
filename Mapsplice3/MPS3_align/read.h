#ifndef __READ_H_INCLUDED__
#define __READ_H_INCLUDED__

#include <string>
#include <string.h>
#include <stdexcept>

using namespace std;

class Read
{
private:

	/*
	 * Member variables
	 */
	string _name;
	string _sequence;
	string _quality;

	/*
	 * Private Methods
	 */
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
		for(int i = 0; i < getLength(); i++)
		{
			rc += reverseComplement(_sequence.at(getLength() - 1 - i));
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

	/*
	 * Constructors
	 */
	Read(string name, string sequence, string quality)
	{
		_name = name;
		_sequence = sequence;
		_quality = quality;
	}

	Read(Read* other)
	{
		_name = other->_name;
		_sequence = other->_sequence;
		_quality = other->_quality;
	}

	/*
	 * Public Methods
	 */
	unsigned int getLength()
	{
		return _sequence.length();
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
	Read* getReverseComplement()
	{
		return  new Read(_name, getReverseComplementString(), getReverseComplementQuaility());
	}

};


#endif
