#ifndef __CHROMOSOME_H_INCLUDED__
#define __CHROMOSOME_H_INCLUDED__

#include <string>

using namespace std;

class Chromosome
{
private:
	string _name;
	string _sequence;
	unsigned int _endPosition;
	int _partNumber;
public:

	Chromosome(string name, int endPosition, int partNumber)
	{
		_name = name;
		_sequence = "";
		_endPosition = endPosition;
		_partNumber = partNumber;
	}

	void setSequence(string value)
	{
		_sequence = value;
	}

	string getName()
	{
		return _name;
	}

	int getLength()
	{
		return _sequence.length();
	}

	unsigned int getEndPosition()
	{
		return _endPosition;
	}

	int getPartNumber()
	{
		return _partNumber;
	}
};

#endif
