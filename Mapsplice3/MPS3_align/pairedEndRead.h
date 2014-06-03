#ifndef __PAIRED_END_READ_H_INCLUDED__
#define __PAIRED_END_READ_H_INCLUDED__

#include <string>
#include <string.h>
#include <stdexcept>

#include "read.h"

using namespace std;

class PairedEndRead
{

public:

	Read* firstPairedEndRead;
	Read* secondPairedEndRead;

	PairedEndRead()
	{}

	~PairedEndRead()
	{
		delete firstPairedEndRead;
		delete secondPairedEndRead;
	}

	PairedEndRead(const string& readName1, const string& readName2,
			const string& readSeq1, const string& readSeq2)
	{
		setReadData(readName1, readName2, readSeq1, readSeq2);
	}

	void setReadData(const string& readName1, const string& readName2,
		const string& readSeq1, const string& readSeq2)
	{
		firstPairedEndRead = new Read(readName1, readSeq1, "*");
		secondPairedEndRead = new Read(readName2, readSeq2, "*");
	}

	Read* getFirstRead()
	{
		return firstPairedEndRead;
	}

	Read* getSecondRead()
	{
		return secondPairedEndRead;
	}

	Read* getFirstReadReverseComplement()
	{
		return firstPairedEndRead->getReverseComplement();
	}

	Read* getSecondReadReverseComplement()
	{
		return secondPairedEndRead->getReverseComplement();
	}
};

#endif
