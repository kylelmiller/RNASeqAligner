#ifndef __SEGMENT_H_INCLUDED__
#define __SEGMENT_H_INCLUDED__

/*
 * Used to store mapped segments of the read
 */
class Segment
{
private:

	static const int LONG_SEGMENT_THRESHOLD = 20;
	unsigned int _length;
	unsigned int _locationInRead;
	unsigned int _alignmentNumber;
	unsigned int _alignmentLocation[CANDALILOC];
public:
	Segment(unsigned int length, unsigned int locationInRead, unsigned int alignmentNumber)
	{
		_length = length;
		_locationInRead = locationInRead;
		_alignmentNumber = alignmentNumber;
	}

	bool isLong()
	{
		return getLength() >= LONG_SEGMENT_THRESHOLD;
	}

	/*
	 * Gets the length of the mapped segment
	 */
	unsigned int getLength()
	{
		return _length;
	}

	/*
	 * Gets the alignment location
	 */
	unsigned int getAlignmentLocation(unsigned int index)
	{
		return _alignmentLocation[index];
	}

	/*
	 * Gets the alignment number of this segment
	 */
	unsigned int getAlignmentNumber()
	{
		return _alignmentNumber;
	}

	/*
	 * Sets the location that the segment is mapped in the read
	 */
	unsigned int getLocationInRead()
	{
		return _locationInRead;
	}

	/*
	 * Sets the length of the mapped segment
	 */
	void setLength(unsigned int value)
	{
		_length = value;
	}

	/*
	 * Sets the location that the segment is mapped in the read
	 */
	void setLocationInRead(unsigned int value)
	{
		_locationInRead = value;
	}

	/*
	 * Sets the alignment location value
	 */
	void setAlignmentLocation(unsigned int value, unsigned int index)
	{
		_alignmentLocation[index] = value;
	}
};

#endif
