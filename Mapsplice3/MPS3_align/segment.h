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
	unsigned int _locationInReference;

public:
	Segment(unsigned int length, unsigned int locationInRead, unsigned int locationInReference)
	{
		_length = length;
		_locationInRead = locationInRead;
		_locationInReference = locationInReference;
	}

	bool isConfident()
	{
		return getLength() >= LONG_SEGMENT_THRESHOLD;
	}

	int getConfidenceValue()
	{
		return _length;
	}

	/*
	 * Gets the length of the mapped segment
	 */
	unsigned int getLength()
	{
		return _length;
	}

	/*
	 * Gets the location in the reference this segment
	 * is mapped to
	 */
	unsigned int getLocationInReference()
	{
		return _locationInReference;
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

	int distanceBetweenSegment(Segment* other)
	{
		unsigned int segmentDistance = getLocationInReference() >= other->getLocationInReference()
			? other->getLocationInReference() - getLocationInReference()
			: getLocationInReference() - other->getLocationInReference();

		return segmentDistance < 300000 ? segmentDistance : 1000000;
	}
};

#endif
