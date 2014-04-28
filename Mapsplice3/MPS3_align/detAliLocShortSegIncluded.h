//#define SEGMENTNUM 20
//#define POSSIBLE_MAP_CASES_MAX 40
//#define MAX_SPLICE_LENGTH 300000

bool checkShortSegWithCandLoc(unsigned int shortSegMapPos, int oldMapLabel, unsigned int* oldSegMapPos,
	Index_Info* indexInfo)
{
	bool checkShortSegWithOldMapLabel = false;

	unsigned int shortSegMapPosInChr, shortSegMapChrNameInt;
	indexInfo->getChrLocation(shortSegMapPos, &shortSegMapChrNameInt, &shortSegMapPosInChr);

	for(int tmpOldMapLabel = 0; tmpOldMapLabel < oldMapLabel; tmpOldMapLabel++)
	{
		unsigned int tmpOldLabelMapPos = oldSegMapPos[tmpOldMapLabel];
		unsigned int oldSegMapPosInChr, oldSegMapChrNameInt;
		indexInfo->getChrLocation(tmpOldLabelMapPos, &oldSegMapChrNameInt, &oldSegMapPosInChr);
		if(shortSegMapChrNameInt != oldSegMapChrNameInt) 	
			continue;
		int mapPosGap = (int)shortSegMapPosInChr - (int)oldSegMapPosInChr;
		if( (mapPosGap < MAX_SPLICE_LENGTH) && (mapPosGap > (-1)*MAX_SPLICE_LENGTH) )
		{
			return true;
		}
	}
	return false;
}

int detAliLocShortSegIncluded(unsigned int segmentNum, unsigned int* segmentLength, unsigned int* segmentLocInRead,
	unsigned int* segmentAlignNum, unsigned int* segmentAlignLoc, unsigned int oldMapLabel, 
	//unsigned int* oldSegMapRangeStart, unsigned int* oldSegMapRangeEnd, 
	unsigned int* oldSegMapPos, 
	//unsigned int* mapLabel, 
	unsigned int* segMapRangeStart, unsigned int* segMapRangeEnd, 
	unsigned int* segMapLoc,
	Index_Info* indexInfo
	)
{
	bool detAliLocShortSegIncluded = false;
	unsigned int segmentNum_max = SEGMENTNUM;
	unsigned int segmentAliNum_max = CANDALILOC;
	map <unsigned int, unsigned int> candAliLocMap;
	map <unsigned int, unsigned int> ::iterator candAliLoc_it;
	unsigned int mapLabel = 0; // normal
	unsigned int candAliLocWeight[segmentNum_max * segmentAliNum_max]; //normal


	#ifdef DEBUG
	cout << "debug -- start second level detAliLoc for alignment..." << endl;
	#endif
	for (unsigned int tmpSegNum = 1; tmpSegNum <= segmentNum; tmpSegNum++)
	{
		if (//(*(norSegmentLength + tmpSegNum - 1) < minValSegLength)||
			(*(segmentAlignNum + tmpSegNum - 1) > POSSIBLE_MAP_CASES_MAX))
		{
			continue;
		}		

		for(unsigned int tmpSegAliLoc = 1; tmpSegAliLoc <= segmentAlignNum[tmpSegNum-1]; tmpSegAliLoc++) 
		{				
			unsigned int tmpAlignLoc = *(segmentAlignLoc + (tmpSegNum-1)*segmentAliNum_max + tmpSegAliLoc - 1) 
			- *(segmentLocInRead + tmpSegNum - 1) + 1;
			if(tmpAlignLoc > indexInfo->indexSize)
			{
				continue;
			}
			candAliLoc_it = candAliLocMap.find(tmpAlignLoc);

			if ((*(segmentLength + tmpSegNum - 1) >= minValSegLength) && (candAliLoc_it == candAliLocMap.end()))  //cannot find, insert
			{				
				candAliLocMap.insert(pair<unsigned int, unsigned int>(tmpAlignLoc, mapLabel));
				candAliLocWeight[mapLabel] = *(segmentLength + tmpSegNum - 1);

				segMapRangeStart[mapLabel] = tmpSegNum;
				segMapRangeEnd[mapLabel] = tmpSegNum;
				segMapLoc[mapLabel] = tmpAlignLoc;
				mapLabel ++;
			}
			else if ((*(segmentLength + tmpSegNum - 1) < minValSegLength) && (candAliLoc_it == candAliLocMap.end()))  //cannot find, insert
			{	
				if(checkShortSegWithCandLoc(tmpAlignLoc, oldMapLabel, oldSegMapPos, indexInfo))
				{			
					candAliLocMap.insert(pair<unsigned int, unsigned int>(tmpAlignLoc, mapLabel));
					candAliLocWeight[mapLabel] = *(segmentLength + tmpSegNum - 1);

					segMapRangeStart[mapLabel] = tmpSegNum;
					segMapRangeEnd[mapLabel] = tmpSegNum;
					segMapLoc[mapLabel] = tmpAlignLoc;
					mapLabel ++;
				}
			}
			else if ((*(segmentLength + tmpSegNum - 1) >= minValSegLength) && (candAliLoc_it != candAliLocMap.end()))// find, add corresponding weight 
			{
				candAliLocWeight[candAliLoc_it->second] = candAliLocWeight[candAliLoc_it->second] + *(segmentLength + tmpSegNum - 1);
				segMapRangeEnd[candAliLoc_it->second] = tmpSegNum;
			}
			else //((*(segmentLength + tmpSegNum - 1) < minValSegLength) && (candAliLoc_it != candAliLocMap.end()))// find, add corresponding weight 
			{
				candAliLocWeight[candAliLoc_it->second] = candAliLocWeight[candAliLoc_it->second] + *(segmentLength + tmpSegNum - 1);
				segMapRangeEnd[candAliLoc_it->second] = tmpSegNum;
			}
		}
	}
	candAliLocMap.clear();	

	return mapLabel;
}


/*
int main(int argc, char** argv)
{


	return 0;
}
*/