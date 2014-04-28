//#define PE_MAP_DISTANCE 300000

/////////////////////////////// filter multiAlignment cadidate cases with pair-end information /////////////////
///Input:  mapLabel_1, segMapRangeStart_1, segMapRangeEnd_1, segMapLoc_1
///		   mapLabel_2, segMapRangeStart_2, segMapRangeEnd_2, segMapLoc_2
///Output: finalMapLabel_1, finalSegMapRangeStart_1, finalMapRangeEnd_1, finalSegMapLoc_1 //////////////////////
///		   finalMapLabel_2, finalSegMapRangeStart_2, finalMapRangeEnd_2, finalSegMapLoc_2 //////////////////////
#include <math.h>

unsigned int filterOutNum_1 = 0;
unsigned int filterOutNum_2 = 0;
unsigned int filterOutNum_1_noLabel = 0;
unsigned int filterOutNum_2_noLabel = 0;
unsigned int NofilterOutNum_long = 0;
unsigned int NofilterOutNum_short = 0;

bool checkTwoIntInRange(int a, int b, int range)
{
	if(range < 0)
	{
		cout << "range < 0" << endl;
		exit;
	}
	else
	{
		return ( (a-b) < range ) && ( (b-a) < range );
	}
}

bool checkTwoUnsigndIntInRange(unsigned int a, unsigned int b, unsigned int range)
{}

bool filterMapLabelCasesWithPairEndInformation(
	int mapLabel_1, unsigned int* segMapRangeStart_1, unsigned int* segMapRangeEnd_1, unsigned int* segMapLoc_1,
	int mapLabel_2, unsigned int* segMapRangeStart_2, unsigned int* segMapRangeEnd_2, unsigned int* segMapLoc_2,
	int* finalMapLabel_1, unsigned int* finalSegMapRangeStart_1, unsigned int* finalSegMapRangeEnd_1, unsigned int* finalSegMapLoc_1,
	int* finalMapLabel_2, unsigned int* finalSegMapRangeStart_2, unsigned int* finalSegMapRangeEnd_2, unsigned int* finalSegMapLoc_2,
	Index_Info* indexInfo/*,
	int* filterOutMapLabel_1, unsigned int* filterOutSegMapRangeStart_1, unsigned int* filterOutMapRangeEnd_1, unsigned int* filterOutSegMapLoc_1,
	int* filterOutMapLabel_2, unsigned int* filterOutSegMapRangeStart_2, unsigned int* filterOutMapRangeEnd_2, unsigned int* filterOutSegMapLoc_2*/)
{
	typedef set<int> MapLabelSet; // <firstMapLabel, Vec<secondMapLabel> >
	MapLabelSet mapLabelSet_1, mapLabelSet_2;
	unsigned int tmpChrInt_1, tmpChrPosInt_1;
	unsigned int tmpChrInt_2, tmpChrPosInt_2;

	for(int tmpLabel_1 = 0; tmpLabel_1 < mapLabel_1; tmpLabel_1++)
	{
		for(int tmpLabel_2 = 0; tmpLabel_2 < mapLabel_2; tmpLabel_2++)
		{
			indexInfo->getChrLocation(segMapLoc_1[tmpLabel_1], &tmpChrInt_1, &tmpChrPosInt_1);
			indexInfo->getChrLocation(segMapLoc_2[tmpLabel_2], &tmpChrInt_2, &tmpChrPosInt_2);
			if(tmpChrInt_1 == tmpChrInt_2)
			{
				if(checkTwoIntInRange((int)tmpChrPosInt_1, (int)tmpChrPosInt_2, PE_MAP_DISTANCE))
				{
					mapLabelSet_1.insert(tmpLabel_1);
					mapLabelSet_2.insert(tmpLabel_2);
				}
				else
				{
					continue;
				}
			}
			else
			{
				continue;
			}	
		}
	}

	int tmpFinalMapLabel = 0;
	for(MapLabelSet::iterator tmpIter = mapLabelSet_1.begin(); tmpIter != mapLabelSet_1.end(); tmpIter++)
	{
		int tmpMapLabel = (*tmpIter);
		finalSegMapRangeStart_1[tmpFinalMapLabel] = segMapRangeStart_1[tmpMapLabel];
		finalSegMapRangeEnd_1[tmpFinalMapLabel] = segMapRangeEnd_1[tmpMapLabel];		
		finalSegMapLoc_1[tmpFinalMapLabel] = segMapLoc_1[tmpMapLabel];
		tmpFinalMapLabel ++;
	}

	tmpFinalMapLabel = 0;
	for(MapLabelSet::iterator tmpIter = mapLabelSet_2.begin(); tmpIter != mapLabelSet_2.end(); tmpIter++)
	{
		int tmpMapLabel = (*tmpIter);
		finalSegMapRangeStart_2[tmpFinalMapLabel] = segMapRangeStart_2[tmpMapLabel];
		finalSegMapRangeEnd_2[tmpFinalMapLabel] = segMapRangeEnd_2[tmpMapLabel];		
		finalSegMapLoc_2[tmpFinalMapLabel] = segMapLoc_2[tmpMapLabel];
		tmpFinalMapLabel ++;
	}

	(*finalMapLabel_1) = mapLabelSet_1.size();

	(*finalMapLabel_2) = mapLabelSet_2.size();

	return true;
}

const int Nor1_Rcm2_FilteredOut = 1;
const int Nor2_Rcm1_FilteredOut = 2;
const int NoPair_FilteredOut = 0;

int filterMapLabelForReducingComplexity(
	int norMapLabel_1, unsigned int* norSegMapRangeStart_1, unsigned int* norSegMapRangeEnd_1, unsigned int* norSegMapLoc_1,
	int rcmMapLabel_1, unsigned int* rcmSegMapRangeStart_1, unsigned int* rcmSegMapRangeEnd_1, unsigned int* rcmSegMapLoc_1,
	int norMapLabel_2, unsigned int* norSegMapRangeStart_2, unsigned int* norSegMapRangeEnd_2, unsigned int* norSegMapLoc_2,
	int rcmMapLabel_2, unsigned int* rcmSegMapRangeStart_2, unsigned int* rcmSegMapRangeEnd_2, unsigned int* rcmSegMapLoc_2,
	unsigned int* norSegmentLength, unsigned int* rcmSegmentLength, 
	unsigned int* norSegmentLength_PE, unsigned int* rcmSegmentLength_PE,
	int readLength_1,
	int readLength_2
	)
{
	if((norMapLabel_1 == 0) || (rcmMapLabel_2 == 0))
	{
		//norMapLabel_1 = 0;
		//rcmMapLabel_2 = 0;
		filterOutNum_1_noLabel ++;
		return Nor1_Rcm2_FilteredOut;
	}
	else if((norMapLabel_2 == 0) || (rcmMapLabel_1 == 0))
	{
		//norMapLabel_2 = 0;
		//rcmMapLabel_1 = 0;
		filterOutNum_2_noLabel ++;
		return Nor2_Rcm1_FilteredOut;
	}

	int complexity_1 = 0, complexity_2 = 0;
	int segMappedLength_1 = 0, segMappedLength_2 = 0;
	int tmpSegNo = 0;

	typedef set<int>::iterator setIter;
	set<int> norMapSeg_1;
	set<int> rcmMapSeg_1;
	set<int> norMapSeg_2;
	set<int> rcmMapSeg_2;

	for(int tmpLabel = 0; tmpLabel < norMapLabel_1; tmpLabel++)
	{
		for(int tmpSeg = norSegMapRangeStart_1[tmpLabel]; tmpSeg <= norSegMapRangeEnd_1[tmpLabel]; tmpSeg++ )
		{
			norMapSeg_1.insert(tmpSeg);
		}
	}

	for(int tmpLabel = 0; tmpLabel < rcmMapLabel_1; tmpLabel++)
	{
		for(int tmpSeg = rcmSegMapRangeStart_1[tmpLabel]; tmpSeg <= rcmSegMapRangeEnd_1[tmpLabel]; tmpSeg++ )
		{
			rcmMapSeg_1.insert(tmpSeg);
		}
	}	

	for(int tmpLabel = 0; tmpLabel < norMapLabel_2; tmpLabel++)
	{
		for(int tmpSeg = norSegMapRangeStart_2[tmpLabel]; tmpSeg <= norSegMapRangeEnd_2[tmpLabel]; tmpSeg++ )
		{
			norMapSeg_2.insert(tmpSeg);
		}
	}

	for(int tmpLabel = 0; tmpLabel < rcmMapLabel_2; tmpLabel++)
	{
		for(int tmpSeg = rcmSegMapRangeStart_2[tmpLabel]; tmpSeg <= rcmSegMapRangeEnd_2[tmpLabel]; tmpSeg++ )
		{
			rcmMapSeg_2.insert(tmpSeg);
		}
	}


	////////////////////////////////// get mapped length /////////////////////////////////////////////////////	
	tmpSegNo = 0;
	for(setIter tmpSeg = norMapSeg_1.begin(); tmpSeg != norMapSeg_1.end(); tmpSeg++)
	{
		tmpSegNo ++;
		segMappedLength_1 += norSegmentLength[tmpSegNo-1];
	}
	tmpSegNo = 0;
	for(setIter tmpSeg = rcmMapSeg_2.begin(); tmpSeg != rcmMapSeg_2.end(); tmpSeg++)
	{
		tmpSegNo ++;
		segMappedLength_1 += rcmSegmentLength_PE[tmpSegNo-1];
	}	
	tmpSegNo = 0;
	for(setIter tmpSeg = norMapSeg_2.begin(); tmpSeg != norMapSeg_2.end(); tmpSeg++)
	{
		tmpSegNo ++;
		segMappedLength_2 += norSegmentLength_PE[tmpSegNo-1];
	}
	tmpSegNo = 0;
	for(setIter tmpSeg = rcmMapSeg_1.begin(); tmpSeg != rcmMapSeg_1.end(); tmpSeg++)
	{
		tmpSegNo ++;
		segMappedLength_2 += rcmSegmentLength[tmpSegNo-1];
	}

	///////////////////////////////// get complexity /////////////////////////////////////////////////////////
	
	unsigned int norValSegStartNo_1, norValSegEndNo_1;
	unsigned int norPossibleMapCaseMax_1;	
	getValSegNoShortFragIncluded(norMapLabel_1, norSegMapRangeStart_1, 
		norSegMapRangeEnd_1, norSegMapLoc_1, &norValSegStartNo_1, 
		&norValSegEndNo_1, &norPossibleMapCaseMax_1);

	for (int tmpLabel = 0; tmpLabel < norPossibleMapCaseMax_1; tmpLabel ++)
	{
		if(norSegMapRangeEnd_1[tmpLabel] == norValSegEndNo_1)
		{
			complexity_1++;
			continue;
		}
		for(int tmpLabel_2 = norPossibleMapCaseMax_1; tmpLabel_2 < norMapLabel_1; tmpLabel_2 ++)
		{
			if(norSegMapRangeEnd_1[tmpLabel_2] == norValSegEndNo_1)
			{
				complexity_1++;
				continue;
			}
			for(int tmpLabel_3 = tmpLabel_2+1; tmpLabel_3 < norMapLabel_1; tmpLabel_3 ++)
			{				
				if(norSegMapRangeEnd_1[tmpLabel_3] == norValSegEndNo_1)
				{
					complexity_1++;
					continue;
				}
			}
		}
	}

	unsigned int norValSegStartNo_2, norValSegEndNo_2;
	unsigned int norPossibleMapCaseMax_2;	
	getValSegNoShortFragIncluded(norMapLabel_2, norSegMapRangeStart_2, 
		norSegMapRangeEnd_2, norSegMapLoc_2, &norValSegStartNo_2, 
		&norValSegEndNo_2, &norPossibleMapCaseMax_2);

	for (int tmpLabel = 0; tmpLabel < norPossibleMapCaseMax_2; tmpLabel ++)
	{
		if(norSegMapRangeEnd_2[tmpLabel] == norValSegEndNo_2)
		{
			complexity_2++;
			continue;
		}
		for(int tmpLabel_2 = norPossibleMapCaseMax_2; tmpLabel_2 < norMapLabel_2; tmpLabel_2 ++)
		{
			if(norSegMapRangeEnd_2[tmpLabel_2] == norValSegEndNo_2)
			{
				complexity_2++;
				continue;
			}
			for(int tmpLabel_3 = tmpLabel_2+1; tmpLabel_3 < norMapLabel_2; tmpLabel_3 ++)
			{				
				if(norSegMapRangeEnd_2[tmpLabel_3] == norValSegEndNo_2)
				{
					complexity_2++;
					continue;
				}
			}
		}
	}

	unsigned int rcmValSegStartNo_1, rcmValSegEndNo_1;
	unsigned int rcmPossibleMapCaseMax_1;	
	getValSegNoShortFragIncluded(rcmMapLabel_1, rcmSegMapRangeStart_1, 
		rcmSegMapRangeEnd_1, rcmSegMapLoc_1, &rcmValSegStartNo_1, 
		&rcmValSegEndNo_1, &rcmPossibleMapCaseMax_1);

	for (int tmpLabel = 0; tmpLabel < rcmPossibleMapCaseMax_1; tmpLabel ++)
	{
		if(rcmSegMapRangeEnd_1[tmpLabel] == rcmValSegEndNo_1)
		{
			complexity_2++;
			continue;
		}
		for(int tmpLabel_2 = rcmPossibleMapCaseMax_1; tmpLabel_2 < rcmMapLabel_1; tmpLabel_2 ++)
		{
			if(rcmSegMapRangeEnd_1[tmpLabel_2] == rcmValSegEndNo_1)
			{
				complexity_2++;
				continue;
			}
			for(int tmpLabel_3 = tmpLabel_2+1; tmpLabel_3 < rcmMapLabel_1; tmpLabel_3 ++)
			{				
				if(rcmSegMapRangeEnd_1[tmpLabel_3] == rcmValSegEndNo_1)
				{
					complexity_2++;
					continue;
				}
			}
		}
	}	

	unsigned int rcmValSegStartNo_2, rcmValSegEndNo_2;
	unsigned int rcmPossibleMapCaseMax_2;	
	getValSegNoShortFragIncluded(rcmMapLabel_2, rcmSegMapRangeStart_2, 
		rcmSegMapRangeEnd_2, rcmSegMapLoc_2, &rcmValSegStartNo_2, 
		&rcmValSegEndNo_2, &rcmPossibleMapCaseMax_2);

	for (int tmpLabel = 0; tmpLabel < rcmPossibleMapCaseMax_2; tmpLabel ++)
	{
		if(rcmSegMapRangeEnd_2[tmpLabel] == rcmValSegEndNo_2)
		{
			complexity_1++;
			continue;
		}
		for(int tmpLabel_2 = rcmPossibleMapCaseMax_2; tmpLabel_2 < rcmMapLabel_2; tmpLabel_2 ++)
		{
			if(rcmSegMapRangeEnd_2[tmpLabel_2] == rcmValSegEndNo_2)
			{
				complexity_1++;
				continue;
			}
			for(int tmpLabel_3 = tmpLabel_2+1; tmpLabel_3 < rcmMapLabel_2; tmpLabel_3 ++)
			{				
				if(rcmSegMapRangeEnd_2[tmpLabel_3] == rcmValSegEndNo_2)
				{
					complexity_1++;
					continue;
				}
			}
		}
	}	
	
	////////////////////////////////// get comparison result /////////////////////////////////////////////////
	double segMappedPercentage_1 = (double)segMappedLength_1/(readLength_1*2);
	double segMappedPercentage_2 = (double)segMappedLength_2/(readLength_2*2);

	//double comparison_1 = log((double)complexity_1)/segMappedPercentage_1;
	//double comparison_2 = log((double)complexity_2)/segMappedPercentage_2;
	/*
	cout << "comparison_1 = " << comparison_1 << " complexity_1 = " << complexity_1 << " segMappedPercentage_1 = " << segMappedPercentage_1 << endl;
	cout << "comparison_2 = " << comparison_2 << " complexity_2 = " << complexity_2 << " segMappedPercentage_2 = " << segMappedPercentage_2 << endl;
	*/
	/*
	if(comparison_1 > comparison_2 * 10)
	{
		return Nor1_Rcm2_FilteredOut;
	}
	else if(comparison_2 > comparison_1 * 10)
	{
		return Nor2_Rcm1_FilteredOut;
	}
	else
	{
		return NoPair_FilteredOut;
	}*/

	if((segMappedPercentage_1 < 0.4) && (segMappedPercentage_2 < 0.4))	
	{
		NofilterOutNum_short++;
		return NoPair_FilteredOut;
	}
	else if((segMappedPercentage_1/segMappedPercentage_2 > 2) //&& (segMappedPercentage_2 < 0.5)
		)
	{
		filterOutNum_2 ++;
		return Nor2_Rcm1_FilteredOut;
	}
	else if((segMappedPercentage_2/segMappedPercentage_1 > 2) //&& (segMappedPercentage_2 > 0.8)
		)
	{
		filterOutNum_1 ++;
		return Nor1_Rcm2_FilteredOut;
	}
	else
	{
		NofilterOutNum_long++;
		return NoPair_FilteredOut;
	}

		
}