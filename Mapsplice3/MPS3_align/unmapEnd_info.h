
class UnmapEnd_Info
{
public:
	//Alignment_Info* otherEndAlignInfo;
	
	string chrNameStr;
	int chrMapPos_start;
	int chrMapPos_end;

	int mapPosIntervalStart;
	int mapPosIntervalEnd;
	int secondLevelIndexNum; // start from 1

	int chrPosStartIn2ndLevelIndex;

	UnmapEnd_Info()
	{}

	bool set_UnmapEnd_Info(//Alignment_Info* alignInfo, bool End1OrEnd2, bool NorOrRcm, Index_Info* indexInfo
		PE_Read_Alignment_Info* peAlignInfo, Index_Info* indexInfo)
	{
		bool setUnmapEndInfo = false;
			Alignment_Info* alignInfo;
			bool End1OrEnd2;
			bool NorOrRcm;
			
			if((peAlignInfo->oriAlignPair_Nor1Rcm2).size() + 
				(peAlignInfo->oriAlignPair_Nor2Rcm1).size() == 0)
			{
				//cout << "no pair !" << endl;
				
				if(((peAlignInfo->norAlignmentInfo_PE_1).size() 
						+ (peAlignInfo->rcmAlignmentInfo_PE_1).size() 
						+ (peAlignInfo->norAlignmentInfo_PE_2).size()
						+ (peAlignInfo->rcmAlignmentInfo_PE_2).size()) == 1)
				{
					//cout << "1 candi alignment !" << endl;
					if((peAlignInfo->norAlignmentInfo_PE_1).size() == 1)
					{
						alignInfo = (peAlignInfo->norAlignmentInfo_PE_1)[0];
						End1OrEnd2 = true; NorOrRcm = true;
					}
					else if((peAlignInfo->rcmAlignmentInfo_PE_1).size() == 1)
					{
						alignInfo = (peAlignInfo->rcmAlignmentInfo_PE_1)[0];
						End1OrEnd2 = true; NorOrRcm = false;
					}
					else if((peAlignInfo->norAlignmentInfo_PE_2).size() == 1)
					{
						alignInfo = (peAlignInfo->norAlignmentInfo_PE_2)[0];
						End1OrEnd2 = false; NorOrRcm = true;
					}
					else
					{
						alignInfo = (peAlignInfo->rcmAlignmentInfo_PE_2)[0];
						End1OrEnd2 = false; NorOrRcm = false;
					}
					//cout << "End1OrEnd2: " << End1OrEnd2 << endl;
					//cout << "NorOrRcm: " << NorOrRcm << endl;
		
					this->setUnmapEndInfo(alignInfo, End1OrEnd2, NorOrRcm, indexInfo);
					setUnmapEndInfo = true;
				}
				else
				{
					//cout << "multiple alignments !" << endl;
				}
			}
			else
			{
				//cout << "pair exits !" << endl;
			}



		return setUnmapEndInfo;	
	}

	void setUnmapEndInfo(Alignment_Info* alignInfo, bool End1OrEnd2, bool NorOrRcm, Index_Info* indexInfo)
	{
		//otherEndAlignInfo = alignInfo;
		chrNameStr = alignInfo->alignChromName;
		chrMapPos_start = alignInfo->alignChromPos;
		chrMapPos_end = alignInfo->getEndMatchedPosInChr();

		bool ForwardOrReverseBool; // true: known -- unknown; false: unknown -- known

		if((End1OrEnd2)&&(NorOrRcm))
		{
			ForwardOrReverseBool = true;
		}
		else if((End1OrEnd2)&&(!NorOrRcm))
		{
			ForwardOrReverseBool = false;
		}
		else if((!End1OrEnd2)&&(NorOrRcm))
		{
			ForwardOrReverseBool = false;
		}
		else
		{
			ForwardOrReverseBool = true;
		}


		if(ForwardOrReverseBool)
		{
			mapPosIntervalStart = chrMapPos_start;
			mapPosIntervalEnd = chrMapPos_end + READ_ALIGN_AREA_LENGTH;
		} 
		else
		{
			mapPosIntervalStart = chrMapPos_start - READ_ALIGN_AREA_LENGTH;
			mapPosIntervalEnd = chrMapPos_end;
		}

		int chrNameInt = indexInfo->convertStringToInt(chrNameStr);

		secondLevelIndexNum = indexInfo->getSecondLevelIndexFromChrAndPos(chrNameInt, chrMapPos_start); // Xinan: need to debug

		//unsigned int wholeGenomePos_start = getWholeGenomeLocation((unsigned int)chrNameInt, (unsigned int)chrMapPos_start);
		//unsigned int wholeGenomePosStartIn2ndLevelIndex = (secondLevelIndexNum - 1) * (indexInfo->secondLevelIndexNormalSize);
		chrPosStartIn2ndLevelIndex = indexInfo->getChrPosFromSecondLevelIndexPos(chrNameInt, secondLevelIndexNum, 1);
	}

};
