#include <string>
#include <string.h>
#include "splice_info.h"

using namespace std;

class FixHeadTailInfo
{
public:

	FixHeadTailInfo()
	{

	}

	void fixHeadTail(PE_Read_Info* peReadInfo, PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJ,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		vector<unsigned int*>& secondLevelLcp,
		vector<unsigned int*>& secondLevelUp,
		vector<unsigned int*>& secondLevelDown,
		vector<unsigned int*>& secondLevelNext,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo)
	{
		int readLength_1 = ((peReadInfo->readInfo_pe1).readSeq).length();
		int readLength_2 = ((peReadInfo->readInfo_pe2).readSeq).length();

		Alignment_Info* tmpAlignmentInfo;

		int readLength = 0;

		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= 4; tmpAlignInfoType++)
		{
			readLength = tmpAlignInfoType <= 2
				? readLength_1
				: readLength_2;

			for(int tmpAlignmentNO = 0; tmpAlignmentNO < (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
				tmpAlignmentNO++)
			{
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				if( (tmpAlignmentInfo->cigarStringJumpCode)[0].type != "S" )
					continue; // no unfixed head

				// start to nor1
				Unfixed_Head unfixedHeadInfo;
				unfixedHeadInfo.getUnfixedHeadInfoFromRecordWithAlignInfoType(peReadInfo, tmpAlignInfoType, tmpAlignmentInfo, indexInfo);
				// finish getUnfixedTailInfoFromRecord

				string readSeqWithDirection;
				if(unfixedHeadInfo.alignDirection == "+")
				{
					readSeqWithDirection = unfixedHeadInfo.readSeqOriginal;
				}
				else
				{
					readSeqWithDirection = covertStringToReverseComplement(unfixedHeadInfo.readSeqOriginal);
				}


				///////////////////////////////////////////////////////////////////////////////////////////
				/////////////////////////  try remapping with splice junction hash ////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////
				//cout << "start SJsearchInSJhash !" << endl;
				bool spliceJunctionFoundInHash;

				if(spliceJunctionHashExists)
				{
					spliceJunctionFoundInHash = unfixedHeadInfo.SJsearchInSJhash(SJ, readSeqWithDirection, indexInfo);
				}
				else
				{
					spliceJunctionFoundInHash = false;
				}

				if(spliceJunctionFoundInHash)
				{

					int tmpDonerEndPosInRead = (unfixedHeadInfo.SJposFromRemapping)[0].first;
					int tmpSpliceJunctionDistance = (unfixedHeadInfo.SJposFromRemapping)[0].second;
					int tmpFirstMatchLength = tmpDonerEndPosInRead;

					Alignment_Info* newTmpAlignInfo = new Alignment_Info();
					newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
						tmpFirstMatchLength, tmpSpliceJunctionDistance, indexInfo);

					//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
					peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpAlignmentNO);

					for(int tmpSJposVec = 1; tmpSJposVec < (unfixedHeadInfo.SJposFromRemapping).size(); tmpSJposVec++)
					{
						tmpDonerEndPosInRead = (unfixedHeadInfo.SJposFromRemapping)[tmpSJposVec].first;
						tmpSpliceJunctionDistance = (unfixedHeadInfo.SJposFromRemapping)[tmpSJposVec].second;
						tmpFirstMatchLength = tmpDonerEndPosInRead;

						Alignment_Info* newTmpAlignInfo = new Alignment_Info();
						newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
							tmpFirstMatchLength, tmpSpliceJunctionDistance, indexInfo);
						//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
						peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
					}
					continue;
				}

				///////////////////////////////////////////////////////////////////////////////////////////
				///////////////////////// remapping with splice junction hash failed //////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////
				//cout << "start to getPossibleSJpos " << endl;
				unfixedHeadInfo.getPossibleSJpos(readSeqWithDirection, indexInfo);
				//cout << "finish getting Possible SJ pos" << endl;
				///////////////////////////////////////////////////////////////////////////////////////////
				///////////////////////////////   check possible SJs    ///////////////////////////////////
				//cout << "unfixedHeadInfo.midPartMapPosInChr: " << unfixedHeadInfo.midPartMapPosInChr << endl;
				int midPartMapPosSecondLevelIndexNO
					= indexInfo->getSecondLevelIndexFromChrAndPos(unfixedHeadInfo.midPartMapChrInt,
						unfixedHeadInfo.midPartMapPosInChr);
				midPartMapPosSecondLevelIndexNO --;

				if(
					((indexInfo->invalidSecondLevelIndexNOset).find(midPartMapPosSecondLevelIndexNO+1))
					!= (indexInfo->invalidSecondLevelIndexNOset).end()
					)
				{
					continue;
				}

				unsigned int midPartMapPosForLongHeadInSecondLevelIndex = unfixedHeadInfo.midPartMapPosInChr -
					((unfixedHeadInfo.midPartMapPosInChr)/(indexInfo->secondLevelIndexNormalSize))*(indexInfo->secondLevelIndexNormalSize);

				//   check GTAG splice junctions
				for(int tmp = 0; tmp < (unfixedHeadInfo.possiGTAGpos).size(); tmp++)
				{
					int tmpSJposInRead = unfixedHeadInfo.possiGTAGpos[tmp];

					if(tmpSJposInRead-1 < min_anchor_length)
						continue;

					string tmpShortAnchorStr = readSeqWithDirection.substr(0, tmpSJposInRead-1);
					string targetMappingStr = tmpShortAnchorStr + "GT";

					char* headChar = const_cast<char*>(targetMappingStr.c_str());
					int targetMappingNum = 0;
					unsigned int targetMappingLoc[100];

					unsigned int finalMidPartMappingPos
						= unfixedHeadInfo.midPartMapPosInWholeGenome + tmpSJposInRead - unfixedHeadInfo.unfixedHeadLength - 1;
					bool headSegMapMain = mapMainSecondLevelForTargetMapping_compressedIndex(headChar,
						secondLevelSa[midPartMapPosSecondLevelIndexNO],
						secondLevelLcpCompress[midPartMapPosSecondLevelIndexNO],
						secondLevelChildTab[midPartMapPosSecondLevelIndexNO],
						secondLevelDetChild[midPartMapPosSecondLevelIndexNO],
						secondLevelChrom[midPartMapPosSecondLevelIndexNO],
						targetMappingStr.length(), 3000000, unfixedHeadInfo.midPartMapPosInWholeGenome,
						midPartMapPosForLongHeadInSecondLevelIndex, &targetMappingNum, targetMappingLoc,
						indexInfo
						);


					if(headSegMapMain)
					{
						int tmpMinSpliceDistance = 300001;
						for(int tmp2 = 0; tmp2 < targetMappingNum; tmp2++)
						{
							int tmpSpliceDistance = finalMidPartMappingPos - (targetMappingLoc[tmp2] + unfixedHeadInfo.possiGTAGpos[tmp] - 2) - 1;
							if((tmpSpliceDistance < 300000)
								&& (tmpSpliceDistance > -4))
							{
								if(tmpSpliceDistance < tmpMinSpliceDistance)
									tmpMinSpliceDistance = tmpSpliceDistance;
							}
						}
						if(tmpMinSpliceDistance < 300000)
						{
							(unfixedHeadInfo.GTAGsjPos).push_back(pair<int,int>(unfixedHeadInfo.possiGTAGpos[tmp], tmpMinSpliceDistance));
						}
					}

				}
				//   check CTAC splice junctions
				//cout << "unfixedHeadInfo.possiCTACpos.size(): " << (unfixedHeadInfo.possiCTACpos).size() << endl;
				for(int tmp = 0; tmp < (unfixedHeadInfo.possiCTACpos).size(); tmp++)
				{
					int tmpSJposInRead = unfixedHeadInfo.possiCTACpos[tmp];

					if(tmpSJposInRead-1 < min_anchor_length)
						continue;

					string tmpShortAnchorStr = readSeqWithDirection.substr(0, tmpSJposInRead-1);
					string targetMappingStr = tmpShortAnchorStr + "CT";
					//cout << "headLength: " << tmpSJposInRead-1 << " targetStr: " << targetMappingStr << endl;
					char* headChar = const_cast<char*>(targetMappingStr.c_str());
					int targetMappingNum = 0;
					unsigned int targetMappingLoc[100];

					unsigned int finalMidPartMappingPos
						= unfixedHeadInfo.midPartMapPosInWholeGenome + tmpSJposInRead - unfixedHeadInfo.unfixedHeadLength - 1;

					bool headSegMapMain = mapMainSecondLevelForTargetMapping_compressedIndex(headChar,
						secondLevelSa[midPartMapPosSecondLevelIndexNO],
						secondLevelLcpCompress[midPartMapPosSecondLevelIndexNO],
						secondLevelChildTab[midPartMapPosSecondLevelIndexNO],
						secondLevelDetChild[midPartMapPosSecondLevelIndexNO],
						secondLevelChrom[midPartMapPosSecondLevelIndexNO],
						targetMappingStr.length(), 3000000, unfixedHeadInfo.midPartMapPosInWholeGenome,
						midPartMapPosForLongHeadInSecondLevelIndex, &targetMappingNum, targetMappingLoc,
						indexInfo);

					if(headSegMapMain)
					{
						int tmpMinSpliceDistance = 300001;
						for(int tmp2 = 0; tmp2 < targetMappingNum; tmp2++)
						{
							int tmpSpliceDistance = finalMidPartMappingPos - (targetMappingLoc[tmp2] + unfixedHeadInfo.possiCTACpos[tmp] - 2) - 1;
							if((tmpSpliceDistance < 300000)
								&& (tmpSpliceDistance > -4))
							{
								if(tmpSpliceDistance < tmpMinSpliceDistance)
									tmpMinSpliceDistance = tmpSpliceDistance;
							}
						}
						if(tmpMinSpliceDistance < 300000)
						{
							(unfixedHeadInfo.CTACsjPos).push_back(pair<int,int>(unfixedHeadInfo.possiCTACpos[tmp], tmpMinSpliceDistance));
						}
					}

				}

				if((unfixedHeadInfo.GTAGsjPos).size() > 0)
				{
					Alignment_Info* newTmpAlignInfo = new Alignment_Info();
					newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
						(unfixedHeadInfo.GTAGsjPos[0]).first - 1, (unfixedHeadInfo.GTAGsjPos[0]).second, indexInfo);

					peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpAlignmentNO);

					for(int tmp = 1; tmp < (unfixedHeadInfo.GTAGsjPos).size(); tmp++)
					{
						Alignment_Info* newTmpAlignInfo
							= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
							(unfixedHeadInfo.GTAGsjPos[tmp]).first - 1, (unfixedHeadInfo.GTAGsjPos[tmp]).second, indexInfo);

						peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
					}

					for(int tmp = 0; tmp < (unfixedHeadInfo.CTACsjPos).size(); tmp++)
					{
						Alignment_Info* newTmpAlignInfo
							= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
							(unfixedHeadInfo.CTACsjPos[tmp]).first - 1, (unfixedHeadInfo.CTACsjPos[tmp]).second, indexInfo);

						peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
					}
				}
				else if((unfixedHeadInfo.CTACsjPos).size() > 0)
				{
					Alignment_Info* newTmpAlignInfo = new Alignment_Info();
					newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
						(unfixedHeadInfo.CTACsjPos[0]).first - 1, (unfixedHeadInfo.CTACsjPos[0]).second, indexInfo);

					peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpAlignmentNO);

					for(int tmp = 1; tmp < (unfixedHeadInfo.CTACsjPos).size(); tmp++)
					{
						Alignment_Info* newTmpAlignInfo
							= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
							(unfixedHeadInfo.CTACsjPos[tmp]).first - 1, (unfixedHeadInfo.CTACsjPos[tmp]).second, indexInfo);
						peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
					}
				}
				else
				{}

			}
		}

		/////////////////// fix tail //////////////////////////////
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= 4; tmpAlignInfoType++)
		{
			if(tmpAlignInfoType <= 2)
			{
				readLength = readLength_1;
			}
			else
			{
				readLength = readLength_2;
			}

			for(int tmpAlignmentNO = 0; tmpAlignmentNO < (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
				tmpAlignmentNO++)
			{
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				int cigarStringJumpCodeSize = (tmpAlignmentInfo->cigarStringJumpCode).size();
				if( (tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].type != "S" )
				{
					continue;
				}

				Unfixed_Tail unfixedTailInfo;
				unfixedTailInfo.getUnfixedTailInfoFromRecordWithAlignInfoType(peReadInfo, tmpAlignInfoType, tmpAlignmentInfo, indexInfo);
				string readSeqWithDirection;
				if(unfixedTailInfo.alignDirection == "+")
				{
					// Do not need to copy because we already know the direction in the 4 cases
					readSeqWithDirection = unfixedTailInfo.readSeqOriginal;
				}
				else
				{
					readSeqWithDirection = covertStringToReverseComplement(unfixedTailInfo.readSeqOriginal);
				}

				///////////////////////////////////////////////////////////////////////////////////////////
				/////////////////////////  try remapping with splice junction hash ////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////
				//cout << "start SJsearchInSJhash " << endl;
				bool spliceJunctionFoundInHash;

				if(spliceJunctionHashExists)
				{
					spliceJunctionFoundInHash = unfixedTailInfo.SJsearchInSJhash(SJ, readSeqWithDirection, indexInfo);
				}
				else
				{
					spliceJunctionFoundInHash = false;
				}

				if(spliceJunctionFoundInHash)
				{
					//if((unfixedTailInfo.SJposFromRemapping).size() == 1)
					int tmpDonerEndPosInRead = (unfixedTailInfo.SJposFromRemapping)[0].first;
					int tmpSpliceJunctionDistance = (unfixedTailInfo.SJposFromRemapping)[0].second;
					int tmpLastMatchLength = readLength - tmpDonerEndPosInRead;

					Alignment_Info* newTmpAlignInfo = new Alignment_Info();
					newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment(
						tmpLastMatchLength, tmpSpliceJunctionDistance, indexInfo);

					//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
					peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpAlignmentNO);
					for(int tmpSJposVec = 1; tmpSJposVec < (unfixedTailInfo.SJposFromRemapping).size(); tmpSJposVec++)
					{
						tmpDonerEndPosInRead = (unfixedTailInfo.SJposFromRemapping)[tmpSJposVec].first;
						tmpSpliceJunctionDistance = (unfixedTailInfo.SJposFromRemapping)[tmpSJposVec].second;
						tmpLastMatchLength = readLength - tmpDonerEndPosInRead;

						Alignment_Info* newTmpAlignInfo = new Alignment_Info();
						newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment(
							tmpLastMatchLength, tmpSpliceJunctionDistance, indexInfo);
						//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
						peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
					}

					continue;
				}

				///////////////////////////////////////////////////////////////////////////////////////////
				///////////////////////// remapping with splice junction hash failed //////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////
				//cout << "start to get possible SJ pos" << endl;
				unfixedTailInfo.getPossibleSJpos(readSeqWithDirection, indexInfo->chromString, indexInfo);
				//cout << "finish getting possible SJ pos " << endl;
				///////////////////////////////////////////////////////////////////////////////////////////
				///////////////////////////////   check possible SJs    ///////////////////////////////////

				int midPartMapPosSecondLevelIndexNO
					= indexInfo->getSecondLevelIndexFromChrAndPos(unfixedTailInfo.midPartMapChrInt,
						unfixedTailInfo.midPartMapPosInChr);
				midPartMapPosSecondLevelIndexNO --;
				//cout << "midPartMapPosSecondLevelIndexNO: " << midPartMapPosSecondLevelIndexNO << endl;
				if(
					((indexInfo->invalidSecondLevelIndexNOset).find(midPartMapPosSecondLevelIndexNO+1))
					!= (indexInfo->invalidSecondLevelIndexNOset).end()
					)
				{
					//delete(seg2ndOriInfo);
					//delete(unmapEndInfo);
					continue;
				}


				unsigned int midPartMapPosForLongTailInSecondLevelIndex = unfixedTailInfo.midPartMapPosInChr -
					((unfixedTailInfo.midPartMapPosInChr)/(indexInfo->secondLevelIndexNormalSize))*(indexInfo->secondLevelIndexNormalSize) ;
				//cout << "SecondLevelIndexNO: " << midPartMapPosSecondLevelIndexNO << endl;

				//   check GTAG splice junctions
				//cout << "SJposInRead: " << endl;
				for(int tmp = 0; tmp < (unfixedTailInfo.possiGTAGpos).size(); tmp++)
				{
					int tmpSJposInRead = unfixedTailInfo.possiGTAGpos[tmp];

					//cout << "SJposInRead: " << tmpSJposInRead << endl;
					if(readLength - tmpSJposInRead + 1 < min_anchor_length)
						continue;
					string tmpShortAnchorStr
						= readSeqWithDirection.substr(tmpSJposInRead-1, readLength - tmpSJposInRead + 1);

					string targetMappingStr = "AG" + tmpShortAnchorStr;
					//cout << "targetStr: " << targetMappingStr << endl;
					char* tailChar = const_cast<char*>(targetMappingStr.c_str());
					int targetMappingNum = 0;
					unsigned int targetMappingLoc[100];

					unsigned int finalMidPartMappingPos
						= unfixedTailInfo.midPartMapPosInWholeGenome + tmpSJposInRead - readLength + unfixedTailInfo.unfixedTailLength - 1;

					bool tailSegMapMain = mapMainSecondLevelForTargetMapping_compressedIndex(tailChar,
							secondLevelSa[midPartMapPosSecondLevelIndexNO],
							secondLevelLcpCompress[midPartMapPosSecondLevelIndexNO],
							secondLevelChildTab[midPartMapPosSecondLevelIndexNO],
							secondLevelDetChild[midPartMapPosSecondLevelIndexNO],
							secondLevelChrom[midPartMapPosSecondLevelIndexNO],
							targetMappingStr.length(), 3000000, unfixedTailInfo.midPartMapPosInWholeGenome,
							midPartMapPosForLongTailInSecondLevelIndex, &targetMappingNum, targetMappingLoc
							, indexInfo);

					if(tailSegMapMain)// select the nearest short anchor location
					{
						int tmpMinSpliceDistance = 300001;
						for(int tmp2 = 0; tmp2 < targetMappingNum; tmp2++)
						{
							int tmpSpliceDistance = targetMappingLoc[tmp2] - finalMidPartMappingPos + 1;
							if((tmpSpliceDistance < 300000)
								&& (tmpSpliceDistance > -4))
							{
								if(tmpSpliceDistance < tmpMinSpliceDistance)
									tmpMinSpliceDistance = tmpSpliceDistance;
							}
						}
						if(tmpMinSpliceDistance < 300000)
						{
							(unfixedTailInfo.GTAGsjPos).push_back(pair<int,int>(unfixedTailInfo.possiGTAGpos[tmp], tmpMinSpliceDistance));
						}
					}
				}

				//   check CTAC splice junctions
				for(int tmp = 0; tmp < (unfixedTailInfo.possiCTACpos).size(); tmp++)
				{
					int tmpSJposInRead = unfixedTailInfo.possiCTACpos[tmp];

					//cout << "SJposInRead: " << tmpSJposInRead << endl;
					if(readLength - tmpSJposInRead + 1 < min_anchor_length)
						continue;
					string tmpShortAnchorStr
						= readSeqWithDirection.substr(tmpSJposInRead-1, readLength - tmpSJposInRead + 1);

					string targetMappingStr = "AC" + tmpShortAnchorStr;
					//cout << "targetStr: " << targetMappingStr << endl;
					char* tailChar = const_cast<char*>(targetMappingStr.c_str());
					int targetMappingNum = 0;
					unsigned int targetMappingLoc[100];

					unsigned int finalMidPartMappingPos
						= unfixedTailInfo.midPartMapPosInWholeGenome + tmpSJposInRead - readLength + unfixedTailInfo.unfixedTailLength - 1;

					bool tailSegMapMain = mapMainSecondLevelForTargetMapping_compressedIndex(tailChar,
							secondLevelSa[midPartMapPosSecondLevelIndexNO],
							secondLevelLcpCompress[midPartMapPosSecondLevelIndexNO],
							secondLevelChildTab[midPartMapPosSecondLevelIndexNO],
							secondLevelDetChild[midPartMapPosSecondLevelIndexNO],
							secondLevelChrom[midPartMapPosSecondLevelIndexNO],
							targetMappingStr.length(), 3000000, unfixedTailInfo.midPartMapPosInWholeGenome,
							midPartMapPosForLongTailInSecondLevelIndex, &targetMappingNum, targetMappingLoc
							, indexInfo);

					if(tailSegMapMain)
					{
						int tmpMinSpliceDistance = 300001;
						for(int tmp2 = 0; tmp2 < targetMappingNum; tmp2++)
						{
							int tmpSpliceDistance = targetMappingLoc[tmp2] - finalMidPartMappingPos + 1;
							if((tmpSpliceDistance < 300000)
								&& (tmpSpliceDistance > -4))
							{
								if(tmpSpliceDistance < tmpMinSpliceDistance)
									tmpMinSpliceDistance = tmpSpliceDistance;
							}
						}
						if(tmpMinSpliceDistance < 300000)
						{
							(unfixedTailInfo.CTACsjPos).push_back(pair<int,int>(unfixedTailInfo.possiCTACpos[tmp], tmpMinSpliceDistance));
						}
					}
				}

				if((unfixedTailInfo.GTAGsjPos).size() > 0)
				{
					Alignment_Info* newTmpAlignInfo = new Alignment_Info();
					newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment(
						readLength - (unfixedTailInfo.GTAGsjPos[0]).first + 1, (unfixedTailInfo.GTAGsjPos[0]).second, indexInfo);

					//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
					peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpAlignmentNO);

					for(int tmp = 1; tmp < (unfixedTailInfo.GTAGsjPos).size(); tmp++)
					{
						Alignment_Info* newTmpAlignInfo
							= tmpAlignmentInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment(
							readLength - (unfixedTailInfo.GTAGsjPos[tmp]).first + 1, (unfixedTailInfo.GTAGsjPos[tmp]).second, indexInfo);
						//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
						peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
					}

					for(int tmp = 0; tmp < (unfixedTailInfo.CTACsjPos).size(); tmp++)
					{
						Alignment_Info* newTmpAlignInfo
							= tmpAlignmentInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment(
							readLength - (unfixedTailInfo.CTACsjPos[tmp]).first + 1, (unfixedTailInfo.CTACsjPos[tmp]).second, indexInfo);
						//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
						peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
					}
				}
				else if((unfixedTailInfo.CTACsjPos).size() > 0)
				{
					Alignment_Info* newTmpAlignInfo = new Alignment_Info();
					newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment(
						readLength - (unfixedTailInfo.CTACsjPos[0]).first + 1, (unfixedTailInfo.CTACsjPos[0]).second, indexInfo);

					//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
					peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpAlignmentNO);

					for(int tmp = 1; tmp < (unfixedTailInfo.CTACsjPos).size(); tmp++)
					{
						Alignment_Info* newTmpAlignInfo
							= tmpAlignmentInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment(
							readLength - (unfixedTailInfo.CTACsjPos[tmp]).first + 1, (unfixedTailInfo.CTACsjPos[tmp]).second, indexInfo);
						//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
						peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
					}
				}
				else
				{
				}

			}
		}

	}

};
