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

	/*void fixHeadTail(PE_Read_Info* peReadInfo, PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJ, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		vector<unsigned int*>& secondLevelLcp,
		vector<unsigned int*>& secondLevelUp,
		vector<unsigned int*>& secondLevelDown,
		vector<unsigned int*>& secondLevelNext, 
		bool load2ndLevelIndexBool_compressedSize,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo)
	{
		int readLength_1 = ((peReadInfo->readInfo_pe1).readSeq).length();
		int readLength_2 = ((peReadInfo->readInfo_pe2).readSeq).length();

		Alignment_Info* tmpAlignmentInfo;

		int readLength = 0;

				//////////////////// fix head //////////////////////////////
				//cout << "start to fix head" << endl;
				for(int tmpAlignInfoType = 1; tmpAlignInfoType <= 4; tmpAlignInfoType++)
				{
					//readLength = readLength_1 * (tmpAlignInfoType <= 2) + readLength_2 * (tmpAlignInfoType > 2);
					
					//cout << "tmpAlignInfoType: " << tmpAlignInfoType << endl;

					if(tmpAlignInfoType <= 2)
					{
						readLength = readLength_1;
					}
					else
					{
						readLength = readLength_2;
					}

					//cout << "readLength: " << readLength << endl;
					for(int tmpAlignmentNO = 0; tmpAlignmentNO < (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType)); 
						tmpAlignmentNO++)
					{
						//cout << "tmpAlignmentNO: " << tmpAlignmentNO << endl;
						//tmpAlignmentInfo = (peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO];
						tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
						if( (tmpAlignmentInfo->cigarStringJumpCode)[0].type != "S" )
						{
							//cout << "no unfixed head !" << endl;
							continue;
						}

						//cout << "start to nor1" << endl;
						//cout << "start to getUnfixedHeadInfoFromRecordWithAlignInfoType ..." << endl;
						Unfixed_Head unfixedHeadInfo;
						unfixedHeadInfo.getUnfixedHeadInfoFromRecordWithAlignInfoType(peReadInfo, tmpAlignInfoType, tmpAlignmentInfo, indexInfo);

						//cout << "finish getUnfixedTailInfoFromRecord ..." << endl;

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
							/////////////////////////////////////////////////////////////////////////////
							////////////////////////     position hash     //////////////////////////////
							/////////////////////////////////////////////////////////////////////////////
							spliceJunctionFoundInHash 
								= unfixedHeadInfo.SJsearchInSJhash(SJ, readSeqWithDirection, indexInfo);
							/////////////////////////////////////////////////////////////////////////////
							////////////////////////      string hash      //////////////////////////////
							/////////////////////////////////////////////////////////////////////////////
							//spliceJunctionFoundInHash
							//	= unfixedHeadInfo.SJsearchInSJhash(SJ, readSeqWithDirection, indexInfo, SJ->areaSize);
						}
						else
						{
							spliceJunctionFoundInHash = false;
						}
						//cout << "spliceJunctionFoundInHash: " << spliceJunctionFoundInHash << endl;

						if(spliceJunctionFoundInHash)
						{

							int tmpDonerEndPosInRead = (unfixedHeadInfo.SJposFromRemappingVec)[0].first;
							int tmpSpliceJunctionDistance = (unfixedHeadInfo.SJposFromRemappingVec)[0].second;
							int tmpFirstMatchLength = tmpDonerEndPosInRead;

							Alignment_Info* newTmpAlignInfo = new Alignment_Info();
							newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
								tmpFirstMatchLength, tmpSpliceJunctionDistance, indexInfo);

							//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
							peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpAlignmentNO);

							for(int tmpSJposVec = 1; tmpSJposVec < (unfixedHeadInfo.SJposFromRemappingVec).size(); tmpSJposVec++)
							{
								tmpDonerEndPosInRead = (unfixedHeadInfo.SJposFromRemappingVec)[tmpSJposVec].first;
								tmpSpliceJunctionDistance = (unfixedHeadInfo.SJposFromRemappingVec)[tmpSJposVec].second;
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

						//cout << "2ndLevel Index: " << midPartMapPosSecondLevelIndexNO << endl;		

						if(
							((indexInfo->invalidSecondLevelIndexNOset).find(midPartMapPosSecondLevelIndexNO+1))
							!= (indexInfo->invalidSecondLevelIndexNOset).end()
							)
						{
							//delete(seg2ndOriInfo);
							//delete(unmapEndInfo);
							continue;
						}

						//cout << "2ndLevel Index: " << midPartMapPosSecondLevelIndexNO << endl;		
			    		unsigned int midPartMapPosForLongHeadInSecondLevelIndex = unfixedHeadInfo.midPartMapPosInChr - 
							((unfixedHeadInfo.midPartMapPosInChr)/(indexInfo->secondLevelIndexNormalSize))*(indexInfo->secondLevelIndexNormalSize);

						//   check GTAG splice junctions
						//cout << "unfixedHeadInfo.possiGTAGpos.size(): " << (unfixedHeadInfo.possiGTAGpos).size() << endl;
						for(int tmp = 0; tmp < (unfixedHeadInfo.possiGTAGpos).size(); tmp++)
						{
							int tmpSJposInRead = unfixedHeadInfo.possiGTAGpos[tmp];

							if(tmpSJposInRead-1 < min_anchor_length)
								continue;

							string tmpShortAnchorStr = readSeqWithDirection.substr(0, tmpSJposInRead-1);
							string targetMappingStr = tmpShortAnchorStr + "GT";
							//cout << "headLength: " << tmpSJposInRead-1 << " targetStr: " << targetMappingStr << endl;
							char* headChar = const_cast<char*>(targetMappingStr.c_str());
							int targetMappingNum = 0;
							unsigned int targetMappingLoc[100];
						
							unsigned int finalMidPartMappingPos
								= unfixedHeadInfo.midPartMapPosInWholeGenome + tmpSJposInRead - unfixedHeadInfo.unfixedHeadLength - 1;
							bool headSegMapMain;
							
							if(!load2ndLevelIndexBool_compressedSize)
							{
								headSegMapMain = mapMainSecondLevelForTargetMapping(headChar, 
									secondLevelSa[midPartMapPosSecondLevelIndexNO], 
									secondLevelLcp[midPartMapPosSecondLevelIndexNO], 
									secondLevelUp[midPartMapPosSecondLevelIndexNO],
									secondLevelDown[midPartMapPosSecondLevelIndexNO],
									secondLevelNext[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedHeadInfo.midPartMapPosInWholeGenome, 
									midPartMapPosForLongHeadInSecondLevelIndex, &targetMappingNum, targetMappingLoc,
									indexInfo
									);
							}
							else
							{
								headSegMapMain = mapMainSecondLevelForTargetMapping_compressedIndex(headChar,
									secondLevelSa[midPartMapPosSecondLevelIndexNO], 
									secondLevelLcpCompress[midPartMapPosSecondLevelIndexNO], 
									secondLevelChildTab[midPartMapPosSecondLevelIndexNO],
									secondLevelDetChild[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedHeadInfo.midPartMapPosInWholeGenome, 
									midPartMapPosForLongHeadInSecondLevelIndex, &targetMappingNum, targetMappingLoc,
									indexInfo
									);
							}

							//cout << "headSegMapMain: " << headSegMapMain << endl;

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
							

							bool headSegMapMain;
							
							if(!load2ndLevelIndexBool_compressedSize)
							{
								headSegMapMain = mapMainSecondLevelForTargetMapping(headChar, 
									secondLevelSa[midPartMapPosSecondLevelIndexNO], 
									secondLevelLcp[midPartMapPosSecondLevelIndexNO], 
									secondLevelUp[midPartMapPosSecondLevelIndexNO],
									secondLevelDown[midPartMapPosSecondLevelIndexNO],
									secondLevelNext[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedHeadInfo.midPartMapPosInWholeGenome, 
									midPartMapPosForLongHeadInSecondLevelIndex, &targetMappingNum, targetMappingLoc,
									indexInfo
									);
							}
							else
							{
								headSegMapMain = mapMainSecondLevelForTargetMapping_compressedIndex(headChar,
									secondLevelSa[midPartMapPosSecondLevelIndexNO], 
									secondLevelLcpCompress[midPartMapPosSecondLevelIndexNO], 
									secondLevelChildTab[midPartMapPosSecondLevelIndexNO],
									secondLevelDetChild[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedHeadInfo.midPartMapPosInWholeGenome, 
									midPartMapPosForLongHeadInSecondLevelIndex, &targetMappingNum, targetMappingLoc,
									indexInfo
									);
							}


							
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
						//cout << "unfixedHeadInfo.GTAGsjPos.size(): " << (unfixedHeadInfo.GTAGsjPos).size() << endl;
						if((unfixedHeadInfo.GTAGsjPos).size() > 0)
						{
							Alignment_Info* newTmpAlignInfo = new Alignment_Info();
							newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
								(unfixedHeadInfo.GTAGsjPos[0]).first - 1, (unfixedHeadInfo.GTAGsjPos[0]).second, indexInfo);

							//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
							peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpAlignmentNO);


							for(int tmp = 1; tmp < (unfixedHeadInfo.GTAGsjPos).size(); tmp++)
							{
								Alignment_Info* newTmpAlignInfo 
									= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
									(unfixedHeadInfo.GTAGsjPos[tmp]).first - 1, (unfixedHeadInfo.GTAGsjPos[tmp]).second, indexInfo);
								//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
								peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
							}
							
							for(int tmp = 0; tmp < (unfixedHeadInfo.CTACsjPos).size(); tmp++)
							{
								Alignment_Info* newTmpAlignInfo 
									= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
									(unfixedHeadInfo.CTACsjPos[tmp]).first - 1, (unfixedHeadInfo.CTACsjPos[tmp]).second, indexInfo);
								//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
								peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
							}
						}
						else if((unfixedHeadInfo.CTACsjPos).size() > 0)
						{
							Alignment_Info* newTmpAlignInfo = new Alignment_Info();
							newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
								(unfixedHeadInfo.CTACsjPos[0]).first - 1, (unfixedHeadInfo.CTACsjPos[0]).second, indexInfo);
							
							//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
							peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpAlignmentNO);

							for(int tmp = 1; tmp < (unfixedHeadInfo.CTACsjPos).size(); tmp++)
							{
								Alignment_Info* newTmpAlignInfo 
									= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
									(unfixedHeadInfo.CTACsjPos[tmp]).first - 1, (unfixedHeadInfo.CTACsjPos[tmp]).second, indexInfo);
								//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
								peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
							}				
						}
						else
						{}

					}	
				}
				//cout << "start to fix tail " << endl;
				/////////////////// fix tail //////////////////////////////
				for(int tmpAlignInfoType = 1; tmpAlignInfoType <= 4; tmpAlignInfoType++)
				{
					//cout << "tmpAlignInfoType: " << tmpAlignInfoType << endl;
					if(tmpAlignInfoType <= 2)
					{
						readLength = readLength_1;
					}
					else
					{
						readLength = readLength_2;
					}
					//cout << "readLength: " << readLength << endl;
					for(int tmpAlignmentNO = 0; tmpAlignmentNO < (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType)); 
						tmpAlignmentNO++)
					{
						//cout << "tmpAlignmentNO: " << tmpAlignmentNO << endl;
						//tmpAlignmentInfo = (peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO];
						tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
						int cigarStringJumpCodeSize = (tmpAlignmentInfo->cigarStringJumpCode).size();
						if( (tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].type != "S" )
						{
							//cout << "no unfixedTail !" << endl; 
							continue;
						}

						//cout << "start to getUnfixedTailInfoFromRecordWithAlignInfoType !" << endl;
						Unfixed_Tail unfixedTailInfo;
						//unfixedTailInfo.getUnfixedTailInfoFromRecord(peReadInfo, true, tmpAlignmentInfo, indexInfo);
						unfixedTailInfo.getUnfixedTailInfoFromRecordWithAlignInfoType(peReadInfo, tmpAlignInfoType, tmpAlignmentInfo, indexInfo);
						//cout << "finish getUnfixedTailInfoFromRecordWithAlignInfoType !" << endl;
						string readSeqWithDirection;
						if(unfixedTailInfo.alignDirection == "+")
						{
							readSeqWithDirection = unfixedTailInfo.readSeqOriginal; // Do not need to copy because we already know the direction in the 4 cases
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
							/////////////////////////////////////////////////////////////////////////////
							////////////////////////     position hash     //////////////////////////////
							/////////////////////////////////////////////////////////////////////////////
							spliceJunctionFoundInHash 
								= unfixedTailInfo.SJsearchInSJhash(SJ, readSeqWithDirection, indexInfo);
							/////////////////////////////////////////////////////////////////////////////
							/////////////////////////     String hash     //////////////////////////////
							/////////////////////////////////////////////////////////////////////////////						
							//spliceJunctionFoundInHash 
							//	= unfixedTailInfo.SJsearchInSJhash_areaStringHash(SJ, readSeqWithDirection, indexInfo, SJ->areaSize);
						}
						else
						{							
							spliceJunctionFoundInHash = false;
						}

						//cout << "spliceJunctionFoundInHash: " << spliceJunctionFoundInHash << endl;
						if(spliceJunctionFoundInHash)
						{
							//if((unfixedTailInfo.SJposFromRemappingVec).size() == 1)
							int tmpDonerEndPosInRead = (unfixedTailInfo.SJposFromRemappingVec)[0].first;
							int tmpSpliceJunctionDistance = (unfixedTailInfo.SJposFromRemappingVec)[0].second;
							int tmpLastMatchLength = readLength - tmpDonerEndPosInRead;

							Alignment_Info* newTmpAlignInfo = new Alignment_Info();
							newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment(
								tmpLastMatchLength, tmpSpliceJunctionDistance, indexInfo);

							//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
							peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpAlignmentNO);
							for(int tmpSJposVec = 1; tmpSJposVec < (unfixedTailInfo.SJposFromRemappingVec).size(); tmpSJposVec++)
							{
								tmpDonerEndPosInRead = (unfixedTailInfo.SJposFromRemappingVec)[tmpSJposVec].first;
								tmpSpliceJunctionDistance = (unfixedTailInfo.SJposFromRemappingVec)[tmpSJposVec].second;
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
							
							bool tailSegMapMain;
							if(!load2ndLevelIndexBool_compressedSize)
							{
								tailSegMapMain = mapMainSecondLevelForTargetMapping(tailChar, 
									secondLevelSa[midPartMapPosSecondLevelIndexNO], 
									secondLevelLcp[midPartMapPosSecondLevelIndexNO], 
									secondLevelUp[midPartMapPosSecondLevelIndexNO],
									secondLevelDown[midPartMapPosSecondLevelIndexNO],
									secondLevelNext[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedTailInfo.midPartMapPosInWholeGenome, 
									midPartMapPosForLongTailInSecondLevelIndex, &targetMappingNum, targetMappingLoc
									, indexInfo);
							}
							else
							{
								tailSegMapMain = mapMainSecondLevelForTargetMapping_compressedIndex(tailChar, 
									secondLevelSa[midPartMapPosSecondLevelIndexNO],
									secondLevelLcpCompress[midPartMapPosSecondLevelIndexNO], 
									secondLevelChildTab[midPartMapPosSecondLevelIndexNO],
									secondLevelDetChild[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedTailInfo.midPartMapPosInWholeGenome, 
									midPartMapPosForLongTailInSecondLevelIndex, &targetMappingNum, targetMappingLoc
									, indexInfo);								
							}

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

							bool tailSegMapMain;
							if(!load2ndLevelIndexBool_compressedSize)
							{
								tailSegMapMain = mapMainSecondLevelForTargetMapping(tailChar, 
									secondLevelSa[midPartMapPosSecondLevelIndexNO], 
									secondLevelLcp[midPartMapPosSecondLevelIndexNO], 
									secondLevelUp[midPartMapPosSecondLevelIndexNO],
									secondLevelDown[midPartMapPosSecondLevelIndexNO],
									secondLevelNext[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedTailInfo.midPartMapPosInWholeGenome, 
									midPartMapPosForLongTailInSecondLevelIndex, &targetMappingNum, targetMappingLoc
									, indexInfo);
							}
							else
							{
								tailSegMapMain = mapMainSecondLevelForTargetMapping_compressedIndex(tailChar, 
									secondLevelSa[midPartMapPosSecondLevelIndexNO],
									secondLevelLcpCompress[midPartMapPosSecondLevelIndexNO], 
									secondLevelChildTab[midPartMapPosSecondLevelIndexNO],
									secondLevelDetChild[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedTailInfo.midPartMapPosInWholeGenome, 
									midPartMapPosForLongTailInSecondLevelIndex, &targetMappingNum, targetMappingLoc
									, indexInfo);								
							}


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

	}*/

	void fixHeadTail_areaAndStringHash(PE_Read_Info* peReadInfo, PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJ, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		vector<unsigned int*>& secondLevelLcp,
		vector<unsigned int*>& secondLevelUp,
		vector<unsigned int*>& secondLevelDown,
		vector<unsigned int*>& secondLevelNext, 
		bool load2ndLevelIndexBool_compressedSize,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo)
	{
		int readLength_1 = ((peReadInfo->readInfo_pe1).readSeq).length();
		int readLength_2 = ((peReadInfo->readInfo_pe2).readSeq).length();

		Alignment_Info* tmpAlignmentInfo;

		int readLength = 0;

				//////////////////// fix head //////////////////////////////
				//cout << "start to fix head" << endl;
				for(int tmpAlignInfoType = 1; tmpAlignInfoType <= 4; tmpAlignInfoType++)
				{
					//readLength = readLength_1 * (tmpAlignInfoType <= 2) + readLength_2 * (tmpAlignInfoType > 2);
					
					//cout << "tmpAlignInfoType: " << tmpAlignInfoType << endl;

					if(tmpAlignInfoType <= 2)
					{
						readLength = readLength_1;
					}
					else
					{
						readLength = readLength_2;
					}

					//cout << "readLength: " << readLength << endl;
					for(int tmpAlignmentNO = 0; tmpAlignmentNO < (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType)); 
						tmpAlignmentNO++)
					{
						//cout << "tmpAlignmentNO: " << tmpAlignmentNO << endl;
						//tmpAlignmentInfo = (peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO];
						tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
						if( (tmpAlignmentInfo->cigarStringJumpCode)[0].type != "S" )
						{
							//cout << "no unfixed head !" << endl;
							continue;
						}

						//cout << "start to nor1" << endl;
						//cout << "start to getUnfixedHeadInfoFromRecordWithAlignInfoType ..." << endl;
						Unfixed_Head unfixedHeadInfo;
						unfixedHeadInfo.getUnfixedHeadInfoFromRecordWithAlignInfoType(peReadInfo, tmpAlignInfoType, tmpAlignmentInfo, indexInfo);

						//cout << "finish getUnfixedTailInfoFromRecord ..." << endl;

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
							/////////////////////////////////////////////////////////////////////////////
							////////////////////////     position hash     //////////////////////////////
							/////////////////////////////////////////////////////////////////////////////
							//spliceJunctionFoundInHash 
							//	= unfixedHeadInfo.SJsearchInSJhash(SJ, readSeqWithDirection, indexInfo);
							/////////////////////////////////////////////////////////////////////////////
							////////////////////////      string hash      //////////////////////////////
							/////////////////////////////////////////////////////////////////////////////
							spliceJunctionFoundInHash
								= unfixedHeadInfo.SJsearchInSJhash_areaStringHash(SJ, readSeqWithDirection, indexInfo, SJ->areaSize);
						}
						else
						{
							spliceJunctionFoundInHash = false;
						}
						//cout << "spliceJunctionFoundInHash: " << spliceJunctionFoundInHash << endl;

						if(spliceJunctionFoundInHash)
						{

							int tmpDonerEndPosInRead = (unfixedHeadInfo.SJposFromRemappingVec)[0].first;
							int tmpSpliceJunctionDistance = (unfixedHeadInfo.SJposFromRemappingVec)[0].second;
							int tmpFirstMatchLength = tmpDonerEndPosInRead;

							int newMismatchToAddInHead = (unfixedHeadInfo.SJposFromRemappingVec_mismatch)[0];

							Alignment_Info* newTmpAlignInfo = new Alignment_Info();
							newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
								tmpFirstMatchLength, tmpSpliceJunctionDistance, indexInfo, newMismatchToAddInHead);

							//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
							peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpAlignmentNO);

							for(int tmpSJposVec = 1; tmpSJposVec < (unfixedHeadInfo.SJposFromRemappingVec).size(); tmpSJposVec++)
							{
								tmpDonerEndPosInRead = (unfixedHeadInfo.SJposFromRemappingVec)[tmpSJposVec].first;
								tmpSpliceJunctionDistance = (unfixedHeadInfo.SJposFromRemappingVec)[tmpSJposVec].second;
								tmpFirstMatchLength = tmpDonerEndPosInRead;	

								int tmpNewMismatchToAddInHead = (unfixedHeadInfo.SJposFromRemappingVec_mismatch)[tmpSJposVec];

								Alignment_Info* newTmpAlignInfo = new Alignment_Info();
								newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
									tmpFirstMatchLength, tmpSpliceJunctionDistance, indexInfo, tmpNewMismatchToAddInHead);
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

						//cout << "2ndLevel Index: " << midPartMapPosSecondLevelIndexNO << endl;		

						if(
							((indexInfo->invalidSecondLevelIndexNOset).find(midPartMapPosSecondLevelIndexNO+1))
							!= (indexInfo->invalidSecondLevelIndexNOset).end()
							)
						{
							//delete(seg2ndOriInfo);
							//delete(unmapEndInfo);
							continue;
						}

						//cout << "2ndLevel Index: " << midPartMapPosSecondLevelIndexNO << endl;		
			    		unsigned int midPartMapPosForLongHeadInSecondLevelIndex = unfixedHeadInfo.midPartMapPosInChr - 
							((unfixedHeadInfo.midPartMapPosInChr)/(indexInfo->secondLevelIndexNormalSize))*(indexInfo->secondLevelIndexNormalSize);

						//   check GTAG splice junctions
						//cout << "unfixedHeadInfo.possiGTAGpos.size(): " << (unfixedHeadInfo.possiGTAGpos).size() << endl;
						for(int tmp = 0; tmp < (unfixedHeadInfo.possiGTAGpos).size(); tmp++)
						{
							int tmpSJposInRead = unfixedHeadInfo.possiGTAGpos[tmp];

							if(tmpSJposInRead-1 < min_anchor_length)
								continue;

							string tmpShortAnchorStr = readSeqWithDirection.substr(0, tmpSJposInRead-1);
							string targetMappingStr = tmpShortAnchorStr + "GT";
							//cout << "headLength: " << tmpSJposInRead-1 << " targetStr: " << targetMappingStr << endl;
							char* headChar = const_cast<char*>(targetMappingStr.c_str());
							int targetMappingNum = 0;
							unsigned int targetMappingLoc[100];
						
							unsigned int finalMidPartMappingPos
								= unfixedHeadInfo.midPartMapPosInWholeGenome + tmpSJposInRead - unfixedHeadInfo.unfixedHeadLength - 1;
							bool headSegMapMain;
							
							if(!load2ndLevelIndexBool_compressedSize)
							{
								headSegMapMain = mapMainSecondLevelForTargetMapping(headChar, 
									secondLevelSa[midPartMapPosSecondLevelIndexNO], 
									secondLevelLcp[midPartMapPosSecondLevelIndexNO], 
									secondLevelUp[midPartMapPosSecondLevelIndexNO],
									secondLevelDown[midPartMapPosSecondLevelIndexNO],
									secondLevelNext[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedHeadInfo.midPartMapPosInWholeGenome, 
									midPartMapPosForLongHeadInSecondLevelIndex, &targetMappingNum, targetMappingLoc,
									indexInfo
									);
							}
							else
							{
								headSegMapMain = mapMainSecondLevelForTargetMapping_compressedIndex(headChar,
									secondLevelSa[midPartMapPosSecondLevelIndexNO], 
									secondLevelLcpCompress[midPartMapPosSecondLevelIndexNO], 
									secondLevelChildTab[midPartMapPosSecondLevelIndexNO],
									secondLevelDetChild[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedHeadInfo.midPartMapPosInWholeGenome, 
									midPartMapPosForLongHeadInSecondLevelIndex, &targetMappingNum, targetMappingLoc,
									indexInfo
									);
							}

							//cout << "headSegMapMain: " << headSegMapMain << endl;

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
									(unfixedHeadInfo.GTAGsjPos_mismatch).push_back(unfixedHeadInfo.possiGTAGpos_mismatch[tmp]);
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
							

							bool headSegMapMain;
							
							if(!load2ndLevelIndexBool_compressedSize)
							{
								headSegMapMain = mapMainSecondLevelForTargetMapping(headChar, 
									secondLevelSa[midPartMapPosSecondLevelIndexNO], 
									secondLevelLcp[midPartMapPosSecondLevelIndexNO], 
									secondLevelUp[midPartMapPosSecondLevelIndexNO],
									secondLevelDown[midPartMapPosSecondLevelIndexNO],
									secondLevelNext[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedHeadInfo.midPartMapPosInWholeGenome, 
									midPartMapPosForLongHeadInSecondLevelIndex, &targetMappingNum, targetMappingLoc,
									indexInfo
									);
							}
							else
							{
								headSegMapMain = mapMainSecondLevelForTargetMapping_compressedIndex(headChar,
									secondLevelSa[midPartMapPosSecondLevelIndexNO], 
									secondLevelLcpCompress[midPartMapPosSecondLevelIndexNO], 
									secondLevelChildTab[midPartMapPosSecondLevelIndexNO],
									secondLevelDetChild[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedHeadInfo.midPartMapPosInWholeGenome, 
									midPartMapPosForLongHeadInSecondLevelIndex, &targetMappingNum, targetMappingLoc,
									indexInfo
									);
							}


							
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
									(unfixedHeadInfo.CTACsjPos_mismatch).push_back(unfixedHeadInfo.possiCTACpos_mismatch[tmp]);
								}			
							}

						}
						//cout << "unfixedHeadInfo.GTAGsjPos.size(): " << (unfixedHeadInfo.GTAGsjPos).size() << endl;
						if((unfixedHeadInfo.GTAGsjPos).size() > 0)
						{
							Alignment_Info* newTmpAlignInfo = new Alignment_Info();

							int tmpNewMismatchToAddInHead = unfixedHeadInfo.GTAGsjPos_mismatch[0];

							newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
								(unfixedHeadInfo.GTAGsjPos[0]).first - 1, (unfixedHeadInfo.GTAGsjPos[0]).second, indexInfo, tmpNewMismatchToAddInHead);

							//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
							peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpAlignmentNO);


							for(int tmp = 1; tmp < (unfixedHeadInfo.GTAGsjPos).size(); tmp++)
							{
								int tmpNewMismatchToAddInHead = unfixedHeadInfo.GTAGsjPos_mismatch[tmp]; //0;

								Alignment_Info* newTmpAlignInfo 
									= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
									(unfixedHeadInfo.GTAGsjPos[tmp]).first - 1, (unfixedHeadInfo.GTAGsjPos[tmp]).second, indexInfo, tmpNewMismatchToAddInHead);
								//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
								peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
							}
							
							for(int tmp = 0; tmp < (unfixedHeadInfo.CTACsjPos).size(); tmp++)
							{
								int tmpNewMismatchToAddInHead = unfixedHeadInfo.CTACsjPos_mismatch[tmp]; //0;

								Alignment_Info* newTmpAlignInfo 
									= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
									(unfixedHeadInfo.CTACsjPos[tmp]).first - 1, (unfixedHeadInfo.CTACsjPos[tmp]).second, indexInfo, tmpNewMismatchToAddInHead);
								//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
								peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
							}
						}
						else if((unfixedHeadInfo.CTACsjPos).size() > 0)
						{
							Alignment_Info* newTmpAlignInfo = new Alignment_Info();

							int tmpNewMismatchToAddInHead = unfixedHeadInfo.CTACsjPos_mismatch[0]; //0;

							newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
								(unfixedHeadInfo.CTACsjPos[0]).first - 1, (unfixedHeadInfo.CTACsjPos[0]).second, indexInfo, tmpNewMismatchToAddInHead);
							
							//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
							peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpAlignmentNO);

							for(int tmp = 1; tmp < (unfixedHeadInfo.CTACsjPos).size(); tmp++)
							{
								int tmpNewMismatchToAddInHead = unfixedHeadInfo.CTACsjPos_mismatch[tmp];

								Alignment_Info* newTmpAlignInfo 
									= tmpAlignmentInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
									(unfixedHeadInfo.CTACsjPos[tmp]).first - 1, (unfixedHeadInfo.CTACsjPos[tmp]).second, indexInfo, tmpNewMismatchToAddInHead);
								//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
								peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
							}				
						}
						else
						{}

					}	
				}
				//cout << "start to fix tail " << endl;
				/////////////////// fix tail //////////////////////////////
				for(int tmpAlignInfoType = 1; tmpAlignInfoType <= 4; tmpAlignInfoType++)
				{
					//cout << "tmpAlignInfoType: " << tmpAlignInfoType << endl;
					if(tmpAlignInfoType <= 2)
					{
						readLength = readLength_1;
					}
					else
					{
						readLength = readLength_2;
					}
					//cout << "readLength: " << readLength << endl;
					for(int tmpAlignmentNO = 0; tmpAlignmentNO < (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType)); 
						tmpAlignmentNO++)
					{
						//cout << "tmpAlignmentNO: " << tmpAlignmentNO << endl;
						//tmpAlignmentInfo = (peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO];
						tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
						int cigarStringJumpCodeSize = (tmpAlignmentInfo->cigarStringJumpCode).size();
						if( (tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].type != "S" )
						{
							//cout << "no unfixedTail !" << endl; 
							continue;
						}

						//cout << "start to getUnfixedTailInfoFromRecordWithAlignInfoType !" << endl;
						Unfixed_Tail unfixedTailInfo;
						//unfixedTailInfo.getUnfixedTailInfoFromRecord(peReadInfo, true, tmpAlignmentInfo, indexInfo);
						unfixedTailInfo.getUnfixedTailInfoFromRecordWithAlignInfoType(peReadInfo, tmpAlignInfoType, tmpAlignmentInfo, indexInfo);
						//cout << "finish getUnfixedTailInfoFromRecordWithAlignInfoType !" << endl;
						string readSeqWithDirection;
						if(unfixedTailInfo.alignDirection == "+")
						{
							readSeqWithDirection = unfixedTailInfo.readSeqOriginal; // Do not need to copy because we already know the direction in the 4 cases
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
							/////////////////////////////////////////////////////////////////////////////
							////////////////////////     position hash     //////////////////////////////
							/////////////////////////////////////////////////////////////////////////////
							//spliceJunctionFoundInHash 
							//	= unfixedTailInfo.SJsearchInSJhash(SJ, readSeqWithDirection, indexInfo);
							/////////////////////////////////////////////////////////////////////////////
							/////////////////////////     String hash     //////////////////////////////
							/////////////////////////////////////////////////////////////////////////////						
							spliceJunctionFoundInHash 
								= unfixedTailInfo.SJsearchInSJhash_areaStringHash(SJ, readSeqWithDirection, indexInfo, SJ->areaSize);
						}
						else
						{							
							spliceJunctionFoundInHash = false;
						}

						//cout << "spliceJunctionFoundInHash: " << spliceJunctionFoundInHash << endl;
						if(spliceJunctionFoundInHash)
						{
							//if((unfixedTailInfo.SJposFromRemapping).size() == 1)
							int tmpDonerEndPosInRead = (unfixedTailInfo.SJposFromRemappingVec)[0].first;
							int tmpSpliceJunctionDistance = (unfixedTailInfo.SJposFromRemappingVec)[0].second;
							int tmpLastMatchLength = readLength - tmpDonerEndPosInRead;

							int tmpNewMismatchToAddInTail = (unfixedTailInfo.SJposFromRemappingVec_mismatch)[0];

							Alignment_Info* newTmpAlignInfo = new Alignment_Info();
							newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment(
								tmpLastMatchLength, tmpSpliceJunctionDistance, indexInfo, tmpNewMismatchToAddInTail);

							//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
							peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpAlignmentNO);
							for(int tmpSJposVec = 1; tmpSJposVec < (unfixedTailInfo.SJposFromRemappingVec).size(); tmpSJposVec++)
							{
								int tmpNewMismatchToAddInTail = (unfixedTailInfo.SJposFromRemappingVec_mismatch)[tmpSJposVec];

								tmpDonerEndPosInRead = (unfixedTailInfo.SJposFromRemappingVec)[tmpSJposVec].first;
								tmpSpliceJunctionDistance = (unfixedTailInfo.SJposFromRemappingVec)[tmpSJposVec].second;
								tmpLastMatchLength = readLength - tmpDonerEndPosInRead;		

								Alignment_Info* newTmpAlignInfo = new Alignment_Info();
								newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment(
									tmpLastMatchLength, tmpSpliceJunctionDistance, indexInfo, tmpNewMismatchToAddInTail);
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
							
							bool tailSegMapMain;
							if(!load2ndLevelIndexBool_compressedSize)
							{
								tailSegMapMain = mapMainSecondLevelForTargetMapping(tailChar, 
									secondLevelSa[midPartMapPosSecondLevelIndexNO], 
									secondLevelLcp[midPartMapPosSecondLevelIndexNO], 
									secondLevelUp[midPartMapPosSecondLevelIndexNO],
									secondLevelDown[midPartMapPosSecondLevelIndexNO],
									secondLevelNext[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedTailInfo.midPartMapPosInWholeGenome, 
									midPartMapPosForLongTailInSecondLevelIndex, &targetMappingNum, targetMappingLoc
									, indexInfo);
							}
							else
							{
								tailSegMapMain = mapMainSecondLevelForTargetMapping_compressedIndex(tailChar, 
									secondLevelSa[midPartMapPosSecondLevelIndexNO],
									secondLevelLcpCompress[midPartMapPosSecondLevelIndexNO], 
									secondLevelChildTab[midPartMapPosSecondLevelIndexNO],
									secondLevelDetChild[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedTailInfo.midPartMapPosInWholeGenome, 
									midPartMapPosForLongTailInSecondLevelIndex, &targetMappingNum, targetMappingLoc
									, indexInfo);								
							}

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
									(unfixedTailInfo.GTAGsjPos_mismatch).push_back(unfixedTailInfo.possiGTAGpos_mismatch[tmp]);
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

							bool tailSegMapMain;
							if(!load2ndLevelIndexBool_compressedSize)
							{
								tailSegMapMain = mapMainSecondLevelForTargetMapping(tailChar, 
									secondLevelSa[midPartMapPosSecondLevelIndexNO], 
									secondLevelLcp[midPartMapPosSecondLevelIndexNO], 
									secondLevelUp[midPartMapPosSecondLevelIndexNO],
									secondLevelDown[midPartMapPosSecondLevelIndexNO],
									secondLevelNext[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedTailInfo.midPartMapPosInWholeGenome, 
									midPartMapPosForLongTailInSecondLevelIndex, &targetMappingNum, targetMappingLoc
									, indexInfo);
							}
							else
							{
								tailSegMapMain = mapMainSecondLevelForTargetMapping_compressedIndex(tailChar, 
									secondLevelSa[midPartMapPosSecondLevelIndexNO],
									secondLevelLcpCompress[midPartMapPosSecondLevelIndexNO], 
									secondLevelChildTab[midPartMapPosSecondLevelIndexNO],
									secondLevelDetChild[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedTailInfo.midPartMapPosInWholeGenome, 
									midPartMapPosForLongTailInSecondLevelIndex, &targetMappingNum, targetMappingLoc
									, indexInfo);								
							}


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
									(unfixedTailInfo.CTACsjPos_mismatch).push_back(unfixedTailInfo.possiCTACpos_mismatch[tmp]);
								}
							}
						}

						if((unfixedTailInfo.GTAGsjPos).size() > 0)
						{
							int tmpNewMismatchToAddInTail = (unfixedTailInfo.GTAGsjPos_mismatch)[0];

							Alignment_Info* newTmpAlignInfo = new Alignment_Info();
							newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment(
								readLength - (unfixedTailInfo.GTAGsjPos[0]).first + 1, (unfixedTailInfo.GTAGsjPos[0]).second, indexInfo, 
									tmpNewMismatchToAddInTail);

							//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
							peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpAlignmentNO);

							for(int tmp = 1; tmp < (unfixedTailInfo.GTAGsjPos).size(); tmp++)
							{
								int tmpNewMismatchToAddInTail = unfixedTailInfo.GTAGsjPos_mismatch[tmp];

								Alignment_Info* newTmpAlignInfo 
									= tmpAlignmentInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment(
									readLength - (unfixedTailInfo.GTAGsjPos[tmp]).first + 1, (unfixedTailInfo.GTAGsjPos[tmp]).second, indexInfo, 
										tmpNewMismatchToAddInTail);
								//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
								peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
							}
							
							for(int tmp = 0; tmp < (unfixedTailInfo.CTACsjPos).size(); tmp++)
							{
								int tmpNewMismatchToAddInTail = unfixedTailInfo.CTACsjPos_mismatch[tmp];//0;

								Alignment_Info* newTmpAlignInfo 
									= tmpAlignmentInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment(
									readLength - (unfixedTailInfo.CTACsjPos[tmp]).first + 1, (unfixedTailInfo.CTACsjPos[tmp]).second, indexInfo,
										tmpNewMismatchToAddInTail);
								//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
								peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
							}
						}
						else if((unfixedTailInfo.CTACsjPos).size() > 0)
						{
							int tmpNewMismatchToAddInTail = unfixedTailInfo.CTACsjPos_mismatch[0];

							Alignment_Info* newTmpAlignInfo = new Alignment_Info();
							newTmpAlignInfo	= tmpAlignmentInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment(
								readLength - (unfixedTailInfo.CTACsjPos[0]).first + 1, (unfixedTailInfo.CTACsjPos[0]).second, indexInfo,
									tmpNewMismatchToAddInTail);
							
							//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
							peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpAlignmentNO);

							for(int tmp = 1; tmp < (unfixedTailInfo.CTACsjPos).size(); tmp++)
							{
								int tmpNewMismatchToAddInTail = unfixedTailInfo.CTACsjPos_mismatch[tmp]; //0;

								Alignment_Info* newTmpAlignInfo 
									= tmpAlignmentInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment(
									readLength - (unfixedTailInfo.CTACsjPos[tmp]).first + 1, (unfixedTailInfo.CTACsjPos[tmp]).second, indexInfo,
										tmpNewMismatchToAddInTail);
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