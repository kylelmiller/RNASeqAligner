#include <string>
#include <string.h>
#include "splice_info.h"

using namespace std;

class FixOneEndUnmappedInfo
{
public:
	FixOneEndUnmappedInfo()
	{}

	void fixOneEndUnmapped(PE_Read_Info* peReadInfo, PE_Read_Alignment_Info* peAlignInfo,
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
		Index_Info* indexInfo
		)
	{
				int readLength;// = 100; // to debug

				bool End1OrEnd2; bool NorOrRcm;

				vector<int> alignmentInfoVecSize_ori;
				alignmentInfoVecSize_ori.push_back(peAlignInfo->norAlignmentInfo_PE_1.size());
				alignmentInfoVecSize_ori.push_back(peAlignInfo->rcmAlignmentInfo_PE_1.size());
				alignmentInfoVecSize_ori.push_back(peAlignInfo->norAlignmentInfo_PE_2.size());
				alignmentInfoVecSize_ori.push_back(peAlignInfo->rcmAlignmentInfo_PE_2.size());

					for(int tmpAlignInfoType = 1; tmpAlignInfoType <= 4; tmpAlignInfoType ++)
					{
						End1OrEnd2 = peReadInfo->checkEnd1OrEnd2WithAlignInfoTypeNo(tmpAlignInfoType);
						NorOrRcm = peReadInfo->checkNorOrRcmWithAlignInfoTypeNo(tmpAlignInfoType);
						readLength = peReadInfo->checkReadLengthWithAlignInfoTypeNo(tmpAlignInfoType);

						int peAlignInfoVectorSize = (peAlignInfo->getAlignInfoVecSize(tmpAlignInfoType));

						for(int tmpAlignmentNO = 0; 
							tmpAlignmentNO < peAlignInfoVectorSize; 
							tmpAlignmentNO ++)
						{

							//cout << "tmpAlignemntNO: " << tmpAlignmentNO << endl;

							Alignment_Info* tmpAlignInfo = peAlignInfo->getAlignInfo(tmpAlignInfoType,
								tmpAlignmentNO);
							UnmapEnd_Info* unmapEndInfo = new UnmapEnd_Info();
							unmapEndInfo->setUnmapEndInfo(tmpAlignInfo, End1OrEnd2, NorOrRcm, indexInfo);

							char* unmapEndReadChar =  const_cast<char*>(peReadInfo->getIncompleteEndReadSeq(End1OrEnd2, NorOrRcm).c_str());

							Seg2ndOri_Info* seg2ndOriInfo = new Seg2ndOri_Info();

							//cout << "start to map in 2nd level index " << endl;
							int secondLevelIndexNO = unmapEndInfo->secondLevelIndexNum - 1;
							//cout << "secondLevelIndexNO: " << secondLevelIndexNO << endl;
							
							if(
								((indexInfo->invalidSecondLevelIndexNOset).find(secondLevelIndexNO+1))
								!= (indexInfo->invalidSecondLevelIndexNOset).end()
								)
							{
								delete(seg2ndOriInfo);
								delete(unmapEndInfo);
								continue;
							}

							bool unmapEndMapBool;
							//cout << "load2ndLevelIndexBool_compressedSize: " << load2ndLevelIndexBool_compressedSize << endl;
							if(!load2ndLevelIndexBool_compressedSize)
							{
								//cout << "start to mapMainSecondLevel !" << endl;
								unmapEndMapBool = seg2ndOriInfo->mapMainSecondLevel(
									unmapEndReadChar, 
									secondLevelSa[secondLevelIndexNO], 
									secondLevelLcp[secondLevelIndexNO], 
									secondLevelUp[secondLevelIndexNO], 
									secondLevelDown[secondLevelIndexNO], 
									secondLevelNext[secondLevelIndexNO], 
									secondLevelChrom[secondLevelIndexNO], 
									readLength, indexInfo->secondLevelIndexNormalSize + 1, indexInfo);
								//cout << "finish doing mapMainSecondLevel !" << endl;
								//cout << "unmapEndMapBool: " << unmapEndMapBool << endl;
							}
							else
							{
								//cout << "start to mapMainSecondLevel_compressedIndex !" << endl;
								unmapEndMapBool = seg2ndOriInfo->mapMainSecondLevel_compressedIndex(
									unmapEndReadChar,
									secondLevelSa[secondLevelIndexNO], 
									secondLevelLcpCompress[secondLevelIndexNO],
									secondLevelChildTab[secondLevelIndexNO],
									secondLevelChrom[secondLevelIndexNO], 
									secondLevelDetChild[secondLevelIndexNO],
									readLength, indexInfo);
								//cout << "finish doing mapMainSecondLevel_compressedIndex !" << endl;
								//cout << "unmapEndMapBool: " << unmapEndMapBool << endl;
							}
							//cout << seg2ndOriInfo->segInfoStr(indexInfo, unmapEndInfo->chrPosStartIn2ndLevelIndex, unmapEndInfo->chrNameStr) << endl;

							if(!unmapEndMapBool)
							{
								delete(seg2ndOriInfo);
								delete(unmapEndInfo);
								continue;								
							}

							Seg_Info* segInfo = new Seg_Info(seg2ndOriInfo, unmapEndInfo->mapPosIntervalStart,
								unmapEndInfo->mapPosIntervalEnd, unmapEndInfo->chrPosStartIn2ndLevelIndex,
								indexInfo, unmapEndInfo->chrNameStr);

							//cout << segInfo->segInfoStr(indexInfo) << endl;

							Path_Info* pathInfo = new Path_Info();
							pathInfo->getPossiPathFromSeg(segInfo);

							//cout << pathInfo->possiPathStr() << endl;
							int pathValidNum = pathInfo->pathValidNumInt();
							if(pathValidNum > 10)
							{
								pathInfo->memoryFree();
								delete(pathInfo);
								delete(segInfo);
								delete(seg2ndOriInfo);
								delete(unmapEndInfo);								
								continue;
							}

							Gap_Info* gapInfo = new Gap_Info();
							gapInfo->fixGapInPath(pathInfo, segInfo, 
								indexInfo, peReadInfo->getIncompleteEndReadSeq(End1OrEnd2, NorOrRcm), readLength);
							
							//cout << pathInfo->fixedPathVecStr(indexInfo, segInfo) << endl;
							//cout << pathInfo->finalFixedPathStr(indexInfo) << endl;
							if(pathInfo->finalPathVec.size() == 1)
							{
								peAlignInfo->pushBackPathInfo2PeAlignInfo(pathInfo, End1OrEnd2, NorOrRcm, indexInfo);
							}

							peAlignInfo->pushBackPathInfo2PeAlignInfo(pathInfo, End1OrEnd2, NorOrRcm, indexInfo);

							//if(pathInfo->finalPathVec.size() > 0)
							//{
							//	inCompleteEndFixedBool = true;
							//}

							delete(gapInfo);
							pathInfo->memoryFree();
							delete(pathInfo);
							delete(segInfo);
							delete(seg2ndOriInfo);
							delete(unmapEndInfo);
						}
					}
	}
};