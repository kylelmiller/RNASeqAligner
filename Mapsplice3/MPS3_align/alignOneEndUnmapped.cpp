#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>

//#include "switch.h"
#include "read_block_test.h"
#include "bwtmap_info.h"
#include "DoubleAnchorScore.h"
#include "sbndm.h"
#include "otherFunc.h"
#include "index_info.h"
#include "constantDefinitions.h"
#include "segmentMapping.h"
//#include "segmentMapping_secondLevel.h"
#include "splice_info.h"
#include "fixGapRelationParameters.h"
#include "read_info.h"
#include "seg_info.h"
#include "gap_info.h"
#include "align_info.h"
#include "spliceJunction_info.h"
#include "unmapEnd_info.h"
#include "unfixedHead.h"
#include "unfixedTail.h"
#include "sam2junc.h"
#include "fixHeadTail.h"
#include "fixOneEndUnmapped.h"
#include "fixPhase1.h"

#define PreIndexSize 268435456

using namespace std;  

////////////////////////////////////////////////////////////////////////////////////////
/////////////////Compare  ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////


int main(int argc, char**argv)
{
    if(argc < 5)
	{
		cout << "Executable <InputRecords> <OutputDir> <SpliceJunction> <threads_num> <wholeGenomeIndex> <2ndLevelIndex>" << endl;
		exit(0);
	}

	clock_t overAll_begin, overAll_end, overAll_cost, 
			loadIndex_begin, loadIndex_end, loadIndex_cost, 
			input_begin, input_end, input_cost,
			fixHead_begin, fixHead_end, fixHead_cost,
			output_begin, output_end, output_cost;

	overAll_begin = clock();
	time_t nowtime;
	nowtime = time(NULL);
	struct tm *local;
	///////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////   switches of seperate processes    ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
	string outputDirStr = argv[2];
   	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	system(mkdirOutputCommand.c_str());
	string testClassFileStr = outputDirStr + "/output.sam";

	string processLogStr = outputDirStr + "/process.log";
   	ofstream log_ofs(processLogStr.c_str());

	bool Do_Phase1_Only = false;
	//Do_Phase1_Only = true;

	bool outputAlignInfoAndSamForAllPairedAlignmentBool = false;
	outputAlignInfoAndSamForAllPairedAlignmentBool = true;

	//bool outputAllSamFileInEveryStep = false;
	//outputAllSamFileInEveryStep = true;

	bool DoSam2JuncBool = false;
	DoSam2JuncBool = true;

	bool load2ndLevelIndexBool = false;
	load2ndLevelIndexBool = true;

	bool load2ndLevelIndexBool_compressedSize = false;
	load2ndLevelIndexBool_compressedSize = true;

	bool DoRemappingOnUnmapEndReadsBool = false;
	DoRemappingOnUnmapEndReadsBool = true;

	bool DoRemappingOnUnfixedHeadTailAlignmentBool = false;
	DoRemappingOnUnfixedHeadTailAlignmentBool = true;



	if(Do_Phase1_Only)
	{
		log_ofs << "Do_Phase1 only!" << endl;
		DoSam2JuncBool = false;
		load2ndLevelIndexBool = false;
		load2ndLevelIndexBool_compressedSize = false;
		DoRemappingOnUnmapEndReadsBool = false;
		DoRemappingOnUnfixedHeadTailAlignmentBool = false;
	}	
	else
	{
		log_ofs << "Do_Phase1_Phase2! " << endl;
		DoSam2JuncBool = true;//false;
		load2ndLevelIndexBool = true;//false;
		load2ndLevelIndexBool_compressedSize = true;//false;
		DoRemappingOnUnmapEndReadsBool = true;//false;
		DoRemappingOnUnfixedHeadTailAlignmentBool = true;//false;
	}

	int normalRecordNum_1stMapping = 1000000;
	int normalRecordNum_fixOneEndUnmapped = 1;//000000;
	int normalRecordNum_fixHeadTail = 1000000;//000000;

	log_ofs << "normalRecordNum_1stMapping: " << normalRecordNum_1stMapping << endl;
	log_ofs << "normalRecordNum_fixOneEndUnmapped: " << normalRecordNum_fixOneEndUnmapped << endl;
	log_ofs << "normalRecordNum_fixHeadTail: " << normalRecordNum_fixHeadTail << endl;


	////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////

	log_ofs << "inputRecords: " << argv[1] << endl;
	log_ofs << "outputSam: " << argv[2] << endl;
	log_ofs << "spliceJunction: " << argv[3] << endl;
	log_ofs << "thread: " << argv[4] << endl;
	//cout << "species: " << argv[5] << endl;

	string threadsNumStr = argv[4];
	//string speciesStr = argv[5];

	int threads_num = atoi(threadsNumStr.c_str());

	omp_set_num_threads(threads_num);
//////////////////////////////////////////////////////////              LOAD INDEX           ////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////              LOAD INDEX           ////////////////////////////////////////////////////////////////

    string headTailSoftClippingFile = argv[1];

    //string OutputSamFile = argv[2];
    //ofstream OutputSamFile_ofs(OutputSamFile.c_str());

    string tmpAlignOneEndUnmapped = headTailSoftClippingFile;//testClassFileStr + ".oneEndUnmapped.alignInfo";
	//ofstream tmpAlignOneEndUnmapped_ofs(tmpAlignOneEndUnmapped.c_str());

	loadIndex_begin = clock();

	string secondLevelIndexStr = argv[6];
	string wholeGenomeIndex = argv[5];
	wholeGenomeIndex += "/";
	string chrom_bit_file = wholeGenomeIndex + "_chrom";
	string parameter_file = wholeGenomeIndex + "_parameter";


    ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
    ifstream parameter_file_ifs(parameter_file.c_str());

    Index_Info* indexInfo = new Index_Info(parameter_file_ifs);
	log_ofs << "indexInfo->MAX: " << indexInfo->indexSize << endl;

	log_ofs << "start to load whole genome" << endl;
	char *chrom;

	chrom = (char*)malloc((indexInfo->indexSize) * sizeof(char));
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->indexSize) * sizeof(char)); 

	indexInfo->chromString = chrom;
	log_ofs << "chromSize = " <<(indexInfo->chromString).size() << endl;
	
	log_ofs << "start to load every chromosome" << endl;
	//log_ofs << "start to load every chromosome" << endl;
	//chromStr[0] = 
	(indexInfo->chromStr).push_back((indexInfo->chromString).substr(0, (indexInfo->chrEndPosInGenome)[0]+1));
	(indexInfo->chromLength).push_back(((indexInfo->chrEndPosInGenome)[0]+1));
	for(int tmp = 1; tmp < indexInfo->chromNum; tmp++)
	{
		//chromStr[tmp] = 
		(indexInfo->chromStr).push_back((indexInfo->chromString).substr((indexInfo->chrEndPosInGenome)[tmp-1]+2, 
			(indexInfo->chrEndPosInGenome)[tmp]-(indexInfo->chrEndPosInGenome)[tmp-1]-1));	
		(indexInfo->chromLength).push_back(((indexInfo->chrEndPosInGenome)[tmp]-(indexInfo->chrEndPosInGenome)[tmp-1]-1));
	}	

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////    	Load Second Level Index      ////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load 2nd level index ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... load 2nd level index starts ......" << endl << endl; 

	vector<char*> secondLevelChrom;
	vector<unsigned int*> secondLevelSa;

	vector<BYTE*> secondLevelLcpCompress;
	vector<unsigned int*> secondLevelChildTab;
	vector<BYTE*> secondLevelDetChild;

	vector<unsigned int*> secondLevelLcp;
	vector<unsigned int*> secondLevelUp;
	vector<unsigned int*> secondLevelDown;
	vector<unsigned int*> secondLevelNext;

	if(load2ndLevelIndexBool)
	{
		log_ofs << "start to load second-level index ..." << endl;
		
		int secondLevelIndexNO = 0;
		for(int tmpChrNO = 0; tmpChrNO < indexInfo->chromNum; tmpChrNO ++)
		{
			for(int tmpSecondLevelIndexNO = 1; tmpSecondLevelIndexNO <= (indexInfo->secondLevelIndexPartsNum)[tmpChrNO]; tmpSecondLevelIndexNO ++)
			{
				char tmpFileNumChar[4];
				sprintf(tmpFileNumChar, "%d", tmpSecondLevelIndexNO);
				string tmpFileNumStr = tmpFileNumChar;
				
				string inputIndexFileStr = secondLevelIndexStr + "/" + indexInfo->chrNameStr[tmpChrNO] + "/" 
					//+ indexInfo->chrNameStr[tmpChrNO] + "_part."
					+ tmpFileNumStr + "/";//"." + "test3_";


				string secondLevelIndexFileChromStr = inputIndexFileStr + "chrom"; 
				ifstream secondLevelChrom_file_ifs(secondLevelIndexFileChromStr.c_str(), ios::binary);
				string secondLevelIndexFileSaStr = inputIndexFileStr + "SA";
				ifstream secondLevelSA_file_ifs(secondLevelIndexFileSaStr.c_str(), ios::binary);

				//if(!load2ndLevelIndexBool_compressedSize)
				//{
				
					string secondLevelIndexFileLcpStr = inputIndexFileStr + "lcp";	
					ifstream secondLevelLcp_file_ifs(secondLevelIndexFileLcpStr.c_str(), ios::binary);	
					string secondLevelIndexFileUpStr = inputIndexFileStr + "up";	
					ifstream secondLevelUp_file_ifs(secondLevelIndexFileUpStr.c_str(), ios::binary);
					string secondLevelIndexFileDownStr = inputIndexFileStr + "down";	
					ifstream secondLevelDown_file_ifs(secondLevelIndexFileDownStr.c_str(), ios::binary);
					string secondLevelIndexFileNextStr = inputIndexFileStr + "next";	
					ifstream secondLevelNext_file_ifs(secondLevelIndexFileNextStr.c_str(), ios::binary);
				
				//}
				//else
				//{
					string secondLevelIndexFileLcpCompressStr = inputIndexFileStr + "_lcpCompress";	
					ifstream secondLevelLcpCompress_file_ifs(secondLevelIndexFileLcpCompressStr.c_str(), ios::binary);	
					string secondLevelIndexFileChildTabStr = inputIndexFileStr + "childTab";	
					ifstream secondLevelChildTab_file_ifs(secondLevelIndexFileChildTabStr.c_str(), ios::binary);
					string secondLevelIndexFileDetChildStr = inputIndexFileStr + "detChild";	
					ifstream secondLevelDetChild_file_ifs(secondLevelIndexFileDetChildStr.c_str(), ios::binary);					
				//}

				int sizeOfIndex = indexInfo->secondLevelIndexNormalSize + 1;
				char* tmpSecondLevelChrom = (char*)malloc(sizeOfIndex * sizeof(char));
				for(int tmpMallocSpace = 0; tmpMallocSpace < sizeOfIndex; tmpMallocSpace++)
				{
					tmpSecondLevelChrom[tmpMallocSpace] = '0';
				}
				secondLevelChrom_file_ifs.read((char*)tmpSecondLevelChrom, sizeOfIndex * sizeof(char));
				if(tmpSecondLevelChrom[sizeOfIndex-1] != 'X')
				{
					(indexInfo->invalidSecondLevelIndexNOset).insert(secondLevelIndexNO + 1);
				}

				bool No_ATGC_Bool = true;
				for(int tmpMallocSpace = 0; tmpMallocSpace < sizeOfIndex; tmpMallocSpace++)
				{
					char ch = tmpSecondLevelChrom[tmpMallocSpace];
					if((ch == 'A')||(ch == 'T')||(ch == 'G')||(ch == 'C'))
					{
						No_ATGC_Bool = false;
						break;
					}
				}				
				if(No_ATGC_Bool)
				{
					(indexInfo->invalidSecondLevelIndexNOset).insert(secondLevelIndexNO + 1);
				}	

				secondLevelChrom.push_back(tmpSecondLevelChrom);
				
				unsigned int* tmpSecondLevelSa = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
				secondLevelSA_file_ifs.read((char*)tmpSecondLevelSa, sizeOfIndex * sizeof(unsigned int));
				secondLevelSa.push_back(tmpSecondLevelSa);

				if(!load2ndLevelIndexBool_compressedSize)
				{
					unsigned int* tmpSecondLevelLcp = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
					secondLevelLcp_file_ifs.read((char*)tmpSecondLevelLcp, sizeOfIndex * sizeof(unsigned int));
					secondLevelLcp.push_back(tmpSecondLevelLcp);
					
					unsigned int* tmpSecondLevelUp = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
					secondLevelUp_file_ifs.read((char*)tmpSecondLevelUp, sizeOfIndex * sizeof(unsigned int));
					secondLevelUp.push_back(tmpSecondLevelUp);
					
					unsigned int* tmpSecondLevelDown = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
					secondLevelDown_file_ifs.read((char*)tmpSecondLevelDown, sizeOfIndex * sizeof(unsigned int));
					secondLevelDown.push_back(tmpSecondLevelDown);
					
					unsigned int* tmpSecondLevelNext = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
					secondLevelNext_file_ifs.read((char*)tmpSecondLevelNext, sizeOfIndex * sizeof(unsigned int));
					secondLevelNext.push_back(tmpSecondLevelNext);
				}
				else
				{
					BYTE* tmpSecondLevelLcpCompress = (BYTE*)malloc(sizeOfIndex * sizeof(BYTE));
					secondLevelLcpCompress_file_ifs.read((char*)tmpSecondLevelLcpCompress, sizeOfIndex * sizeof(BYTE));
					secondLevelLcpCompress.push_back(tmpSecondLevelLcpCompress);
					
					unsigned int* tmpSecondLevelChildTab = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
					secondLevelChildTab_file_ifs.read((char*)tmpSecondLevelChildTab, sizeOfIndex * sizeof(unsigned int));
					secondLevelChildTab.push_back(tmpSecondLevelChildTab);

					BYTE* tmpSecondLevelDetChild = (BYTE*)malloc(sizeOfIndex * sizeof(BYTE));
					secondLevelDetChild_file_ifs.read((char*)tmpSecondLevelDetChild, sizeOfIndex * sizeof(BYTE));
					secondLevelDetChild.push_back(tmpSecondLevelDetChild);
				}


				secondLevelChrom_file_ifs.close();
				secondLevelSA_file_ifs.close();
				//if(!load2ndLevelIndexBool_compressedSize)
				//{
					secondLevelLcp_file_ifs.close();
					secondLevelUp_file_ifs.close();
					secondLevelDown_file_ifs.close();
					secondLevelNext_file_ifs.close();
				//}
				//else
				//{
					secondLevelLcpCompress_file_ifs.close();
					secondLevelChildTab_file_ifs.close();
					secondLevelDetChild_file_ifs.close();					
				//}				

				secondLevelIndexNO ++;

			}
			log_ofs << "finish loading 2nd-level index of " << indexInfo->chrNameStr[tmpChrNO] << endl; 
		}
		log_ofs << "finish loading ALL 2nd-level index !" << endl;
		log_ofs << indexInfo->getInvalidSecondLevelIndexNOstr() << endl;
		//loadIndex_end = clock(); 
	}

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... load 2nd level index ends ......" << endl << endl ; 		
	log_ofs << endl << "[" << asctime(local) << "... load 2nd level index ends ......" << endl << endl ; 	

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////    Do REMAPPING On one end unmapped Reads   //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads starts ......" << endl << endl ; 	 
	log_ofs << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads starts ......" << endl << endl ; 

	string OutputSamFile_oneEndMapped = testClassFileStr + ".oneEndUnmapped.pairedComplete.sam";
	ofstream OutputSamFile_oneEndMapped_ofs(OutputSamFile_oneEndMapped.c_str());

	string OutputSamFile_oneEndMapped_unpairComplete = testClassFileStr + ".oneEndUnmapped.unpairedComplete.sam";
	ofstream OutputSamFile_oneEndMapped_unpairComplete_ofs(OutputSamFile_oneEndMapped_unpairComplete.c_str());

	string OutputSamFile_oneEndMapped_alignInfo = testClassFileStr + ".oneEndUnmapped.pairedComplete.sam_alignInfo";
	ofstream OutputSamFile_oneEndMapped_alignInfo_ofs(OutputSamFile_oneEndMapped_alignInfo.c_str());	
	
	string tmpAlignIncompletePair = testClassFileStr + ".oneEndUnmapped.incompletePair.sam_alignInfo";
	ofstream tmpAlignIncompletePair_ofs(tmpAlignIncompletePair.c_str());
	input_begin = clock();

	unsigned int inputRecordNum = 0;
	unsigned int completeEndNum = 0;
	unsigned int incompleteEndNum = 0;
	unsigned int incompleteEndFixedNum = 0;
	unsigned int incompleteEndUnFixedNum = 0;
	unsigned int bothEndUnmapNum = 0;
	cout << "start doing remapping on long-Head reads" << endl;

	//bool DoRemappingOnUnmapEndReadsBool = false;
	//DoRemappingOnUnmapEndReadsBool = true;


	if(DoRemappingOnUnmapEndReadsBool)
	{
		log_ofs << "start doing remapping on unmapped end reads" << endl;

		string oneEndMappedFileStr = tmpAlignOneEndUnmapped;

		FILE *fp_HeadTail = fopen(oneEndMappedFileStr.c_str(), "r");

		char line1Char[200], line2Char[500], line3Char[500], line4Char[200], line5Char[500],
			 line6Char[500], line7Char[2000], line8Char[2000], line9Char[2000], line10Char[2000],
			 line11Char[200], readNameChar_1[100], readNameChar_2[100], othersChar[100];

		string line1, line2, line3, line4, line5, line6, line7, 
			line8, line9, line10, line11;
		
		int normalRecordNum = normalRecordNum_fixOneEndUnmapped; //1000000;

		bool EndOfRecord = false;

		int tmpTurn = 0;

		int realRecordNum;// = normalRecordNum;

		char line0Char[200];
		fgets(line0Char, sizeof(line0Char), fp_HeadTail);

		vector<string> line1StrVec(normalRecordNum);
		vector<string> line2StrVec(normalRecordNum);
		vector<string> line3StrVec(normalRecordNum);
		vector<string> line4StrVec(normalRecordNum);
		vector<string> line5StrVec(normalRecordNum);
		vector<string> line6StrVec(normalRecordNum);
		vector<string> line7StrVec(normalRecordNum);
		vector<string> line8StrVec(normalRecordNum);
		vector<string> line9StrVec(normalRecordNum);
		vector<string> line10StrVec(normalRecordNum);
		vector<string> peAlignInfoVec(normalRecordNum);
		vector<string> peAlignSamVec(normalRecordNum);
		vector<string> peAlignSamVec_unpair_complete(normalRecordNum);

		vector<string> peAlignInfoVec_pair_complete(normalRecordNum);

		int oneEndUnmappedRecordNum = 0;

		for(tmpTurn = 0; /*tmpTurn < TurnNum*/; tmpTurn++)
		{
			if(EndOfRecord)
				break;

			int recordNum = normalRecordNum;

			//cout << "start to read record" << endl;
			//cout << "start to read record, turn: " << tmpTurn+1 << endl;
			realRecordNum = normalRecordNum;

			for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
			{
				if(feof(fp_HeadTail))
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}

				fgets(line1Char, sizeof(line1Char), fp_HeadTail); line1 = line1Char;
				//line1StrVec.push_back(line1);
				line1StrVec[recordNumTmp] = line1;

				fgets(line2Char, sizeof(line2Char), fp_HeadTail); line2 = line2Char;
				line2StrVec[recordNumTmp] = line2;

				fgets(line3Char, sizeof(line3Char), fp_HeadTail); line3 = line3Char;
				line3StrVec[recordNumTmp] = line3;
				
				fgets(line4Char, sizeof(line4Char), fp_HeadTail); line4 = line4Char;
				line4StrVec[recordNumTmp] = line4;

				fgets(line5Char, sizeof(line5Char), fp_HeadTail); line5 = line5Char;
				line5StrVec[recordNumTmp] = line5;
				fgets(line6Char, sizeof(line6Char), fp_HeadTail); line6 = line6Char;
				line6StrVec[recordNumTmp] = line6;
				fgets(line7Char, sizeof(line7Char), fp_HeadTail); line7 = line7Char;
				line7StrVec[recordNumTmp] = line7;
				fgets(line8Char, sizeof(line8Char), fp_HeadTail); line8 = line8Char;
				line8StrVec[recordNumTmp] = line8;
				fgets(line9Char, sizeof(line9Char), fp_HeadTail); line9 = line9Char;
				line9StrVec[recordNumTmp] = line9;
				fgets(line10Char, sizeof(line10Char), fp_HeadTail); line10 = line10Char;
				line10StrVec[recordNumTmp] = line10;

				fgets(line11Char, sizeof(line11Char), fp_HeadTail); 
			}	
		
			//cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;

			//cout << "finish reading record, turn: " << tmpTurn+1 << endl << endl;

			//cout << "start to fix incomplete end, turn: " << tmpTurn+1 << endl;

			omp_set_num_threads(threads_num);
			//omp_set_num_threads(1);
			#pragma omp parallel for
			for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
			{
				oneEndUnmappedRecordNum ++;
				if(oneEndUnmappedRecordNum < 2202658)
				{
					continue;
				}
				////////////////  parse long head reads record after 1-mapping process  ///////////////////////////////////////
				int Nor1Num = 0, Rcm1Num = 0, Nor2Num = 0, Rcm2Num = 0;
				string readNameStr_1, readNameStr_2;

				int startSearchPos = 0, foundSearchPos;
				foundSearchPos = line1StrVec[tmpOpenMP].find("\t", startSearchPos);
				readNameStr_1 = line1StrVec[tmpOpenMP].substr(startSearchPos, foundSearchPos-1-startSearchPos+1);
				startSearchPos = foundSearchPos + 1;
				foundSearchPos = line1StrVec[tmpOpenMP].find("\t", startSearchPos);
				Nor1Num = atoi((line1StrVec[tmpOpenMP].substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());		
				startSearchPos = foundSearchPos + 1;
				foundSearchPos = line1StrVec[tmpOpenMP].find("\t", startSearchPos);		
				Rcm1Num = atoi((line1StrVec[tmpOpenMP].substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());			

				startSearchPos = 0; 
				foundSearchPos = line4StrVec[tmpOpenMP].find("\t", startSearchPos);
				readNameStr_2 = line4StrVec[tmpOpenMP].substr(startSearchPos, foundSearchPos-1-startSearchPos+1);
				startSearchPos = foundSearchPos + 1;
				foundSearchPos = line4StrVec[tmpOpenMP].find("\t", startSearchPos);
				Nor2Num = atoi((line4StrVec[tmpOpenMP].substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());		
				startSearchPos = foundSearchPos + 1;
				foundSearchPos = line4StrVec[tmpOpenMP].find("\t", startSearchPos);		
				Rcm2Num = atoi((line4StrVec[tmpOpenMP].substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());	

				int readLength_1 = line2StrVec[tmpOpenMP].length() - 1;
				int readLength_2 = line5StrVec[tmpOpenMP].length() - 1;
			

				////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////	
				//inputRecordNum ++;
				//cout << "recordNum: " << inputRecordNum << endl;
				//if(inputRecordNum < 10362)//<= 2000000)//2074795)
				//	continue;
					//break;
					//continue;				

				////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////
				//cout << endl << "readName: " << endl << readNameStr_1 << endl;

				PE_Read_Info* peReadInfo = new PE_Read_Info();

				peReadInfo->getFastaFormatReadInfo(readNameStr_1, readNameStr_2,
					line2StrVec[tmpOpenMP].substr(0, readLength_1), line5StrVec[tmpOpenMP].substr(0, readLength_2));	

				char* read = const_cast<char*>((peReadInfo->readInfo_pe1).readSeq.c_str());
				char* read_PE = const_cast<char*>((peReadInfo->readInfo_pe2).readSeq.c_str());

 		   		char read_RC_PE[(peReadInfo->readInfo_pe1).readSeqLength], read_RC[(peReadInfo->readInfo_pe2).readSeqLength];
    			for(int read_RC_loc = 0; read_RC_loc < (peReadInfo->readInfo_pe1).readSeqLength; read_RC_loc++) // get read_RC
    				*(read_RC + read_RC_loc) = reverseComplement(*(read + (peReadInfo->readInfo_pe1).readSeqLength - 1 - read_RC_loc)); 			
    	 		for(int read_RC_loc = 0; read_RC_loc < (peReadInfo->readInfo_pe2).readSeqLength; read_RC_loc++) // get read_RC_PE
    				*(read_RC_PE + read_RC_loc) = reverseComplement(*(read_PE + (peReadInfo->readInfo_pe2).readSeqLength - 1 - read_RC_loc)); 

	    		string rcmReadSeq_1 = read_RC;
    			(peReadInfo->readInfo_pe1).rcmReadSeq = rcmReadSeq_1.substr(0, (peReadInfo->readInfo_pe1).readSeqLength);

	    		string rcmReadSeq_2 = read_RC_PE;
    			(peReadInfo->readInfo_pe2).rcmReadSeq = rcmReadSeq_2.substr(0, (peReadInfo->readInfo_pe2).readSeqLength);


				PE_Read_Alignment_Info* peAlignInfo = 
					new PE_Read_Alignment_Info(line7StrVec[tmpOpenMP], line8StrVec[tmpOpenMP], line9StrVec[tmpOpenMP], line10StrVec[tmpOpenMP],
					Nor1Num, Rcm1Num, Nor2Num, Rcm2Num);

				//peAlignInfo->generateOtherEndUnmappedBoolVec();

				//peAlignInfo->pairingAlignment();

				
				vector<int> alignmentInfoVecSize_ori;
				alignmentInfoVecSize_ori.push_back(peAlignInfo->norAlignmentInfo_PE_1.size());
				alignmentInfoVecSize_ori.push_back(peAlignInfo->rcmAlignmentInfo_PE_1.size());
				alignmentInfoVecSize_ori.push_back(peAlignInfo->norAlignmentInfo_PE_2.size());
				alignmentInfoVecSize_ori.push_back(peAlignInfo->rcmAlignmentInfo_PE_2.size());
				
				
				int readLength;// = 100; // to debug

				bool End1OrEnd2; bool NorOrRcm;


					//cout << endl << "No Pair Read Name: " << endl << readNameStr_1.substr(0, readNameStr_1.length() - 1) << endl;

					for(int tmpAlignInfoType = 1; tmpAlignInfoType <= 4; tmpAlignInfoType ++)
					{

						End1OrEnd2 = peReadInfo->checkEnd1OrEnd2WithAlignInfoTypeNo(tmpAlignInfoType);
						NorOrRcm = peReadInfo->checkNorOrRcmWithAlignInfoTypeNo(tmpAlignInfoType);
						readLength = peReadInfo->checkReadLengthWithAlignInfoTypeNo(tmpAlignInfoType);

						//int peAlignInfoVectorSize = (peAlignInfo->getAlignInfoVecSize(tmpAlignInfoType));

						int peAlignInfoVecOriSize = alignmentInfoVecSize_ori[tmpAlignInfoType - 1];

						//cout << "peAlignInfoVectorSize: " << peAlignInfoVectorSize << endl;
						//cout << "tmpAlignInfoType: " << tmpAlignInfoType << endl;
						//cout << "peAlignInfoVecOriSize: " << peAlignInfoVecOriSize << endl;
						for(int tmpAlignmentNO = 0; 
							tmpAlignmentNO < peAlignInfoVecOriSize;//peAlignInfoVectorSize; 
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
							
							if(pathInfo->finalPathVec.size() == 1)
							{
								peAlignInfo->pushBackPathInfo2PeAlignInfo(pathInfo, End1OrEnd2, NorOrRcm, indexInfo);
							}

							delete(gapInfo);
							pathInfo->memoryFree();
							delete(pathInfo);
							delete(segInfo);
							delete(seg2ndOriInfo);
							delete(unmapEndInfo);
						}

					}

				peAlignInfo->pairingAlignment();
				peAlignInfo->chooseBestAlignment();

				bool pairExistsBool = peAlignInfo->finalPairExistsBool();
				bool allAlignmentCompleteBool = peAlignInfo->allAlignmentInFinalPairCompleted();
				
				string tmpPeAlignSamStr, tmpPeAlignInfoStr, tmpPeAlignSamStr_unpair_complete, tmpPeAlignInfo_complete_pair;
				if(pairExistsBool && allAlignmentCompleteBool) // some pair exists, all completed, print out paired SAM info
				{
					/*tmpPeAlignSamStr = peAlignInfo->getTmpPEreadAlignInfoInSAMformatForFinalPair(
					 	(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName, 
						(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq);*/
					tmpPeAlignSamStr = peAlignInfo->getSAMformatForFinalPair(
					 	(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName, 
						(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq);
					tmpPeAlignInfoStr = "";
					tmpPeAlignSamStr_unpair_complete = "";
					if(outputAlignInfoAndSamForAllPairedAlignmentBool)
					{
						tmpPeAlignInfo_complete_pair = peAlignInfo->getTmpAlignInfoForFinalPair(
						 	(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName, 
							(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq,
							(peReadInfo->readInfo_pe1).readQual, (peReadInfo->readInfo_pe2).readQual);
					}
				}
				else if(pairExistsBool && (!allAlignmentCompleteBool)) // pair exists, incomplete
				{
					tmpPeAlignSamStr = "";
					tmpPeAlignInfoStr = peAlignInfo->getTmpAlignInfoForFinalPair(
					 	(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName, 
						(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq,
						(peReadInfo->readInfo_pe1).readQual, (peReadInfo->readInfo_pe2).readQual);
					tmpPeAlignSamStr_unpair_complete = "";
					if(outputAlignInfoAndSamForAllPairedAlignmentBool)
					{					
						tmpPeAlignInfo_complete_pair = "";
					}
				}
				else if((!pairExistsBool) && (allAlignmentCompleteBool)) // no pair exists, all complete, print out original SAM info
				{

					/*tmpPeAlignSamStr = peAlignInfo->getTmpPEreadAlignInfoInSAMformat(
						(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName,
						(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq );*/
					tmpPeAlignSamStr = "";
					tmpPeAlignInfoStr = "";
					tmpPeAlignSamStr_unpair_complete = peAlignInfo->getTmpPEreadAlignInfoInSAMformat(
						(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName,
						(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq );
					if(outputAlignInfoAndSamForAllPairedAlignmentBool)
					{					
						tmpPeAlignInfo_complete_pair = "";
					}
				}
				else // no pair exists, incomplete, print out alignInfo
				{
					tmpPeAlignSamStr = "";
					tmpPeAlignInfoStr = peAlignInfo->getTmpAlignInfo(
						(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName,
						(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq, 
						//readQualSeq_1, readQualSeq_2
						"*", "*");// << endl;
					tmpPeAlignSamStr_unpair_complete = "";
					if(outputAlignInfoAndSamForAllPairedAlignmentBool)
					{
						tmpPeAlignInfo_complete_pair = "";
					}
				}

				peAlignSamVec[tmpOpenMP] = tmpPeAlignSamStr;

				peAlignInfoVec[tmpOpenMP] = tmpPeAlignInfoStr;
				peAlignSamVec_unpair_complete[tmpOpenMP] = tmpPeAlignSamStr_unpair_complete;

				if(outputAlignInfoAndSamForAllPairedAlignmentBool)
				{
					peAlignInfoVec_pair_complete[tmpOpenMP] = tmpPeAlignInfo_complete_pair;
				}
				delete(peReadInfo);
				peAlignInfo->memoryFree();
				delete(peAlignInfo);

			}

			//cout << "start to output ... turn: " << tmpTurn+1 << endl;
			for(int tmp = 0; tmp < realRecordNum; tmp++)
			{
				//if((peAlignSamVec[tmp].length() == 0)||( peAlignInfoVec[tmp].length() == 0))
				//	continue;
				if(peAlignSamVec[tmp] != "")
				{
					OutputSamFile_oneEndMapped_ofs << peAlignSamVec[tmp] << endl;
				}
				if(peAlignInfoVec[tmp] != "")
				{
					tmpAlignIncompletePair_ofs << peAlignInfoVec[tmp] << endl;// << endl;
				}
				if(peAlignSamVec_unpair_complete[tmp] != "")
				{
					OutputSamFile_oneEndMapped_unpairComplete_ofs << peAlignSamVec_unpair_complete[tmp] << endl;
				}
				if(outputAlignInfoAndSamForAllPairedAlignmentBool)
				{					
					if(peAlignInfoVec_pair_complete[tmp] != "")
					{
						OutputSamFile_oneEndMapped_alignInfo_ofs << peAlignInfoVec_pair_complete[tmp] << endl;
					}
				}
			}		
			//cout << "finish output, turn: " << tmpTurn+1 << endl << endl;

		}		

	}

	OutputSamFile_oneEndMapped_ofs.close();
	OutputSamFile_oneEndMapped_unpairComplete_ofs.close();
	tmpAlignIncompletePair_ofs.close();
	//tmpAlignInfoForDebugFile_oneEndMapped_ofs.close();

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads ends ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads ends ......" << endl << endl ; 

	cout << "completeEndNum: " << completeEndNum << endl;
 	cout << "incompleteEndNum: " << incompleteEndNum << endl;
 	cout << "incompleteEndFixedNum: " << incompleteEndFixedNum << endl;
	cout << "incompleteEndUnFixedNum: " << incompleteEndUnFixedNum << endl;
	cout << "bothEndUnmapNum: " << bothEndUnmapNum << endl;

	overAll_end = clock();
	loadIndex_cost = loadIndex_end - loadIndex_begin;
	double loadIndexTime = (double)loadIndex_cost/CLOCKS_PER_SEC;

	overAll_cost = overAll_end - overAll_begin;
	double overAllTime = (double)overAll_cost/CLOCKS_PER_SEC;

	cout << "overAllTime = " << overAllTime << endl;
	cout << "loadIndexTime = " << loadIndexTime << endl;
	return 0;
} 