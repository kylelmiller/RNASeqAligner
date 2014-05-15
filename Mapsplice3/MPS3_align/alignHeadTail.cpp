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

using namespace std;  

////////////////////////////////////////////////////////////////////////////////////////
/////////////////Compare  ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////


int main(int argc, char**argv)
{
    if(argc < 7)
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

    string outputDirStr = argv[2];
   	string testClassFileStr = outputDirStr + "/output";


   	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	system(mkdirOutputCommand.c_str());
   	
   	string processLogStr = outputDirStr + "/process.log";
   	ofstream log_ofs(processLogStr.c_str());

	///////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////   switches of seperate processes    ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
	int normalRecordNum_1stMapping = 1000000;
	int normalRecordNum_fixOneEndUnmapped = 1000000;
	int normalRecordNum_fixHeadTail = 1000000;//000000;

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


	cout << "inputRecords: " << argv[1] << endl;
	cout << "outputSam: " << argv[2] << endl;
	cout << "spliceJunction: " << argv[3] << endl;
	cout << "thread: " << argv[4] << endl;
	//cout << "species: " << argv[5] << endl;

	string threadsNumStr = argv[4];
	//string speciesStr = argv[5];
	string wholeGenomeIndexStr = argv[5];
	string secondLevelIndexStr = argv[6];
	secondLevelIndexStr += "/";

	string spliceJunctionFileStr = argv[3];

	int threads_num = atoi(threadsNumStr.c_str());

	int tmpRecordNum = 0;

	omp_set_num_threads(threads_num);
//////////////////////////////////////////////////////////              LOAD INDEX           ////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////              LOAD INDEX           ////////////////////////////////////////////////////////////////

    string headTailSoftClippingFile = argv[1];

	loadIndex_begin = clock();

	//string secondLevelIndexStr = 2ndLevelIndexStr + "/";
	wholeGenomeIndexStr += "/";
	string chrom_bit_file = wholeGenomeIndexStr + "_chrom";
	string parameter_file = wholeGenomeIndexStr + "_parameter";
   	

    ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
    ifstream parameter_file_ifs(parameter_file.c_str());

    Index_Info* indexInfo = new Index_Info(parameter_file_ifs);
	cout << "indexInfo->MAX: " << indexInfo->indexSize << endl;

	cout << "start to load whole genome" << endl;
	char *chrom;

	chrom = (char*)malloc((indexInfo->indexSize) * sizeof(char));
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->indexSize) * sizeof(char)); 

	indexInfo->chromString = chrom;
	cout << "chromSize = " <<(indexInfo->chromString).size() << endl;
	
	cout << "start to load every chromosome" << endl;
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
//////////////////////////////////////   Do REMAPPING On Head/Tail soft-clipping Reads    ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads starts ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads starts ......" << endl << endl ; 

	int junctionNum = 0;
	
	string InputSpliceJunction = spliceJunctionFileStr;
	FILE *fp_spliceJunction = fopen(InputSpliceJunction.c_str(), "r");

	string OutputSamFile_fixHeadTail_complete_pair = testClassFileStr + ".fixHeadTail_complete_pair.sam";
	ofstream OutputSamFile_fixHeadTail_complete_pair_ofs(OutputSamFile_fixHeadTail_complete_pair.c_str());
		
	string OutputSamFile_fixHeadTail_complete_pair_alignInfo = testClassFileStr + ".fixHeadTail_complete_pair.sam_alignInfo";
	ofstream OutputSamFile_fixHeadTail_complete_pair_alignInfo_ofs(OutputSamFile_fixHeadTail_complete_pair_alignInfo.c_str());

	string OutputSamFile_fixHeadTail_incomplete_pair = testClassFileStr + ".fixHeadTail_incomplete_pair.sam";
	ofstream OutputSamFile_fixHeadTail_incomplete_pair_ofs(OutputSamFile_fixHeadTail_incomplete_pair.c_str());
		
	string OutputSamFile_fixHeadTail_incomplete_pair_alignInfo = testClassFileStr + ".fixHeadTail_incomplete_pair.sam_alignInfo";
	ofstream OutputSamFile_fixHeadTail_incomplete_pair_alignInfo_ofs(OutputSamFile_fixHeadTail_incomplete_pair_alignInfo.c_str());	


	string OutputSamFile_fixHeadTail_complete_unpair = testClassFileStr + ".fixHeadTail_complete_unpair.sam";
	ofstream OutputSamFile_fixHeadTail_complete_unpair_ofs(OutputSamFile_fixHeadTail_complete_unpair.c_str());

	string OutputSamFile_fixHeadTail_incomplete_unpair = testClassFileStr + ".fixHeadTail_incomplete_unpair.sam";
	ofstream OutputSamFile_fixHeadTail_incomplete_unpair_ofs(OutputSamFile_fixHeadTail_incomplete_unpair.c_str());	

	if(DoRemappingOnUnfixedHeadTailAlignmentBool)
	{
    	log_ofs << "start to build spliceJunction Hash" << endl;

    	bool spliceJunctionHashExists = true;

		string entryString;
		int tabLocation1, tabLocation2, tabLocation3, tabLocation4, tabLocation5;
		char entry[500];
		int chrInt;
		int spliceStartPos;
		int spliceEndPos;
		string chrIntString;
		string spliceStartPosString;
		string spliceEndPosString;

		SJhash_Info SJ;

		SJ.initiateSJintHash(indexInfo->chromNum);

		fgets(entry, sizeof(entry), fp_spliceJunction);
		while(!feof(fp_spliceJunction))
		{
			fgets(entry, sizeof(entry), fp_spliceJunction);
			if(feof(fp_spliceJunction))
				break;
			junctionNum ++;
			//cout << "entryString: " << entryString << endl;
			entryString = entry;
			tabLocation1 = entryString.find('\t', 0);
			tabLocation2 = entryString.find('\t', tabLocation1+1);
			tabLocation3 = entryString.find('\t', tabLocation2+1);
			chrIntString = entryString.substr(0, tabLocation1);
			spliceStartPosString = entryString.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
			spliceEndPosString = entryString.substr(tabLocation2+1, tabLocation3-tabLocation2-1);
			chrInt = indexInfo->convertStringToInt(chrIntString);
			spliceStartPos = atoi(spliceStartPosString.c_str());
			spliceEndPos = atoi(spliceEndPosString.c_str());

			//cout << "chrInt: " << chrInt << endl;
			//cout << "spliceStartPos: " << spliceStartPos << endl;
			//cout << "spliceEndPos: " << spliceEndPos << endl;

			SJ.insert2SJintHash(chrInt, spliceStartPos, spliceEndPos);
		}
		fclose(fp_spliceJunction);

		if(junctionNum == 0)
		{
			spliceJunctionHashExists = false;
		}

		log_ofs << "junctionNum = " << junctionNum << endl;
		log_ofs << "finish building spliceJunction Hash" << endl;		

		log_ofs << "start doing remapping on unfixed head/tail alignments" << endl;

		//string headTailSoftClippingFile = headTailSoftClippingFile;
		FILE *fp_HeadTail = fopen(headTailSoftClippingFile.c_str(), "r");
		//ifstream headTail_ifs(headTailSoftClippingFile.c_str());

		/*
		char line1Char[200], line2Char[500], line3Char[500], line4Char[200], line5Char[500],
			 line6Char[500], line7Char[2000], line8Char[2000], line9Char[2000], line10Char[2000],
			 line11Char[200], readNameChar_1[100], readNameChar_2[100], othersChar[100];*/

		string line1, line2, line3, line4, line5, line6, line7, 
			line8, line9, line10, line11;

		//int TurnNum = 10;

		int normalRecordNum = normalRecordNum_fixHeadTail; //1000000;

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
		//vector<string> peAlignInfoVec(normalRecordNum);
		vector<string> peAlignSamVec_complete_pair(normalRecordNum);
		vector<string> peAlignSamVec_incomplete_pair(normalRecordNum);
		vector<string> peAlignSamVec_complete_unpair(normalRecordNum);
		vector<string> peAlignSamVec_incomplete_unpair(normalRecordNum);

		vector<string> peAlignSamVec_complete_pair_alignInfo(normalRecordNum);
		vector<string> peAlignSamVec_incomplete_pair_alignInfo(normalRecordNum);		
		//vector<string> peAlignSamVec_incomplete_unpair(normalRecordNum);

		for(tmpTurn = 0; /*tmpTurn < TurnNum*/; tmpTurn++)
		{
			if(EndOfRecord)
				break;

			int recordNum = normalRecordNum;

			//cout << "start to read long-head file record" << endl;
			//cout << "start to read Head/Tail file record, turn: " << tmpTurn+1 << endl;
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
				/*getline(headTail_ifs, line1);
				getline(headTail_ifs, line2);
				getline(headTail_ifs, line3);
				getline(headTail_ifs, line4);
				getline(headTail_ifs, line5);
				getline(headTail_ifs, line6);
				getline(headTail_ifs, line7);			
				getline(headTail_ifs, line8);
				getline(headTail_ifs, line9);
				getline(headTail_ifs, line10);
				getline(headTail_ifs, line11);*/
			}	
		
			//cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;

			//cout << "finish reading Head/Tail records file, turn: " << tmpTurn+1 << endl << endl;

			//cout << "start to fix Head/Tail, turn: " << tmpTurn+1 << endl;

			omp_set_num_threads(threads_num);
			//omp_set_num_threads(1);
			#pragma omp parallel for
			for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
			{

				//tmpRecordNum ++;
				//if(tmpRecordNum < 8761628)
				//	continue;
				//cout << "tmpRecordNum: " << tmpRecordNum << endl;
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
			
				PE_Read_Info* peReadInfo = new PE_Read_Info(readNameStr_1, readNameStr_2,
					line2StrVec[tmpOpenMP].substr(0, readLength_1), line5StrVec[tmpOpenMP].substr(0, readLength_2));

				//cout << endl << "readName_1: " << readNameStr_1 << endl;
				//cout << "readName_2: " << readNameStr_2 << endl;
				//cout << "readSeq_1: " << line2StrVec[tmpOpenMP].substr(0, readLength_1) << endl;
				//cout << "readSeq_2: " << line5StrVec[tmpOpenMP].substr(0, readLength_2) << endl;

				PE_Read_Alignment_Info* peAlignInfo = 
					new PE_Read_Alignment_Info(line7StrVec[tmpOpenMP], line8StrVec[tmpOpenMP], line9StrVec[tmpOpenMP], line10StrVec[tmpOpenMP],
					Nor1Num, Rcm1Num, Nor2Num, Rcm2Num);

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
							spliceJunctionFoundInHash = unfixedHeadInfo.SJsearchInSJhash(SJ, readSeqWithDirection, indexInfo);
						}
						else
						{
							spliceJunctionFoundInHash = false;
						}
						//cout << "spliceJunctionFoundInHash: " << spliceJunctionFoundInHash << endl;

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
							
							/*bool headSegMapMain = mapMainSecondLevelForTargetMapping(headChar, 
								secondLevelSa[midPartMapPosSecondLevelIndexNO], 
								secondLevelLcp[midPartMapPosSecondLevelIndexNO], 
								secondLevelUp[midPartMapPosSecondLevelIndexNO],
								secondLevelDown[midPartMapPosSecondLevelIndexNO],
								secondLevelNext[midPartMapPosSecondLevelIndexNO],
								secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
								targetMappingStr.length(), 3000000, unfixedHeadInfo.midPartMapPosInWholeGenome, 
								midPartMapPosForLongHeadInSecondLevelIndex, &targetMappingNum, targetMappingLoc,
								indexInfo
								);*/

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
						//cout << "unfixedHeadInfo.CTACsjPos.size(): " << (unfixedHeadInfo.CTACsjPos).size() << endl;
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
						{

						}

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
					//	cout << "tmpAlignmentNO: " << tmpAlignmentNO << endl;
						//tmpAlignmentInfo = (peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO];
						tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
						int cigarStringJumpCodeSize = (tmpAlignmentInfo->cigarStringJumpCode).size();
						if( (tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].type != "S" )
						{
							//cout << "no unfixedTail !" << endl; 
							continue;
						}

					//	cout << "start to getUnfixedTailInfoFromRecordWithAlignInfoType !" << endl;
						Unfixed_Tail unfixedTailInfo;
						//unfixedTailInfo.getUnfixedTailInfoFromRecord(peReadInfo, true, tmpAlignmentInfo, indexInfo);
						unfixedTailInfo.getUnfixedTailInfoFromRecordWithAlignInfoType(peReadInfo, tmpAlignInfoType, tmpAlignmentInfo, indexInfo);
					//	cout << "finish getUnfixedTailInfoFromRecordWithAlignInfoType !" << endl;
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
					//	cout << "start SJsearchInSJhash " << endl;
						bool spliceJunctionFoundInHash;

						if(spliceJunctionHashExists)
						{
							spliceJunctionFoundInHash = unfixedTailInfo.SJsearchInSJhash(SJ, readSeqWithDirection, indexInfo);
						}
						else
						{
							spliceJunctionFoundInHash = false;
						}

					//	cout << "spliceJunctionFoundInHash: " << spliceJunctionFoundInHash << endl;
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
							//cout << "targetMappingNum: " << targetMappingNum << endl;
							/*if(tailSegMapMain)// get all possible short anchor locations
							{
								for(int tmp2 = 0; tmp2 < targetMappingNum; tmp2++)
								{
									int tmpSpliceDistance = targetMappingLoc[tmp2] - finalMidPartMappingPos + 1; 
									if((tmpSpliceDistance < 300000) 
										&& (tmpSpliceDistance > -4))
									{
										(unfixedTailInfo.GTAGsjPos).push_back(pair<int,int>(unfixedTailInfo.possiGTAGpos[tmp], tmpSpliceDistance));
									}
								}
							}*/
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
							
							/*bool tailSegMapMain = mapMainSecondLevelForTargetMapping(tailChar, 
								secondLevelSa[midPartMapPosSecondLevelIndexNO], 
									secondLevelLcp[midPartMapPosSecondLevelIndexNO], 
									secondLevelUp[midPartMapPosSecondLevelIndexNO],
								secondLevelDown[midPartMapPosSecondLevelIndexNO],
								secondLevelNext[midPartMapPosSecondLevelIndexNO],
								secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
								targetMappingStr.length(), 3000000, unfixedTailInfo.midPartMapPosInWholeGenome, 
								midPartMapPosForLongTailInSecondLevelIndex, &targetMappingNum, targetMappingLoc, indexInfo);*/

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

				peAlignInfo->pairingAlignment();
				peAlignInfo->chooseBestAlignment();

				bool pairExistsBool = peAlignInfo->finalPairExistsBool();
							

				string tmpPeAlignSamStr_complete_pair, //tmpPeAlignSamStr_unpair, 
					tmpPeAlignSamStr_incomplete_pair, tmpPeAlignSamStr_complete_unpair, tmpPeAlignSamStr_incomplete_unpair,
					tmpPeAlignSamStr_complete_pair_alignInfo, tmpPeAlignSamStr_incomplete_pair_alignInfo;// tmpPeAlignSamStr_incomplete_unpair;

				if(pairExistsBool) // some pair exists, all completed, print out paired SAM info
				{
					bool allFinalPairAlignmentCompleteBool = peAlignInfo->allAlignmentInFinalPairCompleted();	
					if(allFinalPairAlignmentCompleteBool)
					{
						/*tmpPeAlignSamStr_complete_pair = peAlignInfo->getTmpPEreadAlignInfoInSAMformatForFinalPair(
						 	(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName, 
							(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq);*/
						tmpPeAlignSamStr_complete_pair = peAlignInfo->getSAMformatForFinalPair(
						 	(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName, 
							(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq);						
						tmpPeAlignSamStr_incomplete_pair = "";
						
						tmpPeAlignSamStr_complete_unpair = "";
						
						tmpPeAlignSamStr_incomplete_unpair = "";
						if(outputAlignInfoAndSamForAllPairedAlignmentBool)
						{
							tmpPeAlignSamStr_complete_pair_alignInfo = peAlignInfo->getTmpAlignInfoForFinalPair(
						 		(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName, 
								(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq,
								(peReadInfo->readInfo_pe1).readQual, (peReadInfo->readInfo_pe2).readQual);
							tmpPeAlignSamStr_incomplete_pair_alignInfo = "";
						}
					}
					else
					{
						tmpPeAlignSamStr_complete_pair = "";
						
						/*tmpPeAlignSamStr_incomplete_pair = peAlignInfo->getTmpPEreadAlignInfoInSAMformatForFinalPair(
						 	(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName, 
							(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq);*/
						tmpPeAlignSamStr_incomplete_pair = peAlignInfo->getSAMformatForFinalPair(
						 	(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName, 
							(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq);
						tmpPeAlignSamStr_complete_unpair = "";
						
						tmpPeAlignSamStr_incomplete_unpair = "";

						if(outputAlignInfoAndSamForAllPairedAlignmentBool)
						{
							tmpPeAlignSamStr_complete_pair_alignInfo = "";
							tmpPeAlignSamStr_incomplete_pair_alignInfo = peAlignInfo->getTmpAlignInfoForFinalPair(
						 		(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName, 
								(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq,
								(peReadInfo->readInfo_pe1).readQual, (peReadInfo->readInfo_pe2).readQual);
						}						
					}
				}
				else //if((!pairExistsBool) && (allAlignmentCompleteBool)) // no pair exists, all complete, print out original SAM info
				{
					bool allUnpairAlignmentCompleteBool = peAlignInfo->allUnpairAlignmentComplete();
					
					if(allUnpairAlignmentCompleteBool)
					{
						tmpPeAlignSamStr_complete_pair = "";

						tmpPeAlignSamStr_incomplete_pair = "";
						
						tmpPeAlignSamStr_complete_unpair = peAlignInfo->getTmpPEreadAlignInfoInSAMformat(
				 			(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName, 
							(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq);

						tmpPeAlignSamStr_incomplete_unpair = "";	

						if(outputAlignInfoAndSamForAllPairedAlignmentBool)
						{
							tmpPeAlignSamStr_complete_pair_alignInfo = "";
							tmpPeAlignSamStr_incomplete_pair_alignInfo = "";
						}
					}
					else
					{
						tmpPeAlignSamStr_complete_pair = "";

						tmpPeAlignSamStr_incomplete_pair = "";
						
						tmpPeAlignSamStr_complete_unpair = "";

						tmpPeAlignSamStr_incomplete_unpair = peAlignInfo->getTmpPEreadAlignInfoInSAMformat(
				 			(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName, 
							(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq);				
						if(outputAlignInfoAndSamForAllPairedAlignmentBool)
						{
							tmpPeAlignSamStr_complete_pair_alignInfo = "";
							tmpPeAlignSamStr_incomplete_pair_alignInfo = "";
						}
					}


					//tmpPeAlignSamStr_incomplete_unpair = "";
				}

				peAlignSamVec_complete_pair[tmpOpenMP] = tmpPeAlignSamStr_complete_pair;
				
				peAlignSamVec_incomplete_pair[tmpOpenMP] = tmpPeAlignSamStr_incomplete_pair;
				peAlignSamVec_complete_unpair[tmpOpenMP] = tmpPeAlignSamStr_complete_unpair;
				peAlignSamVec_incomplete_unpair[tmpOpenMP] = tmpPeAlignSamStr_incomplete_unpair;

				if(outputAlignInfoAndSamForAllPairedAlignmentBool)
				{
					peAlignSamVec_complete_pair_alignInfo[tmpOpenMP] = tmpPeAlignSamStr_complete_pair_alignInfo;
					
					peAlignSamVec_incomplete_pair_alignInfo[tmpOpenMP] = tmpPeAlignSamStr_incomplete_pair_alignInfo;
				}
				//peAlignSamVec_incomplete_unpair[tmpOpenMP] = tmpPeAlignSamStr_incomplete_unpair;

				peAlignInfo->memoryFree();
				delete(peReadInfo);
				delete(peAlignInfo);

			}
			//cout << "finish fixing Head/Tail, turn: " << tmpTurn+1 << endl << endl;

			//cout << "start to output ... turn: " << tmpTurn+1 << endl;
			for(int tmp = 0; tmp < realRecordNum; tmp++)
			{
				if(peAlignSamVec_complete_pair[tmp] != "")
					OutputSamFile_fixHeadTail_complete_pair_ofs << peAlignSamVec_complete_pair[tmp] << endl;
				if(peAlignSamVec_incomplete_pair[tmp] != "")
					OutputSamFile_fixHeadTail_incomplete_pair_ofs << peAlignSamVec_incomplete_pair[tmp] << endl;
				if(peAlignSamVec_complete_unpair[tmp] != "")
					OutputSamFile_fixHeadTail_complete_unpair_ofs << peAlignSamVec_complete_unpair[tmp] << endl;
				if(peAlignSamVec_incomplete_unpair[tmp] != "")	
					OutputSamFile_fixHeadTail_incomplete_unpair_ofs << peAlignSamVec_incomplete_unpair[tmp] << endl;

				if(outputAlignInfoAndSamForAllPairedAlignmentBool)
				{
					if(peAlignSamVec_complete_pair_alignInfo[tmp] != "")
						OutputSamFile_fixHeadTail_complete_pair_alignInfo_ofs << peAlignSamVec_complete_pair_alignInfo[tmp] << endl;
					if(peAlignSamVec_incomplete_pair_alignInfo[tmp] != "")
						OutputSamFile_fixHeadTail_incomplete_pair_alignInfo_ofs << peAlignSamVec_incomplete_pair_alignInfo[tmp] << endl;					
				}
			}		
			//cout << "finish output, turn: " << tmpTurn+1 << endl << endl;

		}
		//OutputSamFile_fixHeadTail_incomplete_unpair_ofs.close();		

	}

	OutputSamFile_fixHeadTail_complete_pair_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_pair_ofs.close();
	OutputSamFile_fixHeadTail_complete_unpair_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_unpair_ofs.close();
	OutputSamFile_fixHeadTail_complete_pair_alignInfo_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_pair_alignInfo_ofs.close();
	
	overAll_end = clock();
	loadIndex_cost = loadIndex_end - loadIndex_begin;
	double loadIndexTime = (double)loadIndex_cost/CLOCKS_PER_SEC;

	overAll_cost = overAll_end - overAll_begin;
	double overAllTime = (double)overAll_cost/CLOCKS_PER_SEC;

	cout << "overAllTime = " << overAllTime << endl;
	cout << "loadIndexTime = " << loadIndexTime << endl;
	return 0;
} 