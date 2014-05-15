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
#include "readSeqPreProcessing.h"

#define PreIndexSize 268435456

using namespace std;  

clock_t read_file_begin, read_file_end, read_file_end2, read_file_end3, 
		align_begin, align_end, align_cost = 0,
    	overall_begin, overall_end, overall_cost = 0,
    	input_begin, input_end, input_cost = 0,
    	output_begin, output_end, output_cost = 0,
    	segMap_begin, segMap_end, segMap_cost = 0,
    	getPath_begin, getPath_end, getPath_cost = 0,
    	fixGap_begin, fixGap_end, fixGap_cost = 0,
    	getPEalignInfo_begin, getPEalignInfo_end, getPEalignInfo_cost = 0,
    	getSamFormat_begin, getSamFormat_end, getSamFormat_cost = 0,
		insertSamFormat_begin, insertSamFormat_end, insertSamFormat_cost = 0,
		getReadInfo_begin, getReadInfo_end, getReadInfo_cost = 0,
		freeMem_begin, freeMem_end, freeMem_cost = 0,

		generateReadAlignInfo_begin, generateReadAlignInfo_end, generateReadAlignInfo_cost = 0;

unsigned int PairedReadNum = 0, BothUnmappedReadNum = 0, UnpairedReadNum = 0;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char**argv)
{
    if(argc < 9)
	{
		//cout << "Executable <InputReads> <InputReads_PE> <OutputSAM> <threads_num> <Fasta_or_Fastq> <HumanOrMouse>" << endl;
		cout << "Executable <InputReads> <InputReads_PE> <OutputDir> <threads_num> <Fasta_or_Fastq> <wholeGenomeIndexPrefix> <localIndexPrefix> <chromsomeDir>" << endl;

		exit(0);
	}

	#ifdef CAL_TIME
    overall_begin = clock();
    #endif
	//////////////////////////////////////////////////
    string outputDirStr = argv[3];

   	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	system(mkdirOutputCommand.c_str());
   	
   	string processLogStr = outputDirStr + "/process.log";
   	ofstream log_ofs(processLogStr.c_str());

	///////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////   switches of seperate processes    ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////

	bool Do_Phase1_Only = false;
	Do_Phase1_Only = true;

	bool outputAlignInfoAndSamForAllPairedAlignmentBool = false;
	//outputAlignInfoAndSamForAllPairedAlignmentBool = true;

	bool removeAllIntermediateFilesBool = false;
	removeAllIntermediateFilesBool = true;

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
	int normalRecordNum_fixOneEndUnmapped = 1000000;
	int normalRecordNum_fixHeadTail = 1000000;//000000;

	int readTotalNum = 0;

	log_ofs << "normalRecordNum_1stMapping: " << normalRecordNum_1stMapping << endl;
	log_ofs << "normalRecordNum_fixOneEndUnmapped: " << normalRecordNum_fixOneEndUnmapped << endl;
	log_ofs << "normalRecordNum_fixHeadTail: " << normalRecordNum_fixHeadTail << endl;


	////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////

	time_t nowtime;
	nowtime = time(NULL);
	struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... MPS starts ......" << endl << endl;  
	log_ofs << endl << "[" << asctime(local) << "... MPS starts ......" << endl << endl; 
	//////////////////////////////////////////////////////              LOAD INDEX           ////////////////////////////////////////////////////////////////

	log_ofs << "Reads input 1: " << argv[1] << endl;
	log_ofs << "Reads input 2: " << argv[2] << endl;
	log_ofs << "output directory: " << argv[3] << endl;
	//cout << "Annotation: " << argv[4] << endl;
	log_ofs << "threads_num: " << argv[4] << endl;
	log_ofs << "data format: " << argv[5] << endl;
	log_ofs << "wholeGenomeIndex: " << argv[6] << endl;
	log_ofs << "localIndex: " << argv[7] << endl;
	log_ofs << "chromsomeDir: " << argv[8] << endl;

    char *InputReadFile = argv[1];//read sample, exacted from fastq file every time
    char *InputReadFile_PE = argv[2];// another end read for pair-end reads

	string threadsNumStr = argv[4];
	int threads_num = atoi(threadsNumStr.c_str());
	omp_set_num_threads(threads_num);

	string fastqOrFasta = argv[5];

	bool InputAsFastq = false; 
	if((fastqOrFasta == "fastq")||(fastqOrFasta == "Fastq"))
	{
		InputAsFastq = true;
	}
	else if((fastqOrFasta == "fasta")||(fastqOrFasta == "Fasta"))
	{
		InputAsFastq = false;
	}
	else
	{
		cout << "input read format error" << endl;
		log_ofs << "input read format error" << endl;
		exit(0);
	}

	
	//////////////////////////////////////////////////////////              LOAD INDEX           ////////////////////////////////////////////////////////////////
	#ifdef CAL_TIME
	read_file_begin = clock();
	#endif
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 

    //cout << "start to load preIndex ..." << endl;
    log_ofs << "start to load preIndex ..." << endl;
	   

    string preIndexArrayPreStr;
    string indexStr;// = "/data/homes/lxauky/adSA_table/mm9_table/testAll2_table/testAll2Index";
    string chromDirStr;
    string secondLevelIndexStr;

    indexStr = argv[6];
    preIndexArrayPreStr = indexStr;
    chromDirStr = argv[8];
    secondLevelIndexStr = argv[7];

    preIndexArrayPreStr.append("/");
    indexStr.append("/");
    chromDirStr.append("/");
    secondLevelIndexStr.append("/");

	string preIndexMapLengthArrayStr = preIndexArrayPreStr;
	preIndexMapLengthArrayStr.append("_MapLength"); 
	ifstream preIndexMapLengthArray_ifs(preIndexMapLengthArrayStr.c_str(), ios::binary);

	string preIndexIntervalStartArrayStr = preIndexArrayPreStr;
	preIndexIntervalStartArrayStr.append("_IntervalStart");
	ifstream preIndexIntervalStartArray_ifs(preIndexIntervalStartArrayStr.c_str(), ios::binary);

	string preIndexIntervalEndArrayStr = preIndexArrayPreStr;
	preIndexIntervalEndArrayStr.append("_IntervalEnd");
	ifstream preIndexIntervalEndArray_ifs(preIndexIntervalEndArrayStr.c_str(), ios::binary);


	int* preIndexMapLengthArray;
	preIndexMapLengthArray = (int*)malloc(PreIndexSize * sizeof(int));
	preIndexMapLengthArray_ifs.read((char*)preIndexMapLengthArray, PreIndexSize * sizeof(int));

	unsigned int *preIndexIntervalStartArray;
    preIndexIntervalStartArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int));
    preIndexIntervalStartArray_ifs.read((char*)preIndexIntervalStartArray, PreIndexSize * sizeof(int));

	unsigned int *preIndexIntervalEndArray;
    preIndexIntervalEndArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int));
    preIndexIntervalEndArray_ifs.read((char*)preIndexIntervalEndArray, PreIndexSize * sizeof(int));

 	log_ofs << "finish loading preIndex ..." << endl;

	string SA_file = indexStr; SA_file.append("_SA"); ifstream SA_file_ifs(SA_file.c_str(),ios::binary); 
	//CompressIndex
	string lcpCompress_file = indexStr; lcpCompress_file.append("_lcpCompress"); ifstream lcpCompress_file_ifs(lcpCompress_file.c_str(),ios::binary);
	string childTab_file = indexStr; childTab_file.append("_childTab"); ifstream childTab_file_ifs(childTab_file.c_str(),ios::binary);
	string verifyChild_file = indexStr; verifyChild_file.append("_detChild"); ifstream verifyChild_file_ifs(verifyChild_file.c_str(),ios::binary);
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);

	string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);

	Index_Info* indexInfo = new Index_Info(parameter_file_ifs);

	log_ofs << "index: " << indexStr << endl;
	///////////////////////////////////////
 
	log_ofs << "start to load whole genome" << endl;
	char *chrom;

	chrom = (char*)malloc((indexInfo->indexSize) * sizeof(char));
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->indexSize) * sizeof(char)); 

	indexInfo->chromString = chrom;
	log_ofs << "chromSize = " <<(indexInfo->chromString).size() << endl;
	
	log_ofs << "start to load every chromosome" << endl;
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

	log_ofs << "start to load SA" << endl;
    unsigned int *sa;
    sa = (unsigned int*)malloc((indexInfo->indexSize) * sizeof(unsigned int));
    SA_file_ifs.read((char*)sa, (indexInfo->indexSize) * sizeof(unsigned int));
	log_ofs << "start to load lcpCompress" << endl;
	BYTE *lcpCompress;
	lcpCompress = (BYTE*)malloc((indexInfo->indexSize) * sizeof(BYTE));
	lcpCompress_file_ifs.read((char*)lcpCompress, (indexInfo->indexSize) * sizeof(BYTE));	

	//Xinan: Compress Index
	//read childTab file
	log_ofs << "start to load childTab " << endl;
	unsigned int *childTab;
	childTab = (unsigned int*)malloc((indexInfo->indexSize) * sizeof(unsigned int));
	childTab_file_ifs.read((char*)childTab, (indexInfo->indexSize) * sizeof(unsigned int));
	//cout << "childTab loaded" << endl;

	log_ofs << "start to load detChild" << endl;
	BYTE *verifyChild;
	verifyChild = (BYTE*)malloc((indexInfo->indexSize) * sizeof(BYTE));
	verifyChild_file_ifs.read((char*)verifyChild, (indexInfo->indexSize) * sizeof(BYTE));
	
	log_ofs << "All index files loaded" << endl;
	
	#ifdef CAL_TIME
	read_file_end = clock();
	double read_file_time = (double)(read_file_end - read_file_begin)/CLOCKS_PER_SEC;
	log_ofs << "read_file cpu time = " << read_file_time << endl;
	#endif
	//////////////////////////////////////////////////
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... whole genome index loaded ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... whole genome index loaded ......" << endl << endl;

	//////////////////////////////////////////////////       finish LOADing INDEX           ////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////       finish LOADing INDEX           ////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////       1st Mapping Process           ////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////       1st Mapping Process           ////////////////////////////////////////////////////////////////

	/* align main*/



   	string testClassFileStr = outputDirStr + "/output.sam";
   	//testClassFileStr += ".testClass.sam";
   	//ofstream testClassFile_ofs(testClassFileStr.c_str());

	string tmpAlignInfoForDebug = testClassFileStr+ ".alignInfoDebug";

	string tmpAlignUnmap = testClassFileStr + ".unmap";

	string tmpAlignCompleteRead = testClassFileStr + ".completePair.sam";
	ofstream tmpAlignCompleteRead_ofs(tmpAlignCompleteRead.c_str());

	string tmpAlignOneEndUnmapped = testClassFileStr + ".oneEndUnmapped";
	
	if(Do_Phase1_Only)
		tmpAlignOneEndUnmapped += ".sam";	
	else
		tmpAlignOneEndUnmapped += ".alignInfo";

	ofstream tmpAlignOneEndUnmapped_ofs(tmpAlignOneEndUnmapped.c_str());

	string tmpAlignBothEndsUnmapped = testClassFileStr + ".bothEndsUnmapped.sam";
	ofstream tmpAlignBothEndsUnmapped_ofs(tmpAlignBothEndsUnmapped.c_str());

	string tmpAlignIncompletePair = testClassFileStr + ".incomplete.alignInfo"; 
	ofstream tmpAlignIncompletePair_ofs(tmpAlignIncompletePair.c_str());

	string tmpAlignIncompletePair_SAM = testClassFileStr + ".incompletePair.sam"; 
	ofstream tmpAlignIncompletePair_SAM_ofs(tmpAlignIncompletePair_SAM.c_str());	

	string tmpIntermediateJunctionFile = testClassFileStr + ".inter.junc";

	string tmpAlignCompleteRead_alignInfo = testClassFileStr + ".completePair.sam_alignInfo";
	ofstream tmpAlignCompleteRead_alignInfo_ofs(tmpAlignCompleteRead_alignInfo.c_str());

	ifstream inputRead_ifs(InputReadFile);
	ifstream inputRead_PE_ifs(InputReadFile_PE);

    string line1, line2, line3, line4, line1_PE, line2_PE, line3_PE, line4_PE;
    string line2_afterProcess, line2_PE_afterProcess;

	int normalRecordNum = normalRecordNum_1stMapping; //1000000;//1500000;

	#ifdef DEBUG_INFO
	normalRecordNum = 1;
	#endif

	bool EndOfRecord = false;
	int tmpTurn = 0;
	int realRecordNum;

	int readPairNum = 0;

	vector<string> readName1Vec(normalRecordNum);
	vector<string> readSeq1Vec(normalRecordNum);
	//vector<string> readSeq1Vec_RC(normalRecordNum);
	vector<string> readName2Vec(normalRecordNum);
	vector<string> readSeq2Vec(normalRecordNum);

	vector<string> PeAlignSamStrVec_complete(normalRecordNum);

	vector<string> PeAlignInfoStrVec_inCompletePair(normalRecordNum);

	vector<string> PeAlignInfoStrVec_oneEndUnmapped(normalRecordNum);

	vector<string> PeAlignSamStrVec_bothEndsUnmapped(normalRecordNum);

	vector<string> PeAlignSamStrVec_inCompletePair(normalRecordNum);

	//if(outputAlignInfoAndSamForAllPairedAlignmentBool)
	//{
		vector<string> PeAlignInfoStrVec_completePaired(normalRecordNum);
	//}

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... 1st mapping process starts ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... 1st mapping process starts ......" << endl << endl; 

	int readLengthMax_tmp = 0;

	InputReadPreProcess* readPreProcessInfo = new InputReadPreProcess();

	for(tmpTurn = 0; ; tmpTurn++)
	{
		#ifdef CAL_TIME		
		input_begin = clock();
		#endif

		if(EndOfRecord)
			break;

		int recordNum = normalRecordNum;

		//cout << "start to read Fasta file, turn: " << tmpTurn + 1 << endl;

		realRecordNum = normalRecordNum;

		for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
		{
    		if((inputRead_ifs.eof())||(inputRead_PE_ifs.eof()))
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;    			
    		}

    		getline(inputRead_ifs, line1); // readName_1
    		
    		if((inputRead_ifs.eof())||(inputRead_PE_ifs.eof()))
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;    			
    		}
    		
    		readName1Vec[recordNumTmp] = line1.substr(1);
    		getline(inputRead_ifs, line2); // readSeq_1

    		//readSeq1Vec[recordNumTmp] = line2;
    		//int readLength_1 = line2.length();

    		line2_afterProcess = readPreProcessInfo->upperCaseReadSeq(line2);
    		//line2_afterProcess = readPreProcessInfo->readPreProcess_upperCase_trim(line2);
    		readSeq1Vec[recordNumTmp] = line2_afterProcess;
    		int readLength_1 = line2_afterProcess.length();


			if(readLength_1 > readLengthMax_tmp)
			{
				readLengthMax_tmp = readLength_1;
			}
    		if(InputAsFastq)
    		{
    			getline(inputRead_ifs, line3);
    			getline(inputRead_ifs, line4);
    		}

    		getline(inputRead_PE_ifs, line1_PE); // readName_2
    		readName2Vec[recordNumTmp] = line1_PE.substr(1);
    		getline(inputRead_PE_ifs, line2_PE);
    	
    		//readSeq2Vec[recordNumTmp] = line2_PE;
    		//int readLength_2 = line2_PE.length();
    	    		
    		line2_PE_afterProcess = readPreProcessInfo->upperCaseReadSeq(line2_PE);
    		//line2_PE_afterProcess = readPreProcessInfo->readPreProcess_upperCase_trim(line2_PE);
    		readSeq2Vec[recordNumTmp] = line2_PE_afterProcess;
    		int readLength_2 = line2_PE_afterProcess.length();

    		if(readLength_2 > readLengthMax_tmp)
    		{
    			readLengthMax_tmp = readLength_2;
    		}
    		if(InputAsFastq)
    		{
    			getline(inputRead_PE_ifs, line3_PE);
    			getline(inputRead_PE_ifs, line4_PE);
    		}
		}

		readTotalNum += realRecordNum;

		#ifdef CAL_TIME	
		input_end = clock();
		input_cost = input_cost + input_end - input_begin;
		#endif
		//cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;

		//cout << "finish reading Fasta/Fastq file, turn: " << tmpTurn+1 << endl;

		//cout << "start to fix mid, turn: " << tmpTurn+1 << endl;

		#ifdef CAL_TIME	
		align_begin = clock();
		#endif

		omp_set_num_threads(threads_num);
		
		#pragma omp parallel for
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
		{
			//#ifdef DEBUG_INFO
			//readPairNum ++;
			//if((readPairNum != 928)&&(readPairNum < 847283))  
			//	continue;
			//if(readPairNum > 98880)
			//	exit(1);
				//continue;
				//break;
			//#endif
			//cout << "readName: " << readName1Vec[tmpOpenMP] << endl << endl;
			#ifdef CAL_TIME
			getReadInfo_begin = clock();
			#endif

			//cout << "readName_1: " << readName1Vec[tmpOpenMP] << endl;
			//cout << "readSeq_1: " << endl << readSeq1Vec[tmpOpenMP] << endl;
			//cout << "readName_2: " << readName2Vec[tmpOpenMP] << endl;
			//cout << "readSeq_2: " << endl << readSeq2Vec[tmpOpenMP] << endl;

			PE_Read_Info* readInfo = new PE_Read_Info();

			readInfo->getFastaFormatReadInfo(readName1Vec[tmpOpenMP], readName2Vec[tmpOpenMP],
				readSeq1Vec[tmpOpenMP], readSeq2Vec[tmpOpenMP]);

			char* read = const_cast<char*>((readInfo->readInfo_pe1).readSeq.c_str());
			char* read_PE = const_cast<char*>((readInfo->readInfo_pe2).readSeq.c_str());

    		char read_RC[(readInfo->readInfo_pe1).readSeqLength], read_RC_PE[(readInfo->readInfo_pe2).readSeqLength];

    		for(int read_RC_loc = 0; read_RC_loc < (readInfo->readInfo_pe1).readSeqLength; read_RC_loc++) // get read_RC
    			*(read_RC + read_RC_loc) = reverseComplement(*(read + (readInfo->readInfo_pe1).readSeqLength - 1 - read_RC_loc)); 			
     		for(int read_RC_loc = 0; read_RC_loc < (readInfo->readInfo_pe2).readSeqLength; read_RC_loc++) // get read_RC_PE
    			*(read_RC_PE + read_RC_loc) = reverseComplement(*(read_PE + (readInfo->readInfo_pe2).readSeqLength - 1 - read_RC_loc)); 

    		string rcmReadSeq_1 = read_RC;
    		(readInfo->readInfo_pe1).rcmReadSeq = rcmReadSeq_1.substr(0, (readInfo->readInfo_pe1).readSeqLength);

    		string rcmReadSeq_2 = read_RC_PE;
    		(readInfo->readInfo_pe2).rcmReadSeq = rcmReadSeq_2.substr(0, (readInfo->readInfo_pe2).readSeqLength);
			

			#ifdef CAL_TIME 
			getReadInfo_end = clock();
			getReadInfo_cost = getReadInfo_cost + getReadInfo_end - getReadInfo_begin;

    		segMap_begin = clock();
    		#endif

    		
			FixPhase1Info* fixPhase1Info = new FixPhase1Info();

			fixPhase1Info->fixPhase1_segInfo(read, read_RC, read_PE, read_RC_PE, sa, lcpCompress, childTab, 
				chrom, verifyChild, indexInfo, preIndexMapLengthArray, preIndexIntervalStartArray, preIndexIntervalEndArray,
				readInfo);

				//cout << endl << "## do segment mapping for Nor_1 ##" << endl;
				//cout << (fixPhase1Info->segInfo_Nor1)->segInfoStr(indexInfo) << endl;
				//cout << pathInfo_Nor1->possiPathStr() << endl;
				//cout << pathInfo_Nor1->fixedPathVecStr(indexInfo, segInfo_Nor1) << endl;
				//cout << pathInfo_Nor1->finalFixedPathStr(indexInfo) << endl;

		    	//cout << endl << "## do segment mapping for Rcm_1 ##" << endl;
		    	//cout << (fixPhase1Info->segInfo_Rcm1)->segInfoStr(indexInfo) << endl;
				//cout << pathInfo_Rcm1->possiPathStr() << endl;
				//cout << pathInfo_Rcm1->fixedPathVecStr(indexInfo, segInfo_Rcm1) << endl;
				//cout << pathInfo_Rcm1->finalFixedPathStr(indexInfo) << endl;

				//cout << endl << "## do segment mapping for Nor_2 ##" << endl;
				//cout << (fixPhase1Info->segInfo_Nor2)->segInfoStr(indexInfo) << endl;
				//cout << pathInfo_Nor2->possiPathStr() << endl;
				//cout << pathInfo_Nor2->fixedPathVecStr(indexInfo, segInfo_Nor2) << endl;
				//cout << pathInfo_Nor2->finalFixedPathStr(indexInfo) << endl;

				//cout << endl << "## do segment mapping for Rcm_2 ##" << endl;
				//cout << (fixPhase1Info->segInfo_Rcm2)->segInfoStr(indexInfo) << endl;
				//cout << pathInfo_Rcm2->possiPathStr() << endl;
				//cout << pathInfo_Rcm2->fixedPathVecStr(indexInfo, segInfo_Rcm2) << endl;
				//cout << pathInfo_Rcm2->finalFixedPathStr(indexInfo) << endl;		

 
			#ifdef CAL_TIME
			segMap_end = clock();
			segMap_cost = segMap_cost + segMap_end - segMap_begin;

			getPath_begin = clock();
			#endif
			     

			fixPhase1Info->fixPhase1_pathInfo();
			//cout << "finish fixing pathInfo ..." << endl; 
			//cout << fixPhase1Info->pathInfo_Nor1->possiPathStr() << endl;


			#ifdef CAL_TIME
			getPath_end = clock();
			getPath_cost = getPath_cost + getPath_end - getPath_begin;

			fixGap_begin = clock();
			#endif


			fixPhase1Info->fixPhase1_gapInfo(readInfo, indexInfo);


			#ifdef CAL_TIME
			fixGap_end = clock();
			fixGap_cost = fixGap_cost + fixGap_end - fixGap_begin;
			
			getPEalignInfo_begin = clock();
			#endif

			PE_Read_Alignment_Info* peAlignInfo = new PE_Read_Alignment_Info(
				fixPhase1Info->pathInfo_Nor1, fixPhase1Info->pathInfo_Rcm1, 
				fixPhase1Info->pathInfo_Nor2, fixPhase1Info->pathInfo_Rcm2, indexInfo);


			peAlignInfo->pairingAlignment();
			peAlignInfo->chooseBestAlignment();

			#ifdef CAL_TIME
			getPEalignInfo_end = clock();
			getPEalignInfo_cost = getPEalignInfo_cost + getPEalignInfo_end - getPEalignInfo_begin;
			
			getSamFormat_begin = clock();
			#endif		
			
			PeAlignSamStrVec_complete[tmpOpenMP] = "";			
			PeAlignInfoStrVec_inCompletePair[tmpOpenMP] = "";
			PeAlignInfoStrVec_oneEndUnmapped[tmpOpenMP] = "";
			PeAlignSamStrVec_bothEndsUnmapped[tmpOpenMP] = "";
			PeAlignSamStrVec_inCompletePair[tmpOpenMP] = "";

			if(outputAlignInfoAndSamForAllPairedAlignmentBool)
			{
				PeAlignInfoStrVec_completePaired[tmpOpenMP] = "";
			}
			bool pairExistsBool = peAlignInfo->finalPairExistsBool();
			if(pairExistsBool) // some pair exists
			{
				bool allAlignmentCompleteBool = peAlignInfo->allAlignmentInFinalPairCompleted();
				if(allAlignmentCompleteBool)
				{
					/*PeAlignSamStrVec_complete[tmpOpenMP] = peAlignInfo->getTmpPEreadAlignInfoInSAMformatForFinalPair(
					 	(readInfo->readInfo_pe1).readName, (readInfo->readInfo_pe2).readName, 
						(readInfo->readInfo_pe1).readSeq, (readInfo->readInfo_pe2).readSeq);*/
					PeAlignSamStrVec_complete[tmpOpenMP] = peAlignInfo->getSAMformatForFinalPair(
					 	(readInfo->readInfo_pe1).readName, (readInfo->readInfo_pe2).readName, 
						(readInfo->readInfo_pe1).readSeq, (readInfo->readInfo_pe2).readSeq);
					if(outputAlignInfoAndSamForAllPairedAlignmentBool)
					{
						PeAlignInfoStrVec_completePaired[tmpOpenMP] = peAlignInfo->getTmpAlignInfoForFinalPair(
					 		(readInfo->readInfo_pe1).readName, (readInfo->readInfo_pe2).readName, 
							(readInfo->readInfo_pe1).readSeq, (readInfo->readInfo_pe2).readSeq,
							(readInfo->readInfo_pe1).readQual, (readInfo->readInfo_pe2).readQual);
					}
				}
				else
				{
					if(!Do_Phase1_Only)
					{
						PeAlignInfoStrVec_inCompletePair[tmpOpenMP] = peAlignInfo->getTmpAlignInfoForFinalPair(
					 		(readInfo->readInfo_pe1).readName, (readInfo->readInfo_pe2).readName, 
							(readInfo->readInfo_pe1).readSeq, (readInfo->readInfo_pe2).readSeq,
							(readInfo->readInfo_pe1).readQual, (readInfo->readInfo_pe2).readQual);
					}

					/*PeAlignSamStrVec_inCompletePair[tmpOpenMP] = peAlignInfo->getTmpPEreadAlignInfoInSAMformatForFinalPair(
					 	(readInfo->readInfo_pe1).readName, (readInfo->readInfo_pe2).readName, 
						(readInfo->readInfo_pe1).readSeq, (readInfo->readInfo_pe2).readSeq);*/
					PeAlignSamStrVec_inCompletePair[tmpOpenMP] = peAlignInfo->getSAMformatForFinalPair(
					 	(readInfo->readInfo_pe1).readName, (readInfo->readInfo_pe2).readName, 
						(readInfo->readInfo_pe1).readSeq, (readInfo->readInfo_pe2).readSeq);
				}
			}
			else // no pair exists: 1.one end unmapped; 2. both ends unmapped
			{
				bool alignmentExistsBool = peAlignInfo->alignInfoExistsBool();
				if(alignmentExistsBool) // one end unmapped
				{	
					if(Do_Phase1_Only)
					{
						/*PeAlignInfoStrVec_oneEndUnmapped[tmpOpenMP] = peAlignInfo->getTmpPEreadAlignInfoInSAMformat(
			 				(readInfo->readInfo_pe1).readName, (readInfo->readInfo_pe2).readName, 
							(readInfo->readInfo_pe1).readSeq, (readInfo->readInfo_pe2).readSeq);*/ 
						PeAlignInfoStrVec_oneEndUnmapped[tmpOpenMP] = peAlignInfo->getSAMformatForUnpairedAlignments(
							(readInfo->readInfo_pe1).readName, (readInfo->readInfo_pe2).readName, 
							(readInfo->readInfo_pe1).readSeq, (readInfo->readInfo_pe2).readSeq);
					}
					else
					{
						PeAlignInfoStrVec_oneEndUnmapped[tmpOpenMP] = peAlignInfo->getTmpAlignInfo(
						 	(readInfo->readInfo_pe1).readName, (readInfo->readInfo_pe2).readName, 
							(readInfo->readInfo_pe1).readSeq, (readInfo->readInfo_pe2).readSeq,
							(readInfo->readInfo_pe1).readQual, (readInfo->readInfo_pe2).readQual);
					}
				}
				else // both ends unmapped
				{
					/*PeAlignSamStrVec_bothEndsUnmapped[tmpOpenMP] = peAlignInfo->getTmpPEreadAlignInfoInSAMformat(
			 			(readInfo->readInfo_pe1).readName, (readInfo->readInfo_pe2).readName, 
						(readInfo->readInfo_pe1).readSeq, (readInfo->readInfo_pe2).readSeq);*/
					PeAlignSamStrVec_bothEndsUnmapped[tmpOpenMP] = peAlignInfo->getSAMformatForBothEndsUnmapped(
			 			(readInfo->readInfo_pe1).readName, (readInfo->readInfo_pe2).readName, 
						(readInfo->readInfo_pe1).readSeq, (readInfo->readInfo_pe2).readSeq);						
				}
			}
			
			#ifdef CAL_TIME
			getSamFormat_end = clock();
			getSamFormat_cost = getSamFormat_cost + getSamFormat_end - getSamFormat_begin;

			freeMem_begin = clock();
			#endif

			#ifdef DEBUG_INFO
			fixPhase1Info->coutDebugInfo(readInfo, indexInfo);
			#endif

			fixPhase1Info->memoryFree();
			delete (fixPhase1Info);
			delete(readInfo); 
			peAlignInfo->memoryFree(); 
			delete(peAlignInfo);

			#ifdef CAL_TIME
			freeMem_end = clock();
			freeMem_cost = freeMem_cost + freeMem_end - freeMem_begin;
			#endif

		} // read file end
		
		#ifdef CAL_TIME
		align_end = clock();
		align_cost = align_cost + align_end - align_begin;
		
		output_begin = clock();
		#endif

		//cout << "finish fixing mid, turn: " << tmpTurn+1 << endl;

		//cout << "start to output ... turn: " << tmpTurn+1 << endl;
		
		for(int tmp = 0; tmp < realRecordNum; tmp++)
		{
			#ifdef DEBUG_INFO
			if(!((readPairNum == 62223)))  
				continue;
			#endif

			if(PeAlignSamStrVec_complete[tmp] != "")
			{
				tmpAlignCompleteRead_ofs << PeAlignSamStrVec_complete[tmp] << endl;
				PairedReadNum ++;
				//PeAlignSamStrVec_complete[tmp] = "";
			}			

			if(PeAlignInfoStrVec_inCompletePair[tmp] != "")
			{
				tmpAlignIncompletePair_ofs << PeAlignInfoStrVec_inCompletePair[tmp] << endl;
				//PeAlignInfoStrVec_inCompletePair[tmp] = "";
			}
			
			if(PeAlignInfoStrVec_oneEndUnmapped[tmp] != "")
			{
				tmpAlignOneEndUnmapped_ofs << PeAlignInfoStrVec_oneEndUnmapped[tmp] << endl;
				//PeAlignInfoStrVec_oneEndUnmapped[tmp] = "";
			}
			
			if(PeAlignSamStrVec_bothEndsUnmapped[tmp] != "")
			{
				tmpAlignBothEndsUnmapped_ofs << PeAlignSamStrVec_bothEndsUnmapped[tmp] << endl;
				//PeAlignSamStrVec_bothEndsUnmapped[tmp] = "";
				BothUnmappedReadNum ++;
			}

			if(PeAlignSamStrVec_inCompletePair[tmp] != "")
			{
				tmpAlignIncompletePair_SAM_ofs << PeAlignSamStrVec_inCompletePair[tmp] << endl;
			}

			if(outputAlignInfoAndSamForAllPairedAlignmentBool)
			{
				if(PeAlignInfoStrVec_completePaired[tmp] != "")
				{
					tmpAlignCompleteRead_alignInfo_ofs << PeAlignInfoStrVec_completePaired[tmp] << endl;
				}
			}

		}

		#ifdef CAL_TIME
		output_end = clock();
		output_cost = output_cost + output_end - output_begin;
		#endif
		//cout << "finish output, turn: " << tmpTurn+1 << endl << endl;

	}

	log_ofs << "readTotalNum: " << readTotalNum << endl;

	inputRead_ifs.close();
	inputRead_PE_ifs.close();
	tmpAlignCompleteRead_ofs.close();
	//tmpAlignIncompletePair_ofs.close();
	tmpAlignOneEndUnmapped_ofs.close();
	tmpAlignBothEndsUnmapped_ofs.close();
	tmpAlignIncompletePair_SAM_ofs.close();
	if(Do_Phase1_Only)
	{
		tmpAlignIncompletePair_ofs.close();
	}

	//fclose(fp_in);

	free(preIndexMapLengthArray); free(preIndexIntervalStartArray); free(preIndexIntervalEndArray);
	free(sa);free(lcpCompress);//free(child_up);free(child_down);free(child_next);
	free(childTab);free(chrom);
	
	#ifdef CAL_TIME
	overall_end = clock();
	#endif
	
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... 1st mapping process ends ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... 1st mapping process ends ......" << endl << endl; 
	
	log_ofs << endl << "**********************************" << endl << "**********************************";
	
	#ifdef CAL_TIME
	double overall_time = (double)(overall_end - overall_begin)/CLOCKS_PER_SEC;
	double input_time = (double)input_cost/CLOCKS_PER_SEC;
	double align_time = (double)align_cost/CLOCKS_PER_SEC;
	double ouput_time = (double)output_cost/CLOCKS_PER_SEC;

	double getReadInfo_time = (double)getReadInfo_cost/CLOCKS_PER_SEC;
	double segMap_time = (double)segMap_cost/CLOCKS_PER_SEC;
	double getPath_time = (double)getPath_cost/CLOCKS_PER_SEC;
	double fixGap_time = (double)fixGap_cost/CLOCKS_PER_SEC;
	double getPEalignInfo_time = (double)getPEalignInfo_cost/CLOCKS_PER_SEC;
	double getSamFormat_time = (double)getSamFormat_cost/CLOCKS_PER_SEC;
	double insertSamFormat_time = (double)insertSamFormat_cost/CLOCKS_PER_SEC;
	double freeMem_time = (double)freeMem_cost/CLOCKS_PER_SEC;


	cout << endl << "overall_time = " << overall_time << endl;
	cout << endl << "input_time = " << input_time << endl;
	cout << endl << "align_time = " << align_time << endl << endl;
	
	cout << endl << "getReadInfo_time = " << getReadInfo_time << endl;
	cout << endl << "segMap_time = " << segMap_time << endl;
	cout << endl << "getPath_time = " << getPath_time << endl;
	cout << endl << "fixGap_time = " << fixGap_time << endl;
	cout << endl << "getPEalignInfo_time = " << getPEalignInfo_time << endl;
	cout << endl << "freeMem_time = " << freeMem_time << endl;
	cout << endl << "getSamFormat_time = " << getSamFormat_time << endl;
	cout << endl << "insertSamFormat time = " << insertSamFormat_time << endl;

	cout << endl << endl << "ouput_time = " << ouput_time << endl;
	#endif

	log_ofs << endl << "**********************************" << endl << "**********************************" << endl;
	//cout << endl << "totalReadNum = " << read_num << endl;

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
	//////////////////////////////////////   Do REMAPPING On one end unmapped Reads    ///////////////////////////////
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

	int tmpRecordNum_oneEndUnmapped = 0;

	if(DoRemappingOnUnmapEndReadsBool)
	{
		log_ofs << "start doing remapping on unmapped end reads" << endl;

		string oneEndMappedFileStr = tmpAlignOneEndUnmapped;

		ifstream inputRecord_ifs(oneEndMappedFileStr.c_str());

		string line1, line2, line3, line4, line5, line6, line7, 
			line8, line9, line10, line11;
		
		int normalRecordNum = normalRecordNum_fixOneEndUnmapped; //1000000;

		bool EndOfRecord = false;

		int tmpTurn = 0;

		int realRecordNum;// = normalRecordNum;

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

		//getline(inputRecord_ifs, line11);

		for(tmpTurn = 0; /*tmpTurn < TurnNum*/; tmpTurn++)
		{
			if(EndOfRecord)
				break;

			int recordNum = normalRecordNum;

			//cout << "start to read record" << endl;
			//cout << "start to input oneEndUnmapped records, turn: " << tmpTurn+1 << endl;
			realRecordNum = normalRecordNum;

			for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
			{
				if(inputRecord_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}				
				getline(inputRecord_ifs, line11);
				getline(inputRecord_ifs, line1);
				getline(inputRecord_ifs, line2);
				getline(inputRecord_ifs, line3);
				getline(inputRecord_ifs, line4);
				getline(inputRecord_ifs, line5);
				getline(inputRecord_ifs, line6);
				getline(inputRecord_ifs, line7);
				getline(inputRecord_ifs, line8);
				getline(inputRecord_ifs, line9);
				getline(inputRecord_ifs, line10);
				//getline(inputRecord_ifs, line11);

				line1StrVec[recordNumTmp] = line1;
				line2StrVec[recordNumTmp] = line2;
				line3StrVec[recordNumTmp] = line3;
				line4StrVec[recordNumTmp] = line4;
				line5StrVec[recordNumTmp] = line5;
				line6StrVec[recordNumTmp] = line6;
				line7StrVec[recordNumTmp] = line7;
				line8StrVec[recordNumTmp] = line8;
				line9StrVec[recordNumTmp] = line9;
				line10StrVec[recordNumTmp] = line10;

			}	
		
			//cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;

			//cout << "finish reading record, turn: " << tmpTurn+1 << endl;

			//cout << "start to fix oneEndUnmapped, turn: " << tmpTurn+1 << endl;

			omp_set_num_threads(threads_num);
			//omp_set_num_threads(1);
			#pragma omp parallel for
			for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
			{
				////////////////  parse long head reads record after 1-mapping process  ///////////////////////////////////////
				tmpRecordNum_oneEndUnmapped ++;
				//if(tmpRecordNum_oneEndUnmapped <= 2202658)
				//	continue;
				//cout << "recordNum: " << tmpRecordNum_oneEndUnmapped << endl;


				PE_Read_Info* peReadInfo = new PE_Read_Info();
				PE_Read_Alignment_Info* peAlignInfo = new PE_Read_Alignment_Info();
				peAlignInfo->generatePeReadInfoAndPeAlignInfo_Fasta_toFixOneEndUnmapped_getline(line1StrVec[tmpOpenMP], line2StrVec[tmpOpenMP], 
					line4StrVec[tmpOpenMP], line5StrVec[tmpOpenMP],
					line7StrVec[tmpOpenMP], line8StrVec[tmpOpenMP], //line9StrVec[tmpOpenMP],
					line9StrVec[tmpOpenMP], line10StrVec[tmpOpenMP], peReadInfo);		

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


				FixOneEndUnmappedInfo* fixOneEndUnmappedInfo = new FixOneEndUnmappedInfo();
				fixOneEndUnmappedInfo->fixOneEndUnmapped(peReadInfo, peAlignInfo,
					secondLevelChrom,
					secondLevelSa,
					secondLevelLcpCompress,
					secondLevelChildTab,
					secondLevelDetChild,
					secondLevelLcp,
					secondLevelUp,
					secondLevelDown,
					secondLevelNext, 
					load2ndLevelIndexBool_compressedSize,
					indexInfo);

				peAlignInfo->pairingAlignment();
				peAlignInfo->chooseBestAlignment();

				bool pairExistsBool = peAlignInfo->finalPairExistsBool();
				bool allAlignmentCompleteBool = peAlignInfo->allAlignmentInFinalPairCompleted();
				bool allUnpairedAlignmentCompleteBool = peAlignInfo->allUnpairedAlignmentCompleted();

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
				else if((!pairExistsBool) && (allUnpairedAlignmentCompleteBool)) // no pair exists, all complete, print out original SAM info
				{

					/*tmpPeAlignSamStr = peAlignInfo->getTmpPEreadAlignInfoInSAMformat(
						(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName,
						(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq );*/
					tmpPeAlignSamStr = "";
					tmpPeAlignInfoStr = "";
					/*tmpPeAlignSamStr_unpair_complete = peAlignInfo->getTmpPEreadAlignInfoInSAMformat(
						(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName,
						(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq );*/
					tmpPeAlignSamStr_unpair_complete = peAlignInfo->getSAMformatForUnpairedAlignments(
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
			//cout << "finish fixing oneEndUnmapped, turn: " << tmpTurn+1 << endl;// << endl;
			//cout << "start to output ... turn: " << tmpTurn+1 << endl;
			for(int tmp = 0; tmp < realRecordNum; tmp++)
			{
				//if((peAlignSamVec[tmp].length() == 0)||( peAlignInfoVec[tmp].length() == 0))
				//	continue;
				if(peAlignSamVec[tmp] != "")
				{
					OutputSamFile_oneEndMapped_ofs << peAlignSamVec[tmp] << endl;
					PairedReadNum ++;
				}
				if(peAlignInfoVec[tmp] != "")
				{
					tmpAlignIncompletePair_ofs << peAlignInfoVec[tmp] << endl;// << endl;
				}
				if(peAlignSamVec_unpair_complete[tmp] != "")
				{
					OutputSamFile_oneEndMapped_unpairComplete_ofs << peAlignSamVec_unpair_complete[tmp] << endl;
					UnpairedReadNum ++;
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
		inputRecord_ifs.close();
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

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Sam 2 Junc   ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... Sam 2 Junc starts ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... Sam 2 Junc starts ......" << endl << endl; 

	string juncfile;
	string juncInsFile;
	if(DoSam2JuncBool)
	{
		juncfile = tmpIntermediateJunctionFile;
		juncInsFile = tmpIntermediateJunctionFile + "ins";
		int readLength_max = readLengthMax_tmp; //to Debug
		log_ofs << "... readLengthMax = " << readLength_max << endl;
		int read_width = readLength_max;
		string chrom_dir = chromDirStr;
		int min_intron = 50;
		int max_intron = 300000;
		int min_anchor = 1;
		chrom_dir.append("/");
		double max_rank = 0;

		vector<string> comb_mapreads_files;
		comb_mapreads_files.push_back(tmpAlignCompleteRead);
		comb_mapreads_files.push_back(tmpAlignIncompletePair_SAM);
		comb_mapreads_files.push_back(OutputSamFile_oneEndMapped);

		Covert2JuncComb(juncfile.c_str(), comb_mapreads_files, read_width, chrom_dir, max_rank, min_intron, max_intron, min_anchor);

	}
	//string juncInsFile = juncfile.append("ins");

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... Sam 2 Junc ends ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... Sam 2 Junc ends ......" << endl << endl; 

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Do REMAPPING On unfixed head/tail Reads    ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads starts ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads starts ......" << endl << endl ; 

	int junctionNum = 0;
	
	string InputSpliceJunction = tmpIntermediateJunctionFile;
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

	SJhash_Info* SJ = new SJhash_Info();
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

		//SJhash_Info SJ;
		

		SJ->initiateSJintHash(indexInfo->chromNum);

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

			SJ->insert2SJintHash(chrInt, spliceStartPos, spliceEndPos);
		}
		fclose(fp_spliceJunction);

		if(junctionNum == 0)
		{
			spliceJunctionHashExists = false;
		}

		log_ofs << "junctionNum = " << junctionNum << endl;
		log_ofs << "finish building spliceJunction Hash" << endl;		

		log_ofs << "start doing remapping on unfixed head/tail alignments" << endl;

		string headTailSoftClippingFile = tmpAlignIncompletePair;
		//FILE *fp_HeadTail = fopen(headTailSoftClippingFile.c_str(), "r");
		ifstream inputUnfixedHeadTailRecord_ifs(headTailSoftClippingFile.c_str());

		string line1, line2, line3, line4, line5, line6, line7, 
			line8, line9, line10, line11;

		int normalRecordNum = normalRecordNum_fixHeadTail; //1000000;

		bool EndOfRecord = false;

		int tmpTurn = 0;

		int realRecordNum;// = normalRecordNum;

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

		//getline(inputUnfixedHeadTailRecord_ifs, line11);

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

				if(inputUnfixedHeadTailRecord_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}				

				getline(inputUnfixedHeadTailRecord_ifs, line11);
				getline(inputUnfixedHeadTailRecord_ifs, line1);
				getline(inputUnfixedHeadTailRecord_ifs, line2);
				getline(inputUnfixedHeadTailRecord_ifs, line3);
				getline(inputUnfixedHeadTailRecord_ifs, line4);
				getline(inputUnfixedHeadTailRecord_ifs, line5);
				getline(inputUnfixedHeadTailRecord_ifs, line6);
				getline(inputUnfixedHeadTailRecord_ifs, line7);
				getline(inputUnfixedHeadTailRecord_ifs, line8);
				getline(inputUnfixedHeadTailRecord_ifs, line9);
				getline(inputUnfixedHeadTailRecord_ifs, line10);
				//getline(inputUnfixedHeadTailRecord_ifs, line11);

				line1StrVec[recordNumTmp] = line1;
				line2StrVec[recordNumTmp] = line2;
				line3StrVec[recordNumTmp] = line3;
				line4StrVec[recordNumTmp] = line4;
				line5StrVec[recordNumTmp] = line5;
				line6StrVec[recordNumTmp] = line6;
				line7StrVec[recordNumTmp] = line7;
				line8StrVec[recordNumTmp] = line8;
				line9StrVec[recordNumTmp] = line9;
				line10StrVec[recordNumTmp] = line10;
			}	
		
			//cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;

			//cout << "finish reading Head/Tail records file, turn: " << tmpTurn+1 << endl;

			//cout << "start to fix Head/Tail, turn: " << tmpTurn+1 << endl;

			omp_set_num_threads(threads_num);
			//omp_set_num_threads(1);
			#pragma omp parallel for
			for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
			{
				////////////////  parse long head reads record after 1-mapping process  ///////////////////////////////////////

				PE_Read_Info* peReadInfo = new PE_Read_Info();
				PE_Read_Alignment_Info* peAlignInfo = new PE_Read_Alignment_Info();
				peAlignInfo->generatePeReadInfoAndPeAlignInfo_Fasta_toFixIncompleteAlignment_getline(line1StrVec[tmpOpenMP], line2StrVec[tmpOpenMP], 
					line4StrVec[tmpOpenMP], line5StrVec[tmpOpenMP],
					line7StrVec[tmpOpenMP], line8StrVec[tmpOpenMP], //line9StrVec[tmpOpenMP],
					line9StrVec[tmpOpenMP], line10StrVec[tmpOpenMP], peReadInfo);		

				FixHeadTailInfo* fixHeadTailInfo = new FixHeadTailInfo();
				fixHeadTailInfo->fixHeadTail(peReadInfo, peAlignInfo, SJ, 
					secondLevelChrom,
					secondLevelSa,
					secondLevelLcpCompress,
					secondLevelChildTab,
					secondLevelDetChild,
					secondLevelLcp,
					secondLevelUp,
					secondLevelDown,
					secondLevelNext, 
					load2ndLevelIndexBool_compressedSize,
					spliceJunctionHashExists,
					indexInfo);	

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
					bool allUnpairAlignmentCompleteBool = peAlignInfo->allUnpairedAlignmentCompleted();
					
					if(allUnpairAlignmentCompleteBool)
					{
						tmpPeAlignSamStr_complete_pair = "";

						tmpPeAlignSamStr_incomplete_pair = "";
						
						/*tmpPeAlignSamStr_complete_unpair = peAlignInfo->getTmpPEreadAlignInfoInSAMformat(
				 			(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName, 
							(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq);*/

						tmpPeAlignSamStr_complete_unpair = peAlignInfo->getSAMformatForUnpairedAlignments(
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

						/*tmpPeAlignSamStr_incomplete_unpair = peAlignInfo->getTmpPEreadAlignInfoInSAMformat(
				 			(peReadInfo->readInfo_pe1).readName, (peReadInfo->readInfo_pe2).readName, 
							(peReadInfo->readInfo_pe1).readSeq, (peReadInfo->readInfo_pe2).readSeq);*/				

						tmpPeAlignSamStr_incomplete_unpair = peAlignInfo->getSAMformatForUnpairedAlignments(
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
			//cout << "finish fixing Head/Tail, turn: " << tmpTurn+1 << endl;// << endl;

			//cout << "start to output ... turn: " << tmpTurn+1 << endl;
			for(int tmp = 0; tmp < realRecordNum; tmp++)
			{
				if(peAlignSamVec_complete_pair[tmp] != "")
				{
					OutputSamFile_fixHeadTail_complete_pair_ofs << peAlignSamVec_complete_pair[tmp] << endl;
					PairedReadNum ++;
				}
				if(peAlignSamVec_incomplete_pair[tmp] != "")
				{	
					OutputSamFile_fixHeadTail_incomplete_pair_ofs << peAlignSamVec_incomplete_pair[tmp] << endl;
					PairedReadNum ++;
				}
				if(peAlignSamVec_complete_unpair[tmp] != "")
				{
					OutputSamFile_fixHeadTail_complete_unpair_ofs << peAlignSamVec_complete_unpair[tmp] << endl;
					UnpairedReadNum ++;
				}
				if(peAlignSamVec_incomplete_unpair[tmp] != "")	
				{
					OutputSamFile_fixHeadTail_incomplete_unpair_ofs << peAlignSamVec_incomplete_unpair[tmp] << endl;
					UnpairedReadNum ++;
				}

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
		inputUnfixedHeadTailRecord_ifs.close();
	}
	delete(SJ);

	OutputSamFile_fixHeadTail_complete_pair_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_pair_ofs.close();
	OutputSamFile_fixHeadTail_complete_unpair_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_unpair_ofs.close();
	OutputSamFile_fixHeadTail_complete_pair_alignInfo_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_pair_alignInfo_ofs.close();

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);


	cout << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads ends ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads ends ......" << endl << endl ;  

	if(!Do_Phase1_Only)
	{
		cout << "readTotalNum: " << readTotalNum << endl;		
		log_ofs << "readTotalNum: " << readTotalNum << endl;

		//cout << endl << "PairedReadNum: " << PairedReadNum << " " << ((double)(PairedReadNum)/((double)readTotalNum)) * 100 << endl ;  
		//log_ofs << endl << "PairedReadNum: " << PairedReadNum << " " << ((double)(PairedReadNum)/((double)readTotalNum)) * 100 << endl ;  
		//cout << endl << "UnpairedReadNum: " << UnpairedReadNum << " " << ((double)(UnpairedReadNum)/((double)readTotalNum)) * 100 << endl ;  
		//log_ofs << endl << "UnpairedReadNum: " << UnpairedReadNum << " " << ((double)(UnpairedReadNum)/((double)readTotalNum)) * 100 << endl ;  
		//cout << endl << "BothEndsUnmappedReadNum: " << BothUnmappedReadNum << " " << ((double)(BothUnmappedReadNum)/((double)readTotalNum)) * 100 << endl ;  
		//log_ofs << endl << "BothEndsUnmappedReadNum: " << BothUnmappedReadNum << " " << ((double)(BothUnmappedReadNum)/((double)readTotalNum)) * 100 << endl ;  
	}			

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to prepare for final output files ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... start to prepare for final output files ......" << endl << endl ;  

	string finalOutputSam = testClassFileStr;
	if(Do_Phase1_Only)
	{
		string cat_cmd = "cat " + tmpAlignCompleteRead  
			//+ " " + tmpAlignOneEndUnmapped
			+ " " + tmpAlignIncompletePair_SAM 
			+ " " + tmpAlignOneEndUnmapped
			+ " " + tmpAlignBothEndsUnmapped
			+ " > " + finalOutputSam;
		system(cat_cmd.c_str()); 
	}
	else
	{
		string cat_cmd = "cat " + tmpAlignCompleteRead 
			+ " " + OutputSamFile_oneEndMapped 
			+ " " + OutputSamFile_fixHeadTail_complete_pair
			+ " " + OutputSamFile_fixHeadTail_incomplete_pair
			+ " " + OutputSamFile_fixHeadTail_complete_unpair
			+ " " + OutputSamFile_fixHeadTail_incomplete_unpair
			+ " " + OutputSamFile_oneEndMapped_unpairComplete
			+ " " + tmpAlignBothEndsUnmapped 
			+ " > " + finalOutputSam;
		system(cat_cmd.c_str());   		
	}


	if(removeAllIntermediateFilesBool)
	{
		remove(InputSpliceJunction.c_str());
		remove(juncInsFile.c_str());
		remove(tmpAlignIncompletePair.c_str());
		remove(OutputSamFile_fixHeadTail_complete_pair_alignInfo.c_str());
		remove(OutputSamFile_fixHeadTail_incomplete_pair_alignInfo.c_str());
		remove(OutputSamFile_fixHeadTail_complete_pair.c_str());
		remove(OutputSamFile_fixHeadTail_incomplete_pair.c_str());
		remove(OutputSamFile_fixHeadTail_complete_unpair.c_str());
		remove(OutputSamFile_fixHeadTail_incomplete_unpair.c_str());
		remove(OutputSamFile_oneEndMapped.c_str());
		remove(OutputSamFile_oneEndMapped_unpairComplete.c_str());
		remove(tmpAlignOneEndUnmapped.c_str());
		remove(tmpAlignIncompletePair_SAM.c_str());
		remove(tmpAlignBothEndsUnmapped.c_str());
		remove(tmpAlignIncompletePair.c_str());
		remove(tmpAlignCompleteRead_alignInfo.c_str());
		remove(OutputSamFile_oneEndMapped_alignInfo.c_str());
		remove(tmpAlignCompleteRead.c_str());
	}
	delete(indexInfo);

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... all jobs done ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... all jobs done ......" << endl << endl ;  

    return 0;
} //end main
