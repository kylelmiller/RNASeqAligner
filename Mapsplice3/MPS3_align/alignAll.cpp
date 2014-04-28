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
		freeMem_begin, freeMem_end, freeMem_cost = 0;

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
    overall_begin = clock();
    #endif
	//////////////////////////////////////////////////
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
	
	read_file_end = clock();
	double read_file_time = (double)(read_file_end - read_file_begin)/CLOCKS_PER_SEC;
	log_ofs << "read_file_time = " << read_file_time << endl;
	
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
	//ofstream tmpAlignInfoForDebugFile_ofs(tmpAlignInfoForDebug.c_str());

	string tmpAlignUnmap = testClassFileStr + ".unmap";
	//ofstream tmpAlignUnmap_ofs(tmpAlignUnmap.c_str());

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

	string tmpAlignIncompletePair = testClassFileStr + ".incompletePair.alignInfo"; 
	ofstream tmpAlignIncompletePair_ofs(tmpAlignIncompletePair.c_str());

	string tmpAlignIncompletePair_SAM = testClassFileStr + ".incompletePair.sam"; 
	ofstream tmpAlignIncompletePair_SAM_ofs(tmpAlignIncompletePair_SAM.c_str());	

	string tmpIntermediateJunctionFile = testClassFileStr + ".inter.junc";
	//ofstream tmpIntermediateJunctionFile_ofs(tmpIntermediateJunctionFile.c_str());

	string tmpAlignCompleteRead_alignInfo = testClassFileStr + ".completePair.sam_alignInfo";
	ofstream tmpAlignCompleteRead_alignInfo_ofs(tmpAlignCompleteRead_alignInfo.c_str());


    FILE *fp_in = fopen(InputReadFile, "r");    
    FILE *fp_in_PE = fopen(InputReadFile_PE, "r");

    char line1Char[500], line2Char[150], line2Char_RC[150],
    	line1Char_PE[500], line2Char_PE[150], line2Char_PE_RC[150],
    	line0Char[500]; 

    string line1, line2, line2_RC, 
    	line1_PE, line2_PE, line2_PE_RC;

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
			//cout << "recordNumTmp: " << recordNumTmp << endl;
			if((feof(fp_in)) || (feof(fp_in_PE)))
			{
				//cout << "here" << endl;
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;
			}
			//cout << endl << "recordNumTmp: " << recordNumTmp << endl;
			fgets(line1Char, sizeof(line1Char), fp_in); line1 = line1Char;
			line1 = line1.substr(1, line1.length()-2);
			readName1Vec[recordNumTmp] = line1;
			//cout << "readName_1: " << line1 << endl;

			fgets(line2Char, sizeof(line2Char), fp_in); line2 = line2Char;
			int readLength_1 = line2.length() - 1;
			readSeq1Vec[recordNumTmp] = line2.substr(0, readLength_1);
			if(readLength_1 > readLengthMax_tmp)
			{
				readLengthMax_tmp = readLength_1;
			}

    		if(InputAsFastq)
    		{
    			fgets(line0Char, sizeof(line0Char), fp_in);
    			fgets(line0Char, sizeof(line0Char), fp_in);
    		}

			fgets(line1Char_PE, sizeof(line1Char_PE), fp_in_PE); line1_PE = line1Char_PE;
			line1_PE = line1_PE.substr(1, line1_PE.length()-2);
			readName2Vec[recordNumTmp] = line1_PE;

			fgets(line2Char_PE, sizeof(line2Char_PE), fp_in_PE); line2_PE = line2Char_PE;
			int readLength_2 = line2_PE.length() - 1;
			if(readLength_2 > readLengthMax_tmp)
			{
				readLengthMax_tmp = readLength_2;
			}
			readSeq2Vec[recordNumTmp] = line2_PE.substr(0, readLength_2);	

    		if(InputAsFastq)
    		{
    			fgets(line0Char, sizeof(line0Char), fp_in_PE);
    			fgets(line0Char, sizeof(line0Char), fp_in_PE);
    		}
		}

		#ifdef CAL_TIME	
		input_end = clock();
		input_cost = input_cost + input_end - input_begin;
		#endif
		//cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;

		//cout << "finish reading Fasta/Fastq file, turn: " << tmpTurn+1 << endl << endl;

		//cout << "start to fix mid, turn: " << tmpTurn+1 << endl;

		#ifdef CAL_TIME	
		align_begin = clock();
		#endif

		omp_set_num_threads(threads_num);
		
		#pragma omp parallel for
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
		{
			//#ifdef DEBUG_INFO
			/*readPairNum++;
			if(readPairNum < 851606)  
				continue;*/
			//#endif
			//cout << "readName: " << readName1Vec[tmpOpenMP] << endl;
			#ifdef CAL_TIME	
			getReadInfo_begin = clock();
			#endif

			PE_Read_Info* readInfo = new PE_Read_Info();

			readInfo->getFastaFormatReadInfo(readName1Vec[tmpOpenMP], readName2Vec[tmpOpenMP],
				readSeq1Vec[tmpOpenMP], readSeq2Vec[tmpOpenMP]);

			char* read = const_cast<char*>((readInfo->readInfo_pe1).readSeq.c_str());
			char* read_PE = const_cast<char*>((readInfo->readInfo_pe2).readSeq.c_str());

    		char read_RC_PE[(readInfo->readInfo_pe1).readSeqLength], read_RC[(readInfo->readInfo_pe2).readSeqLength];
    		for(int read_RC_loc = 0; read_RC_loc < (readInfo->readInfo_pe1).readSeqLength; read_RC_loc++) // get read_RC
    			*(read_RC + read_RC_loc) = reverseComplement(*(read + (readInfo->readInfo_pe1).readSeqLength - 1 - read_RC_loc)); 			
     		for(int read_RC_loc = 0; read_RC_loc < (readInfo->readInfo_pe2).readSeqLength; read_RC_loc++) // get read_RC_PE
    			*(read_RC_PE + read_RC_loc) = reverseComplement(*(read_PE + (readInfo->readInfo_pe2).readSeqLength - 1 - read_RC_loc)); 

    		//cout << "readName: " << (readInfo->readInfo_pe1).readName << endl;

    		string rcmReadSeq_1 = read_RC;
    		//cout << "rcmReadSeq_1 length: " << rcmReadSeq_1.length() << endl;
    		//cout << "substr length: " << (readInfo->readInfo_pe1).readSeqLength << endl;
    		(readInfo->readInfo_pe1).rcmReadSeq = rcmReadSeq_1.substr(0, (readInfo->readInfo_pe1).readSeqLength);

    		string rcmReadSeq_2 = read_RC_PE;
    		//cout << "rcmReadSeq_2 length: " << rcmReadSeq_2.length() << endl;
    		//cout << "substr length: " << (readInfo->readInfo_pe2).readSeqLength << endl;
    		(readInfo->readInfo_pe2).rcmReadSeq = rcmReadSeq_2.substr(0, (readInfo->readInfo_pe2).readSeqLength);
 
     		//cout << "finish readInfo !" << endl;


			unsigned int norValLength = 0;
			unsigned int norValLength_PE = 0;
			unsigned int rcmValLength = 0;
			unsigned int rcmValLength_PE = 0;
			unsigned int minValLengthToStitch = MIN_LENGTH_TO_STITCH;
			#ifdef CAL_TIME
			getReadInfo_end = clock();
			getReadInfo_cost = getReadInfo_cost + getReadInfo_end - getReadInfo_begin;
			#endif
			////////////////////////////  segment mapping for pair_end_1 reads in both directions  ///////////
    	   	#ifdef CAL_TIME
    	   	segMap_begin = clock(); 	
    		#endif

	    	Seg_Info* segInfo_Nor1 = new Seg_Info();

			bool normalMapMain = segInfo_Nor1->mapMain_SegInfo_preIndex(read, sa, lcpCompress, childTab, 
				chrom, &norValLength, verifyChild, (readInfo->readInfo_pe1).readSeqLength, indexInfo, preIndexMapLengthArray,
				preIndexIntervalStartArray, preIndexIntervalEndArray);

	    	Seg_Info* segInfo_Rcm1 = new Seg_Info();

	    	bool rcmMapMain = segInfo_Rcm1->mapMain_SegInfo_preIndex(read_RC, sa, lcpCompress, childTab,
				chrom, &rcmValLength, verifyChild, (readInfo->readInfo_pe1).readSeqLength, indexInfo, preIndexMapLengthArray,
				preIndexIntervalStartArray, preIndexIntervalEndArray);		

			Seg_Info* segInfo_Nor2 = new Seg_Info();

			bool normalMapMain_PE = segInfo_Nor2->mapMain_SegInfo_preIndex(read_PE, sa, lcpCompress, childTab, 
				chrom, &norValLength_PE, verifyChild, (readInfo->readInfo_pe2).readSeqLength, indexInfo, preIndexMapLengthArray,
				preIndexIntervalStartArray, preIndexIntervalEndArray);

			Seg_Info* segInfo_Rcm2 = new Seg_Info();

			bool rcmMapMain_PE = segInfo_Rcm2->mapMain_SegInfo_preIndex(read_RC_PE, sa, lcpCompress, childTab, 
				chrom, &rcmValLength_PE, verifyChild, (readInfo->readInfo_pe2).readSeqLength, indexInfo, preIndexMapLengthArray,
				preIndexIntervalStartArray, preIndexIntervalEndArray);		
			
			//cout << "segInfo ends!" << endl;// << endl;
			#ifdef CAL_TIME
			segMap_end = clock();
			segMap_cost = segMap_cost + segMap_end - segMap_begin;

			getPath_begin = clock();
			#endif	

			Path_Info* pathInfo_Nor1 = new Path_Info();
			pathInfo_Nor1->getPossiPathFromSeg(segInfo_Nor1);

			Path_Info* pathInfo_Rcm1 = new Path_Info();
			pathInfo_Rcm1->getPossiPathFromSeg(segInfo_Rcm1);

			Path_Info* pathInfo_Nor2 = new Path_Info();
			pathInfo_Nor2->getPossiPathFromSeg(segInfo_Nor2);

			Path_Info* pathInfo_Rcm2 = new Path_Info();
			pathInfo_Rcm2->getPossiPathFromSeg(segInfo_Rcm2);

			//cout << "pathInfo ends !" << endl;
			#ifdef CAL_TIME
			getPath_end = clock();
			getPath_cost = getPath_cost + getPath_end - getPath_begin;

			fixGap_begin = clock();
			#endif

			Gap_Info* gapInfo_Nor1 = new Gap_Info();
			gapInfo_Nor1->fixGapInPath(pathInfo_Nor1, segInfo_Nor1, 
				indexInfo, (readInfo->readInfo_pe1).readSeq, (readInfo->readInfo_pe1).readSeqLength);
			//cout << "gapInfo ends 1 !" << endl;
			//cout << "read_RC:" << endl << (readInfo->readInfo_pe1).rcmReadSeq << endl << "readLength: " 
			//	<< endl << (readInfo->readInfo_pe1).readSeqLength << endl; 
			
		    //	cout << endl << "## do segment mapping for Rcm_1 ##" << endl;
		    //	cout << segInfo_Rcm1->segInfoStr(indexInfo) << endl;
			//	cout << pathInfo_Rcm1->possiPathStr() << endl;

			Gap_Info* gapInfo_Rcm1 = new Gap_Info();
			gapInfo_Rcm1->fixGapInPath(pathInfo_Rcm1, segInfo_Rcm1, 
				indexInfo, (readInfo->readInfo_pe1).rcmReadSeq, (readInfo->readInfo_pe1).readSeqLength);
			//cout << "gapInfo ends 2 !" << endl;

			//cout << endl << "Nor2 segInfo and pathInfo: " << endl;
			//cout << segInfo_Nor2->segInfoStr(indexInfo) << endl;
			//cout << pathInfo_Nor2->possiPathStr() << endl;			
			Gap_Info* gapInfo_Nor2 = new Gap_Info();
			gapInfo_Nor2->fixGapInPath(pathInfo_Nor2, segInfo_Nor2, 
				indexInfo, (readInfo->readInfo_pe2).readSeq, (readInfo->readInfo_pe2).readSeqLength);
			//cout << "gapInfo ends 3 !" << endl;

			
			Gap_Info* gapInfo_Rcm2 = new Gap_Info();
			gapInfo_Rcm2->fixGapInPath(pathInfo_Rcm2, segInfo_Rcm2, 
				indexInfo, (readInfo->readInfo_pe2).rcmReadSeq, (readInfo->readInfo_pe2).readSeqLength);
			//cout << "gapInfo ends 4 !" << endl;
			//cout << "gapInfo ends !" << endl;
			#ifdef CAL_TIME
			fixGap_end = clock();
			fixGap_cost = fixGap_cost + fixGap_end - fixGap_begin;		

			////////////////////////////// get all alignment info //////////////////////////////////
			getPEalignInfo_begin = clock();
			#endif

			PE_Read_Alignment_Info* peAlignInfo = new PE_Read_Alignment_Info(
				pathInfo_Nor1, pathInfo_Rcm1, 
				pathInfo_Nor2, pathInfo_Rcm2, indexInfo);

			//cout << "peAlignInfo ends !" << endl << endl;

			peAlignInfo->pairingAlignment();
			peAlignInfo->chooseBestAlignment();
		
			#ifdef CAL_TIME
			getPEalignInfo_end = clock();
			getPEalignInfo_cost = getPEalignInfo_cost + getPEalignInfo_end - getPEalignInfo_begin;

			insertSamFormat_begin = clock();			
			#endif
			
			PeAlignSamStrVec_complete[tmpOpenMP] = "";			
			//vector<string> PeAlignInfoStrVec_complete(normalRecordNum);

			//vector<string> PeAlignSamStrVec_inCompletePair(normalRecordNum);
			PeAlignInfoStrVec_inCompletePair[tmpOpenMP] = "";
			
			//vector<string> PeAlignSamStrVec_oneEndUnmapped(normalRecordNum);
			PeAlignInfoStrVec_oneEndUnmapped[tmpOpenMP] = "";
			
			PeAlignSamStrVec_bothEndsUnmapped[tmpOpenMP] = "";
			//vector<string> PeAlignInfoStrVec_bothEndsUnmapped(normalRecordNum);

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
						PeAlignInfoStrVec_oneEndUnmapped[tmpOpenMP] = peAlignInfo->getTmpPEreadAlignInfoInSAMformat(
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
					PeAlignSamStrVec_bothEndsUnmapped[tmpOpenMP] = peAlignInfo->getTmpPEreadAlignInfoInSAMformat(
			 			(readInfo->readInfo_pe1).readName, (readInfo->readInfo_pe2).readName, 
						(readInfo->readInfo_pe1).readSeq, (readInfo->readInfo_pe2).readSeq);
				}
			}
			
			#ifdef CAL_TIME
			insertSamFormat_end = clock();
			insertSamFormat_cost = insertSamFormat_cost + insertSamFormat_end - insertSamFormat_begin;
			#endif

			#ifdef DEBUG_INFO
			bool printDebugInfoBool = true; printDebugInfoBool = false;
			if(printDebugInfoBool)
			{
				cout << endl << "##### readName_1: " << (readInfo->readInfo_pe1).readName << " #####" << endl;
				cout << "##### readName_2: " << (readInfo->readInfo_pe2).readName << " #####"<< endl;

				cout << endl << "## do segment mapping for Nor_1 ##" << endl;
				cout << segInfo_Nor1->segInfoStr(indexInfo) << endl;
				cout << pathInfo_Nor1->possiPathStr() << endl;
				cout << pathInfo_Nor1->fixedPathVecStr(indexInfo, segInfo_Nor1) << endl;
				cout << pathInfo_Nor1->finalFixedPathStr(indexInfo) << endl;

		    	cout << endl << "## do segment mapping for Rcm_1 ##" << endl;
		    	cout << segInfo_Rcm1->segInfoStr(indexInfo) << endl;
				cout << pathInfo_Rcm1->possiPathStr() << endl;
				cout << pathInfo_Rcm1->fixedPathVecStr(indexInfo, segInfo_Rcm1) << endl;
				cout << pathInfo_Rcm1->finalFixedPathStr(indexInfo) << endl;

				cout << endl << "## do segment mapping for Nor_2 ##" << endl;
				cout << segInfo_Nor2->segInfoStr(indexInfo) << endl;
				cout << pathInfo_Nor2->possiPathStr() << endl;
				cout << pathInfo_Nor2->fixedPathVecStr(indexInfo, segInfo_Nor2) << endl;
				cout << pathInfo_Nor2->finalFixedPathStr(indexInfo) << endl;

				cout << endl << "## do segment mapping for Rcm_2 ##" << endl;
				cout << segInfo_Rcm2->segInfoStr(indexInfo) << endl;
				cout << pathInfo_Rcm2->possiPathStr() << endl;
				cout << pathInfo_Rcm2->fixedPathVecStr(indexInfo, segInfo_Rcm2) << endl;
				cout << pathInfo_Rcm2->finalFixedPathStr(indexInfo) << endl;
			}
			#endif

			#ifdef CAL_TIME
			freeMem_begin = clock();
			#endif

			//free(read_RC);
			//free(read_RC_PE);
			delete(readInfo);
			delete(gapInfo_Nor1);
			delete(gapInfo_Rcm1);
			delete(gapInfo_Nor2);
			delete(gapInfo_Rcm2);
			pathInfo_Nor1->memoryFree();
			delete(pathInfo_Nor1);
			pathInfo_Rcm1->memoryFree();
			delete(pathInfo_Rcm1);
			pathInfo_Nor2->memoryFree();
			delete(pathInfo_Nor2);
			pathInfo_Rcm2->memoryFree();
			delete(pathInfo_Rcm2);
			delete(segInfo_Nor1);
			delete(segInfo_Nor2);
			delete(segInfo_Rcm1);
			delete(segInfo_Rcm2);
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

		//cout << "finish fixing mid, turn: " << tmpTurn+1 << endl << endl;

		//cout << "start to output ... turn: " << tmpTurn+1 << endl;
		
		output_begin = clock();
		#endif

		for(int tmp = 0; tmp < realRecordNum; tmp++)
		{
			#ifdef DEBUG_INFO
			if(!((readPairNum == 62223)||(readPairNum == 109850)||(readPairNum == 162867)
				||(readPairNum == 212413)||(readPairNum == 280343)||(readPairNum == 359925)
				||(readPairNum == 393988)))  
				continue;
			#endif

			if(PeAlignSamStrVec_complete[tmp] != "")
			{
				tmpAlignCompleteRead_ofs << PeAlignSamStrVec_complete[tmp] << endl;
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

	tmpAlignCompleteRead_ofs.close();
	//tmpAlignIncompletePair_ofs.close();
	tmpAlignOneEndUnmapped_ofs.close();
	tmpAlignBothEndsUnmapped_ofs.close();
	tmpAlignIncompletePair_SAM_ofs.close();
	if(Do_Phase1_Only)
	{
		tmpAlignIncompletePair_ofs.close();
	}

	fclose(fp_in);
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
	//cout << endl << "getSamFormat_time = " << getSamFormat_time << endl;
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
				/*inputRecordNum ++;
				if(inputRecordNum == 2074795)//<= 2000000)//2074795)
					continue;*/
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


				//peAlignInfo->pairingAlignment();

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

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Sam 2 Junc   ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... Sam 2 Junc starts ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... Sam 2 Junc starts ......" << endl << endl; 

	string juncfile;
	if(DoSam2JuncBool)
	{
		juncfile = tmpIntermediateJunctionFile;
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
	string juncInsFile = juncfile.append("ins");

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

		string headTailSoftClippingFile = tmpAlignIncompletePair;
		FILE *fp_HeadTail = fopen(headTailSoftClippingFile.c_str(), "r");

		char line1Char[200], line2Char[500], line3Char[500], line4Char[200], line5Char[500],
			 line6Char[500], line7Char[2000], line8Char[2000], line9Char[2000], line10Char[2000],
			 line11Char[200], readNameChar_1[100], readNameChar_2[100], othersChar[100];

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

			//cout << "finish reading Head/Tail records file, turn: " << tmpTurn+1 << endl << endl;

			//cout << "start to fix Head/Tail, turn: " << tmpTurn+1 << endl;

			omp_set_num_threads(threads_num);
			//omp_set_num_threads(1);
			#pragma omp parallel for
			for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
			{
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
							spliceJunctionFoundInHash = unfixedTailInfo.SJsearchInSJhash(SJ, readSeqWithDirection, indexInfo);
						}
						else
						{
							spliceJunctionFoundInHash = false;
						}

						//cout << "spliceJunctionFoundInHash: " << spliceJunctionFoundInHash << endl;
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

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads ends ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads ends ......" << endl << endl ;  


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
    return 0;
} //end main
