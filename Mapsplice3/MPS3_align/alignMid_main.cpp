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

#include "index_info.h"
#include "constantDefinitions.h"
#include "segmentMapping.h"
#include "read_block_test.h"
#include "bwtmap_info.h"
#include "DoubleAnchorScore.h"
#include "sbndm.h"
#include "splice_info.h"
#include "detAliLocShortSegIncluded.h"
#include "remapping_shortAnchorCharHash.h"
#include "otherFunc.h"
#include "filterMapCasesWithPEinformation.h"
#include "fixGapRelationParameters.h"
//#include "fixSmallExon.h"
#include "smallExon.h"
#include "fixHead.h"
#include "fixTail.h"
#include "fixMid.h"
#include "read_info.h"
#include "fixGap.h"
#include "seg_info.h"
#include "spliceJunction_info.h"
#include "unfixedHead.h"
#include "unfixedTail.h"

#define PreIndexSize 268435456

using namespace std;  

clock_t read_file_begin, read_file_end, read_file_end2, read_file_end3, 
		align_begin, align_end, align_cost = 0,
    	overall_begin, overall_end, overall_cost = 0,
    	input_begin, input_end, input_cost = 0,
    	output_begin, output_end, output_cost = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char**argv)
{
    if(argc < 6)
	{
		cout << "Executable <InputReads> <InputReads_PE> <OutputSAM> <threads_num> <Fasta_or_Fastq> <HumanOrMouse>" << endl;
		exit(0);
	}

//////////////////////////////////////////////////////////              LOAD INDEX           ////////////////////////////////////////////////////////////////
	cout << "# of Arg = " << argc << endl;
	cout << "Command: " << "./run " << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << endl;
	cout << "Reads input 1: " << argv[1] << endl;
	cout << "Reads input 2: " << argv[2] << endl;
	cout << "Sam output: " << argv[3] << endl;
	//cout << "Annotation: " << argv[4] << endl;
	cout << "threads_num: " << argv[4] << endl;
	cout << "data format: " << argv[5] << endl;
	cout << "Species: " << argv[6] << endl;

	string species = argv[6];
	string threadsNumStr = argv[4];
	int threads_num = atoi(threadsNumStr.c_str());

	string fastqOrFasta = argv[5];

	bool InputAsFastq = false; 
	if((fastqOrFasta == "fastq")||(fastqOrFasta == "Fastq"))
	{
		InputAsFastq = true;
	}
	omp_set_num_threads(threads_num);
	//////////////////////////////////////////////////////////              LOAD INDEX           ////////////////////////////////////////////////////////////////

    overall_begin = clock();

    cout << "start to load preIndex ..." << endl;
   
    string preIndexArrayPreStr;
    
    if(species == "Mouse")
    {
		preIndexArrayPreStr = "/data/homes/lxauky/adSA_table/mm9_table/preIndex/preIndexRecord";
	}	
	else if (species == "Human")
	{
		preIndexArrayPreStr = "/data/homes/lxauky/adSA_table/hg19_table/hg19_noRandom/PreIndex/preIndexRecord";
	}
	else
	{
		cout << "species error" << endl;
		exit(0);
	}

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

 	cout << "finish loading preIndex ..." << endl;

    string indexStr;// = "/data/homes/lxauky/adSA_table/mm9_table/testAll2_table/testAll2Index";
    
    if(species == "Mouse")
    {
		indexStr = "/data/homes/lxauky/adSA_table/mm9_table/testAll2_table/testAll2Index";
	}	
	else if (species == "Human")
	{
		indexStr = "/data/homes/lxauky/adSA_table/hg19_table/hg19_noRandom/wholeGenomeIndex";
	}
	else
	{
		cout << "species error" << endl;
		exit(0);
	}

	string SA_file = indexStr; SA_file.append("_SA"); ifstream SA_file_ifs(SA_file.c_str(),ios::binary); 
	//CompressIndex
	string lcpCompress_file = indexStr; lcpCompress_file.append("_lcpCompress"); ifstream lcpCompress_file_ifs(lcpCompress_file.c_str(),ios::binary);
	string childTab_file = indexStr; childTab_file.append("_childTab"); ifstream childTab_file_ifs(childTab_file.c_str(),ios::binary);
	string verifyChild_file = indexStr; verifyChild_file.append("_detChild"); ifstream verifyChild_file_ifs(verifyChild_file.c_str(),ios::binary);
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);

	string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);

	Index_Info* indexInfo = new Index_Info(parameter_file_ifs);

	cout << "index: " << indexStr << endl;
	///////////////////////////////////////
 
	cout << "start to load whole genome" << endl;
	char *chrom;

	chrom = (char*)malloc((indexInfo->indexSize) * sizeof(char));
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->indexSize) * sizeof(char)); 

	chromString = chrom;
	cout << "chromSize = " <<chromString.size() << endl;
	
	cout << "start to load every chromosome" << endl;
	chromStr[0] = chromString.substr(0, (indexInfo->chrEndPosInGenome)[0]+1);
	for(int tmp = 1; tmp < indexInfo->chromNum; tmp++)
	{
		chromStr[tmp] = chromString.substr((indexInfo->chrEndPosInGenome)[tmp-1]+2, 
			(indexInfo->chrEndPosInGenome)[tmp]-(indexInfo->chrEndPosInGenome)[tmp-1]-1);	
	}	

	cout << "start to load SA" << endl;
    unsigned int *sa;
    sa = (unsigned int*)malloc((indexInfo->indexSize) * sizeof(unsigned int));
    SA_file_ifs.read((char*)sa, (indexInfo->indexSize) * sizeof(unsigned int));
	cout << "start to load lcpCompress" << endl;
	BYTE *lcpCompress;
	lcpCompress = (BYTE*)malloc((indexInfo->indexSize) * sizeof(BYTE));
	lcpCompress_file_ifs.read((char*)lcpCompress, (indexInfo->indexSize) * sizeof(BYTE));	

	//Xinan: Compress Index
	//read childTab file
	cout << "start to load childTab " << endl;
	unsigned int *childTab;
	childTab = (unsigned int*)malloc((indexInfo->indexSize) * sizeof(unsigned int));
	childTab_file_ifs.read((char*)childTab, (indexInfo->indexSize) * sizeof(unsigned int));
	//cout << "childTab loaded" << endl;

	cout << "start to load detChild" << endl;
	BYTE *verifyChild;
	verifyChild = (BYTE*)malloc((indexInfo->indexSize) * sizeof(BYTE));
	verifyChild_file_ifs.read((char*)verifyChild, (indexInfo->indexSize) * sizeof(BYTE));
	
	cout << "All index files loaded" << endl;
	
	read_file_end = clock();
	double read_file_time = (double)(read_file_end - read_file_begin)/CLOCKS_PER_SEC;
	cout << "read_file_time = " << read_file_time << endl;
	
	//////////////////////////////////////////////////

	//////////////////////////////////////////////////       finish LOADing INDEX           ////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////       finish LOADing INDEX           ////////////////////////////////////////////////////////////////
	
	/* align main*/
    char *InputReadFile = argv[1];//read sample, exacted from fastq file every time
    char *InputReadFile_PE = argv[2];// another end read for pair-end reads
    char *OutputSamFile = argv[3];

   	string testClassFileStr = argv[3];
   	//testClassFileStr += ".testClass.sam";
   	ofstream testClassFile_ofs(testClassFileStr.c_str());

	string tmpAlignInfoForDebug = testClassFileStr+ ".alignInfoDebug";
	ofstream tmpAlignInfoForDebugFile_ofs(tmpAlignInfoForDebug.c_str());

	string tmpAlignUnmap = testClassFileStr+ ".unmap";
	ofstream tmpAlignUnmap_ofs(tmpAlignUnmap.c_str());

    FILE *fp_in = fopen(InputReadFile, "r");    
    FILE *fp_in_PE = fopen(InputReadFile_PE, "r");

    char line1Char[500], line2Char[150], line2Char_RC[150],
    	line1Char_PE[500], line2Char_PE[150], line2Char_PE_RC[150],
    	line0Char[500]; 

    string line1, line2, line2_RC, 
    	line1_PE, line2_PE, line2_PE_RC;

	int normalRecordNum = 1;//1500000;
	bool EndOfRecord = false;
	int tmpTurn = 0;
	int realRecordNum;

	vector<string> readName1Vec(normalRecordNum);
	vector<string> readSeq1Vec(normalRecordNum);
	vector<string> readSeq1Vec_RC(normalRecordNum);
	vector<string> readName2Vec(normalRecordNum);
	vector<string> readSeq2Vec(normalRecordNum);
	vector<string> readSeq2Vec_RC(normalRecordNum);
	vector<string> PeAlignSamStrVec(normalRecordNum);
	vector<string> PeAlignInfoStrVec(normalRecordNum);

	vector<string> PeAlignSamStrVec_unmap;
	vector<string> PeAlignInfoStrVec_unmap;

	//int readLength_1 = 100;
   	//int readLength_2 = 100;

	for(tmpTurn = 0; ; tmpTurn++)
	{
		input_begin = clock();

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
			//cout << "readSeq_1: " << readSeq1Vec[recordNumTmp] << endl;

    		for(int read_RC_loc = 0; read_RC_loc < readLength_1; read_RC_loc++) // get read_RC
    			*(line2Char_RC + read_RC_loc) = reverseComplement(*(line2Char + readLength_1 - 1 - read_RC_loc)); 			
    		line2_RC = line2Char_RC;
    		readSeq1Vec_RC[recordNumTmp] = line2_RC;
    		//cout << "readSeqRC_1: " << line2_RC << endl;

    		if(InputAsFastq)
    		{
    			fgets(line0Char, sizeof(line0Char), fp_in);
    			fgets(line0Char, sizeof(line0Char), fp_in);
    		}

			fgets(line1Char_PE, sizeof(line1Char_PE), fp_in_PE); line1_PE = line1Char_PE;
			line1_PE = line1_PE.substr(1, line1_PE.length()-2);
			readName2Vec[recordNumTmp] = line1_PE;
			//cout << "readName_2: " << line1_PE << endl;

			fgets(line2Char_PE, sizeof(line2Char_PE), fp_in_PE); line2_PE = line2Char_PE;
			int readLength_2 = line2_PE.length() - 1;
			readSeq2Vec[recordNumTmp] = line2_PE.substr(0, readLength_2);			
			//cout << "readSeq_2: " << readSeq2Vec[recordNumTmp] << endl;

	    	for(int read_RC_loc_PE = 0; read_RC_loc_PE < readLength_2; read_RC_loc_PE++) // get read_RC
    			*(line2Char_PE_RC + read_RC_loc_PE) = reverseComplement(*(line2Char_PE + readLength_2 - 1 - read_RC_loc_PE)); 			
    		line2_PE_RC = line2Char_PE_RC;
    		readSeq2Vec_RC[recordNumTmp] = line2_PE_RC;
    		//cout << "readSeqRC_2: " << line2_PE_RC << endl;

    		if(InputAsFastq)
    		{
    			fgets(line0Char, sizeof(line0Char), fp_in_PE);
    			fgets(line0Char, sizeof(line0Char), fp_in_PE);
    		}
		}

		input_end = clock();
		input_cost = input_cost + input_end - input_begin;

		//cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;

		//cout << "finish reading Fasta file, turn: " << tmpTurn+1 << endl << endl;

		//cout << "start to fix mid, turn: " << tmpTurn+1 << endl;

		align_begin = clock();

		omp_set_num_threads(threads_num);
		
		#pragma omp parallel for
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
		{

			string alignDirectionNor_1 = "+", alignDirectionRcm_1 = "-";
			string alignDirectionNor_2 = "+", alignDirectionRcm_2 = "-";

			string readNameString_1 = readName1Vec[tmpOpenMP];
			string readSeq_1 = readSeq1Vec[tmpOpenMP];//.substr(0, readLength_1);
			int readLength_1 = readSeq_1.length();

			string readStringNor_1 = readSeq_1;
			string readStringRcm_1 = readSeq1Vec_RC[tmpOpenMP];//.substr(0, readLength_1);


			char* read_name = const_cast<char*>(readNameString_1.c_str());
			char* read = const_cast<char*>(readStringNor_1.c_str());
			char* read_RC = const_cast<char*>(readStringRcm_1.c_str());

			string readNameString_2 = readName2Vec[tmpOpenMP];
			string readSeq_2 = readSeq2Vec[tmpOpenMP];//.substr(0, readLength_1);
			int readLength_2 = readSeq_2.length();

			string readStringNor_2 = readSeq_2;
			string readStringRcm_2 = readSeq2Vec_RC[tmpOpenMP];//.substr;//(0, readLength_1);
			
			
			char* read_name_PE = const_cast<char*>(readNameString_2.c_str());
			char* read_PE = const_cast<char*>(readStringNor_2.c_str());
			char* read_RC_PE = const_cast<char*>(readStringRcm_2.c_str());

			string readQualSeq_1 = "*";
			string readQualSeq_2 = "*";

			/*if((readStringNor_1.find("N") != readStringNor_1.npos)||(readStringNor_2.find("N") != readStringNor_2.npos))
			{

				PeAlignSamStrVec[tmpOpenMP] = "Ns";

				PeAlignInfoStrVec[tmpOpenMP] = "Ns";
				continue;
			}*/

	    	bool detAliLocBool, detAliLocBool_RC, fix_gap, fix_gap_RC;
	    	bool detAliLocBool_PE, detAliLocBool_RC_PE, fix_gap_PE, fix_gap_RC_PE;
	  

			unsigned int norValLength = 0;
			unsigned int norValLength_PE = 0;
			unsigned int rcmValLength = 0;
			unsigned int rcmValLength_PE = 0;
			unsigned int minValLengthToStitch = MIN_LENGTH_TO_STITCH;

			//cout << endl << "##### readName_1: " << readNameString_1 << " #####" << endl;
			//cout << "##### readName_2: " << readNameString_2 << " #####"<< endl;

			////////////////////////////  segment mapping for pair_end_1 reads in both directions  ///////////

			//////////////////////////////////////////////////////////////////////////////////////////////////////    	
			//////////////////////////////////segment mapping for normal read of PE_1  /////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////////////////

			//cout << endl << "## do segment mapping for Nor_1 ##" << endl;
	    	Seg_Info* segInfo_Nor1 = new Seg_Info();

			bool normalMapMain = segInfo_Nor1->mapMain_SegInfo_preIndex(read, sa, lcpCompress, childTab, 
				chrom, &norValLength, verifyChild, readLength_1, indexInfo, preIndexMapLengthArray,
				preIndexIntervalStartArray, preIndexIntervalEndArray);


			Frag_Info* fragInfo_Nor1 = new Frag_Info();

			if(norValLength >= minValLengthToStitch) // //debug minium val length for a read
			{
	      		detAliLocBool = fragInfo_Nor1->detAliLoc_FragInfo(segInfo_Nor1, read, chrom, indexInfo);
			}
			else
			{
				fragInfo_Nor1->mapLabelNum = 0;
				detAliLocBool = false;
			}

			/*Path_Info* pathInfo_Nor1 = new Path_Info();
			pathInfo_Nor1->getPossiPathFromSeg(segInfo_Nor1);
			
			cout << segInfo_Nor1->segInfoStr(indexInfo) << endl;
			cout << pathInfo_Nor1->possiPathStr() << endl;
			delete(pathInfo_Nor1);*/
			//////////////////////////////////////////////////////////////////////////////////////////////////////    	
			//////////////////////////////////segment mapping for rcm read of PE_1 ////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////////////////
	    	
	    	//cout << endl << "## do segment mapping for Rcm_1 ##" << endl;
	    	Seg_Info* segInfo_Rcm1 = new Seg_Info();

	    	bool rcmMapMain = segInfo_Rcm1->mapMain_SegInfo_preIndex(read_RC, sa, lcpCompress, childTab,
				chrom, &rcmValLength, verifyChild, readLength_1, indexInfo, preIndexMapLengthArray,
				preIndexIntervalStartArray, preIndexIntervalEndArray);		

			//////////////////////////////  get rcm alignment candidate cases of PE_1  ////////////////////////////////////		
	    	Frag_Info* fragInfo_Rcm1 = new Frag_Info();

			if(rcmValLength >= minValLengthToStitch) // //debug minium val length for a read
			{
				detAliLocBool_RC = fragInfo_Rcm1->detAliLoc_FragInfo(segInfo_Rcm1, read_RC, chrom, indexInfo);	
			}
			else 
			{
				fragInfo_Rcm1->mapLabelNum = 0;
				detAliLocBool_RC = false;
			}

			/*Path_Info* pathInfo_Rcm1 = new Path_Info();
			pathInfo_Rcm1->getPossiPathFromSeg(segInfo_Rcm1);
			
			cout << segInfo_Rcm1->segInfoStr(indexInfo) << endl;
			cout << pathInfo_Rcm1->possiPathStr() << endl;
			delete(pathInfo_Rcm1);*/
			////////////////////////////  segmengt mapping for pair_end_2 reads in both directions /////////////////////////////////////////////

			//////////////////////////////////////////////////////////////////////////////////////////////////////    	
			//////////////////////////////////segment mapping for normal read of PE_2  /////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			//cout << endl << "## do segment mapping for Nor_2 ##" << endl;
			Seg_Info* segInfo_Nor2 = new Seg_Info();

			bool normalMapMain_PE = segInfo_Nor2->mapMain_SegInfo_preIndex(read_PE, sa, lcpCompress, childTab, 
				chrom, &norValLength_PE, verifyChild, readLength_2, indexInfo, preIndexMapLengthArray,
				preIndexIntervalStartArray, preIndexIntervalEndArray);

	        /////////////////////////////////  get normal alignment candidate cases of PE_2 ///////////////////////
			Frag_Info* fragInfo_Nor2 = new Frag_Info();

			if(norValLength_PE >= minValLengthToStitch) // //debug minium val length for a read
			{
				detAliLocBool_PE = fragInfo_Nor2->detAliLoc_FragInfo(segInfo_Nor2, read_PE, chrom, indexInfo);
			}
			else
			{
				fragInfo_Nor2->mapLabelNum = 0;
				detAliLocBool_PE = false;
			}
			/*Path_Info* pathInfo_Nor2 = new Path_Info();
			pathInfo_Nor2->getPossiPathFromSeg(segInfo_Nor2);
			
			cout << segInfo_Nor2->segInfoStr(indexInfo) << endl;
			cout << pathInfo_Nor2->possiPathStr() << endl;
			delete(pathInfo_Nor2);*/
			//////////////////////////////////////////////////////////////////////////////////////////////////////    	
			//////////////////////////////////segment mapping for rcm read of PE_2 ////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			//cout << endl << "## do segment mapping for Rcm_2 ##" << endl;
			Seg_Info* segInfo_Rcm2 = new Seg_Info();

			bool rcmMapMain_PE = segInfo_Rcm2->mapMain_SegInfo_preIndex(read_RC_PE, sa, lcpCompress, childTab, 
				chrom, &rcmValLength_PE, verifyChild, readLength_2, indexInfo, preIndexMapLengthArray,
				preIndexIntervalStartArray, preIndexIntervalEndArray);		

			/////////////////////////////////  get rcm alignment candidate cases of PE_2 ///////////////////////
			Frag_Info* fragInfo_Rcm2 = new Frag_Info();

			if(rcmValLength_PE > minValLengthToStitch) // //debug minium val length for a read
			{
				detAliLocBool_RC_PE = fragInfo_Rcm2->detAliLoc_FragInfo(segInfo_Rcm2, read_RC_PE, chrom, indexInfo);
			}
			else 
			{
				fragInfo_Rcm2->mapLabelNum = 0;
				detAliLocBool_RC_PE = false;
			}
			/*Path_Info* pathInfo_Rcm2 = new Path_Info();
			pathInfo_Rcm2->getPossiPathFromSeg(segInfo_Rcm2);
			
			cout << segInfo_Rcm2->segInfoStr(indexInfo) << endl;
			cout << pathInfo_Rcm2->possiPathStr() << endl;
			delete(pathInfo_Rcm2);*/
			/////////////////////////////// filter multiAlignment cadidate cases with pair-end information /////////////////
			///Input:  mapLabel_1, segMapRangeStart_1, segMapRangeEnd_1, segMapLoc_1
			///		   mapLabel_2, segMapRangeStart_2, segMapRangeEnd_2, segMapLoc_2
			///Output: finalMapLabel_1, finalSegMapRangeStart_1, finalMapRangeEnd_1, finalSegMapLoc_1 //////////////////////
			///		   finalMapLabel_2, finalSegMapRangeStart_2, finalMapRangeEnd_2, finalSegMapLoc_2 //////////////////////

			//////////////////  filter PE_1's normal alignment and PE_2's rcm alignment  ///////////////////////////
			
			int finalNorMapLabel, finalRcmMapLabel_PE;

			unsigned int* finalNorSegMapRangeStart = (unsigned int*)malloc(SEGMENTNUM * CANDALILOC * sizeof(unsigned int));
			unsigned int* finalNorSegMapRangeEnd = (unsigned int*)malloc(SEGMENTNUM * CANDALILOC * sizeof(unsigned int));
			unsigned int* finalNorSegMapLoc = (unsigned int*)malloc(SEGMENTNUM * CANDALILOC * sizeof(unsigned int));
			
			unsigned int* finalRcmSegMapRangeStart_PE = (unsigned int*)malloc(SEGMENTNUM * CANDALILOC * sizeof(unsigned int));
			unsigned int* finalRcmSegMapRangeEnd_PE = (unsigned int*)malloc(SEGMENTNUM * CANDALILOC * sizeof(unsigned int));
			unsigned int* finalRcmSegMapLoc_PE = (unsigned int*)malloc(SEGMENTNUM * CANDALILOC * sizeof(unsigned int));

			filterMapLabelCasesWithPairEndInformation(
				fragInfo_Nor1->mapLabelNum, fragInfo_Nor1->segMapRangeStart, 
				fragInfo_Nor1->segMapRangeEnd, fragInfo_Nor1->segMapLoc,
				fragInfo_Rcm2->mapLabelNum, fragInfo_Rcm2->segMapRangeStart,
				fragInfo_Rcm2->segMapRangeEnd, fragInfo_Rcm2->segMapLoc,

				&finalNorMapLabel, finalNorSegMapRangeStart, finalNorSegMapRangeEnd, finalNorSegMapLoc,
				&finalRcmMapLabel_PE, finalRcmSegMapRangeStart_PE, finalRcmSegMapRangeEnd_PE, finalRcmSegMapLoc_PE, indexInfo);


			//////////////////  filter PE_1's normal alignment and PE_2's rcm alignment  ///////////////////////////
			int finalRcmMapLabel, finalNorMapLabel_PE;
			
			unsigned int* finalRcmSegMapRangeStart = (unsigned int*)malloc(SEGMENTNUM * CANDALILOC * sizeof(unsigned int));
			unsigned int* finalRcmSegMapRangeEnd = (unsigned int*)malloc(SEGMENTNUM * CANDALILOC * sizeof(unsigned int));
			unsigned int* finalRcmSegMapLoc = (unsigned int*)malloc(SEGMENTNUM * CANDALILOC * sizeof(unsigned int));
			
			unsigned int* finalNorSegMapRangeStart_PE = (unsigned int*)malloc(SEGMENTNUM * CANDALILOC * sizeof(unsigned int));
			unsigned int* finalNorSegMapRangeEnd_PE = (unsigned int*)malloc(SEGMENTNUM * CANDALILOC * sizeof(unsigned int));
			unsigned int* finalNorSegMapLoc_PE = (unsigned int*)malloc(SEGMENTNUM * CANDALILOC * sizeof(unsigned int));		
			
			filterMapLabelCasesWithPairEndInformation(
				fragInfo_Rcm1->mapLabelNum, fragInfo_Rcm1->segMapRangeStart,
				fragInfo_Rcm1->segMapRangeEnd, fragInfo_Rcm1->segMapLoc,
				fragInfo_Nor2->mapLabelNum, fragInfo_Nor2->segMapRangeStart, 
				fragInfo_Nor2->segMapRangeEnd, fragInfo_Nor2->segMapLoc,

				&finalRcmMapLabel, finalRcmSegMapRangeStart, finalRcmSegMapRangeEnd, finalRcmSegMapLoc,
				&finalNorMapLabel_PE, finalNorSegMapRangeStart_PE, finalNorSegMapRangeEnd_PE, finalNorSegMapLoc_PE, indexInfo);


			////////////////////////////  generate alignments for pair_end_1 reads /////////////////////////////////////////////
			
			PE_Read_Alignment_Info* peAlignInfo = new PE_Read_Alignment_Info; 

			//////////////////////////////////////  generate normal alignments of PE_1   ////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////////////////		
			//cout << "................generate NORMAL alignments for PE_1 read..................." << endl;

			if(finalNorMapLabel > 0)
			{

				fix_gap = fixGap(segInfo_Nor1->segmentNum, segInfo_Nor1->norSegmentLength, 
					segInfo_Nor1->norSegmentLocInRead, segInfo_Nor1->norSegmentAlignNum,
					finalNorMapLabel, finalNorSegMapRangeStart, finalNorSegMapRangeEnd, finalNorSegMapLoc, 
					read, chrom, readSeq_1, readStringNor_1, alignDirectionNor_1, 
					readNameString_1, readLength_1, peAlignInfo->norAlignmentInfo_PE_1, indexInfo);
			}

			//////////////////////////////////////////////////////////////////////////////////////////////////////    	
			//////////////////////////////////////  generate rcm alignments of PE_1   ////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////////////////	
			//cout << "................generate RCM alignments for PE_1 read..................." << endl;

			if(finalRcmMapLabel > 0)
			{
				fix_gap_RC = fixGap(segInfo_Rcm1->segmentNum, segInfo_Rcm1->norSegmentLength, 
					segInfo_Rcm1->norSegmentLocInRead, segInfo_Rcm1->norSegmentAlignNum,
					finalRcmMapLabel, finalRcmSegMapRangeStart, finalRcmSegMapRangeEnd, finalRcmSegMapLoc, 
					read, chrom, readSeq_1, readStringRcm_1, alignDirectionRcm_1, 
					readNameString_1, readLength_1, peAlignInfo->rcmAlignmentInfo_PE_1, indexInfo);
			}	

			//////////////////////////////////////////////////////////////////////////////////////////////////////   
			////////////////////////////  generate alignmetns for pair_end_2 reads ///////////////////////////////

			//////////////////////////////////////  generate normal alignments of PE_2   /////////////////////////	
			//#ifdef DEBUG
			//cout << "................generate NORMAL alignments for PE_2 read..................." << endl;
			//#endif
			//int norAlignmentNum_2 = 0;
			if(finalNorMapLabel_PE > 0)
			{
				fix_gap_PE = fixGap(segInfo_Nor2->segmentNum, segInfo_Nor2->norSegmentLength,
					segInfo_Nor2->norSegmentLocInRead, segInfo_Nor2->norSegmentAlignNum,
					finalNorMapLabel_PE, finalNorSegMapRangeStart_PE, finalNorSegMapRangeEnd_PE, finalNorSegMapLoc_PE, 
					read_PE, chrom, readSeq_2, readStringNor_2, alignDirectionNor_2, 
					readNameString_2, readLength_2, peAlignInfo->norAlignmentInfo_PE_2, indexInfo);
			}
	  	
			//////////////////////////////////////  generate rcm alignments of PE_2   ////////////////////////////
			//#ifdef DEBUG
			//cout << "................generate RCM alignments for PE_2 read..................." << endl;
			//#endif   
			//int rcmAlignmentNum_2 = 0;
			if(finalRcmMapLabel_PE > 0)
			{
				fix_gap_RC_PE = fixGap(segInfo_Rcm2->segmentNum, segInfo_Rcm2->norSegmentLength,
					segInfo_Rcm2->norSegmentLocInRead, segInfo_Rcm2->norSegmentAlignNum,
					finalRcmMapLabel_PE, finalRcmSegMapRangeStart_PE, finalRcmSegMapRangeEnd_PE, finalRcmSegMapLoc_PE, 
					read_PE, chrom, readSeq_2, readStringRcm_2, alignDirectionRcm_2, 
					readNameString_2, readLength_2, peAlignInfo->rcmAlignmentInfo_PE_2, indexInfo);

			}		

			/////////////  output result //////////////////////////////

			//cout << "start to insert result string to vec ..." << endl;
			string tmpPeAlignSamStr = peAlignInfo->getTmpPEreadAlignInfoInSAMformat(
			 	readNameString_1, readNameString_2, 
				readSeq_1.substr(0,readLength_1), readSeq_2.substr(0,readLength_2));

			PeAlignSamStrVec[tmpOpenMP] = tmpPeAlignSamStr;
			
			string tmpPeAlignInfoStr = peAlignInfo->getTmpAlignInfo(
			 	readNameString_1, readNameString_2, 
				readSeq_1.substr(0,readLength_1), readSeq_2.substr(0,readLength_2),
				"*", "*");

			PeAlignInfoStrVec[tmpOpenMP] = tmpPeAlignInfoStr;

			/*if( ((peAlignInfo->norAlignmentInfo_PE_1).size() > 5) 
				|| ((peAlignInfo->rcmAlignmentInfo_PE_1).size() > 5)
				|| ((peAlignInfo->norAlignmentInfo_PE_2).size() > 5)
				|| ((peAlignInfo->rcmAlignmentInfo_PE_2).size() > 5))
			{
				cout << "Repeat: readName: " << readNameString_1 << endl;
			}*/

			//cout << "finish inserting result string to vec ..." << endl;
			///////////////////////////////////////////////////////////////

			delete(segInfo_Nor1);
			delete(fragInfo_Nor1);
			delete(segInfo_Nor2);
			delete(fragInfo_Nor2);
			delete(segInfo_Rcm1);
			delete(fragInfo_Rcm1);
			delete(segInfo_Rcm2);
			delete(fragInfo_Rcm2);

			delete(peAlignInfo);

			free(finalNorSegMapRangeStart);
			free(finalNorSegMapRangeEnd);
			free(finalNorSegMapLoc);
			free(finalRcmSegMapRangeStart);
			free(finalRcmSegMapRangeEnd);
			free(finalRcmSegMapLoc);	

			free(finalNorSegMapRangeStart_PE);
			free(finalNorSegMapRangeEnd_PE);
			free(finalNorSegMapLoc_PE);
			free(finalRcmSegMapRangeStart_PE);
			free(finalRcmSegMapRangeEnd_PE);
			free(finalRcmSegMapLoc_PE);

		} // read file end
		
		align_end = clock();
		align_cost = align_cost + align_end - align_begin;

		//cout << "finish fixing mid, turn: " << tmpTurn+1 << endl << endl;

		//cout << "start to output ... turn: " << tmpTurn+1 << endl;
		
		output_begin = clock();

		for(int tmp = 0; tmp < realRecordNum; tmp++)
		{
			//if(PeAlignSamStrVec[tmp] == "Ns")
			//	continue;
   			testClassFile_ofs << PeAlignSamStrVec[tmp] << endl;
			tmpAlignInfoForDebugFile_ofs << PeAlignInfoStrVec[tmp] << endl;
		}

		output_end = clock();
		output_cost = output_cost + output_end - output_begin;
		
		//cout << "finish output, turn: " << tmpTurn+1 << endl << endl;

	}

	fclose(fp_in);
	free(sa);free(lcpCompress);//free(child_up);free(child_down);free(child_next);
	free(childTab);free(chrom);
	overall_end = clock();

	cout << endl << "**********************************" << endl << "**********************************";
	double overall_time = (double)(overall_end - overall_begin)/CLOCKS_PER_SEC;
	double input_time = (double)input_cost/CLOCKS_PER_SEC;
	double align_time = (double)align_cost/CLOCKS_PER_SEC;
	double ouput_time = (double)output_cost/CLOCKS_PER_SEC;

	cout << endl << "overall_time = " << overall_time << endl;
	cout << endl << "input_time = " << input_time << endl;
	cout << endl << "align_time = " << align_time << endl;
	cout << endl << "ouput_time = " << ouput_time << endl;

	cout << endl << "**********************************" << endl << "**********************************";
	//cout << endl << "totalReadNum = " << read_num << endl;

    return 0;
} //end main
