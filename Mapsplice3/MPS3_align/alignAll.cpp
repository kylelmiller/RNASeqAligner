#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>

#include "spliceJunctionHash.h"
#include "sam2junc.h"
#include "fixPhase1.h"
#include "pairedEndRead.h"
#include "chromosome.h"
#include "align_info.h"

#define PreIndexSize 268435456

using namespace std;

#ifdef CAL_TIME
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
#endif

unsigned int PairedReadNum = 0, BothUnmappedReadNum = 0, UnpairedReadNum = 0;

void fixHeadTail_areaAndStringHash(PairedEndRead* peReadInfo, PE_Read_Alignment_Info* peAlignInfo,
	SpliceJunctionHash* SJ, SecondLevelChromosomeList* secondLevelChromosomeList, bool spliceJunctionHashExists);

void fixOneEndUnmapped(PairedEndRead* peReadInfo, PE_Read_Alignment_Info* peAlignInfo,
		SecondLevelChromosomeList* secondLevelChromosomeList);

int main(int argc, char**argv)
{
    if(argc != 9)
	{
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

	bool Do_Phase1_Only = true;
	bool outputAlignInfoAndSamForAllPairedAlignmentBool = false;
	bool removeAllIntermediateFilesBool = true;
	bool DoSam2JuncBool = true;
	bool load2ndLevelIndexBool = true;
	bool DoRemappingOnUnmapEndReadsBool = true;
	bool DoRemappingOnUnfixedHeadTailAlignmentBool = true;

	if(Do_Phase1_Only)
	{
		log_ofs << "Do_Phase1 only!" << endl;
		DoSam2JuncBool = false;
		load2ndLevelIndexBool = false;
		DoRemappingOnUnmapEndReadsBool = false;
		DoRemappingOnUnfixedHeadTailAlignmentBool = false;
	}	
	else
	{
		log_ofs << "Do_Phase1_Phase2! " << endl;
		DoSam2JuncBool = true;
		load2ndLevelIndexBool = true;
		DoRemappingOnUnmapEndReadsBool = true;
		DoRemappingOnUnfixedHeadTailAlignmentBool = true;
	}

	int firstMappingRecordNumber = 1000000;
	int normalRecordNum_fixOneEndUnmapped = 1000000;
	int normalRecordNum_fixHeadTail = 1000000;

	log_ofs << "firstMappingRecordNumber: " << firstMappingRecordNumber << endl;
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

	/////////////////  LOAD INDEX  /////////////////
	log_ofs << "Reads input 1: " << argv[1] << endl;
	log_ofs << "Reads input 2: " << argv[2] << endl;
	log_ofs << "output directory: " << argv[3] << endl;
	log_ofs << "threads_num: " << argv[4] << endl;
	log_ofs << "data format: " << argv[5] << endl;
	log_ofs << "wholeGenomeIndex: " << argv[6] << endl;
	log_ofs << "localIndex: " << argv[7] << endl;
	log_ofs << "chromsomeDir: " << argv[8] << endl;

    char *InputReadFile = argv[1];//read sample, exacted from fastq file every time
    char *InputReadFile_PE = argv[2];// another end read for pair-end reads

	omp_set_num_threads(atoi(argv[4]));

	string fastqOrFasta = argv[5];

	bool InputAsFastq;
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
		cout << "Incorrect input read format: Please specify 'Fastq' or 'Fasta'." << endl;
		log_ofs << "Incorrect input read format: Please specify 'Fastq' or 'Fasta'." << endl;
		exit(0);
	}
	
    string indexStr = argv[6];
    string chromDirStr = argv[8];
    string secondLevelIndexStr = argv[7];

    indexStr.append("/");
    chromDirStr.append("/");
    secondLevelIndexStr.append("/");

    string preIndexArrayPreStr = indexStr;

	//// LOAD INDEX ///////
	#ifdef CAL_TIME
	read_file_begin = clock();
	#endif

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 

    log_ofs << "start to load preIndex ..." << endl;

	ifstream preIndexMapLengthArray_ifs((preIndexArrayPreStr + "_MapLength").c_str(), ios::binary);
	int* preIndexMapLengthArray = (int*)malloc(PreIndexSize * sizeof(int));
	preIndexMapLengthArray_ifs.read((char*)preIndexMapLengthArray, PreIndexSize * sizeof(int));

	ifstream preIndexIntervalStartArray_ifs((preIndexArrayPreStr + "_IntervalStart").c_str(), ios::binary);
	unsigned int *preIndexIntervalStartArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int));
    preIndexIntervalStartArray_ifs.read((char*)preIndexIntervalStartArray, PreIndexSize * sizeof(int));

	ifstream preIndexIntervalEndArray_ifs((preIndexArrayPreStr + "_IntervalEnd").c_str(), ios::binary);
	unsigned int *preIndexIntervalEndArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int));
    preIndexIntervalEndArray_ifs.read((char*)preIndexIntervalEndArray, PreIndexSize * sizeof(int));

 	log_ofs << "finish loading preIndex ..." << endl;

 	ifstream SA_file_ifs((indexStr + "_SA").c_str(),ios::binary);
	ifstream lcpCompress_file_ifs((indexStr + "_lcpCompress").c_str(),ios::binary);
	ifstream childTab_file_ifs((indexStr + "_childTab").c_str(),ios::binary);
	ifstream verifyChild_file_ifs((indexStr + "_detChild").c_str(),ios::binary);
	ifstream chrom_bit_file_ifs((indexStr + "_chrom").c_str(),ios::binary);
	ifstream parameter_file_ifs((indexStr + "_parameter").c_str(),ios::binary);

	Index_Info* indexInfo = new Index_Info(parameter_file_ifs);

	log_ofs << "index: " << indexStr << endl;
	///////////////////////////////////////
 
	log_ofs << "start to load whole genome" << endl;
	char *chrom = (char*)malloc((indexInfo->indexSize) * sizeof(char));
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->indexSize) * sizeof(char)); 

	indexInfo->chromString = chrom;
	log_ofs << "chromSize = " <<(indexInfo->chromString).size() << endl;
	
	log_ofs << "start to load every chromosome" << endl;

	(indexInfo->chromStr).push_back((indexInfo->chromString).substr(0, (indexInfo->chrEndPosInGenome)[0]+1));
	(indexInfo->chromLength).push_back(((indexInfo->chrEndPosInGenome)[0]+1));
	for(int tmp = 1; tmp < indexInfo->chromNum; tmp++)
	{
		(indexInfo->chromStr).push_back((indexInfo->chromString).substr((indexInfo->chrEndPosInGenome)[tmp-1]+2, 
			(indexInfo->chrEndPosInGenome)[tmp]-(indexInfo->chrEndPosInGenome)[tmp-1]-1));	
		(indexInfo->chromLength).push_back(((indexInfo->chrEndPosInGenome)[tmp]-(indexInfo->chrEndPosInGenome)[tmp-1]-1));
	}	

	log_ofs << "start to load SA" << endl;
    unsigned int *sa = (unsigned int*)malloc((indexInfo->indexSize) * sizeof(unsigned int));
    SA_file_ifs.read((char*)sa, (indexInfo->indexSize) * sizeof(unsigned int));

	log_ofs << "start to load lcpCompress" << endl;
	BYTE *lcpCompress = (BYTE*)malloc((indexInfo->indexSize) * sizeof(BYTE));
	lcpCompress_file_ifs.read((char*)lcpCompress, (indexInfo->indexSize) * sizeof(BYTE));	

	log_ofs << "start to load childTab " << endl;
	unsigned int *childTab = (unsigned int*)malloc((indexInfo->indexSize) * sizeof(unsigned int));
	childTab_file_ifs.read((char*)childTab, (indexInfo->indexSize) * sizeof(unsigned int));

	log_ofs << "start to load detChild" << endl;
	BYTE *verifyChild = (BYTE*)malloc((indexInfo->indexSize) * sizeof(BYTE));
	verifyChild_file_ifs.read((char*)verifyChild, (indexInfo->indexSize) * sizeof(BYTE));
	
	Chromosome* alignmentChromosome = new Chromosome(
		sa,
		lcpCompress,
		childTab,
		verifyChild,
		chrom,
		indexInfo,
		preIndexMapLengthArray,
		preIndexIntervalStartArray,
		preIndexIntervalEndArray);

	log_ofs << "All index files loaded" << endl;
	
	#ifdef CAL_TIME
	read_file_end = clock();
	double read_file_time = (double)(read_file_end - read_file_begin)/CLOCKS_PER_SEC;
	log_ofs << "read_file cpu time = " << read_file_time << endl;
	#endif
	//////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... whole genome index loaded ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... whole genome index loaded ......" << endl << endl;

	/////// Finished Loading Index //////////

	/////// 1st Mapping Process /////////////

   	string testClassFileStr = outputDirStr + "/output.sam";
	string tmpAlignInfoForDebug = testClassFileStr+ ".alignInfoDebug";
	string tmpAlignUnmap = testClassFileStr + ".unmap";
	string tmpAlignCompleteRead = testClassFileStr + ".completePair.sam";
	string tmpAlignBothEndsUnmapped = testClassFileStr + ".bothEndsUnmapped.sam";
	string tmpAlignIncompletePair = testClassFileStr + ".incomplete.alignInfo";
	string tmpAlignIncompletePair_SAM = testClassFileStr + ".incompletePair.sam";
	string tmpAlignCompleteRead_alignInfo = testClassFileStr + ".completePair.sam_alignInfo";
	string tmpAlignOneEndUnmapped = testClassFileStr + ".oneEndUnmapped";
	
	tmpAlignOneEndUnmapped += Do_Phase1_Only
		? ".sam"
		: ".alignInfo";

	ofstream tmpAlignOneEndUnmapped_ofs(tmpAlignOneEndUnmapped.c_str());
	ofstream tmpAlignBothEndsUnmapped_ofs(tmpAlignBothEndsUnmapped.c_str());
	ofstream tmpAlignIncompletePair_ofs(tmpAlignIncompletePair.c_str());
	ofstream tmpAlignIncompletePair_SAM_ofs(tmpAlignIncompletePair_SAM.c_str());
	ofstream tmpAlignCompleteRead_alignInfo_ofs(tmpAlignCompleteRead_alignInfo.c_str());
	ofstream tmpAlignCompleteRead_ofs(tmpAlignCompleteRead.c_str());

	string tmpIntermediateJunctionFile = testClassFileStr + ".inter.junc";

	ifstream inputRead_ifs(InputReadFile);
	ifstream inputRead_PE_ifs(InputReadFile_PE);

	vector<PairedEndRead*> pairedEndReads;

	vector<string> PeAlignSamStrVec_complete(firstMappingRecordNumber);
	vector<string> PeAlignInfoStrVec_inCompletePair(firstMappingRecordNumber);
	vector<string> PeAlignInfoStrVec_oneEndUnmapped(firstMappingRecordNumber);
	vector<string> PeAlignSamStrVec_bothEndsUnmapped(firstMappingRecordNumber);
	vector<string> PeAlignSamStrVec_inCompletePair(firstMappingRecordNumber);
	vector<string> PeAlignInfoStrVec_completePaired(firstMappingRecordNumber);

	nowtime = time(NULL);

	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... 1st mapping process starts ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... 1st mapping process starts ......" << endl << endl; 

	int readLengthMax_tmp = 0;
	string line;
	int realRecordNum;
	int readTotalNum = 0;
	bool endOfRecord = false;

	while(!endOfRecord)
	{
		#ifdef CAL_TIME		
		input_begin = clock();
		#endif

		for(realRecordNum = 0; realRecordNum < firstMappingRecordNumber; readTotalNum++, realRecordNum++)
		{
    		if((inputRead_ifs.eof())||(inputRead_PE_ifs.eof()))
    		{
				endOfRecord = true;
				break;    			
    		}

    		getline(inputRead_ifs, line); // readName_1
    		
    		if((inputRead_ifs.eof())||(inputRead_PE_ifs.eof()))
    		{
				endOfRecord = true;
				break;    			
    		}
    		
    		string readNameOne = line.substr(1);

    		getline(inputRead_ifs, line); // readSeq_1
    		string readSequenceOne = Utilities::toUpper(line);

    		if(InputAsFastq)
    		{
    			getline(inputRead_ifs, line);
    			getline(inputRead_ifs, line);
    		}

    		getline(inputRead_PE_ifs, line);
    		string readNameTwo = line.substr(1);

    		getline(inputRead_PE_ifs, line);
    		string readSequenceTwo = Utilities::toUpper(line);

    		if(InputAsFastq)
    		{
    			getline(inputRead_PE_ifs, line);
    			getline(inputRead_PE_ifs, line);
    		}

    		pairedEndReads.push_back(new PairedEndRead(
				readNameOne,
				readNameTwo,
				readSequenceOne,
				readSequenceTwo));
		}

		#ifdef CAL_TIME	
		input_end = clock();
		input_cost = input_cost + input_end - input_begin;
		align_begin = clock();
		#endif

		#pragma omp parallel for
		for(int i=0; i<realRecordNum; i++)
		{
			#ifdef CAL_TIME
			getReadInfo_begin = clock();
			#endif

			#ifdef CAL_TIME 
			getReadInfo_end = clock();
			getReadInfo_cost = getReadInfo_cost + getReadInfo_end - getReadInfo_begin;

    		segMap_begin = clock();
    		#endif

			FixPhase1Info* fixPhase1Info = new FixPhase1Info(pairedEndReads[i], alignmentChromosome);

			#ifdef CAL_TIME
			segMap_end = clock();
			segMap_cost = segMap_cost + segMap_end - segMap_begin;

			getPath_begin = clock();
			#endif

			fixPhase1Info->fixPhase1_pathInfo();

			#ifdef CAL_TIME
			getPath_end = clock();
			getPath_cost = getPath_cost + getPath_end - getPath_begin;

			fixGap_begin = clock();
			#endif

			fixPhase1Info->fixPhase1_gapInfo(indexInfo);

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
			
			PeAlignSamStrVec_complete[i] = "";
			PeAlignInfoStrVec_inCompletePair[i] = "";
			PeAlignInfoStrVec_oneEndUnmapped[i] = "";
			PeAlignSamStrVec_bothEndsUnmapped[i] = "";
			PeAlignSamStrVec_inCompletePair[i] = "";

			if(outputAlignInfoAndSamForAllPairedAlignmentBool)
				PeAlignInfoStrVec_completePaired[i] = "";

			if(peAlignInfo->finalPairExistsBool()) // some pair exists
			{
				bool allAlignmentCompleteBool = peAlignInfo->allAlignmentInFinalPairCompleted();
				if(allAlignmentCompleteBool)
				{
					PeAlignSamStrVec_complete[i] =
						peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(pairedEndReads[i]);

					if(outputAlignInfoAndSamForAllPairedAlignmentBool)
					{
						PeAlignInfoStrVec_completePaired[i] =
							peAlignInfo->getTmpAlignInfoForFinalPair(pairedEndReads[i]);
					}
				}
				else
				{
					if(!Do_Phase1_Only)
					{
						PeAlignInfoStrVec_inCompletePair[i] =
							peAlignInfo->getTmpAlignInfoForFinalPair(pairedEndReads[i]);
					}

					PeAlignSamStrVec_inCompletePair[i] =
						peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(pairedEndReads[i]);
				}
			}
			else // no pair exists: 1.one end unmapped; 2. both ends unmapped
			{
				if(peAlignInfo->alignInfoExistsBool()) // one end unmapped
				{	
					PeAlignInfoStrVec_oneEndUnmapped[i] =
						peAlignInfo->getTmpAlignInfo(pairedEndReads[i]);

				}
				else // both ends unmapped
				{
					PeAlignSamStrVec_bothEndsUnmapped[i] =
						peAlignInfo->getSAMformatForBothEndsUnmapped(pairedEndReads[i]);
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

			delete fixPhase1Info;
			delete pairedEndReads[i];
			delete peAlignInfo;

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
			}			

			if(PeAlignInfoStrVec_inCompletePair[tmp] != "")
				tmpAlignIncompletePair_ofs << PeAlignInfoStrVec_inCompletePair[tmp] << endl;
			
			if(PeAlignInfoStrVec_oneEndUnmapped[tmp] != "")
				tmpAlignOneEndUnmapped_ofs << PeAlignInfoStrVec_oneEndUnmapped[tmp] << endl;
			
			if(PeAlignSamStrVec_bothEndsUnmapped[tmp] != "")
			{
				tmpAlignBothEndsUnmapped_ofs << PeAlignSamStrVec_bothEndsUnmapped[tmp] << endl;
				BothUnmappedReadNum ++;
			}

			if(PeAlignSamStrVec_inCompletePair[tmp] != "")
				tmpAlignIncompletePair_SAM_ofs << PeAlignSamStrVec_inCompletePair[tmp] << endl;

			if(outputAlignInfoAndSamForAllPairedAlignmentBool && PeAlignInfoStrVec_completePaired[tmp] != "")
				tmpAlignCompleteRead_alignInfo_ofs << PeAlignInfoStrVec_completePaired[tmp] << endl;
		}

		#ifdef CAL_TIME
		output_end = clock();
		output_cost = output_cost + output_end - output_begin;
		#endif

	}

	log_ofs << "readTotalNum: " << readTotalNum << endl;

	inputRead_ifs.close();
	inputRead_PE_ifs.close();
	tmpAlignCompleteRead_ofs.close();
	tmpAlignOneEndUnmapped_ofs.close();
	tmpAlignBothEndsUnmapped_ofs.close();
	tmpAlignIncompletePair_SAM_ofs.close();
	if(Do_Phase1_Only)
	{
		tmpAlignIncompletePair_ofs.close();
	}

	delete alignmentChromosome;
	
	#ifdef CAL_TIME
	overall_end = clock();
	#endif
	
	nowtime = time(NULL);
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

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////    	Load Second Level Index      ////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);

	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load 2nd level index ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... load 2nd level index starts ......" << endl << endl; 

	vector<SecondLevelChromosome*> secondLevelChromVector;

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
				
				string inputIndexFileStr = secondLevelIndexStr +
					"/" +
					indexInfo->chrNameStr[tmpChrNO] +
					"/"	+
					tmpFileNumStr +
					"/";

				// Name of all of the second level index files
				ifstream secondLevelChrom_file_ifs((inputIndexFileStr + "chrom").c_str(), ios::binary);
				ifstream secondLevelSA_file_ifs((inputIndexFileStr + "SA").c_str(), ios::binary);
				ifstream secondLevelLcpCompress_file_ifs((inputIndexFileStr + "_lcpCompress").c_str(), ios::binary);
				ifstream secondLevelChildTab_file_ifs((inputIndexFileStr + "childTab").c_str(), ios::binary);
				ifstream secondLevelDetChild_file_ifs((inputIndexFileStr + "detChild").c_str(), ios::binary);

				int sizeOfIndex = indexInfo->secondLevelIndexNormalSize + 1;
				char* tmpSecondLevelChrom = (char*)malloc(sizeOfIndex * sizeof(char));
				for(int tmpMallocSpace = 0; tmpMallocSpace < sizeOfIndex; tmpMallocSpace++)
					tmpSecondLevelChrom[tmpMallocSpace] = '0';

				secondLevelChrom_file_ifs.read((char*)tmpSecondLevelChrom, sizeOfIndex * sizeof(char));
				if(tmpSecondLevelChrom[sizeOfIndex-1] != 'X')
					(indexInfo->invalidSecondLevelIndexNOset).insert(secondLevelIndexNO + 1);

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
					(indexInfo->invalidSecondLevelIndexNOset).insert(secondLevelIndexNO + 1);

				unsigned int* tmpSecondLevelSa = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
				secondLevelSA_file_ifs.read((char*)tmpSecondLevelSa, sizeOfIndex * sizeof(unsigned int));

				BYTE* tmpSecondLevelLcpCompress = (BYTE*)malloc(sizeOfIndex * sizeof(BYTE));
				secondLevelLcpCompress_file_ifs.read((char*)tmpSecondLevelLcpCompress, sizeOfIndex * sizeof(BYTE));

				unsigned int* tmpSecondLevelChildTab = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
				secondLevelChildTab_file_ifs.read((char*)tmpSecondLevelChildTab, sizeOfIndex * sizeof(unsigned int));

				BYTE* tmpSecondLevelDetChild = (BYTE*)malloc(sizeOfIndex * sizeof(BYTE));
				secondLevelDetChild_file_ifs.read((char*)tmpSecondLevelDetChild, sizeOfIndex * sizeof(BYTE));

				secondLevelChromVector.push_back(new SecondLevelChromosome(
					tmpSecondLevelChrom,
					tmpSecondLevelSa,
					tmpSecondLevelLcpCompress,
					tmpSecondLevelChildTab,
					tmpSecondLevelDetChild,
					indexInfo));

				secondLevelChrom_file_ifs.close();
				secondLevelSA_file_ifs.close();
				secondLevelLcpCompress_file_ifs.close();
				secondLevelChildTab_file_ifs.close();
				secondLevelDetChild_file_ifs.close();

				secondLevelIndexNO ++;

			} // for(int tmpSecondLevelIndexNO = 1; tmpSecondLevelIndexNO <= (indexInfo->secondLevelIndexPartsNum)[tmpChrNO]; tmpSecondLevelIndexNO ++)

			log_ofs << "finish loading 2nd-level index of " << indexInfo->chrNameStr[tmpChrNO] << endl;

		} // for(int tmpChrNO = 0; tmpChrNO < indexInfo->chromNum; tmpChrNO ++)

		log_ofs << "finish loading ALL 2nd-level index !" << endl;
		log_ofs << indexInfo->getInvalidSecondLevelIndexNOstr() << endl;

	} // if(load2ndLevelIndexBool)

	SecondLevelChromosomeList* secondLevelChroms = new SecondLevelChromosomeList(
		secondLevelChromVector,
		indexInfo);

	nowtime = time(NULL);

	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... load 2nd level index ends ......" << endl << endl ; 		
	log_ofs << endl << "[" << asctime(local) << "... load 2nd level index ends ......" << endl << endl ;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Do REMAPPING On one end unmapped Reads    ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
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

		ifstream inputRecord_ifs(tmpAlignOneEndUnmapped.c_str());

		string line1, line2, line3, line4, line5, line6, line7, 
			line8, line9, line10, line11;
		
		int normalRecordNum = normalRecordNum_fixOneEndUnmapped; //1000000;

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

		bool endOfRecord = false;
		while(!endOfRecord)
		{
			int recordNum = normalRecordNum;

			realRecordNum = normalRecordNum;

			for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
			{
				if(inputRecord_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					endOfRecord = true;
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

			#pragma omp parallel for
			for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
			{
				////////////////  parse long head reads record after 1-mapping process  ///////////////////////////////////////
				tmpRecordNum_oneEndUnmapped ++;

				PairedEndRead* peReadInfo = new PairedEndRead();
				PE_Read_Alignment_Info* peAlignInfo = new PE_Read_Alignment_Info();
				peAlignInfo->generatePeReadInfoAndPeAlignInfo_Fasta_toFixOneEndUnmapped_getline(line1StrVec[tmpOpenMP], line2StrVec[tmpOpenMP], 
					line4StrVec[tmpOpenMP], line5StrVec[tmpOpenMP],
					line7StrVec[tmpOpenMP], line8StrVec[tmpOpenMP], //line9StrVec[tmpOpenMP],
					line9StrVec[tmpOpenMP], line10StrVec[tmpOpenMP], peReadInfo);		

				fixOneEndUnmapped(peReadInfo, peAlignInfo, secondLevelChroms);

				peAlignInfo->pairingAlignment();
				peAlignInfo->chooseBestAlignment();

				bool pairExistsBool = peAlignInfo->finalPairExistsBool();
				bool allAlignmentCompleteBool = peAlignInfo->allAlignmentInFinalPairCompleted();
				bool allUnpairedAlignmentCompleteBool = peAlignInfo->allUnpairedAlignmentCompleted();

				string tmpPeAlignSamStr, tmpPeAlignInfoStr, tmpPeAlignSamStr_unpair_complete, tmpPeAlignInfo_complete_pair;
				if(pairExistsBool && allAlignmentCompleteBool) // some pair exists, all completed, print out paired SAM info
				{
					tmpPeAlignSamStr = peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(peReadInfo);
					tmpPeAlignInfoStr = "";
					tmpPeAlignSamStr_unpair_complete = "";
					if(outputAlignInfoAndSamForAllPairedAlignmentBool)
					{
						tmpPeAlignInfo_complete_pair = peAlignInfo->getTmpAlignInfoForFinalPair(peReadInfo);
					}
				}
				else if(pairExistsBool && (!allAlignmentCompleteBool)) // pair exists, incomplete
				{
					tmpPeAlignSamStr = "";
					tmpPeAlignInfoStr = peAlignInfo->getTmpAlignInfoForFinalPair(peReadInfo);
					tmpPeAlignSamStr_unpair_complete = "";
					if(outputAlignInfoAndSamForAllPairedAlignmentBool)
					{					
						tmpPeAlignInfo_complete_pair = "";
					}
				}
				else if((!pairExistsBool) && (allUnpairedAlignmentCompleteBool)) // no pair exists, all complete, print out original SAM info
				{
					tmpPeAlignSamStr = "";
					tmpPeAlignInfoStr = "";

					tmpPeAlignSamStr_unpair_complete =
						peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot(peReadInfo);
					if(outputAlignInfoAndSamForAllPairedAlignmentBool)
					{					
						tmpPeAlignInfo_complete_pair = "";
					}
				}
				else // no pair exists, incomplete, print out alignInfo
				{
					tmpPeAlignSamStr = "";
					tmpPeAlignInfoStr = peAlignInfo->getTmpAlignInfo(peReadInfo);// << endl;
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
				delete peReadInfo;
				delete peAlignInfo;

			}

			for(int tmp = 0; tmp < realRecordNum; tmp++)
			{
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
		}		
		inputRecord_ifs.close();
	}

	OutputSamFile_oneEndMapped_ofs.close();
	OutputSamFile_oneEndMapped_unpairComplete_ofs.close();
	tmpAlignIncompletePair_ofs.close();

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads ends ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads ends ......" << endl << endl ; 

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Sam 2 Junc   ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
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

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... Sam 2 Junc ends ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... Sam 2 Junc ends ......" << endl << endl; 

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Do REMAPPING On unfixed head/tail Reads    ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
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

	SpliceJunctionHash* SJ = new SpliceJunctionHash();
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

		SJ->initiateSJintHash(indexInfo->chromNum);

		fgets(entry, sizeof(entry), fp_spliceJunction);
		while(!feof(fp_spliceJunction))
		{
			fgets(entry, sizeof(entry), fp_spliceJunction);
			if(feof(fp_spliceJunction))
				break;
			junctionNum ++;

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
		ifstream inputUnfixedHeadTailRecord_ifs(headTailSoftClippingFile.c_str());

		string line1, line2, line3, line4, line5, line6, line7, 
			line8, line9, line10, line11;

		int normalRecordNum = normalRecordNum_fixHeadTail; //1000000;

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
		vector<string> peAlignSamVec_complete_pair(normalRecordNum);
		vector<string> peAlignSamVec_incomplete_pair(normalRecordNum);
		vector<string> peAlignSamVec_complete_unpair(normalRecordNum);
		vector<string> peAlignSamVec_incomplete_unpair(normalRecordNum);

		vector<string> peAlignSamVec_complete_pair_alignInfo(normalRecordNum);
		vector<string> peAlignSamVec_incomplete_pair_alignInfo(normalRecordNum);

		bool endOfRecord = false;
		while(!endOfRecord)
		{
			int recordNum = normalRecordNum;

			//cout << "start to read long-head file record" << endl;
			//cout << "start to read Head/Tail file record, turn: " << tmpTurn+1 << endl;
			realRecordNum = normalRecordNum;

			for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
			{

				if(inputUnfixedHeadTailRecord_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					endOfRecord = true;
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

			#pragma omp parallel for
			for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
			{
				////////////////  parse long head reads record after 1-mapping process  ///////////////////////////////////////

				PairedEndRead* peReadInfo = new PairedEndRead();
				PE_Read_Alignment_Info* peAlignInfo = new PE_Read_Alignment_Info();
				peAlignInfo->generatePeReadInfoAndPeAlignInfo_Fasta_toFixIncompleteAlignment_getline(line1StrVec[tmpOpenMP], line2StrVec[tmpOpenMP], 
					line4StrVec[tmpOpenMP], line5StrVec[tmpOpenMP],
					line7StrVec[tmpOpenMP], line8StrVec[tmpOpenMP], //line9StrVec[tmpOpenMP],
					line9StrVec[tmpOpenMP], line10StrVec[tmpOpenMP], peReadInfo);		

				fixHeadTail_areaAndStringHash(peReadInfo, peAlignInfo, SJ,
					secondLevelChroms, spliceJunctionHashExists);

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
						tmpPeAlignSamStr_complete_pair = peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(peReadInfo);
						tmpPeAlignSamStr_incomplete_pair = "";
						
						tmpPeAlignSamStr_complete_unpair = "";
						
						tmpPeAlignSamStr_incomplete_unpair = "";
						if(outputAlignInfoAndSamForAllPairedAlignmentBool)
						{
							tmpPeAlignSamStr_complete_pair_alignInfo =
								peAlignInfo->getTmpAlignInfoForFinalPair(peReadInfo);
							tmpPeAlignSamStr_incomplete_pair_alignInfo = "";
						}
					}
					else
					{
						tmpPeAlignSamStr_complete_pair = "";
						
						tmpPeAlignSamStr_incomplete_pair = peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(peReadInfo);
						tmpPeAlignSamStr_complete_unpair = "";
						
						tmpPeAlignSamStr_incomplete_unpair = "";

						if(outputAlignInfoAndSamForAllPairedAlignmentBool)
						{
							tmpPeAlignSamStr_complete_pair_alignInfo = "";
							tmpPeAlignSamStr_incomplete_pair_alignInfo =
								peAlignInfo->getTmpAlignInfoForFinalPair(peReadInfo);
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

						tmpPeAlignSamStr_complete_unpair =
							peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot(peReadInfo);

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

						tmpPeAlignSamStr_incomplete_unpair =
							peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot(peReadInfo);
						if(outputAlignInfoAndSamForAllPairedAlignmentBool)
						{
							tmpPeAlignSamStr_complete_pair_alignInfo = "";
							tmpPeAlignSamStr_incomplete_pair_alignInfo = "";
						}
					}
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

				delete peReadInfo;
				delete peAlignInfo;
			}

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
		}
		inputUnfixedHeadTailRecord_ifs.close();
	}
	delete SJ;

	OutputSamFile_fixHeadTail_complete_pair_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_pair_ofs.close();
	OutputSamFile_fixHeadTail_complete_unpair_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_unpair_ofs.close();
	OutputSamFile_fixHeadTail_complete_pair_alignInfo_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_pair_alignInfo_ofs.close();

	nowtime = time(NULL);
	local = localtime(&nowtime);


	cout << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads ends ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads ends ......" << endl << endl ;  

	if(!Do_Phase1_Only)
	{
		cout << "readTotalNum: " << readTotalNum << endl;		
		log_ofs << "readTotalNum: " << readTotalNum << endl;
	}			

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to prepare for final output files ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... start to prepare for final output files ......" << endl << endl ;  

	string finalOutputSam = testClassFileStr;
	if(Do_Phase1_Only)
	{
		string cat_cmd = "cat " + tmpAlignCompleteRead
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

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... all jobs done ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... all jobs done ......" << endl << endl ;  

    return 0;
} //end main

void fixHeadTail_areaAndStringHash(PairedEndRead* peReadInfo, PE_Read_Alignment_Info* peAlignInfo,
	SpliceJunctionHash* SJ, SecondLevelChromosomeList* secondLevelChromosomeList, bool spliceJunctionHashExists)
{
	Alignment_Info* tmpAlignmentInfo;

	//Read currentRead;

	/* FIXME - THIS NEEDS TO BE REWORKED
	 * 5/28/14 KLM
	//////////////////// fix head //////////////////////////////
	//cout << "start to fix head" << endl;
	for(int i = 1; i <= 4; i++)
	{
		currentRead = i <= 2
			? peReadInfo->firstPairedEndRead
			: peReadInfo->secondPairedEndRead;

		for(int tmpAlignmentNO = 0; tmpAlignmentNO < (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
			tmpAlignmentNO++)
		{
			tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
			if( (tmpAlignmentInfo->cigarStringJumpCode)[0].type != "S" )
			{
				//cout << "no unfixed head !" << endl;
				continue;
			}

			Unfixed_Head unfixedHeadInfo;
			unfixedHeadInfo.getUnfixedHeadInfoFromRecordWithAlignInfoType(peReadInfo, tmpAlignInfoType, tmpAlignmentInfo, indexInfo);

			//cout << "finish getUnfixedTailInfoFromRecord ..." << endl;

			string readSeqWithDirection = unfixedHeadInfo.alignDirection == "+"
				? unfixedHeadInfo.readSeqOriginal
				: covertStringToReverseComplement(unfixedHeadInfo.readSeqOriginal);

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
				bool headSegMapMain = mapMainSecondLevelForTargetMapping_compressedIndex(headChar,
					secondLevelSa[midPartMapPosSecondLevelIndexNO],
					secondLevelLcpCompress[midPartMapPosSecondLevelIndexNO],
					secondLevelChildTab[midPartMapPosSecondLevelIndexNO],
					secondLevelDetChild[midPartMapPosSecondLevelIndexNO],
					secondLevelChrom[midPartMapPosSecondLevelIndexNO],
					targetMappingStr.length(), 3000000, unfixedHeadInfo.midPartMapPosInWholeGenome,
					midPartMapPosForLongHeadInSecondLevelIndex, &targetMappingNum, targetMappingLoc,
					indexInfo);

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
*/
}

void fixOneEndReadForward(Read* read, Read* incompleteEndRead, PE_Read_Alignment_Info* peAlignInfo,
	vector<Alignment_Info*> alignmentVector, SecondLevelChromosomeList* secondLevelChromosomeList)
{
	string chrNameStr;
	int chrMapPos_start;
	int chrMapPos_end;

	int mapPosIntervalStart;
	int mapPosIntervalEnd;
	int secondLevelIndexNum; // start from 1

	int chrPosStartIn2ndLevelIndex;

	for(int i = 0; i < alignmentVector.size();i++)
	{
		Alignment_Info* alignmentInfo = alignmentVector[i];

		mapPosIntervalStart = alignmentInfo->alignChromPos;
		mapPosIntervalEnd = alignmentInfo->getEndMatchedPosInChr() + READ_ALIGN_AREA_LENGTH;

		SecondLevelChromosome* secondLevelChromosome =
			secondLevelChromosomeList->getSecondLevelChromosome(
					alignmentInfo->alignChromName,
					alignmentInfo->alignChromPos);

		if(secondLevelChromosome == NULL)
			continue;

		MappedRead* mappedRead = new MappedRead(read, secondLevelChromosome);

		Path* pathInfo = new Path();
		pathInfo->getPossiPathFromSeg(mappedRead);

		int pathValidNum = pathInfo->pathValidNumInt();
		if(pathValidNum > 10)
		{
			delete pathInfo;
			delete mappedRead;
			continue;
		}

		Gap* gapInfo = new Gap();
		gapInfo->fixGapInPath(pathInfo, mappedRead, secondLevelChromosome->getIndexInfo());

		/* FIXME - FIX THIS LATER 5/28/14 KLM
		if(pathInfo->finalPathVec.size() == 1)
			peAlignInfo->pushBackPathInfo2PeAlignInfo(pathInfo, End1OrEnd2, NorOrRcm, indexInfo);

		peAlignInfo->pushBackPathInfo2PeAlignInfo(pathInfo, End1OrEnd2, NorOrRcm, indexInfo);
		 */
		delete gapInfo;
		delete pathInfo;
		delete mappedRead;
	}
}

void fixOneEndUnmapped(PairedEndRead* peReadInfo, PE_Read_Alignment_Info* peAlignInfo,
		SecondLevelChromosomeList* secondLevelChromosomeList)
{
	// Fix the first read's ends
	fixOneEndReadForward(peReadInfo->getFirstRead(), peReadInfo->getSecondReadReverseComplement(),
		peAlignInfo, peAlignInfo->norAlignmentInfo_PE_1, secondLevelChromosomeList); // End1OrEnd2=true, NorOrRcm=true

	fixOneEndReadForward(peReadInfo->getSecondReadReverseComplement(), peReadInfo->getFirstRead(),
		peAlignInfo, peAlignInfo->rcmAlignmentInfo_PE_2, secondLevelChromosomeList); // End1OrEnd2=false, NorOrRcm=false
}

