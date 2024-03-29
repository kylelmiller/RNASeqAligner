#include <stdlib.h>
#include <string>
#include <string.h>
#include <map>

using namespace std;

#define LONGLCP 255

typedef unsigned char BYTE;

class BuildIndex_Info
{
public:
	
	unsigned int genomeLength;

	int chromNum;

	vector<string> chrNameStr; // size = chromNum

	vector<unsigned int> chrEndPosInGenome;

	map<string, int> chrNameMap;
	
	int secondLevelIndexNormalSize;
	vector<int> secondLevelIndexPartsNum;
	int secondLevelIndexPartsNumSum;

	unsigned int null_num; // 2654911540 for mm9_noRandom genome
	unsigned int indexSize; //2654911539  //sequence length + 1, the length of sa-lcp-down-next 
	//int NULL_NUM; // 2654911540 for mm9_noRandom genome
	//int MAX; //2654911539  //sequence length + 1, the length of sa-lcp-down-next 

	BuildIndex_Info()
	{}

	void getIndexInfoFromParameterFile(ifstream& inputIndexInfoFile)
	{
		string s;
		string chromNumLine;
		string chromNameLine;
		string chromEndPosInGenomeLine;
		string secondLevelIndexSizeLine;
		string chrom2ndLevelIndexNumLine;
		getline(inputIndexInfoFile, s);
		getline(inputIndexInfoFile, chromNumLine);
		getline(inputIndexInfoFile, s);
		getline(inputIndexInfoFile, chromNameLine);
		getline(inputIndexInfoFile, s);
		getline(inputIndexInfoFile, chromEndPosInGenomeLine);
		getline(inputIndexInfoFile, s);
		getline(inputIndexInfoFile, secondLevelIndexSizeLine);
		getline(inputIndexInfoFile, s);
		getline(inputIndexInfoFile, chrom2ndLevelIndexNumLine);

		chromNum = atoi( (chromNumLine.substr(0, chromNumLine.length())).c_str() );
		cout << "chromNum: " << chromNum << endl;

		int startSearchPos = 0;
		int foundSearchPos;
		string tmpChromNameStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromNameLine.find(",", startSearchPos);
			tmpChromNameStr = chromNameLine.substr(startSearchPos+1, foundSearchPos - 2 - startSearchPos - 1 + 1);
			chrNameStr.push_back(tmpChromNameStr);
			//cout << tmp+1 << " tmpChromNameStr: " << tmpChromNameStr << " strLen: " << tmpChromNameStr.length() << endl;
			startSearchPos = foundSearchPos + 1;
		}

		startSearchPos = 0;
		string tmpChromEndPosStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromEndPosInGenomeLine.find(",", startSearchPos);
			tmpChromEndPosStr = chromEndPosInGenomeLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			unsigned int tmpChromEndPos = strtoul(tmpChromEndPosStr.c_str(), NULL, 10);
			chrEndPosInGenome.push_back(tmpChromEndPos);
			//cout << tmp+1 << " tmpChromEndPos: " << tmpChromEndPos << endl;
			startSearchPos = foundSearchPos + 1;
		}		

		secondLevelIndexNormalSize = atoi( (secondLevelIndexSizeLine.substr(0, secondLevelIndexSizeLine.length())).c_str() );
		cout << "secondLevelIndexNormalSize: " << secondLevelIndexNormalSize << endl;

		startSearchPos = 0;
		string tmpChrom2ndLevelIndexNumStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chrom2ndLevelIndexNumLine.find(",", startSearchPos);
			tmpChrom2ndLevelIndexNumStr = chrom2ndLevelIndexNumLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			int tmpChrom2ndLevelIndexNum = atoi(tmpChrom2ndLevelIndexNumStr.c_str());
			secondLevelIndexPartsNum.push_back(tmpChrom2ndLevelIndexNum);
			//cout << tmp+1 << "tmp2ndLevelIndexNum: " << tmpChrom2ndLevelIndexNum << endl;
			startSearchPos = foundSearchPos + 1;
		}

		indexSize = chrEndPosInGenome[chrEndPosInGenome.size()-1] + 1;
		//cout << "MAX: " << indexSize << endl;
		null_num = indexSize + 1;
		//cout << "NULL_NUM: " << null_num << endl;
		this->buildChrNameMap();
	}

	void buildChrNameMap()
	{
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			chrNameMap.insert(pair <string, int> (chrNameStr[tmp], tmp));
		}
	}

	void compressLcp2Lcpcompress(unsigned int *lcp, 
		BYTE *lcpCompress, unsigned int IndexLength)
	{
		for(unsigned int Num = 0; Num < IndexLength; Num++)
		{
			if(lcp[Num] < LONGLCP+1)
				lcpCompress[Num] = lcp[Num];
			else
			{
				lcpCompress[Num] = LONGLCP;
			}
		}
		return;
	}

	void compressUpDownNext2ChildtabVerifyChild(
		unsigned int *up, unsigned int *down, unsigned int *next,
		unsigned int *child, BYTE *verifyChild, 
		unsigned int IndexLength)	
	{
		/////////////////////  build ChildTable ////////////////
		//next -> child
		for(unsigned int Num = 0; Num < IndexLength; Num++)
		{
			child[Num] = next[Num];
		}
		//down -> child
		for(unsigned int Num = 0; Num < IndexLength; Num++)
		{
			if((down[Num] > 0) && /*(child[Num] == 0)*/(down[Num] != up[next[Num]]))
			{
				if(child[Num]!=0)
				{
					cout << "error child[Num]!= 0 " << endl;
					exit;
				}
				child[Num] = down[Num];
			}
		}
		//up -> child
		for(unsigned int Num = 0; Num < IndexLength; Num++)
		{
			if((up[Num] > 0) /*&&(child[Num] > 0)*/)
			{
				if(child[Num-1] != 0)
				{
					cout << "error child[Num-1] > 0" << endl;
					exit;
				}
				child[Num-1] = up[Num]; 
			}
		}

		/////////////////////  build VerifyChild ////////////////
		for(unsigned int index = 0; index < IndexLength; index++)
		{
   	    	verifyChild[index] = 0;
    		if(up[index] > 0) //up
    		{
    			verifyChild[index-1] = 1;
    		}
    		if((down[index] > 0)&&(next[index] == 0))  //down 
    		{	
    			verifyChild[index] = 2;
    		}  	  	
    		if(next[index] > 0)                //next
    		{
    			verifyChild[index] = 3;
    		}
			if((down[index] > 0)&&(next[index] > 0)) 
    		{
    			verifyChild[index] = 4;
    		}
		}

		return;
	}

	bool compressUpDownNext2ChildtabVerifyChild_bool(
		unsigned int *up, unsigned int *down, unsigned int *next,
		unsigned int *child, BYTE *verifyChild, 
		unsigned int IndexLength)	
	{
		bool compressCorrect = true;
		/////////////////////  build ChildTable ////////////////
		//next -> child
		cout << "next->child" << endl;
		for(unsigned int Num = 0; Num < IndexLength; Num++)
		{
			child[Num] = next[Num];
		}
		//down -> child
		cout << "down->child" << endl;
		for(unsigned int Num = 0; Num < IndexLength; Num++)
		{
			if((down[Num] > 0) && ((next[Num] < 0)||((next[Num] > IndexLength-1))))
			{
				compressCorrect = false;
				continue;
			}
			if((down[Num] > 0) && /*(child[Num] == 0)*/(down[Num] != up[next[Num]]))
			{
				if(child[Num]!=0)
				{
					//cout << "error child[Num]!= 0 " << endl;
					compressCorrect = false;
					//exit;
				}
				child[Num] = down[Num];
			}
		}
		//up -> child
		cout << "up->child" << endl;
		for(unsigned int Num = 0; Num < IndexLength; Num++)
		{
			if((up[Num] > 0) /*&&(child[Num] > 0)*/)
			{
				if(Num == 0)
				{
					compressCorrect = false;
					continue;
				}
				if(child[Num-1] != 0)
				{
					compressCorrect = false;
					//cout << "error child[Num-1] > 0" << endl;
					//exit;
				}
				child[Num-1] = up[Num]; 
			}
		}
		cout << "build VerifyChild" << endl;
		/////////////////////  build VerifyChild ////////////////
		for(unsigned int index = 0; index < IndexLength; index++)
		{
   	    	verifyChild[index] = 0;
    		if(up[index] > 0) //up
    		{
				if(index == 0)
				{
					compressCorrect = false;
					//continue;
				}
				else
				{
	    			verifyChild[index-1] = 1;
    			}
    		}
    		if((down[index] > 0)&&(next[index] == 0))  //down 
    		{	
    			verifyChild[index] = 2;
    		}  	  	
    		if(next[index] > 0)                //next
    		{
    			verifyChild[index] = 3;
    		}
			if((down[index] > 0)&&(next[index] > 0)) 
    		{
    			verifyChild[index] = 4;
    		}
		}

		return compressCorrect;
	}

	int convertStringToInt(const string& chrName)
	{
		int chrNameInt;
		map<string, int>::iterator chrNameMapIter;
		chrNameMapIter = chrNameMap.find(chrName);

		if(chrNameMapIter != chrNameMap.end())
		{
			chrNameInt = chrNameMapIter->second;
		}
		else
		{
			cout << "...... chrom name error! ...... " << endl;
		}
		return chrNameInt;
	}

};

class BuildIndexParameter_info
{
private:

	string InputChromFolderStr;
	string OutputIndexFolderStr;
	string wholeGenomeStr;

	set<string> chrNameStrSet;
	vector<string> chrNameStrVec;

	vector<unsigned int> chromLengthVec;
	vector<unsigned int> chrEndPosInGenomeVec;

	unsigned int secondLevelIndexSize;

public:

	// Setting up an iterator
	typedef std::vector<string>::iterator iterator;
	typedef std::vector<string>::const_iterator const_iterator;

	iterator begin() { return chrNameStrVec.begin(); }
	const_iterator begin() const { return chrNameStrVec.begin(); }
    iterator end() { return chrNameStrVec.end(); }
    const_iterator end() const { return chrNameStrVec.end(); }
    // End of iterator

    // Constructor
    // Input - input chromosome folder
    // Output - indexes generated from the chromosome
	BuildIndexParameter_info(char* inputChromFolder, char* outputIndexFolder)
	{
		secondLevelIndexSize = 3000000;
		InputChromFolderStr = inputChromFolder;
		OutputIndexFolderStr = ((string)outputIndexFolder) + "/";
	}

	void printChromLength(ofstream& OutputStream)
	{
		//string chromLengthStr;
		for(unsigned int tmpChr = 0; tmpChr < chromLengthVec.size(); tmpChr++)
		{
			OutputStream << chrNameStrVec[tmpChr] << " length: " 
				<< chromLengthVec[tmpChr] << endl; 
		}
	}

	void getChrEndPosVec()
	{
		unsigned int chrEndPosBase = 0;
		unsigned int tmpChr = 0;
		for (tmpChr = 0; tmpChr < chromLengthVec.size()-1; tmpChr++)
		{
			chrEndPosBase += chromLengthVec[tmpChr] + 1;
			chrEndPosInGenomeVec.push_back(chrEndPosBase - 2);
			//cout << "chrEndPos " << tmpChr+1 << ": " << chrEndPosBase - 2 << endl;
		}
		tmpChr = chromLengthVec.size()-1;
		chrEndPosBase += chromLengthVec[tmpChr] + 1;
		chrEndPosInGenomeVec.push_back(chrEndPosBase - 1);
		//cout << "chrEndPos " << tmpChr+1 << ": " << chrEndPosBase - 1 << endl;
	}

	void initiateAllOutputIndexFile()
	{
		string outputIndexFileStr = OutputIndexFolderStr + "/WholeGenomeIndex2";

		string log_file = outputIndexFileStr; log_file.append("_LOG"); ofstream log_file_ofs(log_file.c_str());

		string parameter_file = outputIndexFileStr; parameter_file.append("_parameter"); ofstream parameter_file_ofs(parameter_file.c_str());
		string chrom_bit_file = outputIndexFileStr; chrom_bit_file.append("_chromBit"); ofstream chrom_bit_file_ofs(chrom_bit_file.c_str(),ios::binary);

		string SA_file = outputIndexFileStr; SA_file.append("_SA"); ofstream SA_file_ofs(SA_file.c_str(),ios::binary); 
		string lcp_file = outputIndexFileStr; lcp_file.append("_lcp"); ofstream lcp_file_ofs(lcp_file.c_str(),ios::binary);
		string up_file = outputIndexFileStr; up_file.append("_up"); ofstream up_file_ofs(up_file.c_str(),ios::binary);
		string down_file = outputIndexFileStr; down_file.append("_down"); ofstream down_file_ofs(down_file.c_str(),ios::binary);
		string next_file = outputIndexFileStr; next_file.append("_next"); ofstream next_file_ofs(next_file.c_str(),ios::binary);

		string lcpCompress_file = outputIndexFileStr; lcpCompress_file.append("_lcpCompress"); ofstream lcpCompress_file_ofs(lcpCompress_file.c_str(),ios::binary);
		string longLcpIndex_file = outputIndexFileStr; longLcpIndex_file.append("_longLcpIndex"); ofstream longLcpIndex_file_ofs(longLcpIndex_file.c_str(),ios::binary);
		string longLcpValue_file = outputIndexFileStr; longLcpValue_file.append("_longLcpValue"); ofstream longLcpValue_file_ofs(longLcpValue_file.c_str(),ios::binary);
		string childTab_file = outputIndexFileStr; childTab_file.append("_childTab"); ofstream childTab_file_ofs(childTab_file.c_str(),ios::binary);
		string detChild_file = outputIndexFileStr; detChild_file.append("_detChild"); ofstream detChild_file_ofs(detChild_file.c_str(), ios::binary);
	}

	void outputChromSeq(ofstream& chromFileOfs)
	{
		//int chromBaseNum;
		unsigned int wholeGenomeChromBaseNum = 0;
		//int wholeGenomeCharacterNum = 0;
		//int wholeGenomeValuableCharacterNum = 0;
		int chromosomeNumber = chrNameStrVec.size();
		unsigned int wholeGenomeIndexSize = chrEndPosInGenomeVec[chromosomeNumber-1] + 1;
		//cout << endl << "wholeGenomeIndexSize: " << wholeGenomeIndexSize << endl;
		char* wholeGenomeSeqChar = (char*)malloc(wholeGenomeIndexSize * sizeof(char));

		for(unsigned int tmp = 0; tmp < chrNameStrVec.size(); tmp++)
		{
			int chromBaseNum = 0;
			int characterNum = 0;
			int valuableCharacterNum = 0;
			
			char *chrSeqChar = (char*)malloc((chromLengthVec[tmp]+1) * sizeof(char));

			string tmpChromFileNameStr = InputChromFolderStr + "/" + chrNameStrVec[tmp];
			FILE *fp_chr = fopen(tmpChromFileNameStr.c_str(), "r");

			char firstLineCharArray[500];
			char ch;
			fgets(firstLineCharArray, sizeof(firstLineCharArray), fp_chr);
			
			while((ch = fgetc(fp_chr)) != EOF)
			{
				if((ch == 'A')||(ch == 'a')) 
				{ 
					chrSeqChar[chromBaseNum] = 'A'; chromBaseNum ++; 
					wholeGenomeSeqChar[wholeGenomeChromBaseNum] = 'A'; wholeGenomeChromBaseNum ++;
					characterNum ++; valuableCharacterNum++;
				}
				else if((ch == 'C')||(ch == 'c')) 
				{ 
					chrSeqChar[chromBaseNum] = 'C';
					chromBaseNum++; 
					wholeGenomeSeqChar[wholeGenomeChromBaseNum] = 'C'; wholeGenomeChromBaseNum ++;
					characterNum ++; valuableCharacterNum++;
				}
				else if((ch == 'G')||(ch == 'g')) 
				{ 
					chrSeqChar[chromBaseNum] = 'G';
					chromBaseNum++; 
					wholeGenomeSeqChar[wholeGenomeChromBaseNum] = 'G'; wholeGenomeChromBaseNum ++;
					characterNum ++; valuableCharacterNum++;
				}
				else if((ch == 'T')||(ch == 't')) 
				{ 
					chrSeqChar[chromBaseNum] = 'T';
					chromBaseNum++; 
					wholeGenomeSeqChar[wholeGenomeChromBaseNum] = 'T'; wholeGenomeChromBaseNum ++;
					characterNum ++; valuableCharacterNum++;
				}
				else if((ch == 'N')||(ch == 'n')) 
				{ 
					chrSeqChar[chromBaseNum] = 'N';
					chromBaseNum++; 
					wholeGenomeSeqChar[wholeGenomeChromBaseNum] = 'N'; wholeGenomeChromBaseNum ++;
					characterNum ++; 
				}
				else if((ch == '\t')||(ch == '\n')) 
				{
					characterNum ++; 
					//continue;
				}
				else 
				{
					cout << "illegal input is " << ch << endl; characterNum ++; 
					break;
				}
			}
			chrSeqChar[chromBaseNum] = 'X';
			wholeGenomeSeqChar[wholeGenomeChromBaseNum] = 'X'; wholeGenomeChromBaseNum ++;
			cout << endl << "chromName: " << chrNameStrVec[tmp] << " chromBaseNum: " << chromBaseNum << endl;
			//cout << "characterNum: " << characterNum << endl
			//	<< "valuableCharNum: " << valuableCharacterNum << endl 
			//	<< "Num: " << chromLengthVec[tmp] << endl;
			cout << "chrSeqChar[last-1]:" << chrSeqChar[(chromBaseNum-1)] << " chrSeqChar[last]:" 
				<< chrSeqChar[chromBaseNum] << endl;
			fclose(fp_chr);
			free(chrSeqChar);
		}
		//wholeGenomeSeqChar[wholeGenomeChromBaseNum] = 'X';
		chromFileOfs.write((const char*) wholeGenomeSeqChar, wholeGenomeChromBaseNum * sizeof(char));

		cout << "wholeGenomeSeqChar[last-2]: " << wholeGenomeSeqChar[wholeGenomeChromBaseNum-3] << endl;
		cout << "wholeGenomeSeqChar[last-1]: " << wholeGenomeSeqChar[wholeGenomeChromBaseNum-2] << endl;
		cout << "wholeGenomeSeqChar[last]: " << wholeGenomeSeqChar[wholeGenomeChromBaseNum-1] << endl;
		cout << "wholeGenomeChromBaseNum: " << wholeGenomeChromBaseNum << endl;
		free(wholeGenomeSeqChar);
	}

	void outputIndexParameter(ofstream& parameterFileOfs)
	{
		int chromosomeNum = chrNameStrVec.size();
		parameterFileOfs << "chromNum:" << endl << chrNameStrVec.size() << endl << "chromName:" << endl;
		for(int tmp = 0; tmp < chromosomeNum; tmp++)
		{
			parameterFileOfs << "\"" << chrNameStrVec[tmp].substr(0, chrNameStrVec[tmp].length()-3) << "\"" << ","; 
		}
		parameterFileOfs << endl << "chromEndPosInGenome:" << endl;
		for(int tmp = 0; tmp < chromosomeNum; tmp++)
		{
			parameterFileOfs << chrEndPosInGenomeVec[tmp] << ",";
		}
		parameterFileOfs << endl << "2ndLevelIndexSize:" << endl << secondLevelIndexSize << endl
		<< "chrom2ndLevelIndexNum:" << endl;
		for(int tmp = 0; tmp < chromosomeNum; tmp++)
		{
			unsigned int tmpChrPartsNum = chromLengthVec[tmp]/secondLevelIndexSize;
			if(tmpChrPartsNum * secondLevelIndexSize == chromLengthVec[tmp])
			{}
			else if( (tmpChrPartsNum * secondLevelIndexSize < chromLengthVec[tmp]) 
				&& ( (tmpChrPartsNum+1) * secondLevelIndexSize > chromLengthVec[tmp] ))
			{
				tmpChrPartsNum ++;
			}
			else
			{
				cout << "error in outputIndexParameter from buildIndexParameter.h" << endl;
				cout << "chromLength " << tmp+1 << ": " << chromLengthVec[tmp] << endl;
				cout << "chrPartsNum: " << tmpChrPartsNum << endl; 
			}

			parameterFileOfs << tmpChrPartsNum << ",";
		}
		parameterFileOfs << endl << "chromosome length: " << endl;

		for(int tmp = 0; tmp < chromosomeNum; tmp++)
		{
			parameterFileOfs << chromLengthVec[tmp] << ",";
		}		
		parameterFileOfs << endl;
	}

	// Getters
	string const GetInputChromFolder()
	{
		return InputChromFolderStr;
	}

	string const GetOutputIndexFolder()
	{
		return OutputIndexFolderStr;
	}

	unsigned int GetMax()
	{
		return chrEndPosInGenomeVec[chrNameStrVec.size()-1] + 1 + 1;
	}

	// Setters
	void AddChromName(string newName)
	{
		chrNameStrSet.insert(newName);
		chrNameStrVec.push_back(newName);
	}

	void AddChromLength(int length)
	{
		chromLengthVec.push_back(length);
	}
};
