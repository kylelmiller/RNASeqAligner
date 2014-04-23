#include <stdlib.h>
#include <string>
#include <string.h>
#include <map>

using namespace std;

//#define MAX 249250622  //string length. when getting r[] array, end with the maximum number Range. the chrome number is MAX-1
#define Range 6 //the max range of the element in the string. a paremeter of sorting in function da
#define NULL_NUM 4294967290 //shoud be -1 in the previous algorithm, changed to NULL_NUM because of unsigned int.

void build_up_down(unsigned int *lcptab, unsigned int *up, unsigned int *down, unsigned int n)
{
    unsigned int lastIndex = NULL_NUM;	
	stack<unsigned int> up_down;
	up_down.push(0);
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		while (lcptab[i] < lcptab[up_down.top()])
		{	lastIndex = up_down.top();
			up_down.pop();
			if((lcptab[i] <= lcptab[up_down.top()]) && (lcptab[up_down.top()] != lcptab[lastIndex]))
				down[up_down.top()] = lastIndex;
		}
		// now lcptab[i] >= lcptab[up_down.top()] holds
		if(lastIndex != NULL_NUM)
		{
			up[i] = lastIndex;
			lastIndex = NULL_NUM;
		}
		up_down.push(i);	
	}
	return;
}

void build_next(unsigned int *lcptab, unsigned int *next, unsigned int n)
{
	stack<unsigned int> nextIndex;
	unsigned int lastIndex;
	unsigned int j;
	nextIndex.push(0);
	for(j = 0; j < n; j++)
	{
		while(lcptab[j] < lcptab[nextIndex.top()])
			nextIndex.pop();
		if(lcptab[j] == lcptab[nextIndex.top()])
		{
			lastIndex = nextIndex.top();
			nextIndex.pop();
			next[lastIndex] = j;
		}
		nextIndex.push(j);
	}	
	return;
}

void build_lcp(unsigned int *r, unsigned int *sa, unsigned int *lcp, unsigned int *rank, unsigned int n)
{

	unsigned int i, j; 
	int k=0;
	for (i = 0; i < n; i++) rank[sa[i]] = i;
		//cout << "stop10_new" << endl;
	for (i = 0; i < n; lcp[rank[i++]] = k) 
	for (k?k--:0, (rank[i] == 0)?(j=0):(j=sa[rank[i]-1]); r[i+k] == r[j+k]; k++);
	lcp[0] = 0;
		//cout << "stop11" << endl; 
	return;
}

unsigned int cmp(unsigned int *r,unsigned int a,unsigned int b,unsigned int l)
{return r[a]==r[b]&&r[a+l]==r[b+l];}  

void da(unsigned int *r,unsigned int *sa,unsigned int n,unsigned int m)
{   
	//cout << "stop3" << endl;
	unsigned int maxn = n;
	//unsigned int wa[maxn],wb[maxn],wv[maxn],ws[maxn];
    unsigned int *wa = (unsigned int*)malloc(maxn * sizeof(unsigned int));
    unsigned int *wb = (unsigned int*)malloc(maxn * sizeof(unsigned int));
    unsigned int *wv = (unsigned int*)malloc(maxn * sizeof(unsigned int));
    unsigned int *ws = (unsigned int*)malloc(maxn * sizeof(unsigned int));
    unsigned int i,j,p,*x=wa,*y=wb,*t;
	//cout << "stop4" << endl;    
    for(i=0;i<m;i++) ws[i]=0;
	//cout << "stop5" << endl;
    for(i=0;i<n;i++) ws[x[i]=r[i]]++; 
	//cout << "stop6" << endl;
	for(i=1;i<m;i++) ws[i]+=ws[i-1];
	//cout << "stop7" << endl;
    for(i=n-1;i>=0;i--) 
    {
    	sa[--ws[x[i]]]=i; 
		if(i == 0)
			break;
	}
	//cout << "stop8" << endl;
    for(j=1,p=1;p<n;j*=2,m=p)
    {
        for(p=0,i=n-j;i<n;i++) y[p++]=i;  
        for(i=0;i<n;i++) if(sa[i]>=j) y[p++]=sa[i]-j; 
        for(i=0;i<n;i++) wv[i]=x[y[i]];  
        for(i=0;i<m;i++) ws[i]=0;
        for(i=0;i<n;i++) ws[wv[i]]++;
        for(i=1;i<m;i++) ws[i]+=ws[i-1];
        for(i=n-1;i>=0;i--) 
        {
        	sa[--ws[wv[i]]]=y[i];  
			if(i == 0)
				break;
		}
		for(t=x,x=y,y=t,p=1,x[sa[0]]=0,i=1;i<n;i++)
        x[sa[i]]=cmp(y,sa[i-1],sa[i],j)?p-1:p++; 
    }
    return;
}

int INDEX_KMER_LENGTH = 14;

int baseChar2intArray[26] = {0, 100, 1, 100, 100, 100, 2,
			100, 100, 100, 100, 100, 100, 100,
			100, 100, 100, 100, 100, 3, 
			100, 100, 100, 100, 100, 100};

int baseCharCount2intArray[14][26] = {0};


class SecondLevelIndex_Info
{
public:
	string chromString;
	unsigned int genomeLength;

	int chromNum;
	vector<string> chrNameStr; // size = chromNum
	vector<int> chromLength; // size = chromNum
	vector<string> chromStr;
	vector<unsigned int> chrEndPosInGenome;

	map<string, int> chrNameMap;
	//map<string, int>::iterator chrNameMapIter;

	int secondLevelIndexNormalSize;// = 3000000;
	vector<int> secondLevelIndexPartsNum;
	int secondLevelIndexPartsNumSum;

	unsigned int null_num; // 2654911540 for mm9_noRandom genome
	unsigned int indexSize; //2654911539  //sequence length + 1, the length of sa-lcp-down-next 
	//omp_lock_t lock;

	SecondLevelIndex_Info()
	{}

	SecondLevelIndex_Info(ifstream& inputIndexInfoFile)
	{
		for(int tmp1 = 0; tmp1 < INDEX_KMER_LENGTH; tmp1++)
		{
			for(int tmp2 = 0; tmp2 < 26/*# of letters alphabet*/; tmp2++)
			{
				baseCharCount2intArray[tmp1][tmp2] = 100;
			}
		}
		int tmpBaseCount = 1;
		for(int tmp3 = 0; tmp3 < INDEX_KMER_LENGTH; tmp3++)
		{
			baseCharCount2intArray[tmp3][0] = 0*tmpBaseCount;
			baseCharCount2intArray[tmp3][2] = 1*tmpBaseCount;
			baseCharCount2intArray[tmp3][6] = 2*tmpBaseCount;
			baseCharCount2intArray[tmp3][19] = 3*tmpBaseCount;
			tmpBaseCount = 4*tmpBaseCount;
		}


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

		indexSize = chrEndPosInGenome[chrEndPosInGenome.size()-1] + 2;
		//cout << "MAX: " << indexSize << endl;
		null_num = indexSize + 1;
		//cout << "NULL_NUM: " << null_num << endl;
		this->buildChrNameMap();
	}

	void Genome2ChromString(const string& genomeString)
	{
		string tmpChromString = genomeString.substr(0, (chrEndPosInGenome)[0]+1);
		chromStr.push_back(tmpChromString);
		cout << "chrom chr.1: " << tmpChromString.length() << endl; 
		for(int tmp = 1; tmp < chromNum; tmp++)
		{	
			tmpChromString = genomeString.substr((chrEndPosInGenome)[tmp-1]+2, 
				(chrEndPosInGenome)[tmp]-(chrEndPosInGenome)[tmp-1]-1);	
			chromStr.push_back(tmpChromString);
			cout << "chrom " << chrNameStr[tmp] << ": " << tmpChromString.length() << endl; 
		}
	}

	void creatAllSubFolder(string output2ndLevelIndexPrefix)
	{
		string all2ndLevelIndexFolder = output2ndLevelIndexPrefix;								
		string mkAll2ndLevelIndexFolderCmd = "mkdir " + all2ndLevelIndexFolder;
		system(mkAll2ndLevelIndexFolderCmd.c_str());
		cout << endl << "creat 2ndlevel index folder: " << all2ndLevelIndexFolder << endl << endl;

		for(int tmpChromNum = 0; tmpChromNum < chromNum; tmpChromNum++)
		{
			string tmp2ndLevelIndexChromFolder = output2ndLevelIndexPrefix 
				+ "/" + chrNameStr[tmpChromNum];
			string mkTmp2ndLevelIndexChromFolderCmd = "mkdir " + tmp2ndLevelIndexChromFolder;
			system(mkTmp2ndLevelIndexChromFolderCmd.c_str());
			cout << endl << "... creat 2ndlevel index subChrfolder " << chrNameStr[tmpChromNum]  << ": " << tmp2ndLevelIndexChromFolder << endl;
			
			int tmp2ndLevelIndexPartsNum = secondLevelIndexPartsNum[tmpChromNum];
			for(int tmpPartNum = 0; 
				tmpPartNum < tmp2ndLevelIndexPartsNum; tmpPartNum++)
			{
				char tmpFileNumChar[4];
				sprintf(tmpFileNumChar, "%d", tmpPartNum+1);
				string tmpFileNumStr = tmpFileNumChar;
				string tmp2ndLevelIndexChromPartFolder = tmp2ndLevelIndexChromFolder
					+ "/" + tmpFileNumStr;
				string mkTmp2ndLevelIndexChromPartFolderCmd = "mkdir " 
					+ tmp2ndLevelIndexChromPartFolder;
				system(mkTmp2ndLevelIndexChromPartFolderCmd.c_str());
				cout << endl << "...... creat subChrPartfolder: " << tmp2ndLevelIndexChromPartFolder << endl;
			}			
			cout << endl;
		}

	}

	void generate2ndLevelIndexChrom(string output2ndLevelIndexPrefix)
	{
		for(int tmpChromNum = 0; tmpChromNum < chromNum; tmpChromNum++)
		{
			cout << endl << endl << "generating chrom file, name: " << chrNameStr[tmpChromNum] << endl;
			int tmpChromLength = chromStr[tmpChromNum].length();
			int tmp2ndLevelIndexPartsNum = secondLevelIndexPartsNum[tmpChromNum];
			string tmp2ndLevelIndexPrefixChrom = output2ndLevelIndexPrefix 
				+ "/" + chrNameStr[tmpChromNum] + "/";  
			for(int tmpPartNum = 0; 
				tmpPartNum < tmp2ndLevelIndexPartsNum-1; tmpPartNum++)
			{
				char tmpFileNumChar[4];
				sprintf(tmpFileNumChar, "%d", tmpPartNum+1);
				string tmpFileNumStr = tmpFileNumChar;
				string tmpChromFileStr = tmp2ndLevelIndexPrefixChrom
					+ tmpFileNumStr + "/chrom";

				string tmp2ndLevelIndexChromStr 
					= chromStr[tmpChromNum].substr(tmpPartNum*secondLevelIndexNormalSize, secondLevelIndexNormalSize);
				tmp2ndLevelIndexChromStr += "X";
				ofstream tmp2ndLevelIndexChromFile_ofs(tmpChromFileStr.c_str());
				tmp2ndLevelIndexChromFile_ofs << tmp2ndLevelIndexChromStr << endl;
				tmp2ndLevelIndexChromFile_ofs.close();
			}	
			char fileMaxNumChar[4];
			sprintf(fileMaxNumChar, "%d", tmp2ndLevelIndexPartsNum);
			string fileMaxNumStr = fileMaxNumChar;

			string tmpChromFileStr = tmp2ndLevelIndexPrefixChrom
				+ fileMaxNumStr + "/chrom";
			string tmp2ndLevelIndexChromStr 
				= chromStr[tmpChromNum].substr((tmp2ndLevelIndexPartsNum-1)*secondLevelIndexNormalSize);
			tmp2ndLevelIndexChromStr += "X";
			ofstream tmp2ndLevelIndexChromFile_ofs(tmpChromFileStr.c_str());
			tmp2ndLevelIndexChromFile_ofs << tmp2ndLevelIndexChromStr << endl;
			tmp2ndLevelIndexChromFile_ofs.close();

		}

	}

	void generate2ndLevelIndexOriginalSize(string output2ndLevelIndexPrefix)
	{
		for(int tmpChromNum = 0; tmpChromNum < chromNum; tmpChromNum ++)
		{
			string tmp2ndLevelIndexPrefixChr = output2ndLevelIndexPrefix + "/" 
				+ chrNameStr[tmpChromNum] + "/"; 
			int tmp2ndLevelIndexPartsNum = secondLevelIndexPartsNum[tmpChromNum];
			
			cout << endl << "... start to generate original size 2ndlevel chr index: " << tmp2ndLevelIndexPrefixChr << endl;
			for(int tmpChrPartNO = 0; tmpChrPartNO < tmp2ndLevelIndexPartsNum; tmpChrPartNO++)
			{
				char tmpFileNumChar[4];
				sprintf(tmpFileNumChar, "%d", tmpChrPartNO+1);
				string tmpFileNumStr = tmpFileNumChar;
				string tmp2ndLevelIndexPrefixChrPart = tmp2ndLevelIndexPrefixChr + tmpFileNumStr + "/";			
			
				cout << endl << "...... start to generate original size 2ndlevel chrPart index: " 
					<< tmp2ndLevelIndexPrefixChrPart << endl;
			
				string tmp2ndLevelIndexPrefixChrPart_chrom = tmp2ndLevelIndexPrefixChrPart + "chrom";
				ifstream tmp2ndLevelIndexPrefixChrPart_chrom_ifs(tmp2ndLevelIndexPrefixChrPart_chrom.c_str(),ios::binary);

				string tmp2ndLevelIndexPrefixChrPart_SA = tmp2ndLevelIndexPrefixChrPart + "SA";	
				ofstream tmp2ndLevelIndexPrefixChrPart_SA_ofs(tmp2ndLevelIndexPrefixChrPart_SA.c_str(), ios::binary);
				
				string tmp2ndLevelIndexPrefixChrPart_lcp = tmp2ndLevelIndexPrefixChrPart + "lcp";
				ofstream tmp2ndLevelIndexPrefixChrPart_lcp_ofs(tmp2ndLevelIndexPrefixChrPart_lcp.c_str(), ios::binary);

				string tmp2ndLevelIndexPrefixChrPart_up = tmp2ndLevelIndexPrefixChrPart + "up";
				ofstream tmp2ndLevelIndexPrefixChrPart_up_ofs(tmp2ndLevelIndexPrefixChrPart_up.c_str(), ios::binary);

				string tmp2ndLevelIndexPrefixChrPart_down = tmp2ndLevelIndexPrefixChrPart + "down";
				ofstream tmp2ndLevelIndexPrefixChrPart_down_ofs(tmp2ndLevelIndexPrefixChrPart_down.c_str(), ios::binary);

				string tmp2ndLevelIndexPrefixChrPart_next = tmp2ndLevelIndexPrefixChrPart + "next";
				ofstream tmp2ndLevelIndexPrefixChrPart_next_ofs(tmp2ndLevelIndexPrefixChrPart_next.c_str(), ios::binary);

				unsigned int MAX = secondLevelIndexNormalSize + 1;
				if(tmpChrPartNO == tmp2ndLevelIndexPartsNum - 1)
				{
					MAX = chromStr[tmpChromNum].length() - (tmp2ndLevelIndexPartsNum-1)*secondLevelIndexNormalSize + 1;
				}

				char *chrom;
				chrom = (char*)malloc(MAX * sizeof(char));				
				tmp2ndLevelIndexPrefixChrPart_chrom_ifs.read((char*)chrom, MAX * sizeof(char));
				string tmpChromString = chrom;
				tmpChromString = tmpChromString.substr(0,MAX);
				if(tmpChromString.find_first_of("ATGC") == string::npos)
				{
					cout << "... no ACGT ..." << endl;
					continue;
				}
				else
				{
					cout << "first location of ACGT is: " 
						<< tmpChromString.find_first_of("ATGC")
						<< endl;
					cout << "stringLength: " << tmpChromString.length() << endl;
					cout << "MAX: " << MAX << endl;
					cout << "char: " << tmpChromString.at(tmpChromString.find_first_of("ATGC")) << endl;
				}
				//chrom[secondLevelIndexNormalSize] = 'X';
				unsigned int *r = (unsigned int*)malloc(MAX * sizeof(unsigned int));
				for(int tmpBaseNum = 0; tmpBaseNum < MAX; tmpBaseNum++)
				{
					if(chrom[tmpBaseNum] == 'A') 
						{r[tmpBaseNum] = 1;}
					else if(chrom[tmpBaseNum] == 'C') 
						{r[tmpBaseNum] = 2;}
					else if(chrom[tmpBaseNum] == 'G') 
						{r[tmpBaseNum] = 3;}
					else if(chrom[tmpBaseNum] == 'T') 
						{r[tmpBaseNum] = 4;}
					else if(chrom[tmpBaseNum] == 'N') 
						{r[tmpBaseNum] = 5;}
					else if(chrom[tmpBaseNum] == 'X') 
						{r[tmpBaseNum] = 6;}
					else 
					{
						printf("\n illegal input is '%c'",chrom[tmpBaseNum]); 
						break;
					}
				}

				r[MAX-1] = Range;
				chrom[MAX-1] = 'X';

				unsigned int *sa = (unsigned int*)malloc(MAX * sizeof(unsigned int));
				da(r,sa,MAX,Range+1);

			    unsigned int *rank = (unsigned int*)malloc(MAX * sizeof(unsigned int));
			    unsigned int *lcp = (unsigned int*)malloc(MAX * sizeof(unsigned int));
			    unsigned int *up = (unsigned int*)malloc(MAX * sizeof(unsigned int));
			    unsigned int *down = (unsigned int*)malloc(MAX * sizeof(unsigned int));
			    unsigned int *next = (unsigned int*)malloc(MAX * sizeof(unsigned int));
			    build_lcp(r, sa, lcp, rank, MAX); 
			    build_up_down(lcp, up, down, MAX); 
			    build_next(lcp, next, MAX);

				tmp2ndLevelIndexPrefixChrPart_SA_ofs.write((const char*) sa, MAX * sizeof(unsigned int));
				tmp2ndLevelIndexPrefixChrPart_lcp_ofs.write((const char*) lcp, MAX * sizeof(unsigned int));
				tmp2ndLevelIndexPrefixChrPart_up_ofs.write((const char*) up, MAX * sizeof(unsigned int));
				tmp2ndLevelIndexPrefixChrPart_down_ofs.write((const char*) down, MAX * sizeof(unsigned int));
				tmp2ndLevelIndexPrefixChrPart_next_ofs.write((const char*) next, MAX * sizeof(unsigned int));
				free(r); free(rank);
				free(chrom); free(lcp); free(sa); free(up); free(down); free(next);
				tmp2ndLevelIndexPrefixChrPart_SA_ofs.close();
				tmp2ndLevelIndexPrefixChrPart_lcp_ofs.close();
				tmp2ndLevelIndexPrefixChrPart_up_ofs.close();
				tmp2ndLevelIndexPrefixChrPart_down_ofs.close();
				tmp2ndLevelIndexPrefixChrPart_next_ofs.close();				
				tmp2ndLevelIndexPrefixChrPart_chrom_ifs.close();
				//chrom_bit_file_ofs.write((const char*) chrom, MAX * sizeof(char));
			}
		}		
	}

	void generate2ndLevelIndexCompressedSize(string output2ndLevelIndexPrefix)
	{
		for(int tmpChromNum = 0; tmpChromNum < chromNum; tmpChromNum ++)
		{
			string tmp2ndLevelIndexPrefixChr = output2ndLevelIndexPrefix + "/" 
				+ chrNameStr[tmpChromNum] + "/"; 
			int tmp2ndLevelIndexPartsNum = secondLevelIndexPartsNum[tmpChromNum];

			cout << endl << "... start to generate compressed size chr index: " << tmp2ndLevelIndexPrefixChr << endl;
			for(int tmpChrPartNO = 0; tmpChrPartNO < tmp2ndLevelIndexPartsNum; tmpChrPartNO++)
			{
				char tmpFileNumChar[4];
				sprintf(tmpFileNumChar, "%d", tmpChrPartNO+1);
				string tmpFileNumStr = tmpFileNumChar;
				string tmp2ndLevelIndexPrefixChrPart = tmp2ndLevelIndexPrefixChr + tmpFileNumStr + "/";			
				
				cout << endl << "...... start to generate compressed size chr index: " 
					<< tmp2ndLevelIndexPrefixChrPart << endl;

				string tmp2ndLevelIndexPrefixChrPart_chrom = tmp2ndLevelIndexPrefixChrPart + "chrom";
				ifstream tmp2ndLevelIndexPrefixChrPart_chrom_ifs(tmp2ndLevelIndexPrefixChrPart_chrom.c_str(),ios::binary);
				//cout << "finish generating chrom" << endl;	
				string tmp2ndLevelIndexPrefixChrPart_lcp = tmp2ndLevelIndexPrefixChrPart + "lcp";
				ifstream lcp_file_ifs(tmp2ndLevelIndexPrefixChrPart_lcp.c_str(), ios::binary);
				//cout << "finish generating lcp" << endl;	
				string tmp2ndLevelIndexPrefixChrPart_up = tmp2ndLevelIndexPrefixChrPart + "up";
				ifstream up_file_ifs(tmp2ndLevelIndexPrefixChrPart_up.c_str(), ios::binary);
				//cout << "finish generating up" << endl;	
				string tmp2ndLevelIndexPrefixChrPart_down = tmp2ndLevelIndexPrefixChrPart + "down";
				ifstream down_file_ifs(tmp2ndLevelIndexPrefixChrPart_down.c_str(), ios::binary);
				//cout << "finish generating down" << endl;	
				string tmp2ndLevelIndexPrefixChrPart_next = tmp2ndLevelIndexPrefixChrPart + "next";
				ifstream next_file_ifs(tmp2ndLevelIndexPrefixChrPart_next.c_str(), ios::binary);
				//cout << "finish generating next" << endl;	
				string childTab_file = tmp2ndLevelIndexPrefixChrPart; childTab_file.append("childTab"); 
				ofstream childTab_file_ofs(childTab_file.c_str(),ios::binary);
				
				string detChild_file = tmp2ndLevelIndexPrefixChrPart; detChild_file.append("detChild"); 
				ofstream detChild_file_ofs(detChild_file.c_str(), ios::binary);

				string lcpCompress_file = tmp2ndLevelIndexPrefixChrPart; lcpCompress_file.append("_lcpCompress"); 
				ofstream lcpCompress_file_ofs(lcpCompress_file.c_str(),ios::binary);

				cout << "finish creating index files ... " << endl;
				unsigned int indexSize = secondLevelIndexNormalSize + 1;
				if(tmpChrPartNO == tmp2ndLevelIndexPartsNum - 1)
				{
					indexSize = chromStr[tmpChromNum].length() - (tmp2ndLevelIndexPartsNum-1)*secondLevelIndexNormalSize + 1;
				}
				//cout << "indexSize: " << indexSize << endl;
				char *chrom;
				chrom = (char*)malloc(indexSize * sizeof(char));
				//cout << "after mallocing chrom ..." << endl;				
				tmp2ndLevelIndexPrefixChrPart_chrom_ifs.read((char*)chrom, indexSize * sizeof(char));
				string tmpChromString = chrom;
				tmpChromString = tmpChromString.substr(0,indexSize);
				//cout << "indexSize: " << indexSize << endl;
				// ACGT_pos = tmpChromString.find_first_of("ATGC")
				if(tmpChromString.find_first_of("ATGC") == string::npos)
				{
					cout << "... no ACGT ..." << endl;
					continue;
				}
				else
				{
					cout << "first location of ACGT is: " 
						<< tmpChromString.find_first_of("ATGC")
						<< endl;
					//cout << "stringLength: " << tmpChromString.length() << endl;
					cout << "indexSize: " << indexSize << endl;
					//cout << "char: " << tmpChromString.at(tmpChromString.find_first_of("ATGC")) << endl;
				}

				//cout << "start to load lcp" << endl;
			    unsigned int *lcp;
			    lcp = (unsigned int*)malloc(indexSize * sizeof(unsigned int));
			    lcp_file_ifs.read((char*)lcp, (indexSize * sizeof(unsigned int)));

				//cout << "start to load up" << endl;
			    unsigned int *up;
			    up = (unsigned int*)malloc(indexSize * sizeof(unsigned int));
			    up_file_ifs.read((char*)up, (indexSize * sizeof(unsigned int)));

				//cout << "start to load down" << endl;
			    unsigned int *down;
			    down = (unsigned int*)malloc(indexSize * sizeof(unsigned int));
			    down_file_ifs.read((char*)down, (indexSize * sizeof(unsigned int)));

				//cout << "start to load next" << endl;
			    unsigned int *next;
			    next = (unsigned int*)malloc(indexSize * sizeof(unsigned int));
			    next_file_ifs.read((char*)next, (indexSize * sizeof(unsigned int)));

			    //cout << "All index files loaded ..." << endl;

			    unsigned int *childTab;
			    childTab = (unsigned int*)malloc(indexSize * sizeof(unsigned int));

			   	BYTE *verifyChild;
				verifyChild = (BYTE*)malloc(indexSize * sizeof(BYTE));

				BYTE *lcpCompress = (BYTE*)malloc(indexSize * sizeof(BYTE));

				BuildIndex_Info* tmpBuildIndexInfo = new BuildIndex_Info();

				cout << "start to compress child ..." << endl;
				bool compressChildBool = tmpBuildIndexInfo->compressUpDownNext2ChildtabVerifyChild_bool(
					up, down, next, childTab, verifyChild, indexSize);

				if(!compressChildBool)
				{
					cout << "...... CompressChild error ! " << endl;
				}
				//cout << "finish compressing child" << endl;
				cout << "start to compress Lcp" << endl;				
				tmpBuildIndexInfo->compressLcp2Lcpcompress(lcp, 
					lcpCompress, indexSize);
				cout << "finish compressing child & lcp" << endl;
				lcpCompress_file_ofs.write((const char*) lcpCompress, indexSize * sizeof(BYTE));	
				childTab_file_ofs.write((const char*) childTab, indexSize * sizeof(unsigned int));
				detChild_file_ofs.write((const char*) verifyChild, indexSize * sizeof(BYTE));	
			
				lcp_file_ifs.close();
				up_file_ifs.close();
				down_file_ifs.close();
				next_file_ifs.close();
				lcpCompress_file_ofs.close();
				childTab_file_ofs.close();
				detChild_file_ofs.close();
			
				free(chrom);
				free(lcpCompress);
				free(lcp);
				free(up); free(down); free(next); free(childTab); free(verifyChild);
			}	
		}	
	}

	void buildChrNameMap()
	{
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			chrNameMap.insert(pair <string, int> (chrNameStr[tmp], tmp));
		}

	}

	int getSecondLevelIndexFromChrAndPos(int chrNameInt, int chrMapPos)
	{
		int tmpTimes = chrMapPos/secondLevelIndexNormalSize;
		int partsTimeBase = 0;
		for(int tmp = 0; tmp < chrNameInt; tmp++)
		{
			partsTimeBase += secondLevelIndexPartsNum[tmp];
		}
		return	(partsTimeBase + tmpTimes + 1); 
	}

	int getChrPosFromSecondLevelIndexPos(int chrNameInt, int secondLevelIndexNum, int secondLevelIndexPos)
	{
		int partsTimeBase = 0;
		for(int tmp = 0; tmp < chrNameInt; tmp++)
		{
			partsTimeBase += secondLevelIndexPartsNum[tmp];
		}

		int tmpSecondLevelIndexNO = secondLevelIndexNum - partsTimeBase;
		return ( (tmpSecondLevelIndexNO-1) * secondLevelIndexNormalSize + secondLevelIndexPos);	
	}

	void getChrLocation(unsigned int locationInWholeGenome, unsigned int *chr_name_int, unsigned int *chr_local_location)
	{
		if(locationInWholeGenome <= chrEndPosInGenome[0])
		{
			(*chr_name_int) = 0;
			*chr_local_location = locationInWholeGenome;
		}
		else
		{
			for(int tmp = 1; tmp < chrEndPosInGenome.size(); tmp++)
			{
				if( (locationInWholeGenome >= chrEndPosInGenome[tmp-1] + 2) 
					&& (locationInWholeGenome <= chrEndPosInGenome[tmp]) )
				{
					*chr_name_int = tmp;
					*chr_local_location = locationInWholeGenome - chrEndPosInGenome[tmp-1] - 2;
				}
				else
				{
					continue;
				}
			}
		}
	}


	unsigned int getWholeGenomeLocation(unsigned int chr_name_int, unsigned int locationInWholeGenome)
	{
		//cout << "in function chr_name_int: " << chr_name_int << endl;
		unsigned int chr_local_location;
		if(chr_name_int == 0)
		{
			chr_local_location = locationInWholeGenome;
		}
		else if(chr_name_int < chromNum)
		{
			//cout << "< chromNum chr_name_int: " << chr_name_int << endl;
			chr_local_location = locationInWholeGenome 
				+ chrEndPosInGenome[chr_name_int-1] + 2; 
			//cout << "chr_local_location: " << chr_local_location << endl;
		}
		else
		{
			cout << "chr_name_int error: " << chr_name_int << endl;
		}
		return chr_local_location;
	}

	int getChr(unsigned int locationInWholeGenome)
	{
		int chrInt;
		if(locationInWholeGenome <= chrEndPosInGenome[0])
		{
			chrInt = 0;
		}
		else
		{
			for(int tmp = 1; tmp < chromNum; tmp++)
			{
				if( (locationInWholeGenome >= chrEndPosInGenome[tmp-1] + 2) 
					&& (locationInWholeGenome <= chrEndPosInGenome[tmp]) )
				{
					chrInt = tmp;
					break;
				}
				else
				{
					continue;
				}				
			}
		}
		return chrInt;
	}


	int convertStringToInt(const string& chrName)
	{
		
		//omp_init_lock(&lock);
		map<string, int>::iterator chrNameMapIter;
		//omp_set_lock(&lock);
		int chrNameInt = 1000;
		//cout << "chrName = " << chrName << endl;
		chrNameMapIter = chrNameMap.find(chrName);
		if(chrNameMapIter != chrNameMap.end())
		{
			chrNameInt = chrNameMapIter->second;
		}
		else
		{
			cout << "...... chrom name error! ...... " << endl;
		}
		//cout << "chrNameInt: " << chrNameInt;
		//omp_unset_lock(&lock);
		//omp_destroy_lock(&lock);
		return chrNameInt;
	}

};