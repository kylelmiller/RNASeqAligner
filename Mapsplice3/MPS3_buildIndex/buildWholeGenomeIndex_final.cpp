/**********************************
 * This file:
 *
 * 1 - Builds a suffix array on a given chromosome
 * 	The suffix array building code uses the DC3
 * 	algorithm which constructs the suffix array
 * 	in linear time. The algorithm is explained
 * 	in the paper "Linear Work Suffix Array Construction".
 *
 * 2 - Builds a longest common prefix array for
 * 	our suffix array. This allows for fastest query times
 * 	by allowing the suffix array to be used like a suffix
 * 	tree.
 *
 * 3 - Builds up, down and nextIndex tables.
 * 	The reason for this is outlined in the paper "Optimal
 * 	Exact String Matching Based on Suffix Arrays". These
 * 	tables,	along with the lcp array, allows for suffix array
 * 	queries to be done in O(k) time, where k is the query string.
 *********************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <sys/types.h>    
#include <dirent.h>    
#include <stdio.h>    
#include <errno.h>
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
#include <set>
#include "buildIndexParameter.h"

// Constants
#define Range 6
#define LONGLCP 255
#define NULL_NUM 4294967290 
#define PREINDEX_STRINGLENGTH 14
#define CANDALILOC 100
#define SEGMENTNUM 20
#define minValSegLength 20
#define PreIndexSize 268435456

#define EXPECTED_ARGUMENT_COUNT 3

// End Constants

using namespace std;
typedef unsigned char BYTE;
int read_num = 0; // to calculate the read num;

int baseChar2intArray[26] =
		{ 0, 100, 1, 100, 100, 100, 2, 100, 100, 100, 100, 100, 100, 100, 100,
				100, 100, 100, 100, 3, 100, 100, 100, 100, 100, 100 };

// used to catch exceptions for custom processing
void handler(int sig) {
	void *array[10];
	size_t size;

	// get void's for all entries on the stack
	size = backtrace(array, 10);

	// print out all the frames to stderr
	fprintf(stderr, "Error: %s\n", strsignal(sig));
	backtrace_symbols_fd(array, size, STDERR_FILENO);
	raise(SIGABRT);
}

inline string int_to_str(int numerical) {
	char c[100];
	sprintf(c, "%d", numerical);
	string str(c);
	return str;
}

unsigned int getPreIndexNO(const string& readPreStr) {
	int preIndexStrSize = readPreStr.length();

	unsigned int preIndexNO = 0;

	int baseForCount = 1;

	for (int tmp = preIndexStrSize - 1; tmp >= 0; tmp--) {
		preIndexNO = preIndexNO
				+ baseChar2intArray[(readPreStr.at(tmp) - 'A')] * baseForCount;
		baseForCount = baseForCount * 4;
	}
	return preIndexNO;
}

unsigned int getChildUpValue(unsigned int *child, BYTE *verifyChild,
		unsigned int index) {
	return ((verifyChild[index - 1] == 1) * child[index - 1]);
}

unsigned int getChildDownValue(unsigned int *child, BYTE *verifyChild,
		unsigned int index) {
	if (verifyChild[index] == 4)
		return child[child[index] - 1];
	else if (verifyChild[index] == 2)
		return child[index];
	else
		return 0;
}

unsigned int getChildNextValue(unsigned int *child, BYTE *verifyChild,
		unsigned int index) {
	return ((verifyChild[index] > 2) * child[index]);
}

void getFirstInterval(char ch, unsigned int *interval_begin,
		unsigned int *interval_end, unsigned int* child, BYTE* verifyChild) {
	unsigned int child_next_tmp = ((verifyChild[0] > 2) * child[0]);
	switch (ch) {
	case 'C':
		*interval_begin = child_next_tmp;
		*interval_end = ((verifyChild[*interval_begin] > 2)
				* child[*interval_begin]) - 1;
		break;
	case 'A':
		*interval_begin = 0;
		*interval_end = child_next_tmp - 1;
		break;
	case 'G':
		*interval_begin = ((verifyChild[child_next_tmp] > 2)
				* child[child_next_tmp]);
		*interval_end = ((verifyChild[*interval_begin] > 2)
				* child[*interval_begin]) - 1;
		break;
	default: // case 'T'
		unsigned int child_next_tmp2 = ((verifyChild[child_next_tmp] > 2)
				* child[child_next_tmp]);
		*interval_begin = ((verifyChild[child_next_tmp2] > 2)
				* child[child_next_tmp2]);
		*interval_end = ((verifyChild[*interval_begin] > 2)
				* child[*interval_begin]) - 1;
		break;
	}
}

void getInterval(unsigned int start, unsigned int end, unsigned int position,
		char ch, unsigned int *interval_begin, unsigned int *interval_end,
		unsigned int* sa, unsigned int* child, char* chrom, BYTE* verifyChild) {
	unsigned int index_begin;
	unsigned int pos;
	*interval_end = 0;
	*interval_begin = 1;

	if (ch == 'C') {
		unsigned int child_up_value_tmp = ((verifyChild[end] == 1) * child[end]);
		if ((start < child_up_value_tmp) && (end >= child_up_value_tmp))
			index_begin = child_up_value_tmp; //getChildUpValue(child, verifyChild, end+1);
		else
			index_begin = getChildDownValue(child, verifyChild, start);

		if (chrom[sa[index_begin] + position] == 'C') {
			*interval_begin = index_begin;
			*interval_end =
					((verifyChild[index_begin] > 2) * child[index_begin]) - 1;
			if ((*interval_end < start) || (*interval_end > end))
				*interval_end = end;
		} else if (chrom[sa[start] + position] == 'C') {
			*interval_begin = start;
			*interval_end = index_begin - 1;
			if ((*interval_end < start) || (*interval_end > end))
				*interval_end = end;
		}
	} else if (ch == 'A') // ch == 'A'
			{
		unsigned int child_up_value_tmp = ((verifyChild[end] == 1) * child[end]);
		if ((start < child_up_value_tmp) && (end >= child_up_value_tmp))
			index_begin = child_up_value_tmp; //getChildUpValue(child, verifyChild, end+1);
		else
			index_begin = getChildDownValue(child, verifyChild, start);

		if (chrom[sa[start] + position] == 'A') {
			*interval_begin = start;
			*interval_end = index_begin - 1;
			if ((*interval_end < start) || (*interval_end > end))
				*interval_end = end;
		}
	} else // ch == 'G' or ch == 'T'
	{
		if (ch == 'G') {
			unsigned int child_up_value_tmp = ((verifyChild[end] == 1)
					* child[end]);
			if ((start < child_up_value_tmp) && (end >= child_up_value_tmp))
				index_begin = child_up_value_tmp; //getChildUpValue(child, verifyChild, end+1);
			else
				index_begin = getChildDownValue(child, verifyChild, start);

			pos = ((verifyChild[index_begin] > 2) * child[index_begin]);
			if (chrom[sa[pos] + position] == 'G') {
				*interval_begin = pos;
				*interval_end = ((verifyChild[pos] > 2) * child[pos]) - 1;
				if ((*interval_end < start) || (*interval_end > end))
					*interval_end = end;
			} else if (chrom[sa[start] + position] == 'G') {
				*interval_begin = start;
				*interval_end = index_begin - 1;
				if ((*interval_end < start) || (*interval_end > end))
					*interval_end = end;
			} else if (chrom[sa[index_begin] + position] == 'G') {
				*interval_begin = index_begin;
				*interval_end = ((verifyChild[index_begin] > 2)
						* child[index_begin]) - 1;
				if ((*interval_end < start) || (*interval_end > end))
					*interval_end = end;
			}
		} else //(ch == 'T')
		{

			unsigned int child_up_value_tmp = ((verifyChild[end] == 1)
					* child[end]);
			if ((start < child_up_value_tmp) && (end >= child_up_value_tmp))
				index_begin = child_up_value_tmp; //getChildUpValue(child, verifyChild, end+1);
			else
				index_begin = getChildDownValue(child, verifyChild, start);

			if (chrom[sa[start] + position] == 'T') {
				*interval_begin = start;
				*interval_end = index_begin - 1;
				if ((*interval_end < start) || (*interval_end > end))
					*interval_end = end;
			} else if (chrom[sa[index_begin] + position] == 'T') {
				*interval_begin = index_begin;
				//*interval_end = child_next[index_begin]-1;
				*interval_end = ((verifyChild[index_begin] > 2)
						* child[index_begin]) - 1;
				if ((*interval_end < start) || (*interval_end > end))
					*interval_end = end;
			} else {
				//pos = child_next[index_begin];
				pos = ((verifyChild[index_begin] > 2) * child[index_begin]);
				if (chrom[sa[pos] + position] == 'T') {
					*interval_begin = pos;
					//*interval_end = child_next[pos] - 1;
					*interval_end = ((verifyChild[pos] > 2) * child[pos]) - 1;
					if ((*interval_end < start) || (*interval_end > end))
						*interval_end = end;
				} else {
					//pos = child_next[pos];
					pos = ((verifyChild[pos] > 2) * child[pos]);
					if (chrom[sa[pos] + position] == 'T') {
						*interval_begin = pos;
						//*interval_end = child_next[pos] - 1;
						*interval_end = ((verifyChild[pos] > 2) * child[pos])
								- 1;
						if ((*interval_end < start) || (*interval_end > end))
							*interval_end = end;
					}
				}
			}
		}
	}
}

unsigned int getlcp(unsigned int start, unsigned int end, BYTE* lcpCompress,
		unsigned int* child, BYTE* verifyChild) {
	unsigned int tmpIndex;

	unsigned int child_up_value_tmp = ((verifyChild[end] == 1) * child[end]);
	if ((start < child_up_value_tmp) && (end >= child_up_value_tmp))
		tmpIndex = child_up_value_tmp;
	else
		tmpIndex = getChildDownValue(child, verifyChild, start);

	return lcpCompress[tmpIndex];
}

//input : read, sa, up, down, next, chrom;
//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc
bool mapMainForIndexStringHash(char *read, unsigned int* sa, BYTE* lcpCompress,
		unsigned int* child, char* chrom, int* mappedLength,
		unsigned int* indexIntervalStart, unsigned int* indexIntervalEnd,
		BYTE* verifyChild, int readLength, int MAX) {
	bool mapMain = false;
	unsigned int stop_loc = 0; // location in one segment for iterations
	unsigned int read_length = readLength; //READ_LENGTH;
	unsigned int interval_begin, interval_end;

	unsigned int n = MAX; //size of SA
	char* read_local = read;
	bool queryFound = false;

	if ((*read_local != 'A') && (*read_local != 'C') && (*read_local != 'G')
			&& (*read_local != 'T')) {
		(*mappedLength) = 0;
		return true;
	}
	unsigned int lcp_length = 0;
	unsigned int start = 0, end = n - 1;
	unsigned int Min;
	unsigned int c = 0;
	//cout << "read " << (*read_local) << endl;
	getFirstInterval(*read_local, &interval_begin, &interval_end, child,
			verifyChild);
	//cout << "firstInterval: " << interval_begin << " ~ " << interval_end << endl;

	unsigned int iterateNum = 0; //debug;
	//cout << "interval_begin = " << interval_begin << endl << "interval_end = " << interval_end << endl;
	while (c < read_length) {
		iterateNum++;
		//cout << "iterateNum: " << iterateNum << endl;
		if (iterateNum > read_length)
			return false;

		unsigned int c_old = c;

		if (interval_begin != interval_end) {
			lcp_length = getlcp(interval_begin, interval_end, lcpCompress,
					child, verifyChild);
			Min = min(lcp_length, read_length);

			unsigned int loc_pos = 0;
			for (loc_pos = 0; loc_pos < Min - c_old && !queryFound; loc_pos++)
				queryFound = (*(read_local + c_old + loc_pos)
						== *(chrom + sa[interval_begin] + c_old + loc_pos));

			//cout << "queryFound: " << queryFound << endl;
			if (!queryFound) {
				stop_loc = c_old + loc_pos;
				(*mappedLength) = stop_loc;
				(*indexIntervalStart) = interval_begin;
				(*indexIntervalEnd) = interval_end;
				return true;
			}

			c = Min;
			if (*(read_local + c) == 'N') {
				stop_loc = c;
				(*mappedLength) = stop_loc;
				(*indexIntervalStart) = interval_begin;
				(*indexIntervalEnd) = interval_end;
				return true;
			}
			start = interval_begin;
			end = interval_end;
			if (c == read_length) {
				(*mappedLength) = read_length;
				(*indexIntervalStart) = interval_begin;
				(*indexIntervalEnd) = interval_end;
				return true;
			}
			unsigned int interval_begin_ori = interval_begin;
			unsigned int interval_end_ori = interval_end;
			getInterval(start, end, c, *(read_local + c), &interval_begin,
					&interval_end, sa, child, chrom, verifyChild);

			//cout << "firstInterval: " << interval_begin << " ~ " << interval_end << endl;
			if (interval_begin > interval_end) {
				queryFound = false;
				stop_loc = c - 1;
				//cout << "interval_begin > interval_end" << endl; cout << "stop_loc = " << stop_loc << endl;
				(*mappedLength) = stop_loc;
				(*indexIntervalStart) = interval_begin_ori;
				(*indexIntervalEnd) = interval_end_ori;
				return true;
			}
		} //end if(interval_begin != interval_end)
		else {
			unsigned int loc_pos = 0;
			for (loc_pos = 0; loc_pos < read_length - c && !queryFound;
					loc_pos++)
				queryFound = (*(read_local + c + loc_pos)
						== *(chrom + sa[interval_begin] + c + loc_pos));

			if (queryFound) {
				(*mappedLength) = read_length;
				(*indexIntervalStart) = interval_begin;
				(*indexIntervalEnd) = interval_end;
				return true;
				//align_length[101] ++;
			} else {
				stop_loc = c + loc_pos;
				(*mappedLength) = stop_loc;
				(*indexIntervalStart) = interval_begin;
				(*indexIntervalEnd) = interval_end;
				return true;
			}
		}
	} //end while((c < read_length))

	return true;
}

void build_up_down(unsigned int *lcptab, unsigned int *up, unsigned int *down,
		unsigned int n) {
	unsigned int lastIndex = NULL_NUM;
	stack<unsigned int> up_down;
	up_down.push(0);
	unsigned int i;
	for (i = 0; i < n; i++) {
		while (lcptab[i] < lcptab[up_down.top()]) {
			lastIndex = up_down.top();
			up_down.pop();
			if ((lcptab[i] <= lcptab[up_down.top()])
					&& (lcptab[up_down.top()] != lcptab[lastIndex]))
				down[up_down.top()] = lastIndex;
		}
		// now lcptab[i] >= lcptab[up_down.top()] holds
		if (lastIndex != NULL_NUM) {
			up[i] = lastIndex;
			lastIndex = NULL_NUM;
		}
		up_down.push(i);
	}
	return;
}

void build_next(unsigned int *lcptab, unsigned int *next, unsigned int n) {
	stack<unsigned int> nextIndex;
	unsigned int lastIndex;
	unsigned int j;
	nextIndex.push(0);
	for (j = 0; j < n; j++) {
		while (lcptab[j] < lcptab[nextIndex.top()])
			nextIndex.pop();
		if (lcptab[j] == lcptab[nextIndex.top()]) {
			lastIndex = nextIndex.top();
			nextIndex.pop();
			next[lastIndex] = j;
		}
		nextIndex.push(j);
	}
	return;
}

///////
// Builds the longest common prefix array for the chromosome
// A LCP augments the suffix array and allows the LCP array to efficiently
// simulate top-down and bottom-up traversals of the suffix tree, speeds up pattern
// matching on the suffix array and is a prerequisite for compressed suffix trees.
///////
void build_LongestCommonPrefix(const unsigned int *r, const unsigned int *sa,
		unsigned int *lcp, unsigned int *rank, const unsigned int n) {

	unsigned int i, j;
	int k = 0;
	for (i = 0; i < n; i++)
		rank[sa[i]] = i;
	cout << "stop10_new" << endl;

	for (i = 0; i < n; i++) {
		if (k)
			k--;

		if (rank[i] == 0)
			j = 0;
		else
			j = sa[rank[i] - 1];

		while (r[i + k] == r[j + k])
			k++;

		lcp[rank[i]] = k;
	}

	lcp[0] = 0;
	cout << "stop11" << endl;
	return;
}

const unsigned int compare(unsigned int *r, unsigned int a, unsigned int b,
		unsigned int l) {
	return r[a] == r[b] && r[a + l] == r[b + l];
}

inline bool leq(int a1, int a2, int b1, int b2) // lexicographic order
{
	return (a1 < b1 || (a1 == b1 && a2 <= b2));
} // for pairs

inline bool leq(int a1, int a2, int a3, int b1, int b2, int b3)
{
	return (a1 < b1 || (a1 == b1 && leq(a2, a3, b2, b3)));
} // and triples

// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
static void radixPass(unsigned int* a, unsigned int* b, unsigned int* r, int n, int K)
{
	// count occurrences
	int* c = new int[K + 1]; // counter array

	for (int i = 0; i <= K; i++)
		c[i] = 0; // reset counters

	for (int i = 0; i < n; i++)
		c[r[a[i]]]++; // count occurrences

	for (int i = 0, sum = 0; i <= K; i++) // exclusive prefix sums
	{
		int t = c[i];
		c[i] = sum;
		sum += t;
	}

	for (int i = 0; i < n; i++)
		b[c[r[a[i]]]++] = a[i]; // sort

	delete[] c;
}

// find the suffix array SA of T[0..n-1] in {1..K}^n
// require T[n]=T[n+1]=T[n+2]=0, n>=2
void suffixArray(unsigned int* T, unsigned int* SA, int n, int K) {
	int n0 = (n + 2) / 3,
		n1 = (n + 1) / 3,
		n2 = n / 3,
		n02 = n0 + n2;

	unsigned int* R = new unsigned int[n02 + 3];
	R[n02] = R[n02 + 1] = R[n02 + 2] = 0;
	unsigned int* SA12 = new unsigned int[n02 + 3];
	SA12[n02] = SA12[n02 + 1] = SA12[n02 + 2] = 0;
	unsigned int* R0 = new unsigned int[n0];
	unsigned int* SA0 = new unsigned int[n0];

	//******* Step 0: Construct sample ********
	// generate positions of mod 1 and mod 2 suffixes
	// the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
	for (int i = 0, j = 0; i < n + (n0 - n1); i++)
		if (i % 3 != 0)
			R[j++] = i;

	//******* Step 1: Sort sample suffixes ********
	// lsb radix sort the mod 1 and mod 2 triples
	radixPass(R, SA12, T + 2, n02, K);
	radixPass(SA12, R, T + 1, n02, K);
	radixPass(R, SA12, T, n02, K);

	// find lexicographic names of triples and
	// write them to correct places in R
	int name = 0, c0 = -1, c1 = -1, c2 = -1;
	for (int i = 0; i < n02; i++)
	{
		if (T[SA12[i]] != c0 || T[SA12[i] + 1] != c1 || T[SA12[i] + 2] != c2)
		{
			name++;
			c0 = T[SA12[i]];
			c1 = T[SA12[i] + 1];
			c2 = T[SA12[i] + 2];
		}

		if (SA12[i] % 3 == 1)
		{
			// write to R1
			R[SA12[i] / 3] = name;
		}
		else
		{
			// write to R2
			R[SA12[i] / 3 + n0] = name;
		}
	}

	// recurse if names are not yet unique
	if (name < n02)
	{
		suffixArray(R, SA12, n02, name);

		// store unique names in R using the suffix array
		for (int i = 0; i < n02; i++)
			R[SA12[i]] = i + 1;
	}
	else // generate the suffix array of R directly
		for (int i = 0; i < n02; i++)
			SA12[R[i] - 1] = i;

	//******* Step 2: Sort nonsample suffixes ********
	// stably sort the mod 0 suffixes from SA12 by their first character
	for (int i = 0, j = 0; i < n02; i++)
		if (SA12[i] < n0)
			R0[j++] = 3 * SA12[i];
	radixPass(R0, SA0, T, n0, K);

	//******* Step 3: Merge ********
	// merge sorted SA0 suffixes and sorted SA12 suffixes
	for (int p = 0, t = n0 - n1, k = 0; k < n; k++)
	{
		#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
		int i = GetI(); // pos of current offset 12 suffix
		int j = SA0[p]; // pos of current offset 0 suffix
		if (SA12[t] < n0 // different compares for mod 1 and mod 2 suffixes
				? leq(T[i], R[SA12[t] + n0], T[j], R[j / 3])
				: leq(T[i], T[i + 1], R[SA12[t] - n0 + 1], T[j], T[j + 1], R[j / 3 + n0]))
		{ 	// suffix from SA12 is smaller

			SA[k] = i;
			t++;
			if (t == n02) // done --- only SA0 suffixes left
				for (k++; p < n0; p++, k++)
					SA[k] = SA0[p];
		}
		else // suffix from SA0 is smaller
		{
			SA[k] = j;
			p++;
			if (p == n0) // done --- only SA12 suffixes left
				for (k++; t < n02; t++, k++)
					SA[k] = GetI();
		}
	}
	delete[] R;
	delete[] SA12;
	delete[] SA0;
	delete[] R0;
}

int main(int argc, char** argv) {
	// Installs our error handler
	signal(SIGSEGV, handler);

	// Validating Input
	if (argc != EXPECTED_ARGUMENT_COUNT) {
		cout
				<< "Correct execution is: Executable <InputChromosomesFolder> <outputIndexFolder>"
				<< endl;
		exit(0);
	}

	BuildIndexParameter_info* buildIndexInfo = new BuildIndexParameter_info(
			argv[1], argv[2]);

	string chrom_file_str = buildIndexInfo->GetOutputIndexFolder()
			+ "_chromMerge";

	ofstream chrom_file_ofs(chrom_file_str.c_str(), ios::binary);
	ofstream parameter_file_ofs(
			(buildIndexInfo->GetOutputIndexFolder() + "_parameter").c_str());

	/////////////////////////////////////////////////////////////////////////////////
	////////////// scan all chromosomes, generate index_parameter, index_chrom //////
	/////////////////////////////////////////////////////////////////////////////////

	DIR *dp;
	struct dirent *dirp;
	if ((dp = opendir(buildIndexInfo->GetInputChromFolder().c_str())) == NULL) {
		printf("can't open %s", buildIndexInfo->GetInputChromFolder().c_str());
	}

	while ((dirp = readdir(dp)) != NULL)
		if (dirp->d_type == 8)
			buildIndexInfo->AddChromName(dirp->d_name);

	closedir(dp);

	char chrFileLine[100];
	string chrFileLineStr;
	int lineNum = 0;
	int tmpChromEndPosInGenome = 0;

	for (vector<string>::iterator myIterator = buildIndexInfo->begin();
			myIterator != buildIndexInfo->end(); myIterator++) {
		lineNum = 0;
		string tmpChromFileName = buildIndexInfo->GetInputChromFolder() + "/"
				+ *myIterator;

		FILE *fp_tmpChr = fopen(tmpChromFileName.c_str(), "r");
		//cout << endl << "readFileName " << tmpChromNum + 1 << ": " << tmpChromFileName << endl;

		fgets(chrFileLine, sizeof(chrFileLine), fp_tmpChr);
		chrFileLineStr = chrFileLine;
		fgets(chrFileLine, sizeof(chrFileLine), fp_tmpChr);
		chrFileLineStr = chrFileLine;
		lineNum++;
		int chrBaseNum = 0;
		while (!feof(fp_tmpChr)) {
			lineNum++;
			fgets(chrFileLine, sizeof(chrFileLine), fp_tmpChr);

			chrFileLineStr = chrFileLine;

			if (feof(fp_tmpChr)) {
				if (chrBaseNum == 0)
					chrBaseNum = 50;

				buildIndexInfo->AddChromLength((lineNum - 2) * 50 + chrBaseNum);
				break;
			}

			if (chrFileLineStr.length() != 51)
				chrBaseNum = chrFileLineStr.length() - 1;

		}
		fclose(fp_tmpChr);

	}

	buildIndexInfo->getChrEndPosVec();
	buildIndexInfo->outputChromSeq(chrom_file_ofs);
	buildIndexInfo->outputIndexParameter(parameter_file_ofs);

	/////////////////////////////////////////////////
	//////generate original size index
	/////////////////////////////////////////////////

	int MAX = buildIndexInfo->GetMax();
	cout << "MAX: " << MAX << endl;

	ofstream SA_file_ofs(
			(buildIndexInfo->GetOutputIndexFolder() + "_SA").c_str(),
			ios::binary);
	ofstream lcp_file_ofs(
			(buildIndexInfo->GetOutputIndexFolder() + "_lcp").c_str(),
			ios::binary);
	ofstream up_file_ofs(
			(buildIndexInfo->GetOutputIndexFolder() + "_up").c_str(),
			ios::binary);
	ofstream down_file_ofs(
			(buildIndexInfo->GetOutputIndexFolder() + "_down").c_str(),
			ios::binary);
	ofstream next_file_ofs(
			(buildIndexInfo->GetOutputIndexFolder() + "_next").c_str(),
			ios::binary);
	ofstream chrom_bit_file_ofs(
			(buildIndexInfo->GetOutputIndexFolder() + "_chrom").c_str(),
			ios::binary);

	FILE *fp = fopen(chrom_file_str.c_str(), "r");
	unsigned int *r = (unsigned int*) malloc(MAX * sizeof(unsigned int));
	char ch;
	char head[100];
	char base[Range] = { 'X', 'A', 'C', 'G', 'T', 'N' };
	char *chrom = (char*) malloc(MAX * sizeof(char));

	unsigned int chrom_base_num = 0;

	// Read from the chromosome file and build two chromosome arrays
	while ((ch = fgetc(fp)) != EOF) {
		switch (ch) {
		case 'A':
		case 'a':
			chrom[chrom_base_num] = 'A';
			r[chrom_base_num] = 1;
			break;
		case 'C':
		case 'c':
			chrom[chrom_base_num] = 'C';
			r[chrom_base_num] = 2;
			break;
		case 'G':
		case 'g':
			chrom[chrom_base_num] = 'G';
			r[chrom_base_num] = 3;
			break;
		case 'T':
		case 't':
			chrom[chrom_base_num] = 'T';
			r[chrom_base_num] = 4;
			break;
		case 'N':
		case 'n':
			chrom[chrom_base_num] = 'N';
			r[chrom_base_num] = 5;
			break;
		case 'X':
			chrom[chrom_base_num] = 'X';
			r[chrom_base_num] = 6;
			break;
		case '\t':
		case '\n':
			continue;
		default:
			printf("\n illegal input is '%c'", ch);
			exit(1);
			break;

		}
		chrom_base_num++;
	}

	cout << "the number of bases in Chromo is " << chrom_base_num << endl;
	cout << "chrom is ready" << endl;
	r[MAX - 1] = Range;
	chrom[MAX - 1] = 'X';

	cout << "chrom[MAX-3]: " << chrom[MAX - 3] << endl;
	cout << "chrom[MAX-2]: " << chrom[MAX - 2] << endl;
	cout << "chrom[MAX-1]: " << chrom[MAX - 1] << endl;

	chrom_bit_file_ofs.write((const char*) chrom, MAX * sizeof(char));

	unsigned int *sa = (unsigned int*) malloc(MAX * sizeof(unsigned int));

	cout << "start to build SA array" << endl;
	suffixArray(r, sa, MAX, Range + 1);

	cout << "SA is ready" << endl;
	//////////////////////////////////////////////////
	//unsigned int rank[MAX]={0}, lcp[MAX]={0}, up[MAX] = {0}, down[MAX] = {0}, next[MAX] = {0};
	unsigned int *rank = (unsigned int*) malloc(MAX * sizeof(unsigned int));
	unsigned int *lcp = (unsigned int*) malloc(MAX * sizeof(unsigned int));
	unsigned int *up = (unsigned int*) malloc(MAX * sizeof(unsigned int));
	unsigned int *down = (unsigned int*) malloc(MAX * sizeof(unsigned int));
	unsigned int *next = (unsigned int*) malloc(MAX * sizeof(unsigned int));

	cout << "start to build LCP array" << endl;
	build_LongestCommonPrefix(r, sa, lcp, rank, MAX); //build lcp array
	free(r);
	free(rank);

	cout << "lcp is ready" << endl;
	cout << "start to build up & down array" << endl;
	build_up_down(lcp, up, down, MAX); // build up and down array
	cout << "up & down are ready" << endl;
	cout << "start to build next array" << endl;
	build_next(lcp, next, MAX); //build next array
	cout << "next is ready" << endl;

	cout << "start to output Index to file " << endl;
	SA_file_ofs.write((const char*) sa, MAX * sizeof(unsigned int));
	lcp_file_ofs.write((const char*) lcp, MAX * sizeof(unsigned int));
	up_file_ofs.write((const char*) up, MAX * sizeof(unsigned int));
	down_file_ofs.write((const char*) down, MAX * sizeof(unsigned int));
	next_file_ofs.write((const char*) next, MAX * sizeof(unsigned int));
	//free(sa);
	cout << "finish writing index to files " << endl;
	/////////////////////////////////////////////////////////////////////
	/////  compress index for whole genome original size index
	/////////////////////////////////////////////////////////////////////

	unsigned int indexSize = MAX;

	ofstream childTab_file_ofs(
			(buildIndexInfo->GetOutputIndexFolder() + "_childTab").c_str(),
			ios::binary);
	ofstream detChild_file_ofs(
			(buildIndexInfo->GetOutputIndexFolder() + "_detChild").c_str(),
			ios::binary);

	unsigned int *childTab;
	childTab = (unsigned int*) malloc(indexSize * sizeof(unsigned int));

	BYTE *verifyChild;
	verifyChild = (BYTE*) malloc(indexSize * sizeof(BYTE));

	cout << "start to compress original size index " << endl;
	BuildIndex_Info* tmpBuildIndexInfo = new BuildIndex_Info();

	tmpBuildIndexInfo->compressUpDownNext2ChildtabVerifyChild(up, down, next,
			childTab, verifyChild, indexSize);
	//cout << "finish compressing index" << endl;
	childTab_file_ofs.write((const char*) childTab,
			indexSize * sizeof(unsigned int));
	detChild_file_ofs.write((const char*) verifyChild,
			indexSize * sizeof(BYTE));

	childTab_file_ofs.close();
	detChild_file_ofs.close();

	cout << "finish compressing index" << endl;
	free(up);
	free(down);
	free(next);

	ofstream lcpCompress_file_ofs(
			(buildIndexInfo->GetOutputIndexFolder() + "_lcpCompress").c_str(),
			ios::binary);

	//BuildIndex_Info* tmpBuildIndexInfo = new BuildIndex_Info();
	cout << "start to compress Lcp array" << endl;
	BYTE *lcpCompress = (BYTE*) malloc(indexSize * sizeof(BYTE));
	tmpBuildIndexInfo->compressLcp2Lcpcompress(lcp, lcpCompress, indexSize);
	cout << "finish compressing Lcp array " << endl;
	lcpCompress_file_ofs.write((const char*) lcpCompress,
			indexSize * sizeof(BYTE));

	free(lcp);
	//free(lcpCompress);
	lcpCompress_file_ofs.close();

	////////////////////////////////////////////////////////
	/////  start to generate preIndex String 
	////////////////////////////////////////////////////////

	int preIndexStringLength = PREINDEX_STRINGLENGTH;
	string generateStringFilePrefix = buildIndexInfo->GetOutputIndexFolder()
			+ "_preIndexString";
	int stringLengthInHash = preIndexStringLength;
	cout << "start to write preIndex string files" << endl;

	string baseStr[4] = { "A", "C", "G", "T" };
	string generateStringFile[stringLengthInHash];

	generateStringFile[0] = generateStringFilePrefix + ".1";

	ofstream GenerateStringFile_ofs(generateStringFile[0].c_str());
	GenerateStringFile_ofs << "A" << endl << "C" << endl << "G" << endl << "T"; // need to debug ...
	GenerateStringFile_ofs.close();
	string readString;
	char readChar[20];
	for (int tmp = 1; tmp < stringLengthInHash; tmp++) {
		char tmpChar[2];
		sprintf(tmpChar, "%d", tmp + 1);
		string tmpString = tmpChar;
		generateStringFile[tmp] = generateStringFilePrefix + "." + tmpString;

		FILE *fp_in = fopen(generateStringFile[tmp - 1].c_str(), "r");
		ofstream GenerateStringFile_ofs(generateStringFile[tmp].c_str());

		while (!feof(fp_in)) {
			fgets(readChar, sizeof(readChar), fp_in);
			readString = readChar;
			readString = readString.substr(0, tmp);
			for (int tmpBase = 0; tmpBase < 4; tmpBase++) {
				string resultString = readString + baseStr[tmpBase];
				GenerateStringFile_ofs << resultString << endl;
			}
		}
		fclose(fp_in);
		GenerateStringFile_ofs.close();

	}
	cout << "finish writing preIndex string files" << endl;

	cout << "start to generate preIndex string mapping record" << endl;
	string preIndexStringFileStr = generateStringFilePrefix + "."
			+ int_to_str(preIndexStringLength);
	FILE *fp_in_2 = fopen(preIndexStringFileStr.c_str(), "r");

	ofstream StringHashRecordFile_ofs(
			(buildIndexInfo->GetOutputIndexFolder() + "_preIndexRecord").c_str());

	char read[20];
	int readLength = stringLengthInHash;
	int preIndexStringNum = 0;
	while (!feof(fp_in_2)) {
		preIndexStringNum++;
		//cout << "preIndexStringNum: " << preIndexStringNum << endl;
		fgets(read, sizeof(read), fp_in_2);

		int mappedLength;
		unsigned int indexIntervalStart, indexIntervalEnd;
		bool mapMain = mapMainForIndexStringHash(read, sa, lcpCompress,
				childTab, chrom, &mappedLength, &indexIntervalStart,
				&indexIntervalEnd, verifyChild, readLength, MAX);
		string readString = read;
		readString = readString.substr(0, readLength);
		if (mapMain) {
			StringHashRecordFile_ofs << readString << "\t" << mappedLength
					<< "\t" << indexIntervalStart << "\t" << indexIntervalEnd
					<< endl;
		} else {
			cout << "mapMain error " << endl;
		}
	}

	fclose(fp_in_2);
	StringHashRecordFile_ofs.close();
	free(sa);
	free(lcpCompress); //free(child_up);free(child_down);free(child_next);
	free(childTab);
	free(chrom);
	free(verifyChild);
	cout << "finish writing preIndex string mapping record !" << endl;

	cout << "start to generate preIndex array" << endl;

	FILE* fp_in_preIndex =
			fopen(
					(buildIndexInfo->GetOutputIndexFolder() + "_preIndexRecord").c_str(),
					"r");

	cout << "start to creat mapLength, intervalBegin, intervalEnd files"
			<< endl;

	ofstream preIndexMapLengthArray_ofs(
			(buildIndexInfo->GetOutputIndexFolder() + "_MapLength").c_str(),
			ios::binary);
	ofstream preIndexIntervalStartArray_ofs(
			(buildIndexInfo->GetOutputIndexFolder() + "_IntervalStart").c_str(),
			ios::binary);
	ofstream preIndexIntervalEndArray_ofs(
			(buildIndexInfo->GetOutputIndexFolder() + "_IntervalEnd").c_str(),
			ios::binary);

	cout << "finish creating mapLength, intervalBegin, intervalEnd files"
			<< endl;

	int* preIndexMapLengthArray;
	preIndexMapLengthArray = (int*) malloc(PreIndexSize * sizeof(int));

	unsigned int *preIndexIntervalStartArray;
	preIndexIntervalStartArray = (unsigned int*) malloc(
			PreIndexSize * sizeof(unsigned int));

	unsigned int *preIndexIntervalEndArray;
	preIndexIntervalEndArray = (unsigned int*) malloc(
			PreIndexSize * sizeof(unsigned int));

	char preIndexRecordChar[100];
	char preIndexStrChar[20];
	char preIndexMapLengthChar[10];
	char preIndexIntervalStartChar[20];
	char preIndexIntervalEndChar[20];

	int preIndexMapLengthInt;
	unsigned int preIndexIntervalStartInt;
	unsigned int preIndexIntervalEndInt;

	unsigned int preIndexNO = 0;

	cout << "start to extract record from record file" << endl;
	while (!feof(fp_in_preIndex)) {
		//cout << "preIndexNO: " << preIndexNO;// << endl;
		fgets(preIndexRecordChar, sizeof(preIndexRecordChar), fp_in_preIndex);
		//cout << "record: " << preIndexRecordChar << endl;

		sscanf(preIndexRecordChar, "%s\t%s\t%s\t%s", preIndexStrChar,
				preIndexMapLengthChar, preIndexIntervalStartChar,
				preIndexIntervalEndChar);
		preIndexMapLengthInt = atoi(preIndexMapLengthChar);
		preIndexIntervalStartInt = strtoul(preIndexIntervalStartChar, NULL, 10);
		preIndexIntervalEndInt = strtoul(preIndexIntervalEndChar, NULL, 10);

		string tmpPreIndexStr = preIndexStrChar;
		preIndexNO = getPreIndexNO(tmpPreIndexStr);

		preIndexMapLengthArray[preIndexNO] = preIndexMapLengthInt;
		preIndexIntervalStartArray[preIndexNO] = preIndexIntervalStartInt;
		preIndexIntervalEndArray[preIndexNO] = preIndexIntervalEndInt;

	}
	cout << "finish extracting record from record file" << endl;
	//delete(preIndexSegInfo);
	fclose(fp_in_preIndex);

	preIndexMapLengthArray_ofs.write((const char*) preIndexMapLengthArray,
			PreIndexSize * sizeof(int));
	preIndexIntervalStartArray_ofs.write(
			(const char*) preIndexIntervalStartArray,
			PreIndexSize * sizeof(unsigned int));
	preIndexIntervalEndArray_ofs.write((const char*) preIndexIntervalEndArray,
			PreIndexSize * sizeof(unsigned int));

	//cout << "mapLengthArray[0] = " << preIndexMapLengthArray[0] << endl;
	//cout << "mapLengthArray[PreIndexSize-2] = " << preIndexMapLengthArray[PreIndexSize-2] << endl;
	//cout << "mapLengthArray[PreIndexSize-1] = " << preIndexMapLengthArray[PreIndexSize-1] << endl;

	//cout << "preIndexIntervalStartArray[0] = " << preIndexIntervalStartArray[0] << endl;
	//cout << "preIndexIntervalStartArray[PreIndexSize-2] = " << preIndexIntervalStartArray[PreIndexSize-2] << endl;
	//cout << "preIndexIntervalStartArray[PreIndexSize-1] = " << preIndexIntervalStartArray[PreIndexSize-1] << endl;

	//cout << "preIndexIntervalEndArray[0] = " << preIndexIntervalEndArray[0] << endl;
	//cout << "preIndexIntervalEndArray[PreIndexSize-2] = " << preIndexIntervalEndArray[PreIndexSize-2] << endl;
	//cout << "preIndexIntervalEndArray[PreIndexSize-1] = " << preIndexIntervalEndArray[PreIndexSize-1] << endl;

	free(preIndexMapLengthArray);
	free(preIndexIntervalStartArray);
	free(preIndexIntervalEndArray);

	cout << "finish generating preIndex array" << endl;
	cout << endl << "all index building jobs done" << endl;
	return 0;
}
