/*    
 *    PE_match_junction.cpp		
 *    MapSplice
 *
 *    Copyright (C) 2010 University of Kentucky and
 *                       Yin Hu
 *
 *    Authors: Yin Hu
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define UNIX

#ifdef UNIX
#include <fstream>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cstdlib> 
#include <cmath>
#include <ctime>
#else
#include <fstream>
#include <stdio.h>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <assert.h>
#endif


using namespace std;

const long FUSION_CUTOFF = 50000;

const long MAX = 2000000000;
const long MAX_CHR_COMB_NUM = 1000;
const int MAX_CHOICE_NUM = 1000000;
const int MAX_COMB_LIMIT = 10;


string targetChromosome[MAX_CHR_COMB_NUM + 1];
unsigned long targetChrNum;
ofstream outputfile[MAX_CHR_COMB_NUM + 1];
ofstream outputfile_chrname;
ofstream dropped_reads;

string dirPrefix;

class reads
{
public:
	string name;
	string chromosome;
	string strand;
	long startpoint;
	long endpoint;
	string end;

	void compute_endpoint();
};

reads* startlist[MAX_CHOICE_NUM];
reads* endlist[MAX_CHOICE_NUM];
int startlist_count;
int endlist_count;

ofstream allDistanceToFile;


void reads::compute_endpoint()
{
	long tmp = 0;

	endpoint = startpoint - 1;

	for (unsigned long tmpLoop = 0; tmpLoop < end.size(); ++tmpLoop)
	{
		if (end[tmpLoop] == 'M')
		{
			endpoint += tmp;
			tmp = 0;
		} 
		else if (end[tmpLoop] == 'N')
		{
			endpoint += tmp;
			tmp = 0;
		}
		else if (end[tmpLoop] == 'I')
		{
			tmp = 0;
		}
		else if (end[tmpLoop] == 'D' || end[tmpLoop] == 'S')
		{
			endpoint += tmp;
			tmp = 0;
		}
		else if (end[tmpLoop] >= '0' && end[tmpLoop] <= '9')
		{
			tmp = tmp * 10 + end[tmpLoop] - 48;
		}
		else
		{
			tmp = 0;
		}
	}

	return;
}



long compute_distance(reads* start_read, reads* end_read)
{
	//compute the min distance between the start_read and the end_read

//	return MAX;
	return end_read->startpoint - start_read->endpoint;
}


void dec2bin(long decimal, char *temp)
{
	int  k = 0, n = 0;
	int  neg_flag = 0;
	int  remain;
	int  old_decimal;  // for test

	// take care of negative input
	if (decimal < 0)
	{      
		decimal = -decimal;
		neg_flag = 1;
	}
	do 
	{
		old_decimal = decimal;   // for test
		remain    = decimal % 2;
		// whittle down the decimal number
		decimal   = decimal / 2;
		// this is a test to show the action
		//printf("%d/2 = %d  remainder = %d\n", old_decimal, decimal, remain);
		// converts digit 0 or 1 to character '0' or '1'
		temp[k++] = remain + '0';
	} while (decimal > 0);

	if (neg_flag)
		temp[k++] = '-';       // add - sign

	// reverse the spelling
//	while (k > 0)
//		binary[n++] = temp[--k];

//	binary[n] = '\0';         // end with NULL
}

void outputPEreads(reads* start_read, reads* end_read)
{
	//output a single read
	string outpurfilename;
	string chromosome;
	unsigned long targetChrIndex, tmp;
	
	chromosome = start_read->chromosome + "~" + end_read->chromosome;

	for (tmp = 1; tmp <= targetChrNum; ++tmp)
		if (targetChromosome[tmp].compare(chromosome) == 0)
			break;

	if (tmp > targetChrNum)
	{
		if (tmp > MAX_CHR_COMB_NUM)
			//reach the maximum combination number 
			return;

		//chromosome not found
		outputfile_chrname << chromosome << endl;

		outpurfilename = dirPrefix + "data/parsedPER/" + chromosome + ".txt";
		targetChromosome[++targetChrNum] = chromosome;
		outputfile[targetChrNum].open(outpurfilename.c_str());
	}

	targetChrIndex = tmp;

	//CHANGE: change the output format, merge + - and - + into + -, also merge + + and - - into + +
//	outputfile[targetChrIndex] << start_read->name << "\t" << start_read->chromosome << "\t" << start_read->strand << "\t" << start_read->startpoint << "\t" << start_read->end << "\t" << end_read->chromosome << "\t" << end_read->strand << "\t" << end_read->startpoint << "\t" << end_read->end << endl; 

	if (start_read->strand.compare("0") == 0 && end_read->strand.compare("16") == 0 || start_read->strand.compare("16") == 0 && end_read->strand.compare("0") == 0 || start_read->strand.compare("+") == 0 && end_read->strand.compare("-") == 0 || start_read->strand.compare("-") == 0 && end_read->strand.compare("+") == 0)
	{
		outputfile[targetChrIndex] << start_read->name << "\t" << start_read->chromosome << "\t" << "0" << "\t" << start_read->startpoint << "\t" << start_read->end << "\t" << end_read->chromosome << "\t" << "16" << "\t" << end_read->startpoint << "\t" << end_read->end << endl; 
	}
	else// if (strcmp(start_read->strand, "0") == 0 && strcmp(end_read->strand, "0") == 0 || strcmp(start_read->strand, "16") == 0 && strcmp(end_read->strand, "16") == 0 || strcmp(start_read->strand, "+") == 0 && strcmp(end_read->strand, "+") == 0 || strcmp(start_read->strand, "-") == 0 && strcmp(end_read->strand, "-") == 0)
	{
		outputfile[targetChrIndex] << start_read->name << "\t" << start_read->chromosome << "\t" << "0" << "\t" << start_read->startpoint << "\t" << start_read->end << "\t" << end_read->chromosome << "\t" << "0" << "\t" << end_read->startpoint << "\t" << end_read->end << endl; 
	}

// 	long dist;
// 	dist = end_read->startpoint - start_read->endpoint - 1;
// 
// 	if (dist >= 0 && dist < 50000)
// 	{
// 		allDistanceToFile << dist << endl;
// 	}
	return;
}

// void outputPEfile()
// {
// 	//unique alignment
// 
// 	int i, j, best_i, best_j;
// 	long distance_i_j, best_distance;
// 
// 	best_distance = MAX;
// 	best_i = 1;
// 	best_j = 1;
// 	
// 	if (startlist_count * endlist_count > MAX_COMB_LIMIT)
// 	{
// 		for (i = 1; i <= startlist_count; i++)
// 		{
// 			dropped_reads << startlist[i]->name << "\t1\t" << startlist[i]->chromosome << '\t' << startlist[i]->startpoint << '\t' << startlist[i]->end << endl;
// 		}
// 		for (j = 1; j <= endlist_count; j++)
// 		{
// 			dropped_reads << endlist[j]->name << "\t2\t" << endlist[j]->chromosome << '\t' << endlist[j]->startpoint << '\t' << endlist[j]->end << endl;
// 		}
// 
// 		return;
// 	}
// 	
// 	for (i = 1; i <= startlist_count; i++)
// 	{
// 		for (j = 1; j <= endlist_count; j++)
// 		{
// // 			if (strcmp(startlist[i]->chromosome, endlist[j]->chromosome) != 0)
// // 			{
// // 				//not on the same chromosome
// // 				distance_i_j = MAX;
// // 			} 
// // 			else
// // 			{
// // 				//on the same chromosome
// // 				if (startlist[i]->endpoint < endlist[j]->startpoint)
// // 				{
// // 					//i is start and j is end
// // 					distance_i_j = compute_distance(startlist[i], endlist[j]);
// // 				}
// // 				else if (startlist[i]->startpoint >endlist[j]->endpoint)
// // 				{
// // 					//i is end and j is start
// // 					distance_i_j = compute_distance(endlist[j], startlist[i]);
// // 				}
// // 				else
// // 				{
// // 					//overlap
// // 					distance_i_j = 0;
// // 				}
// // 			}
// // 
// // 			if (distance_i_j < best_distance)
// // 			{
// // 				best_i = i;
// // 				best_j = j;
// // 			}
// 
// 			if (strcmp(startlist[i]->chromosome, endlist[j]->chromosome) != 0)
// 			{
// 				//not on the same chromosome
// 				outputPEreads(startlist[i], endlist[j]);
// 			} 
// 			else
// 			{
// 				//on the same chromosome
// 				if (startlist[i]->startpoint <= endlist[j]->startpoint)
// 				{
// 					//i is start and j is end
// 					outputPEreads(startlist[i], endlist[j]);
// 				}
// 				else
// 				{
// 					//i is end and j is start
// 					outputPEreads(endlist[j], startlist[i]);
// 				}
// 			}
// 		}		
// 	}
// 
// 	//output reads best_i and best_j
// 	//outputPEreads(startlist[best_i], endlist[best_j])
// 
// 	return;
// }

// void outputPEfile()
// {
// 	//multiple alignment, keep alignments on the same chr
// 
// 	int i, j, ii, jj;
// 
// 	if (startlist_count * endlist_count <= MAX_COMB_LIMIT)
// 	{
// 		for (i = 1; i <= startlist_count; i++)
// 		{
// 			for (j = 1; j <= endlist_count; j++)
// 			{
// 
// 				if (strcmp(startlist[i]->chromosome, endlist[j]->chromosome) != 0)
// 				{
// 					//not on the same chromosome
// 					//outputPEreads(startlist[i], endlist[j]);
// 					if (startlist[i]->chrNum < endlist[j]->chrNum)
// 					{
// 						outputPEreads(startlist[i], endlist[j]);
// 					} 
// 					else
// 					{
// 						outputPEreads(endlist[j], startlist[i]);
// 					}
// 				} 
// 				else
// 				{
// 					//on the same chromosome
// 					if (startlist[i]->endpoint <= endlist[j]->startpoint && endlist[j]->startpoint - startlist[i]->endpoint >= 50000)
// 					{
// 						//i is start and j is end
// 						outputPEreads(startlist[i], endlist[j]);
// 					}
// 					else if (endlist[j]->endpoint <= startlist[i]->startpoint && startlist[i]->startpoint - endlist[j]->endpoint >= 50000)
// 					{
// 						//i is end and j is start
// 						outputPEreads(endlist[j], startlist[i]);
// 					}
// 				}
// 			}
// 		}
// 
// 		return;
// 	} 
// 	else
// 	{
// 		for (i = 1; i <= startlist_count; i++)
// 		{
// 			for (j = 1; j <= endlist_count; j++)
// 			{
// 				if (strcmp(startlist[i]->chromosome, endlist[j]->chromosome) != 0)
// 				{
// 					//not on the same chromosome
// // 					for (ii = 1; ii <= startlist_count; ii++)
// // 					{
// // 						dropped_reads << startlist[ii]->name << "\t1\t" << startlist[ii]->chromosome << '\t' << startlist[ii]->startpoint << '\t' << startlist[ii]->end << endl;
// // 					}
// // 					for (jj = 1; jj <= endlist_count; jj++)
// // 					{
// // 						dropped_reads << endlist[jj]->name << "\t2\t" << endlist[jj]->chromosome << '\t' << endlist[jj]->startpoint << '\t' << endlist[jj]->end << endl;
// // 					}
// 				} 
// 				else
// 				{
// 					//on the same chromosome
// 					if (startlist[i]->endpoint <= endlist[j]->startpoint && endlist[j]->startpoint - startlist[i]->endpoint >= 50000)
// 					{
// 						//i is start and j is end
// 						outputPEreads(startlist[i], endlist[j]);
// 					}
// 					else if (endlist[j]->endpoint <= startlist[i]->startpoint && startlist[i]->startpoint - endlist[j]->endpoint >= 50000)
// 					{
// 						//i is end and j is start
// 						outputPEreads(endlist[j], startlist[i]);
// 					}
// 				}
// 			}
// 		}
// 	}
// 
// 	return;
// }

void outputPEfile()
{
	/************************************************************************/
	/* Specifically for fusion clusters!!!                                  */
	/************************************************************************/
	int i, j;
	bool interesting = true;

	for (i = 1; i <= startlist_count; i++)
	{
		for (j = 1; j <= endlist_count; j++)
		{
			if (startlist[i]->chromosome.compare(endlist[j]->chromosome) == 0 && startlist[i]->strand.compare(endlist[j]->strand) != 0 && (startlist[i]->endpoint <= endlist[j]->startpoint && endlist[j]->startpoint - startlist[i]->endpoint < FUSION_CUTOFF || endlist[j]->endpoint <= startlist[i]->startpoint && startlist[i]->startpoint - endlist[j]->endpoint < FUSION_CUTOFF))
			{
				interesting = false;
			} 
		}
	}

	if (interesting == true && startlist_count * endlist_count == 1)
	{
		for (i = 1; i <= startlist_count; i++)
		{
			for (j = 1; j <= endlist_count; j++)
			{
				if (startlist[i]->chromosome.compare(endlist[j]->chromosome) != 0)
				{
					//not on the same chromosome
					//outputPEreads(startlist[i], endlist[j]);
//					if (startlist[i]->chrNum < endlist[j]->chrNum)
//					{
						outputPEreads(startlist[i], endlist[j]);
//					} 
//					else
//					{
//						outputPEreads(endlist[j], startlist[i]);
//					}
				} 
				else
				{
					//on the same chromosome
//					if (startlist[i]->endpoint <= endlist[j]->startpoint && endlist[j]->startpoint - startlist[i]->endpoint >= 50000)
//					{
						//i is start and j is end
						outputPEreads(startlist[i], endlist[j]);
//					}
//					else if (endlist[j]->endpoint <= startlist[i]->startpoint && startlist[i]->startpoint - endlist[j]->endpoint >= 50000)
//					{
//						//i is end and j is start
//						outputPEreads(endlist[j], startlist[i]);
//					}
				}
			}
		}

		return;
	} 

	return;
}

void cleanup()
{
	reads* ptr;

	for (int i = 1; i <= startlist_count; i++)
	{
		ptr = startlist[i];
		delete ptr;
	}
	for (int j = 1; j <= endlist_count; j++)
	{
		ptr = endlist[j];
		delete ptr;
	}

	startlist_count = 0;
	endlist_count = 0;

	return;
}


void parse(string inputfilename)
{
	string outpurfilename;
	outpurfilename = dirPrefix + "data/parsedPER/ChromosomeName.txt";
	outputfile_chrname.open(outpurfilename.c_str());
//	sprintf(outpurfilename, "%sdata/parsedPER/DroppedReads.txt", dirPrefix);
//	dropped_reads.open(outpurfilename);

	reads *ptr;
	ptr = NULL;

	startlist_count = 0;
	endlist_count = 0;

	ifstream inputfile;
	inputfile.open(inputfilename.c_str());
		
	string name, prevName = "xxx";
	int start_end_switch;
	unsigned long tmp;
	bool writeswitch = false;
	string info;

	char field2[10000], field3[10000], field4[10000], field5[10000], field6[10000], field7[10000], field8[10000];
	char binaryflagbits[500];
	long flag;

	while (inputfile >> name)
	{
		//find start_end_switch
		if (name[name.size()-1] == '2')
			start_end_switch = 2;
		else
			start_end_switch = 1;
		
		name.erase(name.size()-2, 2);

		inputfile >> flag;
		inputfile >> field2;
		inputfile >> field3;
		inputfile >> field4;
		inputfile >> field5;
		inputfile >> field6;
		inputfile >> field7;
		inputfile >> field8;
		getline(inputfile, info);

		dec2bin(flag, binaryflagbits);

		if (prevName.compare(name) != 0)
		{
			if (startlist_count * endlist_count > 0)
			{
				outputPEfile();
			}
			prevName = name;
			cleanup();
			writeswitch = false;
		}

		if (start_end_switch == 1)
		{
			//start point
			if (field6[0] == '*' || field6[0] == '=')
			{
				//not fusion read				
				startlist_count++;
				ptr = new reads;

				ptr->name = name;

				if (binaryflagbits[4] == '0')
					ptr->strand = "+";
				else
					ptr->strand = "-";
				ptr->chromosome = field2;
				ptr->startpoint = strtol(field3, NULL, 10);
				ptr->end = field5;

				ptr->compute_endpoint();

				startlist[startlist_count] = ptr;
			} 
			else
			{
				//fusion read, discard
			}
		} 
		else if (start_end_switch == 2)
		{
			//end point
			if (field6[0] == '*' || field6[0] == '=')
			{
				//not fusion read				
				endlist_count++;
				ptr = new reads;

				ptr->name = name;
				
				if (binaryflagbits[4] == '0')
					ptr->strand = "+";
				else
					ptr->strand = "-";
				ptr->chromosome = field2;
				ptr->startpoint = strtol(field3, NULL, 10);
				ptr->end = field5;

				ptr->compute_endpoint();

				endlist[endlist_count] = ptr;
			} 
			else
			{
				//fusion read, discard
			}

			writeswitch = true;
		}
	}

	if (name.empty() == false)
	{
		if (startlist_count * endlist_count > 0)
		{
			outputPEfile();
		}
		cleanup();
	}


	inputfile.close();
	outputfile_chrname.close();
// 	dropped_reads.close();

	for (tmp = 1; tmp <= targetChrNum; tmp++)
	{
		outputfile[tmp].close();
	}

	return;
}



int main(int argc, char* argv[])
{
	targetChrNum = 0;

	if (argc != 3)
	{
		cout << argv[0] << "\t<filename>\t<targetPath>" << endl;
		return 1;
	}

	string inputfilename;

	inputfilename = argv[1];
	dirPrefix = argv[2];
//	allDistanceToFile.open("C:\\Users\\yin\\Desktop\\PEmatch\\human35\\allDistance.txt");

	parse(inputfilename);

//	allDistanceToFile.close();

	return 0;
}


