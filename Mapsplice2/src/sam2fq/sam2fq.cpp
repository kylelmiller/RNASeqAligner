/*    
 *    sam2fq.cpp		
 *    MapSplice
 *
 *    Copyright (C) 2010 University of Kentucky and
 *                       Kai Wang
 *
 *    Authors: Kai Wang
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


#include <iostream>
#include <vector>

#include <string>

#include <fstream>
#include <sstream>

#include <algorithm>
#include <dirent.h>
#include <iomanip>
#include <map>
#include <queue>
#include <list>

#include <cmath>
#include <errno.h>
#include <time.h>
#include <string.h>

using namespace std;


int main(int argc, char* argv[])
{
	if (argc < 4)
	{
		cout << "fq_file flag sam_file..."<<endl;
		fprintf(stderr,"error: too few arguments\n");exit(1);
	}

	string fq_file = argv[1];

	int flag = atoi(argv[2]);

	ofstream ofs(fq_file.c_str());

	size_t max_len = 0, min_len = -1, count_line = 0, count_total = 0;

	for (size_t i = 3; i < argc; ++i)
	{
		string samfile = argv[i];

		ifstream ifs(samfile.c_str());

		//VEC_COVER_BLOCK v_coverage_block;

		if (ifs.is_open())
		{
			while (!ifs.eof() )
			{
				string line;
				getline(ifs, line);
				if (line == "" || line[0] == '@')
					continue;

				char chromname[1000], readname[1000], chromseq[1000], qualseq[1000], spliceway[2000];
				char strand = '+';
				size_t /*prim, */ prefixst, strand_t, incorrect, mate_offest;

				int mate_diff;

				unsigned short mis_match;

				char mate_match;

				//61KTMAAXX:2:9:9288:4384#0	4	*	0	0	*	*	0	0	TGTTGTGCTTCTTGGCAAAGCGCATGTTCCTCAGGAACTTGGGCTTTACGAGGGCCTTGATAGCCTCGGCACGTG	ffffffffffffffffffffffffffffffffffedfffdfffdfddfedddddccdccddddfbbcaacbac`b	RG:Z:61KTM.2

				//%dM\t%dN\t%dM
				sscanf(line.c_str(), "%s\t%llu\t%s\t%llu\t%llu\t%s\t%c\t%llu\t%d\t%s\t%s", 
					readname, &strand_t, chromname, &prefixst, &incorrect, spliceway, &mate_match, &mate_offest, &mate_diff, chromseq, qualseq);

				if (flag == 0 || spliceway[0] == '*')
				{
					ofs << '@' << readname << endl << chromseq << endl << '+' << endl << qualseq << endl;

					size_t read_len = strlen(chromseq);

					if (read_len < min_len)
						min_len = read_len;

					if (read_len > max_len)
						max_len = read_len;

					++count_line;
				}

				++count_total;
			}
		}
	}

	cout << "min_len: "<< min_len << endl << "max_len: "<< max_len << endl << "read_count: "<<count_line <<endl << "total_count: "<< count_total << endl;
}
