/*    
 *    filter_1hits.cpp		
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

#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
#include <string.h>
using namespace std;

int main(int argc, char* argv[])
{
	ifstream ifs(argv[1]);

	cout <<argv[1]<<endl;

	unsigned short min_mismatch = (unsigned short) atoi(argv[2]);

	//double min_ave_mis = (double) atof(argv[3]);

	double min_lpq = (double) atof(argv[3]);

	//double min_entropy = (double) atof(argv[5]);

	string canon_junc = argv[4]; ofstream ofs_canon_junc(canon_junc.c_str());

	string noncanon_junc = argv[5]; ofstream ofs_noncanon_junc(noncanon_junc.c_str());

	size_t canon_count = 0, non_canon_count = 0, entropy_mis = 0, del_entropy = 0, pcr = 10, mismatch_count = 0, noncanon_mismatch_count = 0, non_canon_entropy = 0;

	if (ifs.is_open())
	{
		string firstline;
		getline(ifs, firstline);

		ofs_canon_junc <<firstline<<endl;
		ofs_noncanon_junc <<firstline<<endl;

		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);
			if (line == "")
				continue;

			char chromname[100], juncname[100], strand[10], rgb[100], blocks[1000], blocksoffset[1000], il[100], flankchr[100];
			int juncst, juncend, prefixend, suffixst, kinds, hits, flankcase;

			char fusionseq[1000];

			double lpq, entropy;

			unsigned short min_mis, max_mis;

			double ave_mis;

			//chr10~chr10	100150	100170	JUNC_6	2	-+	255,0,0	2	31,20	0,20	0	3	CTGC	:I?CI;DI/IIG6'0,3293?04%+(0(5*)+'%(&6&$*%)%%&&&&$%&	4	4	4	0	20	11
			//chr10~chr10	176036	94374495	JUNC_9	1	-+	255,0,0	2	25,26	0,94198459	0	0	ATTT	I+CGI'*<'3:':5-33:1.'('&&*$&+&&'%%%%%'&##%%$#&$&%#%	2	2	2	25	0	1

			sscanf(line.c_str(), "%s\t%d\t%d\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%lf\t%d\t%s\t%s\t%hu\t%hu\t%lf", chromname,  &prefixend, &suffixst, 
				juncname, &hits, strand, rgb, &kinds, blocks, blocksoffset, &entropy, &flankcase, flankchr, fusionseq, &min_mis, &max_mis, &ave_mis);
			
			if (/*lpq < min_lpq ||*/ min_mis >= min_mismatch)// || ave_mis >= min_ave_mis || entropy < min_entropy)
			{
				++non_canon_count;
				++mismatch_count;
				ofs_noncanon_junc << line << "\tmin_mis:"<<min_mis << endl;				
			}
			else if (min_mis >= 1 && flankcase == 0 && entropy < 1)
			{
				++non_canon_count;
				++noncanon_mismatch_count;
				ofs_noncanon_junc << line <<"\tflankcase:"<<flankcase<<"\tmin_mis:"<<min_mis <<"\tentropy:"<<entropy <<  endl;
			}
			else if (flankcase == 0 && entropy == 0)
			{
				++non_canon_count;
				++non_canon_entropy;
				ofs_noncanon_junc << line <<"\tflankcase:"<<flankcase <<"\tentropy:"<<entropy <<  endl;
			}
			//else if (min_mis >= 1 && entropy <= 0.00001)
			//{
			//	++non_canon_count;
			//	++entropy_mis;
			//	ofs_noncanon_junc << line <<"\tmin_mis:"<<min_mis <<"\tentropy:"<<entropy <<  endl;
			//}
			//else if (min_mis > 0 && entropy == 0)
			//{
			//	++entropy_mis;
			//	ofs_noncanon_junc << line << endl;
			//}
			//else if (flankcase == 0 && (suffixst - prefixend) < 70 && entropy < 1)
			//{
			//	++del_entropy;
			//	ofs_noncanon_junc << line << endl;
			//}
			//else if (hits > 10 && entropy == 0)
			//{
			//	++pcr;
			//	ofs_noncanon_junc << line << endl;
			//}
			else
			{
				++canon_count;
				ofs_canon_junc << line << endl;

				//unsigned short prefix_len, suffix_len;

				//sscanf(blocks, "%hu,%hu", &prefix_len, &suffix_len);

				//if (prefix_len + suffix_len < read_len || prefix_len < min_anchor || suffix_len < min_anchor)
				//{
				//	
				//}
				//else
				//{
				//	++canon_count;
				//	ofs_canon_junc << line << endl;
				//}
			}

			//if (hits > 1)
			//{
			//	++canon_count;
			//	ofs_canon_junc << line << endl;
			//}
			//else
			//{
			//	unsigned short prefix_len, suffix_len;

			//	sscanf(blocks, "%hu,%hu", &prefix_len, &suffix_len);

			//	if (prefix_len + suffix_len < read_len || prefix_len < min_anchor || suffix_len < min_anchor)
			//	{
			//		++non_canon_count;
			//		ofs_noncanon_junc << line << endl;
			//	}
			//	else
			//	{
			//		++canon_count;
			//		ofs_canon_junc << line << endl;
			//	}
			//}
		}
	}

	cout <<"in filtered 1hits junction: "<<canon_count<<endl;

	cout <<"not in filtered 1hits junction: "<<non_canon_count<<endl;

	cout <<"mismatch_count: "<<mismatch_count<<endl;

	cout <<"noncanon_mismatch_count: "<<noncanon_mismatch_count<<endl;

	cout <<"entropy_mis: "<<entropy_mis<<endl;

	cout <<"non_canon_entropy: "<<non_canon_entropy<<endl;
	
	if (canon_count == 0)
		return 100;
}

