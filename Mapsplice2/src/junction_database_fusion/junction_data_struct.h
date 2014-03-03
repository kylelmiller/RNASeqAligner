/*    
 *    junction_data_struct.h		
 *    MapSplice
 *
 *    Copyright (C) 2010 University of Kentucky and
 *                       Zeng Zheng
 *
 *    Authors: Zeng Zheng
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

#ifndef JUNCTION_H
#define JUNCTION_H

#include <string>

using namespace std;

#pragma warning(disable:4996)

enum Junc_Type{NORMAL, FUSION};
enum Fusion_Type{START, END};


struct Junction
{
public:
	Junction() 
	{}

	~Junction()
	{}

	void setValue(string line)
	{
		clear();
		char chrom_tmp[1000];
		sscanf(line.c_str(), "%s\t%d\t%d",chrom_tmp, &start, &end);        
		chrom=chrom_tmp;
		junc_type=NORMAL;
	}

	void setValue_fusion(string line, Fusion_Type f_type)
	{
		clear();
		char chrom_tmp[1000], strand_tmp[1000];
		if(f_type==START)
		{
			sscanf(line.c_str(), "%[^~]%*[~]%*s\t%d\t%*d\t%*s%*d%1s", chrom_tmp, &start, strand_tmp);
			end=start;
		}
		else
		{
			sscanf(line.c_str(), "%*[^~]%*[~]%s\t%*d\t%d\t%*s%*d%*1s%1s", chrom_tmp, &end, strand_tmp);
			start=end;
		}
		fusion_type=f_type;
		chrom=chrom_tmp;
		strand=strand_tmp;
		junc_type=FUSION;
	}

	void clear() 
	{ 
		chrom.clear();
		strand.clear();
	}

	string chrom;
	int start; 
	int end;
	Junc_Type junc_type;
	Fusion_Type fusion_type;
	string strand;
	int junc_name;
};


#endif

