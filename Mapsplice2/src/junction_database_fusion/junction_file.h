/*    
 *    junction_file.h		
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

#ifndef JUNCTION_FILE_H
#define JUNCTION_FILE_H

#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <map>
#include "junction_data_struct.h"

using namespace std;

class Junction_File
{
private:
	ifstream junc_fs;
	string line;
	int current_junc_id;

public:
	Junction_File(char* filename, bool has_header)
	{
		junc_fs.open(filename);
		if( !junc_fs ) 
		{
			fprintf(stderr,"error: open junction file error\n");
			exit(1);
		}  
		if(has_header)
			getline(junc_fs, line);// skip header
		current_junc_id = 0;
	}

	Junction_File(char* filename, bool has_header, int junc_id_start)
	{
		junc_fs.open(filename);
		if( !junc_fs ) 
		{
			fprintf(stderr,"error: open junction file error\n");
			exit(1);
		}  
		if(has_header)
			getline(junc_fs, line);// skip header
		current_junc_id = junc_id_start;
	}

	~Junction_File()
	{

	}

	int load_all_normal_junctions(vector<Junction>& new_junction_set)    // load all the rest of junctions
	{
		while(getline(junc_fs, line))
		{
			if(line=="")
				break;
			Junction new_junc;
			new_junc.setValue(line);
			new_junc.junc_name = current_junc_id;
			current_junc_id ++;
			new_junction_set.push_back(new_junc);
		}
		return current_junc_id;
	}

	int load_all_fusion_junctions(vector<Junction>& new_junction_set)    // load all the rest of junctions
	{
		while(getline(junc_fs, line))
		{
			if(line=="")
				break;
			Junction new_junc_start;   
			Junction new_junc_end;  

			new_junc_start.setValue_fusion(line, START);
			new_junc_start.junc_name = current_junc_id;
			new_junction_set.push_back(new_junc_start);

			new_junc_end.setValue_fusion(line, END);
			new_junc_end.junc_name = current_junc_id;
			new_junction_set.push_back(new_junc_end);
			current_junc_id ++;
		}
		return current_junc_id;
	}

	void close()
	{
		junc_fs.close();
	}
};

#endif


