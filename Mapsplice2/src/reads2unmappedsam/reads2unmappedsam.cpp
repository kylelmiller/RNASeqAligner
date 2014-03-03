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
		cout << "sam_file flag fq_file1 ..."<<endl;
		fprintf(stderr,"error: too few arguments\n");exit(1);
	}

	string sam_file = argv[1];

	int flag = atoi(argv[2]);

	ofstream ofs(sam_file.c_str());

	size_t max_len = 0, min_len = -1, count_line = 0, count_total = 0;

	for (size_t i = 3; i < argc; ++i)
	{
		string fq_file = argv[i];

		ifstream ifs(fq_file.c_str());

		if (ifs.is_open())
		{
			while (!ifs.eof() )
			{
				string line1, line2, line3, line4;

				getline(ifs, line1);

				getline(ifs, line2);

				if (line1 == "")
					continue;

				ofs << line1.c_str() + 1 << "\t4\t*\t0\t0\t*\t*\t0\t0\t"<<line2 << '\t';

				if (flag)
				{
					getline(ifs, line3);

					getline(ifs, line4);

					ofs << line4 <<"\tIH:i:0\tHI:i:0" << endl;
				}
				else
				{
					string Is(line2.length(), 'I');

					ofs << Is <<"\tIH:i:0\tHI:i:0" <<endl;
				}

				if (line2.length() < min_len)
					min_len = line2.length();

				if (line2.length() > max_len)
					max_len = line2.length();

				++count_total;
			}
		}
	}

	cout << "min_len: "<< min_len << endl << "max_len: "<< max_len << endl << "read_count: "<<count_line <<endl << "total_count: "<< count_total << endl;
}
