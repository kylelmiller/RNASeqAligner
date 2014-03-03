/*    
 *    SetUnmappedBitFlag.cpp	
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

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <stdlib.h>

using namespace std;

#define IS_PAIRED 0x0001
#define IS_PAIRED_MAPPED 0x0002
#define IS_UNMAPPED 0x0004
#define MATE_UNMAPPED 0x0008
#define IS_REVERSE 0x0010
#define IS_MATE_REVERSE 0x0020
#define IS_FIRST_END 0x040
#define IS_SECOND_END 0x0080
#define IS_PRIMARY 0x0100
#define IS_FAILED_QUAL_CHECK 0x0200
#define IS_PCR_DUP 0x0400

void
RemoveDupMapreads(const char* infile, const char* outfile)
{
	ifstream ifs(infile);

	ofstream ofs(outfile);

	string prev_tagname = "";

	map<int, string> mapped_reads;

	size_t count = 0;

	if (ifs.is_open())
	{
		string line;

		while (getline(ifs,line))
		{
			if (line == "")
				continue;

			char tagname[1000];

			unsigned short bitflag;

			//TRAN00000027662:59:252  16      chr14   56062712        255     8M157375N92M    *       0       0       TATTATTTTCCGCTTTCCCTGGGCTTACAGAGAATCCTTGCCCTTCTTGTACTGTGTCACTTTATGGGGTTGGTGCTTGCCACACTTCTTACAGAAAGTC 
			//#######$#%'*,*++*,/121222012122233455766655666555666677989:;<<=>>>>>==>??>>>>>>>>>>>>>>>>>>>=>=>>>>> NM:i:5  8:G>T,24:C>T,33:T>A,63:G>A,81:T>A

			size_t read_count = sscanf(line.c_str(), "%s\t%hu", tagname, &bitflag);

			string tagbasename = tagname;

			int tagid = tagbasename[tagbasename.length() - 1] - '0';

			if (tagid != 1 && tagid != 2)
				cerr << "unknown tagid:"<< tagid << endl << line << endl;

			tagbasename = tagbasename.substr(0, tagbasename.length() - 2);

			if (prev_tagname.empty() || prev_tagname == tagbasename)
			{
				mapped_reads[tagid] = line;

				if (prev_tagname.empty())
					prev_tagname = tagbasename;
			}
			else
			{
				map<int, string>::iterator msi_iter;

				for (msi_iter = mapped_reads.begin(); msi_iter != mapped_reads.end(); ++msi_iter)
				{
					string& curline = msi_iter->second;

					char curtagname[1000];

					unsigned short curbitflag;

					size_t read_count = sscanf(curline.c_str(), "%s\t%hu", curtagname, &curbitflag);

					if (mapped_reads.size() == 2)
					{
						if (msi_iter->first == 1)
						{
							curbitflag = curbitflag | IS_PAIRED;

							curbitflag = curbitflag | MATE_UNMAPPED;

							curbitflag = curbitflag | IS_FIRST_END;
						}
						else if (msi_iter->first == 2)
						{
							curbitflag = curbitflag | IS_PAIRED;

							curbitflag = curbitflag | MATE_UNMAPPED;

							curbitflag = curbitflag | IS_SECOND_END;
						}
						else
						{
							cerr << "unknown tagid:"<< tagid << endl << curline << endl;
						}

						ofs << curtagname <<'\t' << curbitflag<<'\t' << curline.substr(curline.find("*	0	0"))<<endl;

						++count;
					}
					else if (mapped_reads.size() == 1)
					{
						if (msi_iter->first == 1)
						{
							curbitflag = curbitflag | IS_PAIRED;

							curbitflag = curbitflag | IS_FIRST_END;
						}
						else if (msi_iter->first == 2)
						{
							curbitflag = curbitflag | IS_PAIRED;

							curbitflag = curbitflag | IS_SECOND_END;
						}
						else
						{
							cerr << "unknown tagid:"<< tagid << endl<< curline << endl;
						}

						ofs << curtagname <<'\t' << curbitflag<<'\t' << curline.substr(curline.find("*	0	0"))<<endl;

						++count;

					}
					else
					{
						cerr <<"unmapped both ends size:" <<mapped_reads.size()<<endl;
					}
				}

				mapped_reads.clear();

				mapped_reads[tagid] = line;

				prev_tagname = tagbasename;
			}			
		}
		ifs.close();

	}
	else
	{
		cout << "can't open file "<< infile<<endl; exit(1);
	}

	map<int, string>::iterator msi_iter;

	for (msi_iter = mapped_reads.begin(); msi_iter != mapped_reads.end(); ++msi_iter)
	{
		string& curline = msi_iter->second;

		char curtagname[1000];

		unsigned short curbitflag;

		size_t read_count = sscanf(curline.c_str(), "%s\t%hu", curtagname, &curbitflag);

		if (mapped_reads.size() == 2)
		{
			if (msi_iter->first == 1)
			{
				curbitflag = curbitflag | IS_PAIRED;

				curbitflag = curbitflag | MATE_UNMAPPED;

				curbitflag = curbitflag | IS_FIRST_END;
			}
			else if (msi_iter->first == 2)
			{
				curbitflag = curbitflag | IS_PAIRED;

				curbitflag = curbitflag | MATE_UNMAPPED;

				curbitflag = curbitflag | IS_SECOND_END;
			}
			else
			{
				cerr << "unknown tagid:"<< curbitflag << endl<< curline << endl;
			}

			ofs << curtagname <<'\t' << curbitflag<<'\t' << curline.substr(curline.find("*	0	0"))<<endl;

			++count;
		}
		else if (mapped_reads.size() == 1)
		{
			if (msi_iter->first == 1)
			{
				curbitflag = curbitflag | IS_PAIRED;

				curbitflag = curbitflag | IS_FIRST_END;
			}
			else if (msi_iter->first == 2)
			{
				curbitflag = curbitflag | IS_PAIRED;

				curbitflag = curbitflag | IS_SECOND_END;
			}
			else
			{
				cerr << "unknown tagid:"<< curbitflag << endl<< curline << endl;
			}

			ofs << curtagname <<'\t' << curbitflag<<'\t' << curline.substr(curline.find("*	0	0"))<<endl;

			++count;

		}
		else
		{
			cerr <<"unmapped both ends size:" <<mapped_reads.size()<<endl;
		}
	}

	cout << "unmapped_read\t"<<count << endl;
}

int main(int argc, char** argv)
{
	if (argc < 3)
	{
		cout << "infile outfile" <<endl;
		exit(0);
	}
	const char* infile = argv[1];

	const char* outfile = argv[2];

	RemoveDupMapreads(infile, outfile);
}
