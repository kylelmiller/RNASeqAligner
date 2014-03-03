#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
using namespace std;



int main(int argc, char** argv)
{
	ifstream input_fs(argv[1]);
	ofstream output_fs(argv[2]);
	if(!input_fs)
	{
			fprintf(stderr,"error: open input file error\n");
			exit(1);
	}
	if(!output_fs)
	{
			fprintf(stderr,"error: open output file error\n");
			exit(1);
	}
	string line;
	while(getline(input_fs, line))
	{
		size_t tab_pos = line.find('\t');
		if(tab_pos > 2 && line[tab_pos - 2] == '/' && (line[tab_pos - 1] == '1' || line[tab_pos - 1] == '2') )
		{
			output_fs << line.substr(0, tab_pos - 2) << line.substr(tab_pos) << endl;
			continue;
		}
		output_fs << line << endl;
	}
	input_fs.close();
	output_fs.close();
	return 0;
}