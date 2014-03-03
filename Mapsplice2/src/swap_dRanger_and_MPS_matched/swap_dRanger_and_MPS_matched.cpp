#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
using namespace std;

int main(int argc, char** argv)
{
	if (argc < 3)
	{
		cout << "in_matched_fusions output_swapped_fusions"<<endl;
		exit(0);
	}

	ifstream ifs(argv[1]);

	ofstream ofs(argv[2]);

	if (ifs.is_open())
	{
		string line;

		while (getline(ifs,line))
		{

			if (line.empty())
				continue;

			string dRanger_line, mps_line;

			dRanger_line = line.substr(0, line.find("	----	----	"));

			mps_line = line.substr(line.find("	----	----	") + 11);

			while (mps_line[mps_line.length() - 1] == ' ' || mps_line[mps_line.length() - 1] == '\t')
				mps_line = mps_line.substr(0, mps_line.length() - 1);

			ofs << mps_line << "	----	----	"<<dRanger_line << endl;
		}

		ifs.close();
	}
	else
	{
		cout << "can't open file "<< argv[1]<<endl; exit(1);
	}
}