#include "JunctionHandler.h"

int main(int argc, char** argv)
{
	if (argc < 6)
		cout << "min_anchor max_anchor max_repeat junc_file chrom_dir syn_seq"<<endl;

	vector<string> junction_file1;

	JunctionHandler junc_handler1;

	vector<string> dump_files;

	junc_handler1.Init(dump_files, 76, 3, 10, argv[5], 4, 0, 15);

	junction_file1.push_back(argv[4]);

	junc_handler1.ReadJunction(junction_file1);

	junc_handler1.LoadJuncToSortVec();

	junc_handler1.GetMinimumExon();

	junc_handler1.GenerateJunctionSequence(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), argv[6]);
	
}