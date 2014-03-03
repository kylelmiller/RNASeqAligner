//#include "AlignmentHandler.h"
#include "JunctionHandler.h"

int main(int argc, char** argv)
{
	if (argc < 9)
		cout << "converted_junction max_read_len max_insert chrom_dir sam_file1 sam_file2 ..."<<endl;

	char* converted_junction = argv[1];

	size_t max_read_len = (size_t) atoi(argv[2]);

	int max_insert = atoi(argv[3]);

	int max_deletion = atoi(argv[4]);

	char* chrom_dir = argv[5];

	size_t min_anchor = (size_t) atoi(argv[6]);

	size_t min_junc_anchor = (size_t) atoi(argv[7]);

	size_t paired = (size_t) atoi(argv[8]);

	size_t do_filter = (size_t) atoi(argv[9]);

	size_t min_mismatch = (size_t) atoi(argv[10]);

	string converted_junction_ins = converted_junction; converted_junction_ins.append(".ins");

	string converted_junction_del_str = converted_junction; converted_junction_del_str.append(".del");

	string converted_junction_fusion_str = converted_junction; converted_junction_fusion_str.append(".fusion");

	string filtered_junction = converted_junction; filtered_junction.append(".filtered");

	vector<string> sam_files;

	for (int i = 11; i < argc; ++i)
		sam_files.push_back(argv[i]);

	JunctionHandler junc_handler;

	junc_handler.Init(sam_files, max_read_len, max_insert, max_deletion, chrom_dir, min_anchor, min_mismatch, min_junc_anchor);

	cout << "read sam "<<endl;

	junc_handler.ReadSam();

	if (do_filter)
		junc_handler.MarkFiltered(paired > 0);

	cout << "write junction"<<endl;

	junc_handler.CollectStats();

	string junc_stat = converted_junction; 

	junc_stat.append(".stat");

	junc_handler.WriteStats(junc_stat, "final junction stats");

	junc_handler.WriteJunction(converted_junction, converted_junction_ins.c_str(), converted_junction_del_str.c_str(), converted_junction_fusion_str.c_str(), filtered_junction.c_str());

	
}