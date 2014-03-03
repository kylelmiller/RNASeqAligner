//#include "AlignmentHandler.h"
#include "JunctionHandler.h"

int main(int argc, char** argv)
{
	if (argc < 3)
		cout << "merged_junction junction_file1 junction_file2 junction_file3 ..."<<endl;

	char* merged_junction = argv[1];

	string merged_junc_str = merged_junction;

	string merged_junc_ins_str = merged_junction; merged_junc_ins_str.append(".ins");

	string merged_junc_del_str = merged_junction; merged_junc_del_str.append(".del");

	string merged_junc_fusion_str = merged_junction; merged_junc_fusion_str.append(".fusion");

	string merged_junc_filtered_str = merged_junction; merged_junc_filtered_str.append(".filtered");

	vector<string> junction_files;
	for (int i = 2; i < argc; ++i)
		junction_files.push_back(argv[i]);

	JunctionHandler junc_handler;

	vector<string> sam_files;

	//(sam_files, max_read_len, max_insert, max_deletion, chrom_dir, min_anchor, min_junc_anchor);

	junc_handler.Init(sam_files, 1000, 3, 10, "", 1, 1, 1);

	cout << "read junction "<<endl;
	cout << junc_handler.ReadJunction(junction_files)<<endl;

	cout << "write junction"<<endl;
	junc_handler.WriteJunction(merged_junc_str.c_str(), merged_junc_ins_str.c_str(), merged_junc_del_str.c_str(), merged_junc_fusion_str.c_str(), merged_junc_filtered_str.c_str());
}