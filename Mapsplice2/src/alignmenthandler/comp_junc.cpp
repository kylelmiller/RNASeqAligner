//#include "AlignmentHandler.h"
#include "JunctionHandler.h"

int main(int argc, char** argv)
{

	vector<string> junction_file1, junction_file2;

	JunctionHandler junc_handler1, junc_handler2;

	junction_file1.push_back(argv[1]);

	junction_file2.push_back(argv[2]);

	junc_handler1.ReadJunction(junction_file1);

	junc_handler1.LoadJuncToSortVec();

	junc_handler1.GetMinimumExon();

	junc_handler2.ReadJunction(junction_file2);

	junc_handler2.LoadJuncToSortVec();

	junc_handler1.ConfirmJunction(junc_handler2, atoi(argv[3]), false);
	
}