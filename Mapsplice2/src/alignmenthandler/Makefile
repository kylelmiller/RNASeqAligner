All:
	g++ -O3 -o sam2junc sam2junc.cpp sharedlib.cpp JunctionHandler.cpp SamRec.cpp JunctionSeed.cpp 2>log
	g++ -O3 -o test_handler test_handler.cpp JunctionHandler.cpp SamRec.cpp sharedlib.cpp AlignmentHandler.cpp JunctionSeed.cpp FusionSamRec.cpp UnionExpressedRegions.cpp disjointset.cpp 2>log
	g++ -O3 -o merge_junction merge_junction.cpp JunctionHandler.cpp SamRec.cpp sharedlib.cpp JunctionSeed.cpp 2>log
	g++ -O3 -o comp_junc comp_junc.cpp JunctionHandler.cpp SamRec.cpp sharedlib.cpp JunctionSeed.cpp 2>log
	g++ -O3 -o junc_db junc_db.cpp JunctionHandler.cpp SamRec.cpp sharedlib.cpp JunctionSeed.cpp 2>log
	g++ -O3 -o SepSam SepSam.cpp 2>log
