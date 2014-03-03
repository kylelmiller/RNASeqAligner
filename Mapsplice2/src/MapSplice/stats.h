#ifndef STATS_H
#define STATS_H

#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>

class Mapping_Stats
{
	public:
		Mapping_Stats()
		{
			read_processed = 0;
			read_aligned = 0;
			read_spliced = 0;
			read_fusion = 0;
		}
		
		int read_processed;
		int read_aligned;
		int read_spliced;
		int read_fusion;
};

#endif


