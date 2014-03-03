#ifndef SPLICEWAY_H
#define SPLICEWAY_H

#include "sharedlib.h"
//chr1    11420   133803  JUNC_1  3       -       11420   133803  255,0,0 2       10,66,  0,122450,       0       6       CTAC    0.99965 0.44492 2       2       2

struct SpliceWay {
	string* chrom_name_ptr;

	size_t start;
	
	vector<pair<size_t, int> >* spliceway_vec_ptr;

	SpliceWay();

	SpliceWay(string* cnp, size_t st, vector<pair<size_t, int> >* svp);
};

struct SpliceWayTrue {
	string chrom_name;

	size_t start;
	
	vector<pair<size_t, int> > spliceway_vec;

	SpliceWayTrue();

	SpliceWayTrue(const SpliceWay& sw);

	SpliceWayTrue(string* cnp, size_t st, vector<pair<size_t, int> >* svp);
};

#endif