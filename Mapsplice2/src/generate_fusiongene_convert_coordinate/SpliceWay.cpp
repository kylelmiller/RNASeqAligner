#include "SpliceWay.h"


SpliceWay::SpliceWay() : chrom_name_ptr(NULL), start(0), spliceway_vec_ptr(NULL)
{	
}

SpliceWay::SpliceWay(string* cnp, size_t st, vector<pair<size_t, int> >* svp) : chrom_name_ptr(cnp), start(st), spliceway_vec_ptr(svp)
{
}

SpliceWayTrue::SpliceWayTrue() : start(0)
{
}

SpliceWayTrue::SpliceWayTrue(const SpliceWay& sw) : chrom_name(*(sw.chrom_name_ptr)), start(sw.start), spliceway_vec(*(sw.spliceway_vec_ptr))
{
	++spliceway_true_count;

	*spliceway_true_ofs << chrom_name << '\t' << start << '\t';

	for (size_t i = 0; i < spliceway_vec.size(); ++i)
	{
		*spliceway_true_ofs <<spliceway_vec[i].first <<'\t' << spliceway_vec[i].second << '\t';
	}

	*spliceway_true_ofs << endl;
}

SpliceWayTrue::SpliceWayTrue(string* cnp, size_t st, vector<pair<size_t, int> >* svp) : chrom_name(*cnp), start(st), spliceway_vec(*svp)
{
	++spliceway_true_count;

	*spliceway_true_ofs << chrom_name << '\t' << start << '\t';

	for (size_t i = 0; i < spliceway_vec.size(); ++i)
	{
		*spliceway_true_ofs <<spliceway_vec[i].first <<'\t' << spliceway_vec[i].second << '\t';
	}

	*spliceway_true_ofs << endl;
}

