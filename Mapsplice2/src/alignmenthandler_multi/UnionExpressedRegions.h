#ifndef UNIONEXPRESSEDREGIONS_H
#define UNIONEXPRESSEDREGIONS_H

#include "sharedlib.h"

struct UnionExpressedRegions {
	vector<size_t> m_unioned_regions_ids;

	vector<pair<size_t, size_t> > m_unioned_regions;

	size_t m_unioned_region_length;

	UnionExpressedRegions();
};

#endif