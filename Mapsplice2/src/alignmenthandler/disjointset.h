#ifndef DISJOINTSET_H
#define DISJOINTSET_H

#include "sharedlib.h"

struct DSet{
	int p;
	int rank;

	DSet(int parent, int rk);

};

class DisjointSet{

private:
	//vector<DSet> m_disjoint_set;
	vector<size_t> p;

	vector<int> rank;	

public:

	DisjointSet();

	DisjointSet(size_t size);

	void Resize(size_t size);

	void MakeSet(size_t x);

	void Union(size_t x, size_t y);

	void Link(size_t x, size_t y);

	size_t FindSet(size_t x);

};

#endif