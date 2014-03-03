#include "disjointset.h"

DSet::DSet(int parent, int rk) : p(parent), rank(rk)
{

}

DisjointSet::DisjointSet(size_t size)
{
	p.resize(size);

	rank.resize(size, 0);

	for (size_t i = 0; i < size; ++i)
		p[i] = i;
}

DisjointSet::DisjointSet()
{
}

void
DisjointSet::Resize(size_t size)
{
	p.resize(size);

	rank.resize(size, 0);

	for (size_t i = 0; i < size; ++i)
		p[i] = i;
}

void
DisjointSet::MakeSet(size_t x)
{
	p[x] = x;

	rank[x] = 0;
}

void
DisjointSet::Union(size_t x, size_t y)
{
	Link(FindSet(x), FindSet(y));
}

void
DisjointSet::Link(size_t x, size_t y)
{
	if (rank[x] > rank[y])
		p[y] = x;
	else
	{
		p[x] = y;

		if (rank[x] == rank[y])
			rank[y] = rank[y] + 1;
	}
}

size_t
DisjointSet::FindSet(size_t x)
{
	if (x != p[x])
		p[x] = FindSet(p[x]);

	return p[x];
}
