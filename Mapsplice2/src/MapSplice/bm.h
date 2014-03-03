#ifndef BM_H
#define BM_H


#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm> 
#include <map>

char alphabet[5] = {'G','T','A','C','N'};

class BM_Search
{
public:
	BM_Search(int _pattern_size)
	{
		pattern_size = _pattern_size;
		alphabet_size = 5;
		suffix_table = new int[pattern_size];
		good_table_fw = new int[pattern_size];
		good_table_rc = new int[pattern_size];
		for(int i = 0; i < alphabet_size; i++)
		{
			bad_table_fw[alphabet[i]] = pattern_size;
			bad_table_rc[alphabet[i]] = pattern_size;
		}
	}

	~BM_Search()
	{
		delete[] suffix_table;
		delete[] good_table_fw;
		delete[] good_table_rc;
	}

 	void build_bad_table(string &x, map<char, int>& bad_table)
 	{
 	  map<char, int>::iterator it;
 	  for (it = bad_table.begin(); it != bad_table.end(); it++)
    	it->second = pattern_size;
    for (int i = 0; i < pattern_size - 1; ++i)
      bad_table[x[i]] = pattern_size - i - 1;
 	}

 	void suffixes(string& x, int *suffix_table) 
 	{
   int f = 0, g = 0, i = 0;
   suffix_table[pattern_size - 1] = pattern_size;
   g = pattern_size - 1;
   for (i = pattern_size - 2; i >= 0; --i) 
   {
      if (i > g && suffix_table[i + pattern_size - 1 - f] < i - g)
         suffix_table[i] = suffix_table[i + pattern_size - 1 - f];
      else 
      {
         if (i < g)
            g = i;
         f = i;
         while (g >= 0 && x[g] == x[g + pattern_size - 1 - f])
            --g;
         suffix_table[i] = f - g;
      }
    }
	}

	void build_good_table(string& x, int* good_table) 
	{
   int i, j;
   suffixes(x, suffix_table);
   for (i = 0; i < pattern_size; ++i)
      good_table[i] = pattern_size;
   j = 0;
   for (i = pattern_size - 1; i >= 0; --i)
      if (suffix_table[i] == i + 1)
         for (; j < pattern_size - 1 - i; ++j)
            if (good_table[j] == pattern_size)
               good_table[j] = pattern_size - 1 - i;
   for (i = 0; i <= pattern_size - 2; ++i)
      good_table[pattern_size - 1 - suffix_table[i]] = pattern_size - 1 - i;
	}

	bool search_bm(string& x, string& y, int& start_pos, int& end_pos, int* good_table, map<char, int> bad_table, int &pos)
	{   
      int i = 0;
      int j = start_pos;
      int end_bound = end_pos - pattern_size;
   		while (j <= end_bound) 
   		{
      	for (i = pattern_size - 1; i >= 0 && x[i] == y[i + j]; --i);
     	 		if (i < 0) 
     	 		{
         		pos = j;
         		start_pos = j + good_table[0];
         		return true;
      		}
     		  else
						j += max(good_table[i], bad_table[y[i + j]] - pattern_size + 1 + i);
			}
			return false;
	}

	bool search_turbo_bm(string& x, string& y, int& start_pos, int& end_pos, int* good_table, map<char, int> bad_table, int &pos)
	{
   int bcShift = 0, i = 0, j = start_pos, shift = 0, u = 0, v = 0, turboShift = 0;
   shift = pattern_size;
   int end_bound = end_pos - pattern_size;
   while (j <= end_bound) 
   {
      i = pattern_size - 1;
      while (i >= 0 && x[i] == y[i + j]) {
         --i;
         if (u != 0 && i == pattern_size - 1 - shift)
            i -= u;
      }
      if (i < 0) {
         pos = j;
				 return true;
      }
      else {
         v = pattern_size - 1 - i;
         turboShift = u - v;
         bcShift = bad_table[y[i + j]] - pattern_size + 1 + i;
         shift = max(turboShift, bcShift);
         shift = max(shift, good_table[i]);
         if (shift == good_table[i])
            u = min(pattern_size - shift, v);
         else {
           if (turboShift < bcShift)
              shift = max(shift, u + 1);
           u = 0;
         }
      }
      j += shift;
   	}
   	return false;
	}

	bool search_next(string& text, int& start_pos, int& end_pos, int &pos, bool& map_fw)
	{
		if(search_fw && found_fw && update_fw)
		{
			next_start_fw = max(next_start_fw, start_pos);
			found_fw = search_bm(anchor_forward, text, next_start_fw, end_pos, good_table_fw, bad_table_fw, current_forward_anchor);
			next_start_fw = current_forward_anchor + 1;
			update_fw = false;
		}
		if(search_rc && found_rc && update_rc)
		{
			next_start_rc = max(next_start_rc + 1, start_pos);
			found_rc = search_bm(anchor_reverse, text, next_start_rc, end_pos, good_table_rc, bad_table_rc, current_reverse_anchor);
			next_start_rc = current_reverse_anchor + 1;
			update_rc = false;
		}
		int report = 0;
		if(found_fw && found_rc)
		{
			if(current_forward_anchor <= current_reverse_anchor)	
				report = 1;
			else
				report = -1;
		}
		else if(found_fw)
			report = 1;
		else if(found_rc)
			report = -1;
		if(report == 1)
		{
			pos = current_forward_anchor;
			map_fw = true;
			update_fw = true;
			return true;
		}
		else if(report == -1)	
		{
				pos = current_reverse_anchor;
				map_fw = false;
				update_rc = true;
				return true;
		}	
		else
			return false;
	}
	
	void set(string& m, bool _search_fw, bool _search_rc)
	{
		search_fw = _search_fw;
		search_rc = _search_rc;
		found_fw = search_fw;
	  found_rc = search_rc;
		update_fw = true;
		update_rc = true;
		current_forward_anchor = -1;
		current_reverse_anchor = -1;
		next_start_fw = -1;
		next_start_rc = -1;
		if(search_fw)
		{
			anchor_forward = m;
 			build_bad_table(anchor_forward, bad_table_fw);
 			build_good_table(anchor_forward, good_table_fw);
 		}
 		if(search_rc)
 		{
 			anchor_reverse = revcomp(m);
 			build_bad_table(anchor_reverse, bad_table_rc);
 			build_good_table(anchor_reverse, good_table_rc);
 		}
	}

	void reset()
	{
		found_fw = search_fw;
	  found_rc = search_rc;
		update_fw = true;
		update_rc = true;
		current_forward_anchor = -1;
		current_reverse_anchor = -1;
		next_start_fw = -1;
		next_start_rc = -1;
	}

private:
	int pattern_size;
	int alphabet_size;
	string anchor_forward;
	string anchor_reverse;
	int* suffix_table;
	int* good_table_fw;
	map<char, int> bad_table_fw;
	int* good_table_rc;
	map<char, int> bad_table_rc;
	int current_forward_anchor;
	int current_reverse_anchor;
	int next_start_fw;
	int next_start_rc;
	bool search_fw;
	bool search_rc;
	bool found_fw;
	bool found_rc;
	bool update_fw;
	bool update_rc;
};

#endif

