#ifndef KMP_H
#define KMP_H


#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm> 
#include <map>


class KMP_Search
{
public:
	KMP_Search(int _pattern_size)
	{
		pattern_size = _pattern_size;
		next_forward = new int[pattern_size];
		next_reverse = new int[pattern_size];
	}

	~KMP_Search()
	{
		delete[] next_forward;
		delete[] next_reverse;
	}

  void build_next_table(int* next, string& m)
  {
  	next[0] = 0;
  	int temp; 
    for(int i = 1; i < pattern_size; i++)
    {     
     temp = next[i - 1];
     while(temp > 0 && m[i] != m[temp])
     {  
				temp = next[temp - 1];
     }  
     if(m[i] == m[temp])
        next[i]= temp + 1;
     else next[i]=0;
    }	
  }

	bool search_kmp(string& text, const int& start_pos, const int& end_pos, string& m, int* next, int &pos)
	{   
    int tp = 0;
 		int mp = 0;  // text pointer and match string pointer;
		for(tp = start_pos; tp < end_pos; tp++)
  	{  
   		while(text[tp] != m[mp] && mp)
   		{
    		mp = next[mp - 1];
    	}
   		if(text[tp] == m[mp])
    		mp++;
   		if(mp == pattern_size)
   		{  
   			pos = tp - mp + 1;
      	return true;
   		}
  	}
   	return false;
	}

	bool search_next(string& text, const int& start_pos, const int& end_pos, int &pos, bool& map_fw)
	{
		if(search_fw && found_fw && update_fw)
		{
			found_fw = search_kmp(text, max(current_forward_anchor + 1, start_pos), end_pos, anchor_forward, next_forward, current_forward_anchor);
			update_fw = false;
		}
		if(search_rc && found_rc && update_rc)
		{
			found_rc = search_kmp(text, max(current_reverse_anchor + 1, start_pos), end_pos, anchor_reverse, next_reverse, current_reverse_anchor);
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
		if(search_fw)
		{
			anchor_forward = m;
 			build_next_table(next_forward, anchor_forward);
 		}
 		if(search_rc)
 		{
 			anchor_reverse = revcomp(m);
 			build_next_table(next_reverse, anchor_reverse);
 		}
	}

	void reset()
	{
		found_fw = search_fw;
	  found_rc = search_rc;
		update_fw = true;
		update_rc = true;
	}

private:
	int pattern_size;
	string anchor_forward;
	string anchor_reverse;
	int* next_forward;
	int* next_reverse;
	int current_forward_anchor;
	int current_reverse_anchor;
	bool search_fw;
	bool search_rc;
	bool found_fw;
	bool found_rc;
	bool update_fw;
	bool update_rc;
};

#endif