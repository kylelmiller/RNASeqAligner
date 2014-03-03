#include "SamRec.h"
#include "ReadNextTagAlignHandler.h"

size_t convert_point_coordinate(size_t cur_point, size_t x, size_t y, size_t x_prime, size_t y_prime, char strand)
{
	size_t cur_point_prime;

	if (strand == '+')
	{
		cur_point_prime = cur_point - x + x_prime;
	}
	else if (strand == '-')
	{
		cur_point_prime = y - cur_point + x_prime;
	}
	else
	{
		cerr << "strand?"<<endl;

		exit(0);
	}

	return cur_point_prime;
}


struct Fusion {

	string chr1;

	string chr2;

	size_t start;

	size_t end;

	char strand1;

	char strand2;

	set<string> gene1;

	set<string> gene2;

	string fusionline;

	Fusion (string line)
	{
		char chromname[1000]/*, strand1, strand2*/;

		char from[1000], type[1000], genestrand[100];

		char gene1[1000], gene2[1000];

		char skip[1000];

		//from_fusion	normal		POGK,	POGK,	chr1~chr1	166823269	166823352	JUNC_779	3	--	255,0,0	2	34,32,161,210,	0,116,	1.098612	
		//5	GTAG	0	0	0	21	16	8	3	0	3	1	2	0	3	0	2	166823303	166823320	166823269,62M48P50M|	166823144,84M53P72M|	
		//0	0	0.506667	0.528889	112	112	156	156	not_matched	not_matched	GTAACTAGGACAGATGAAAAATTAG	GGGAGCTTAAAAGATTTTACAAGAC

		//if (line.find("JUNC") != string::npos)

		sscanf(line.c_str(), "%s\t%s", from, type);

		if (string(type) == "intergenic")
			sscanf(line.c_str(), "%s\t%s\t%s\t%s\t%s\t%llu\t%llu\t%s\t%s\t%c%c", from, type, gene1, gene2, chromname, &start, &end, skip, skip, &strand1, &strand2);
		else
			sscanf(line.c_str(), "%s\t%s\t%s\t%s\t%s\t%s\t%llu\t%llu\t%s\t%s\t%c%c", from, type, genestrand, gene1, gene2, chromname, &start, &end, skip, skip, &strand1, &strand2);

		string chrnamestr = chromname;

		size_t idx = chrnamestr.find("chr", 4);

		idx = idx;

		chr1 = chrnamestr.substr(0, idx - 1);

		chr2 = chrnamestr.substr(idx, chrnamestr.length() - idx);

		string gene1str = gene1, gene2str = gene2;

		//CTSD,RP11-295K3.1,CTSD,

		{
			size_t startidx = 0, nextidx = 0;

			while (true)
			{
				nextidx = gene1str.find(',', startidx);

				string gene = gene1str.substr(startidx, nextidx - startidx);

				this->gene1.insert(gene);

				startidx = nextidx + 1;

				if (startidx >= gene1str.length())
					break;
			}
		}

		
		{
			size_t startidx = 0, nextidx = 0;

			while (true)
			{
				nextidx = gene2str.find(',', startidx);

				string gene = gene2str.substr(startidx, nextidx - startidx);

				this->gene2.insert(gene);

				startidx = nextidx + 1;

				if (startidx >= gene2str.length())
					break;
			}
		}

		fusionline = line;

		//else
		//	sscanf(line.c_str(), "%s\t%s\t%s\t%s\t%s\t%s\t%llu\t%llu\t%c%c", from, type, genestrand, gene1, gene2, chromname, &start, &end, &strand1, &strand2);

	}

	string tostring();
};

struct Exon {

	string name;

	string chr;

	size_t start;

	size_t end;

	char strand;

	string geneid;

	string transcriptid;

	size_t exonid;

	string genename;

	string transcriptname;
	//

	string line;

	//char protein_type[1000], exontype[1000];

	////size_t start, end;

	//char skip1[1000], skip2[1000];

	//char gene_id[1000], geneidname[1000]; 

	//char transcript_id[1000], transcriptidname[1000]; 

	//char exon_id[1000], exonidname[1000]; 

	//char gname[1000], gene_name[1000]; 

	//char tname[1000], transcript_name[1000]; 

	//char protein_id[1000], protein_id_name[1000];

	string protein_type_exontype;

	string after_end1;

	string after_end2;

	Exon(string line)
	{
		char chr[1000];//, protein_type[1000], exontype[1000];

		////size_t start, end;

		char protein_type[1000], exontype[1000];

		char skip1[1000], skip2[1000];

		char gene_id[1000], geneidname[1000]; 

		char transcript_id[1000], transcriptidname[1000]; 

		char exon_id[1000], exonidname[1000]; 

		char gname[1000], gene_name[1000]; 

		char tname[1000], transcript_name[1000]; 

		char protein_id[1000], protein_id_name[1000];

		//char skip1[1000], skip2[1000];

		//char gene_id[1000], geneidname[1000]; 

		//char transcript_id[1000], transcriptidname[1000]; 

		//char exon_id[1000], exonidname[1000]; 

		//char gname[1000], gene_name[1000]; 

		//char tname[1000], transcript_name[1000]; 

		//char protein_id[1000], protein_id_name[1000];

		//16	pseudogene	exon	61691	61763	.	+	.	 
		//gene_id "ENSG00000242566"; transcript_id "ENST00000430178"; 
		//exon_number "1"; gene_name "DDX11L10"; transcript_name "DDX11L10-001";

		sscanf(line.c_str(), "%s\t%s\t%s\t%llu\t%llu\t%s\t%c\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", chr, protein_type, exontype,
			&start, &end, skip1, &strand, skip2, gene_id, geneidname, transcript_id, transcriptidname, exon_id, exonidname,
			gname, gene_name, tname, transcript_name/*, protein_id, protein_id_name*/);

		this->chr = string("chr") + string(chr);

		char protein_type_exontype_chr[1000];

		sprintf(protein_type_exontype_chr, "%s\t%s", protein_type, exontype);

		protein_type_exontype = protein_type_exontype_chr;

		char after_end_chr1[5000];

		sprintf(after_end_chr1, "%s", skip1);

		after_end1 = after_end_chr1;

		char after_end_chr2[5000];

		sprintf(after_end_chr2, "%s\t%s %s %s %s %s %s %s %s %s %s", skip2, gene_id, geneidname, transcript_id, transcriptidname, exon_id, exonidname,
			gname, gene_name, tname, transcript_name/*, protein_id, protein_id_name*/);

		after_end2 = after_end_chr2;

		//if (string(gname) == "gene_biotype")
		//	continue;

		//if (string(exontype) != "exon")
		//	continue;

		geneid = geneidname; 

		geneid = geneid.substr(1, geneid.length() - 3);

		//

		transcriptid = transcriptidname; 

		transcriptid = transcriptid.substr(1, transcriptid.length() - 3);

		//
		string exonidstr = exonidname; 

		exonidstr = exonidstr.substr(1, exonidstr.length() - 3);

		exonid = atoi(exonidstr.c_str());

		//
		genename = gene_name; 

		genename = genename.substr(1, genename.length() - 3);

		//
		transcriptname = transcript_name;

		transcriptname = transcriptname.substr(1, transcriptname.length() - 3);

		this->line = line;

	}

	string tostring(string newchrname, size_t newstart, size_t newend, char strand)
	{
		char exonbuf[5000];

		if (this->strand == strand)
			strand = '+';
		else if (this->strand == '+')
			strand = '-';
		else if (this->strand == '-')
			strand = '+';

		sprintf(exonbuf, "%s\t%s\t%llu\t%llu\t%s\t%c\t%s", newchrname.c_str(), protein_type_exontype.c_str(),
			newstart, newend, after_end1.c_str(), strand, after_end2.c_str()/*, protein_id, protein_id_name*/);

		//ofs << corchrname <<'\t' << protein_type<<'\t' << exontype<<'\t' <<start<<'\t' <<end<<'\t' <<skip1<<'\t' <<strand<<'\t' <<skip2
		//	<<'\t' <<gene_id<<' ' << geneidname<<' ' << transcript_id<<' ' << transcriptidname<<' ' << exon_id<<' ' << exonidname
		//	<<' ' <<gname<<' ' << gene_name<<' ' << tname<<' ' << transcript_name<</*'\t' << protein_id<<'\t' << protein_id_name<<*/endl;
		return exonbuf;
	}
};

struct Transcript {

	string name;

	string chr;

	size_t start;

	size_t end;

	char strand;

	vector<Exon> exons;

	Transcript(string& name_, char strand_, vector<Exon>& exons_) : name(name_), strand(strand_), exons(exons_)
	{
	}

	Transcript()
	{
	}
};

struct Gene {
	string name;

	string chr;

	size_t start;

	size_t end;

	char strand;

	vector<Transcript> transcripts;

	Gene (const string& name_) : name(name_)
	{
	}

	Gene ()
	{
	}
};

struct FusionGene {

	Gene* gene1ptr;
	
	Gene* gene2ptr;

	vector<Fusion*> fusions;

	char strand1, strand2;

	size_t gene1_left_ext, gene1_right_ext;

	size_t gene2_left_ext, gene2_right_ext;

	size_t a, b, c, d;

	size_t gene1_start_ext, gene1_end_ext;

	size_t gene2_start_ext, gene2_end_ext;

	string fusiongenename;

	FusionGene(Gene* gene_1_ptr, Gene* gene_2_ptr, string& fusiongenename, Fusion* fusion_ptr, size_t gene1_left_extension, size_t gene1_right_extension, size_t gene2_left_extension, size_t gene2_right_extension) : 
		gene1ptr(gene_1_ptr), gene2ptr(gene_2_ptr), gene1_left_ext(gene1_left_extension), gene1_right_ext(gene1_right_extension), gene2_left_ext(gene2_left_extension), gene2_right_ext(gene2_right_extension)
	{
		strand1 = fusion_ptr->strand1;//gene_1_ptr->strand;

		strand2 = fusion_ptr->strand2;//gene_2_ptr->strand;

		this->fusiongenename = fusiongenename;

		a = 0;

		b = gene_1_ptr->end - gene_1_ptr->start + gene1_left_ext + gene1_right_ext;

		c = b + 1;

		d = gene_2_ptr->end - gene_2_ptr->start + c + gene2_left_ext + gene2_right_ext;

		gene1_start_ext = gene_1_ptr->start - gene1_left_ext;

		gene1_end_ext = gene_1_ptr->end + gene1_right_ext;

		gene2_start_ext = gene_2_ptr->start - gene2_left_ext;

		gene2_end_ext = gene_2_ptr->end + gene2_right_ext;

		fusions.push_back(fusion_ptr);		
	}

	FusionGene()
	{
	}

	void set(size_t gene1_left_extension, size_t gene1_right_extension, size_t gene2_left_extension, size_t gene2_right_extension)
	{
		if (gene1_left_extension > gene1_left_ext)
			gene1_left_ext = gene1_left_extension;
		
		if (gene1_right_extension > gene1_right_ext)
			gene1_right_ext = gene1_right_extension;
		
		if (gene2_left_extension > gene2_left_ext)
			gene2_left_ext = gene2_left_extension;
		
		if (gene2_right_extension > gene2_right_ext)
			gene2_right_ext = gene2_right_extension;

		a = 0;

		b = gene1ptr->end - gene1ptr->start + gene1_left_ext + gene1_right_ext;

		c = b + 1;

		d = gene2ptr->end - gene2ptr->start + c + gene2_left_ext + gene2_right_ext;

		gene1_start_ext = gene1ptr->start - gene1_left_ext;

		gene1_end_ext = gene1ptr->end + gene1_right_ext;

		gene2_start_ext = gene2ptr->start - gene2_left_ext;

		gene2_end_ext = gene2ptr->end + gene2_right_ext;
	}
};

void read_fusion(char* fusion_file, vector<Fusion>& fusions)
{
	ifstream ifs(fusion_file);

	if (ifs.is_open())
	{
		string line;
		while (getline(ifs,line))
		{
			if (line == "" || line.find("track") != string::npos)
				continue;

			fusions.push_back(Fusion(line));
		}

		ifs.close();
	}
	else
	{
		cout << "can't open file "<< fusion_file<<endl; exit(1);
	}
}

void read_annotation(char* annotation_file, vector<Exon>& exons, string outoutgene)
{
	ifstream ifs(annotation_file);

	string nonchrgene = outoutgene; nonchrgene.append(".notchrgenes"); ofstream ofs(nonchrgene.c_str());

	string notwantedgenetype = outoutgene; notwantedgenetype.append(".strange_gene_type"); ofstream ofsstrange(notwantedgenetype.c_str());

	if (ifs.is_open())
	{
		string line;
		while (getline(ifs,line))
		{
			if (line == "" || line.find("track") != string::npos)
				continue;

			//cout << line << endl;

			char chr[1000];

			char protein_type[1000];

			sscanf(line.c_str(), "%s\t%s", chr, protein_type);

			string protein_type_str = protein_type;

			if (((chr[0] >= '0' && chr[0] <= '9') || chr[0] == 'X' || chr[0] == 'Y' || chr[0] == 'M') &&
				(protein_type_str == "protein_coding" || 
				protein_type_str == "processed_transcript" || 
				protein_type_str == "IG_C_gene" || 
				protein_type_str == "IG_D_gene" || 
				protein_type_str == "IG_J_gene" || 
				protein_type_str == "IG_V_gene"))
			{
				exons.push_back(Exon(line));
			}
			else
			{
				ofs << line << endl;
			}
		}

		ifs.close();
	}
	else
	{
		cout << "can't open file "<< annotation_file<<endl; exit(1);
	}
}

bool check_same_gene_name(map<string, vector<Exon> >& transcripts)
{
	map<string, vector<Exon> >::iterator transcript_iter;

	set<string> gene_ids;

	for (transcript_iter = transcripts.begin(); transcript_iter != transcripts.end(); ++transcript_iter)
	{
		vector<Exon>::iterator exon_iter;

		for (exon_iter = transcript_iter->second.begin(); exon_iter != transcript_iter->second.end(); ++exon_iter)
		{
			gene_ids.insert(exon_iter->geneid);
		}
	}

	if (gene_ids.size() == 1)
		return false;
	else
		return true;
}

bool comp_gene(Gene* lhs, Gene* rhs)
{
	if (lhs->start == rhs->start)
		return lhs->end < rhs->end;
	else
		return lhs->start < rhs->start;
}

size_t exon2gene(vector<Exon>& exons, map<string, Gene>& genes, size_t extend_length, string same_genme, map<string, vector<Gene*> >& sorted_genes)
{
	string same_gene_file = same_genme; same_gene_file.append(".same_name");

	ofstream ofs(same_gene_file.c_str());

	size_t max_gene_len = 0;

	string max_gene;

	map<string, map<string, vector<Exon> > > stored_genes;

	vector<Exon>::iterator exon_iter;

	for (exon_iter = exons.begin(); exon_iter != exons.end(); ++exon_iter)
	{
		stored_genes[exon_iter->genename][exon_iter->transcriptid].push_back(*exon_iter);
	}

	map<string, map<string, vector<Exon> > >::iterator gene_iter;

	map<string, vector<Exon> >::iterator transcript_iter;

	for (gene_iter = stored_genes.begin(); gene_iter != stored_genes.end(); ++gene_iter)
	{
		if (check_same_gene_name(gene_iter->second))
		{
			for (transcript_iter = gene_iter->second.begin(); transcript_iter != gene_iter->second.end(); ++transcript_iter)
			{
				vector<Exon>::iterator exon_iter;

				for (exon_iter = transcript_iter->second.begin(); exon_iter != transcript_iter->second.end(); ++exon_iter)
				{
					ofs << exon_iter->line<<endl;
				}
			}
			continue;
		}

		Gene gene(gene_iter->first);

		//Gene gene(gene_iter->second.begin()->second.begin()->genename);

		for (transcript_iter = gene_iter->second.begin(); transcript_iter != gene_iter->second.end(); ++transcript_iter)
		{
			Transcript transcript(transcript_iter->second.front().transcriptid, transcript_iter->second.front().strand, transcript_iter->second);

			vector<Exon>::iterator exon_iter;

			size_t min_start = -1, max_end = 0;

			for (exon_iter = transcript_iter->second.begin(); exon_iter != transcript_iter->second.end(); ++exon_iter)
			{
				if (min_start > exon_iter->start)
					min_start = exon_iter->start;

				if (max_end < exon_iter->end)
					max_end = exon_iter->end;
			}

			transcript.chr = transcript_iter->second.front().chr;

			transcript.start = min_start;

			transcript.end = max_end;

			gene.transcripts.push_back(transcript);
		}

		vector<Transcript>::iterator transcript_iter;

		size_t min_start = -1, max_end = 0;

		for (transcript_iter = gene.transcripts.begin(); transcript_iter != gene.transcripts.end(); ++transcript_iter)
		{
			if (min_start > transcript_iter->start)
				min_start = transcript_iter->start;

			if (max_end < transcript_iter->end)
				max_end = transcript_iter->end;
		}

		if (min_start > extend_length)
			min_start = min_start - extend_length;
		else
			min_start = 1;

		max_end = max_end + extend_length;

		gene.chr = gene.transcripts.front().chr;

		gene.start = min_start;

		gene.end = max_end;

		gene.strand = gene.transcripts.front().strand;

		genes[gene.name] = gene;

		if (max_gene_len < max_end - min_start)
		{
			max_gene = gene_iter->first;
			max_gene_len = max_end - min_start;
		}
	}

	map<string, Gene>::iterator gene_iter2;

	for (gene_iter2 = genes.begin(); gene_iter2 != genes.end(); ++gene_iter2)
	{
		sorted_genes[gene_iter2->second.chr].push_back(&(gene_iter2->second));
	}

	map<string, vector<Gene*> >::iterator sorted_gene_iter;

	for (sorted_gene_iter = sorted_genes.begin(); sorted_gene_iter != sorted_genes.end(); ++sorted_gene_iter)
	{
		sort(sorted_gene_iter->second.begin(), sorted_gene_iter->second.end(), comp_gene);
	}

	cout << "max_gene:"<<max_gene<<endl;

	cout << "max_gene_len:"<<max_gene_len<<endl;

	return max_gene_len;
}

Gene* find_closest_gene(string& chr, size_t pos, char strand, map<string, vector<Gene* > >& sorted_genes, int& left_or_right)
{
	if (sorted_genes.find(chr) == sorted_genes.end())
		return 0;

	vector<Gene* >& chr_sorted_genes = sorted_genes[chr];

	Gene dump_gene;

	dump_gene.start = pos;

	dump_gene.end = pos;

	vector<Gene* >::iterator gene_lb_iter = lower_bound(chr_sorted_genes.begin(), chr_sorted_genes.end(), &dump_gene, comp_gene);

	if (gene_lb_iter == chr_sorted_genes.begin())
	{
		//if ((*gene_lb_iter)->strand == strand)
		//{
			left_or_right = 1;

			return *gene_lb_iter;
		//}
		//else
		//	return 0;
	}
	else if (gene_lb_iter == chr_sorted_genes.end())
	{
		//if ((*(gene_lb_iter - 1))->strand == strand)
		//{
			left_or_right = -1;

			return *(gene_lb_iter - 1);
		//}
		//else
		//	return 0;
	}
	else
	{
		if (((*gene_lb_iter)->start - pos > pos - (*(gene_lb_iter - 1))->end)/* && (*(gene_lb_iter - 1))->strand == strand*/)
		{
			left_or_right = -1;

			return *(gene_lb_iter - 1);
		}
		else /*if ((*gene_lb_iter)->strand == strand)*/
		{
			left_or_right = 1;

			return *gene_lb_iter;
		}
	}
}

void locate_corresponding_gene(vector<Fusion>& fusions, map<string, Gene>& genes, map<string, FusionGene>& fusiongenes, vector<pair<string, size_t> >& syn_chrom_name_length, map<string, vector<Gene* > >& sorted_genes, size_t extend_length)
{
	vector<Fusion>::iterator fusion_iter;

	for (fusion_iter = fusions.begin(); fusion_iter != fusions.end(); ++fusion_iter)
	{
		//resolve multiple genes

		string gene1name, gene2name;

		if (fusion_iter->gene1.size() > 1)
		{
			size_t max_gene = 0;

			string max_gene_name;

			set<string>::iterator genename_iter;

			for (genename_iter = fusion_iter->gene1.begin(); genename_iter != fusion_iter->gene1.end(); ++genename_iter)
			{
				if (genes.find(*genename_iter) == genes.end())
				{
					cerr << "can't find multiple gene:" << *genename_iter << endl;

					exit(0);
				}

				if (max_gene < genes[*genename_iter].end - genes[*genename_iter].start)
				{
					max_gene = genes[*genename_iter].end - genes[*genename_iter].start;

					max_gene_name = *genename_iter;
				}
			}

			gene1name = max_gene_name;
		}
		else
		{
			gene1name = (*(fusion_iter->gene1.begin()));
		}

		if (fusion_iter->gene2.size() > 1)
		{
			size_t max_gene = 0;

			string max_gene_name;

			set<string>::iterator genename_iter;

			for (genename_iter = fusion_iter->gene2.begin(); genename_iter != fusion_iter->gene2.end(); ++genename_iter)
			{
				if (genes.find(*genename_iter) == genes.end())
				{
					cerr << "can't find multiple gene:" << *genename_iter << endl;

					exit(0);
				}

				if (max_gene < genes[*genename_iter].end - genes[*genename_iter].start)
				{
					max_gene = genes[*genename_iter].end - genes[*genename_iter].start;

					max_gene_name = *genename_iter;
				}
			}

			gene2name = max_gene_name;
		}
		else
		{
			gene2name = (*(fusion_iter->gene2.begin()));
		}

		size_t gene1_left_extension = 0;

		size_t gene1_right_extension = 0;

		size_t gene2_left_extension = 0;

		size_t gene2_right_extension = 0;

		if (gene1name == "-")
		{
			int left_or_right = 0;

			Gene* closest_gene = find_closest_gene(fusion_iter->chr1, fusion_iter->start, fusion_iter->strand1, sorted_genes, left_or_right);

			if (closest_gene == 0)
				continue;
			else if (left_or_right == -1)
			{
				gene1_right_extension = static_cast<int>(fusion_iter->start - closest_gene->end + extend_length);

				gene1name = closest_gene->name;
			}
			else if (left_or_right == 1)
			{
				gene1_left_extension = static_cast<int>(closest_gene->start - fusion_iter->start + extend_length);

				gene1name = closest_gene->name;
			}
			else
			{
				cout << "what ?"<<endl;
			}
		}
		else if (genes.find(gene1name) == genes.end())
		{
			continue;
		}
		else
		{
		}

		if (gene2name == "-")
		{
			int left_or_right = 0;

			Gene* closest_gene = find_closest_gene(fusion_iter->chr2, fusion_iter->end, fusion_iter->strand2, sorted_genes, left_or_right);

			if (closest_gene == 0)
				continue;
			else if (left_or_right == -1)
			{
				gene2_right_extension = static_cast<int>(fusion_iter->end - closest_gene->end + extend_length);

				gene2name = closest_gene->name;
			}
			else if (left_or_right == 1)
			{
				gene2_left_extension = static_cast<int>(closest_gene->start - fusion_iter->end + extend_length);

				gene2name = closest_gene->name;
			}
			else
			{
				cout << "what ?"<<endl;
			}
		}
		else if (genes.find(gene2name) == genes.end())
		{
			continue;
		}
		else
		{
		}

		string combgene = gene1name + "~" + gene2name;

		if (gene1_left_extension || gene1_right_extension || gene2_left_extension || gene2_right_extension)
		{
			char combgenechr[5000];

			sprintf(combgenechr, "%llu~%s~%llu~%llu~%s~%llu", gene1_left_extension, gene1name.c_str(), gene1_right_extension, gene2_left_extension, gene2name.c_str(), gene2_right_extension);

			combgene = combgenechr;
		}

		if (fusiongenes.find(combgene) == fusiongenes.end())
		{
			FusionGene fusiongene(&(genes[gene1name]), &(genes[gene2name]), combgene, &(*fusion_iter), gene1_left_extension, gene1_right_extension, gene2_left_extension, gene2_right_extension);

			fusiongenes[combgene] = fusiongene;

			size_t syn_chrom_length = genes[gene1name].end - genes[gene1name].start + genes[gene2name].end - genes[gene2name].start + gene1_left_extension + gene1_right_extension + gene2_left_extension + gene2_right_extension;

			syn_chrom_name_length.push_back(make_pair(combgene, syn_chrom_length));
		}
		else
		{
			fusiongenes[combgene].fusions.push_back(&(*fusion_iter));

			fusiongenes[combgene].set(gene1_left_extension, gene1_right_extension, gene2_left_extension, gene2_right_extension);
		}
	}	
}

struct GeneGenomeLoc {

	string chrname;

	size_t start;

	size_t end;	
	
	bool start_or_end;

	Gene* gene_ptr;

	FusionGene* fusiongene_ptr;

	GeneGenomeLoc (string& chrname_, size_t start_, size_t end_, bool start_or_end_, Gene* gene_ptr_, FusionGene* fusiongene_ptr_) : chrname(chrname_), 
		start(start_), end(end_),  start_or_end(start_or_end_), gene_ptr(gene_ptr_), fusiongene_ptr(fusiongene_ptr_) {}

	GeneGenomeLoc (Gene* gene_ptr_, FusionGene* fusiongene_ptr_, bool start_or_end_) : gene_ptr(gene_ptr_), fusiongene_ptr(fusiongene_ptr_), start_or_end(start_or_end_)
	{
		chrname = gene_ptr->chr;

		if (start_or_end_)
		{
			start = fusiongene_ptr->gene1_start_ext;

			end = fusiongene_ptr->gene1_end_ext;
		}
		else
		{
			start = fusiongene_ptr->gene2_start_ext;

			end = fusiongene_ptr->gene2_end_ext;
		}			
	}
};

bool comp_genome_loc(const GeneGenomeLoc& lhs, const GeneGenomeLoc& rhs)
{
	if (lhs.start == rhs.start)
		return lhs.end < rhs.end;
	else
		return lhs.start < rhs.start;
}

bool comp_genome_loc2(const GeneGenomeLoc& lhs, const GeneGenomeLoc& rhs)
{
	//if (lhs.start == rhs.start)
	//return lhs.end < rhs.end;
	//else
	return lhs.start < rhs.start;
}

void gene2genome_coordinate(map<string, FusionGene>& fusiongenes, map<string, vector<GeneGenomeLoc> >& genes_genome_locs)
{
	map<string, FusionGene>::iterator gene_iter;

	for (gene_iter = fusiongenes.begin(); gene_iter != fusiongenes.end(); ++gene_iter)
	{
		genes_genome_locs[gene_iter->second.gene1ptr->chr].push_back(GeneGenomeLoc(gene_iter->second.gene1ptr, &((gene_iter->second)), true));

		genes_genome_locs[gene_iter->second.gene2ptr->chr].push_back(GeneGenomeLoc(gene_iter->second.gene2ptr, &((gene_iter->second)), false));
	}

	map<string, vector<GeneGenomeLoc> >::iterator chromosome_iter;

	for (chromosome_iter = genes_genome_locs.begin(); chromosome_iter != genes_genome_locs.end(); ++chromosome_iter)
	{
		sort(chromosome_iter->second.begin(), chromosome_iter->second.end(), comp_genome_loc);
	}
}

bool MergeStdFusion(vector<SamRec>::iterator sam_rec_ptr1, vector<SamRec>::iterator sam_rec_ptr2)
{
	//fusion 1
	sam_rec_ptr1->splice_way2 = sam_rec_ptr2->splice_way;

	sam_rec_ptr1->ori_splice_way2 = sam_rec_ptr2->ori_splice_way;

	sam_rec_ptr1->chrom_name2 = sam_rec_ptr2->chrom_name;

	sam_rec_ptr1->start2 = sam_rec_ptr2->start;

	sam_rec_ptr1->end2 = sam_rec_ptr2->end;

	sam_rec_ptr1->strand_t2 = sam_rec_ptr2->strand_t;

	sam_rec_ptr1->mappedlen1 = sam_rec_ptr1->mappedlen;

	sam_rec_ptr1->mappedlen2 = sam_rec_ptr2->mappedlen;

	sam_rec_ptr1->mappedlen = sam_rec_ptr1->mappedlen1 + sam_rec_ptr1->mappedlen2;

	sam_rec_ptr1->fusion_prefix_len = sam_rec_ptr1->mappedlen1;

	sam_rec_ptr1->fusion_suffix_len = sam_rec_ptr1->mappedlen2;

	//cout << "sam_rec_ptr1->fusion_prefix_len:"<<sam_rec_ptr1->fusion_prefix_len<<endl;

	//cout << "sam_rec_ptr1->fusion_suffix_len:"<<sam_rec_ptr1->fusion_suffix_len<<endl;

	sam_rec_ptr1->intron_size = sam_rec_ptr1->intron_size + sam_rec_ptr2->intron_size;

	sam_rec_ptr1->spliceway_vec2 = sam_rec_ptr2->spliceway_vec;

	sam_rec_ptr1->is_fusion = true;

	if (sam_rec_ptr1->strand_t & IS_REVERSE)
	{
		sam_rec_ptr1->strand1 = '-';
	}
	else
	{
		sam_rec_ptr1->strand1 = '+';
	}

	if (sam_rec_ptr1->strand_t2 & IS_REVERSE)
	{
		sam_rec_ptr1->strand2 = '-';
	}
	else
	{
		sam_rec_ptr1->strand2 = '+';
	}

	if (sam_rec_ptr1->strand1 == '+')
	{
		sam_rec_ptr1->fusion_prefix_st = sam_rec_ptr1->start;

		sam_rec_ptr1->fusion_prefix_end = sam_rec_ptr1->end;
	}
	else
	{
		sam_rec_ptr1->fusion_prefix_st = sam_rec_ptr1->end;

		sam_rec_ptr1->fusion_prefix_end = sam_rec_ptr1->start;
	}

	if (sam_rec_ptr1->strand2 == '+')
	{
		sam_rec_ptr1->fusion_suffix_st = sam_rec_ptr1->start2;

		sam_rec_ptr1->fusion_suffix_end = sam_rec_ptr1->end2;
	}
	else
	{
		sam_rec_ptr1->fusion_suffix_st = sam_rec_ptr1->end2;

		sam_rec_ptr1->fusion_suffix_end = sam_rec_ptr1->start2;
	}

	if (sam_rec_ptr1->min_anchor > sam_rec_ptr2->min_anchor)
	{
		sam_rec_ptr1->min_anchor = sam_rec_ptr2->min_anchor;
	}

	if (sam_rec_ptr1->max_anchor < sam_rec_ptr2->max_anchor)
	{
		sam_rec_ptr1->max_anchor = sam_rec_ptr2->max_anchor;
	}

	if (sam_rec_ptr1->chrom_name < sam_rec_ptr1->chrom_name2 || (sam_rec_ptr1->chrom_name == sam_rec_ptr1->chrom_name2 && sam_rec_ptr1->start <= sam_rec_ptr1->start2))
	{
	}
	else
	{
		sam_rec_ptr1->Swap();
	}

	return true;
}

void fusiostd2oneline(pair<vector<SamRec>, vector<SamRec> >& stored_tag_aligns_pe, pair<vector<SamRec*>, vector<SamRec*> >& merged_tag_aligns_pe)
{
	vector<SamRec>::iterator align_iter;

	for (align_iter = stored_tag_aligns_pe.first.begin(); align_iter != stored_tag_aligns_pe.first.end(); ++align_iter)
	{
		if (align_iter->is_fusion_newfmt)
		{
			if (align_iter + 1 == stored_tag_aligns_pe.first.end() || (align_iter + 1)->is_fusion_newfmt == false)
			{
				cout << "new format fusion alignment not together:"<<endl <<align_iter->tostring(0, 0 ) << endl;

				continue;
			}

			MergeStdFusion(align_iter, align_iter + 1);

			merged_tag_aligns_pe.first.push_back(&(*align_iter));

			++align_iter;
		}
		else
		{
			merged_tag_aligns_pe.first.push_back(&(*align_iter));
		}
	}

	for (align_iter = stored_tag_aligns_pe.second.begin(); align_iter != stored_tag_aligns_pe.second.end(); ++align_iter)
	{
		if (align_iter->is_fusion_newfmt)
		{
			if (align_iter + 1 == stored_tag_aligns_pe.second.end() || (align_iter + 1)->is_fusion_newfmt == false)
			{
				cout << "new format fusion alignment not together:"<<endl <<align_iter->tostring(0, 0 ) << endl;

				continue;
			}

			MergeStdFusion(align_iter, align_iter + 1);

			merged_tag_aligns_pe.second.push_back(&(*align_iter));

			++align_iter;
		}
		else
		{
			merged_tag_aligns_pe.second.push_back(&(*align_iter));
		}
	}
}

bool find_corresponding_fusion_gene(SamRec* sam_rec_ptr, map<string, vector<GeneGenomeLoc> >& genes_genome_locs, size_t max_gene_length, vector<GeneGenomeLoc*>& fitted_genome_locs)
{
	if (sam_rec_ptr->isunmapped)
		return false;

	SamRec& samrec = *sam_rec_ptr;

	map<string, vector<GeneGenomeLoc > >::iterator chr_iter = genes_genome_locs.find(samrec.chrom_name);

	if (chr_iter == genes_genome_locs.end())
		return false;

	if (sam_rec_ptr->tag_name == "HWUSI-EAS1531_0001:4:100:1066:348#0/1")
	{
		cout << "HWUSI-EAS1531_0001:4:100:1066:348#0/1:"<<chr_iter->first<<endl;
	}

	if (sam_rec_ptr->tag_name == "HWUSI-EAS1531_0001:4:100:1066:348#0/1")
	{
		cout << "HWUSI-EAS1531_0001:4:100:1066:348#0/1:"<<sam_rec_ptr->start<<endl;
	}

	GeneGenomeLoc cur_loc(samrec.chrom_name, 0, 0, true, 0, 0);

	if (samrec.start > max_gene_length)
		cur_loc = GeneGenomeLoc(samrec.chrom_name, samrec.start - max_gene_length, samrec.start - max_gene_length, true, 0, 0);

	if (sam_rec_ptr->tag_name == "HWUSI-EAS1531_0001:4:100:1066:348#0/1")
	{
		cout << "HWUSI-EAS1531_0001:4:100:1066:348#0/1:"<<cur_loc.start<<endl;
	}

	vector<GeneGenomeLoc >::iterator genomeloc_lb_iter = lower_bound(chr_iter->second.begin(), chr_iter->second.end(), cur_loc, comp_genome_loc2);

	vector<GeneGenomeLoc >::iterator genomeloc_iter;;

	//size_t ori_start = samrec.start;

	for (genomeloc_iter = genomeloc_lb_iter; genomeloc_iter != chr_iter->second.end(); ++genomeloc_iter)
	{
		if (sam_rec_ptr->tag_name == "HWUSI-EAS1531_0001:4:100:1066:348#0/1")
		{
			cout << "HWUSI-EAS1531_0001:4:100:1066:348#0/1:"<<genomeloc_iter->chrname << '\t' <<  genomeloc_iter->start << '\t' << genomeloc_iter->end << endl;
		}

		if (genomeloc_iter->start > samrec.end)//genomeloc_iter->end > samrec.start
		{
			if (sam_rec_ptr->tag_name == "HWUSI-EAS1531_0001:4:100:1066:348#0/1")
			{
				cout << "genomeloc_iter->end:"<<genomeloc_iter->end<<endl;

				cout << "samrec.start:"<<samrec.start<<endl;

				cout << "HWUSI-EAS1531_0001:4:100:1066:348#0/1:"<<"break"<<endl;
			}
			break;
		}

//HWUSI-EAS1531_0001:4:100:1066:348#0/1:unique
//HWUSI-EAS1531_0001:4:100:1066:348#0/1:chr17
//HWUSI-EAS1531_0001:4:100:1066:348#0/1:59445791
//HWUSI-EAS1531_0001:4:100:1066:348#0/1:0
//HWUSI-EAS1531_0001:4:100:1066:348#0/1:chr17	19348246	19351169
//genomeloc_iter->end:19351169
//samrec.start:59445791
//HWUSI-EAS1531_0001:4:100:1066:348#0/1:break
//HWUSI-EAS1531_0001:4:100:1066:348#0/1:0	0

		if (genomeloc_iter->end >= samrec.end &&  genomeloc_iter->start <= samrec.start)// fitted
		{
			if (samrec.is_fusion == false)
				fitted_genome_locs.push_back(&(*genomeloc_iter));
			else
			{
				if (genomeloc_iter->start_or_end)
				{
					Gene* the_other_gene_ptr = genomeloc_iter->fusiongene_ptr->gene2ptr;

					if (the_other_gene_ptr->chr == samrec.chrom_name2 &&
						genomeloc_iter->fusiongene_ptr->gene2_end_ext >= samrec.end2 &&
						genomeloc_iter->fusiongene_ptr->gene2_start_ext <= samrec.start2)
					{
						fitted_genome_locs.push_back(&(*genomeloc_iter));
					}
				}
				else
				{
					Gene* the_other_gene_ptr = genomeloc_iter->fusiongene_ptr->gene1ptr;

					if (the_other_gene_ptr->chr == samrec.chrom_name2 &&
						genomeloc_iter->fusiongene_ptr->gene1_end_ext >= samrec.end2 &&
						genomeloc_iter->fusiongene_ptr->gene1_start_ext <= samrec.start2)
					{
						fitted_genome_locs.push_back(&(*genomeloc_iter));
					}
				}
			}
		}
	}

	if (fitted_genome_locs.size())
		return true;
	else
		return false;
}

//size_t convert_point_coordinate(size_t cur_point, size_t x, size_t y, size_t x_prime, size_t y_prime, char strand)
//{
//	size_t cur_point_prime;
//
//	if (strand == '+')
//	{
//		cur_point_prime = cur_point - x + x_prime;
//	}
//	else if (strand == '-')
//	{
//		cur_point_prime = y - cur_point + x_prime;
//	}
//	else
//	{
//		cerr << "strand?"<<endl;
//
//		exit(0);
//	}
//
//	return cur_point_prime;
//}

string reverse_jump_code(vector<pair<size_t, int> >& spliceway_vec)
{
	string jump_code;

	for (int i = spliceway_vec.size() - 1; i >= 0; --i)
	{
		char jump_code_chr[100];

		if (spliceway_vec[i].second > 0)
		{
			if ((i != spliceway_vec.size() - 1) && spliceway_vec[i + 1].second > 0)
			{
				int intron = spliceway_vec[i + 1].first - spliceway_vec[i].first - spliceway_vec[i].second;

				sprintf(jump_code_chr, "%dN%dM", intron, spliceway_vec[i].second);
			}
			else
			{
				sprintf(jump_code_chr, "%dM", spliceway_vec[i].second);
			}
		}
		else
		{
			int insertlen = -spliceway_vec[i].second;

			sprintf(jump_code_chr, "%dI", insertlen);
		}

		jump_code.append(jump_code_chr);
	}

	return jump_code;
}

void convert_alignment_coordinate(GeneGenomeLoc* genome_loc_ptr, SamRec& samrec)
{
	size_t x_prime, y_prime, x_prime2, y_prime2;

	Gene* the_other_gene_loc_ptr;

	char strand1, strand2;

	size_t the_other_start, the_other_end;

	if (genome_loc_ptr->start_or_end)
	{
		x_prime = genome_loc_ptr->fusiongene_ptr->a;

		y_prime = genome_loc_ptr->fusiongene_ptr->b;

		x_prime2 = genome_loc_ptr->fusiongene_ptr->c;

		y_prime2 = genome_loc_ptr->fusiongene_ptr->d;

		strand1 = genome_loc_ptr->fusiongene_ptr->strand1;

		strand2 = genome_loc_ptr->fusiongene_ptr->strand2;

		the_other_start = genome_loc_ptr->fusiongene_ptr->gene2_start_ext;

		the_other_end = genome_loc_ptr->fusiongene_ptr->gene2_end_ext;

		the_other_gene_loc_ptr = genome_loc_ptr->fusiongene_ptr->gene2ptr;
		
	}
	else
	{
		x_prime = genome_loc_ptr->fusiongene_ptr->c;

		y_prime = genome_loc_ptr->fusiongene_ptr->d;

		x_prime2 = genome_loc_ptr->fusiongene_ptr->a;

		y_prime2 = genome_loc_ptr->fusiongene_ptr->b;

		strand1 = genome_loc_ptr->fusiongene_ptr->strand2;

		strand2 = genome_loc_ptr->fusiongene_ptr->strand1;

		the_other_start = genome_loc_ptr->fusiongene_ptr->gene1_start_ext;

		the_other_end = genome_loc_ptr->fusiongene_ptr->gene1_end_ext;

		the_other_gene_loc_ptr = genome_loc_ptr->fusiongene_ptr->gene1ptr;
	}

	if (samrec.is_fusion == false)
	{
		samrec.start = convert_point_coordinate(samrec.start, genome_loc_ptr->start, genome_loc_ptr->end, x_prime, y_prime, strand1);

		samrec.end = convert_point_coordinate(samrec.end, genome_loc_ptr->start, genome_loc_ptr->end, x_prime, y_prime, strand1);

		if (strand1 == '-')
		{
			samrec.ori_splice_way = reverse_jump_code(samrec.spliceway_vec);

			samrec.mapped_seq = revcomp(samrec.mapped_seq);

			reverse(samrec.qual_str.begin(), samrec.qual_str.end());

			swap(samrec.start, samrec.end);
		}

		samrec.chrom_name = string("chr") + genome_loc_ptr->fusiongene_ptr->fusiongenename;
	}
	else
	{
		string mapped_seq1 = samrec.mapped_seq.substr(0, samrec.mappedlen1);
		
		string mapped_seq2 = samrec.mapped_seq.substr(samrec.mapped_seq.length() - samrec.mappedlen2, samrec.mappedlen2);

		if (strand1 == '-')
		{
			samrec.splice_way = reverse_jump_code(samrec.spliceway_vec);

			mapped_seq1 = revcomp(mapped_seq1);

			if (samrec.strand_t & IS_REVERSE)
				samrec.strand_t = samrec.strand_t - IS_REVERSE;
			else
				samrec.strand_t |= IS_REVERSE;
		}

		if (strand2 == '-')
		{
			samrec.splice_way2 = reverse_jump_code(samrec.spliceway_vec2);

			mapped_seq2 = revcomp(mapped_seq2);

			if (samrec.strand_t2 & IS_REVERSE)
				samrec.strand_t2 = samrec.strand_t2 - IS_REVERSE;
			else
				samrec.strand_t2 |= IS_REVERSE;
		}

		samrec.start = convert_point_coordinate(samrec.start, genome_loc_ptr->start, genome_loc_ptr->end, x_prime, y_prime, strand1);

		samrec.end = convert_point_coordinate(samrec.end, genome_loc_ptr->start, genome_loc_ptr->end, x_prime, y_prime, strand1);
		
		if (strand1 == '-')
			swap(samrec.start, samrec.end);

		samrec.start2 = convert_point_coordinate(samrec.start2, the_other_start, the_other_end, x_prime2, y_prime2, strand2);

		samrec.end2 = convert_point_coordinate(samrec.end2, the_other_start, the_other_end, x_prime2, y_prime2, strand2);

		if (strand2 == '-')
			swap(samrec.start2, samrec.end2);

		if (samrec.start2 > samrec.start)
		{
			if (samrec.start2 < samrec.end)
			{
				cout << "samrec.start2 < samrec.end" <<endl;

				cout << samrec.tag_base_name << endl;
			}

			size_t intron = samrec.start2 - samrec.end - 1;

			char intronchr[1000];

			sprintf(intronchr, "%lluN", intron);

			samrec.ori_splice_way = samrec.splice_way + string(intronchr) + samrec.splice_way2;

			samrec.chrom_name = string("chr") + genome_loc_ptr->fusiongene_ptr->fusiongenename;

			samrec.is_fusion = false;

			samrec.is_fusion_newfmt = false;

			samrec.end = samrec.end2;

			samrec.mapped_seq = mapped_seq1 + mapped_seq2;
		}
		else
		{
			if (samrec.start < samrec.end2)
			{
				cout << "samrec.start < samrec.end2" <<endl;

				cout << samrec.tag_base_name << endl;
			}

			size_t intron = samrec.start - samrec.end2 - 1;

			char intronchr[1000];

			sprintf(intronchr, "%lluN", intron);

			samrec.ori_splice_way = samrec.splice_way2 + string(intronchr) + samrec.splice_way;

			samrec.chrom_name = string("chr") + genome_loc_ptr->fusiongene_ptr->fusiongenename;

			samrec.is_fusion = false;

			samrec.is_fusion_newfmt = false;

			samrec.start = samrec.start2;

			samrec.strand_t = samrec.strand_t2;
			//samrec.end = samrec.end2;

			samrec.mapped_seq = mapped_seq2 + mapped_seq1;
		}
	}
}

void set_mate_pair_distance(SamRec& samrec1, SamRec& samrec2)
{
	int st1 = samrec1.start, end1 = samrec1.end, st2 = samrec2.start, end2 = samrec2.end;
	
	int mate_dist1, mate_dist2;

	mate_dist1 = ((st1 - end2));

	mate_dist2 = ((end1 - st2));

	int outter_dist;

	if (abs(mate_dist1) < abs(mate_dist2))
	{
		if (mate_dist2 > 0)
			outter_dist = mate_dist2 + 1;
		else
			outter_dist = mate_dist2 - 1;
	}
	else
	{
		if (mate_dist1 > 0)
			outter_dist = mate_dist1 + 1;
		else
			outter_dist = mate_dist1 - 1;
	}
	
	samrec1.mate_offset = samrec2.start;

	samrec1.mate_diff = -outter_dist;

	samrec1.mate_match = "=";

	if (samrec2.strand_t & IS_REVERSE)
		samrec1.strand_t |= IS_MATE_REVERSE;

	samrec1.strand_t |= IS_PAIRED_MAPPED;

	samrec1.strand_t |= IS_PAIRED;

	samrec1.strand_t |= IS_FIRST_END;
	//

	samrec2.mate_offset = samrec1.start;

	samrec2.mate_diff = outter_dist;

	samrec2.mate_match = "=";

	if (samrec1.strand_t & IS_REVERSE)
		samrec2.strand_t |= IS_MATE_REVERSE;

	samrec2.strand_t |= IS_PAIRED_MAPPED;

	samrec2.strand_t |= IS_PAIRED;

	samrec2.strand_t |= IS_FIRST_END;
}

void extract_and_convert_sam(map<string, vector<GeneGenomeLoc> >& genes_genome_locs, char* sam_file, char* output_sam_file, size_t max_gene_len, vector<pair<string, size_t> >& syn_chrom_name_length)
{
	ifstream ifs(sam_file);

	ofstream ofs(output_sam_file);

	//write header of sam file

	//@SQ     SN:chr20        LN:63025520

	vector<pair<string, size_t> >::iterator chriter;

	for (chriter = syn_chrom_name_length.begin(); chriter != syn_chrom_name_length.end(); ++chriter)
	{
		ofs <<"@SQ\tSN:chr"<< chriter->first<<"\tLN:" << chriter->second<<endl;
	}

	ReadNextTagAlignHandler<SamRec> readnextalignmenthandler(sam_file, 5);

	pair<vector<SamRec>, vector<SamRec> > stored_tag_aligns_pe;

	size_t cur_read_id;

	while(readnextalignmenthandler.ReadNextTagAlignPE(stored_tag_aligns_pe, cur_read_id))
	{
		pair<vector<SamRec*>, vector<SamRec*> > merged_tag_aligns_pe;

		fusiostd2oneline(stored_tag_aligns_pe, merged_tag_aligns_pe);

		if (merged_tag_aligns_pe.first.size() > 1 || merged_tag_aligns_pe.second.size() > 1)//multiple alignments
		{
			if (merged_tag_aligns_pe.first.size() && merged_tag_aligns_pe.first.front()->tag_name == "HWUSI-EAS1531_0001:4:100:1066:348#0/1")
			{
				cout << "HWUSI-EAS1531_0001:4:100:1066:348#0/1:multiple"<<endl;
			}

			//cout << "multiple "<<endl;

			//if (merged_tag_aligns_pe.first.size() > 1)
			//	cout << merged_tag_aligns_pe.first.front()->tag_name << endl;
			//else if (merged_tag_aligns_pe.second.size() > 1)
			//	cout << merged_tag_aligns_pe.second.front()->tag_name << endl;
		}
		else if (merged_tag_aligns_pe.first.size() == 1 && merged_tag_aligns_pe.second.size() && 1)
		{
			if (merged_tag_aligns_pe.first.size() && merged_tag_aligns_pe.first.front()->tag_name == "HWUSI-EAS1531_0001:4:100:1066:348#0/1")
			{
				cout << "HWUSI-EAS1531_0001:4:100:1066:348#0/1:unique"<<endl;
			}

			vector<GeneGenomeLoc*> fitted_genome_locs1;

			vector<GeneGenomeLoc*> fitted_genome_locs2;
			
			bool find_loc1 = find_corresponding_fusion_gene(merged_tag_aligns_pe.first.front(), genes_genome_locs, max_gene_len, fitted_genome_locs1);

			bool find_loc2 = find_corresponding_fusion_gene(merged_tag_aligns_pe.second.front(), genes_genome_locs, max_gene_len, fitted_genome_locs2);

			vector<GeneGenomeLoc*>::iterator genome_loc_iter1, genome_loc_iter2;

			if (merged_tag_aligns_pe.first.size() && merged_tag_aligns_pe.first.front()->tag_name == "HWUSI-EAS1531_0001:4:100:1066:348#0/1")
			{
				cout << "HWUSI-EAS1531_0001:4:100:1066:348#0/1:"<<find_loc1 <<'\t' <<  find_loc2<<endl;
			}

			if (find_loc1 && find_loc2)
			{
				vector<pair<GeneGenomeLoc*, GeneGenomeLoc*> > valid_genome_loc_pairs;

				for (genome_loc_iter1 = fitted_genome_locs1.begin(); genome_loc_iter1 != fitted_genome_locs1.end(); ++genome_loc_iter1)
				{
					for (genome_loc_iter2 = fitted_genome_locs2.begin(); genome_loc_iter2 != fitted_genome_locs2.end(); ++genome_loc_iter2)
					{
						if ((*genome_loc_iter1)->fusiongene_ptr == (*genome_loc_iter2)->fusiongene_ptr)//both ends fall in same fusion gene!
						{
							if (merged_tag_aligns_pe.first.size() && merged_tag_aligns_pe.first.front()->tag_name == "HWUSI-EAS1531_0001:4:100:1066:348#0/1")
							{
								cout << "HWUSI-EAS1531_0001:4:100:1066:348#0/1:paired matched:"<<find_loc1 <<'\t' <<  find_loc2<<endl;
							}

							valid_genome_loc_pairs.push_back(make_pair(*genome_loc_iter1, *genome_loc_iter2));
						}
					}
				}

				if (!valid_genome_loc_pairs.empty())
				{
					vector<pair<GeneGenomeLoc*, GeneGenomeLoc*> >::iterator genome_loc_pair_iter;

					for (genome_loc_pair_iter = valid_genome_loc_pairs.begin(); genome_loc_pair_iter != valid_genome_loc_pairs.end(); ++genome_loc_pair_iter)
					{
						SamRec samrecdup1 = *(merged_tag_aligns_pe.first.front());

						convert_alignment_coordinate(genome_loc_pair_iter->first, samrecdup1);

						SamRec samrecdup2 = *(merged_tag_aligns_pe.second.front());

						convert_alignment_coordinate(genome_loc_pair_iter->second, samrecdup2);

						set_mate_pair_distance(samrecdup1, samrecdup2);

						ofs << samrecdup1.tostring(1, 1) << endl;

						ofs << samrecdup2.tostring(1, 1) << endl;
					}
				}
				else
				{
				}
			}
			else if (find_loc1)
			{
				for (genome_loc_iter1 = fitted_genome_locs1.begin(); genome_loc_iter1 != fitted_genome_locs1.end(); ++genome_loc_iter1)
				{
					SamRec samrecdup1 = *(merged_tag_aligns_pe.first.front());

					convert_alignment_coordinate((*genome_loc_iter1), samrecdup1);

					ofs << samrecdup1.tostring(1, 1) << endl;
				}				
			}
			else if (find_loc2)
			{
				for (genome_loc_iter1 = fitted_genome_locs2.begin(); genome_loc_iter1 != fitted_genome_locs2.end(); ++genome_loc_iter1)
				{
					SamRec samrecdup1 = *(merged_tag_aligns_pe.second.front());

					convert_alignment_coordinate((*genome_loc_iter1), samrecdup1);

					ofs << samrecdup1.tostring(1, 1) << endl;
				}
			}
			else
			{
			}
		}
		else if (merged_tag_aligns_pe.first.size() == 1)
		{
			vector<GeneGenomeLoc*> fitted_genome_locs1;
			
			bool find_loc1 = find_corresponding_fusion_gene(merged_tag_aligns_pe.first.front(), genes_genome_locs, max_gene_len, fitted_genome_locs1);

			vector<GeneGenomeLoc*>::iterator genome_loc_iter1;

			if (find_loc1)
			{
				for (genome_loc_iter1 = fitted_genome_locs1.begin(); genome_loc_iter1 != fitted_genome_locs1.end(); ++genome_loc_iter1)
				{
					SamRec samrecdup1 = *(merged_tag_aligns_pe.first.front());

					convert_alignment_coordinate((*genome_loc_iter1), samrecdup1);

					ofs << samrecdup1.tostring(1, 1) << endl;
				}
			}
		}
		else if (merged_tag_aligns_pe.second.size() == 1)
		{
			vector<GeneGenomeLoc*> fitted_genome_locs1;
			
			bool find_loc1 = find_corresponding_fusion_gene(merged_tag_aligns_pe.second.front(), genes_genome_locs, max_gene_len, fitted_genome_locs1);

			vector<GeneGenomeLoc*>::iterator genome_loc_iter1;

			if (find_loc1)
			{
				for (genome_loc_iter1 = fitted_genome_locs1.begin(); genome_loc_iter1 != fitted_genome_locs1.end(); ++genome_loc_iter1)
				{
					SamRec samrecdup1 = *(merged_tag_aligns_pe.second.front());

					convert_alignment_coordinate((*genome_loc_iter1), samrecdup1);

					ofs << samrecdup1.tostring(1, 1) << endl;
				}
			}
		}
		else
		{
			cout << "what"<<endl;
		}
	}
}

void fusiongene2chromseq(map<string, FusionGene>& fusiongenes, string chrom_dir, char* outout_junc_seq_file)
{
	ofstream ofs_junc(outout_junc_seq_file);

	map<string, FusionGene>::iterator fusiongene_iter;

	map<string, string> loaded_chromos;

	for (fusiongene_iter = fusiongenes.begin(); fusiongene_iter != fusiongenes.end(); ++fusiongene_iter)
	{
		string chr1, chr2;

		chr1 = fusiongene_iter->second.gene1ptr->chr;

		chr2 = fusiongene_iter->second.gene2ptr->chr;

		//if (chr1 == "")
		//	cout << fusiongene_iter->second.gene1ptr->transcripts.front().exons.front().line<<endl;

		//if (chr2 == "")
		//	cout << fusiongene_iter->second.gene2ptr->transcripts.front().exons.front().line<<endl;

		if (loaded_chromos.find(chr1) == loaded_chromos.end())
		{
			string chromfile = chrom_dir + "/" + chr1 + ".fa";

			readchrom(chromfile.c_str(), loaded_chromos[chr1]);
		}

		if (loaded_chromos.find(chr2) == loaded_chromos.end())
		{
			string chromfile = chrom_dir + "/" + chr2 + ".fa";

			readchrom(chromfile.c_str(), loaded_chromos[chr2]);
		}

		string& chr1seq = loaded_chromos[chr1];

		string& chr2seq = loaded_chromos[chr2];

		char strand1 = fusiongene_iter->second.strand1;

		char strand2 = fusiongene_iter->second.strand2;

		//
		size_t gene1st = fusiongene_iter->second.gene1ptr->start - fusiongene_iter->second.gene1_left_ext;

		size_t gene1end = fusiongene_iter->second.gene1ptr->end + fusiongene_iter->second.gene1_right_ext;

		size_t length1 = gene1end - gene1st + 1;

		//
		size_t gene2st = fusiongene_iter->second.gene2ptr->start - fusiongene_iter->second.gene2_left_ext;

		size_t gene2end = fusiongene_iter->second.gene2ptr->end + fusiongene_iter->second.gene2_right_ext;

		size_t length2 = gene2end - gene2st + 1;

		string chr1substr, chr2substr;

		if (gene1st + length1 > chr1seq.length())
		{
			cerr << "doner anchor exceed chromosome boundary:" << gene1st + length1 << "\tVS\t"<<chr1seq.length() << endl;

			continue;
		}

		chr1substr = chr1seq.substr(gene1st, length1);

		if (strand1 == '-')
		{
			chr1substr = chr1seq.substr(gene1st - 2, length1);

			chr1substr = revcomp(chr1substr);
		}
		else
		{
		}

		if (gene2st + length2 > chr2seq.length())
		{
			cerr << "acceptor anchor exceed chromosome boundary:" << gene2st + length2 << "\tVS\t"<<chr2seq.length() << endl;

			continue;
		}

		chr2substr = chr2seq.substr(gene2st, length2);

		if (strand2 == '-')
		{
			chr2substr = chr2seq.substr(gene2st - 2, length2);

			chr2substr = revcomp(chr2substr);
		}

		string combstr = chr1substr + chr2substr;

		ofs_junc <<">"<< fusiongene_iter->second.fusiongenename <<endl; 

		size_t chunks = combstr.length() / 50;

		for (size_t i = 0; i < chunks + 1; ++i)
		{
			ofs_junc << combstr.substr(i * 50, 50) << endl;

			if (i * 50 + 50 >= combstr.length())
				break;
		}
	}
}

void gene2fusiongene(map<string, FusionGene>& fusiongenes, char* output_gene_file)
{
	ofstream ofs_gene(output_gene_file);

	map<string, FusionGene>::iterator fusiongene_iter;

	for (fusiongene_iter = fusiongenes.begin(); fusiongene_iter != fusiongenes.end(); ++fusiongene_iter)
	{
		vector<Transcript>::iterator transcript_iter;

		for (transcript_iter = fusiongene_iter->second.gene1ptr->transcripts.begin(); transcript_iter != fusiongene_iter->second.gene1ptr->transcripts.end(); ++transcript_iter)
		{
			vector<Exon>::iterator exon_iter;

			for (exon_iter = transcript_iter->exons.begin(); exon_iter != transcript_iter->exons.end(); ++exon_iter)
			{
				size_t newstart = convert_point_coordinate(exon_iter->start, fusiongene_iter->second.gene1_start_ext, fusiongene_iter->second.gene1_end_ext, 
					fusiongene_iter->second.a, fusiongene_iter->second.b, fusiongene_iter->second.strand1);

				size_t newend = convert_point_coordinate(exon_iter->end, fusiongene_iter->second.gene1_start_ext, fusiongene_iter->second.gene1_end_ext, 
					fusiongene_iter->second.a, fusiongene_iter->second.b, fusiongene_iter->second.strand1);

				if (fusiongene_iter->second.strand1 == '+')
					ofs_gene << exon_iter->tostring(fusiongene_iter->second.fusiongenename, newstart, newend, fusiongene_iter->second.strand1) << endl;
				else
					ofs_gene << exon_iter->tostring(fusiongene_iter->second.fusiongenename, newend, newstart, fusiongene_iter->second.strand1) << endl;
			}
		}

		for (transcript_iter = fusiongene_iter->second.gene2ptr->transcripts.begin(); transcript_iter != fusiongene_iter->second.gene2ptr->transcripts.end(); ++transcript_iter)
		{
			vector<Exon>::iterator exon_iter;

			for (exon_iter = transcript_iter->exons.begin(); exon_iter != transcript_iter->exons.end(); ++exon_iter)
			{
				size_t newstart = convert_point_coordinate(exon_iter->start, fusiongene_iter->second.gene2_start_ext, fusiongene_iter->second.gene2_end_ext, 
					fusiongene_iter->second.c, fusiongene_iter->second.d, fusiongene_iter->second.strand2);

				size_t newend = convert_point_coordinate(exon_iter->end, fusiongene_iter->second.gene2_start_ext, fusiongene_iter->second.gene2_end_ext, 
					fusiongene_iter->second.c, fusiongene_iter->second.d, fusiongene_iter->second.strand2);

				if (fusiongene_iter->second.strand2 == '+')
					ofs_gene << exon_iter->tostring(fusiongene_iter->second.fusiongenename, newstart, newend, fusiongene_iter->second.strand2) << endl;
				else
					ofs_gene << exon_iter->tostring(fusiongene_iter->second.fusiongenename, newend, newstart, fusiongene_iter->second.strand2) << endl;
			}
		}
	}
}

void fusiongene2profile(map<string, FusionGene>& fusiongenes, char* fusion_gene_profile)
{
	ofstream ofs_gene(fusion_gene_profile);

	map<string, FusionGene>::iterator fusiongene_iter;

	for (fusiongene_iter = fusiongenes.begin(); fusiongene_iter != fusiongenes.end(); ++fusiongene_iter)
	{
		vector<Fusion*>::iterator fusion_iter;
		
		ofs_gene << fusiongene_iter->second.fusiongenename<<endl;

		size_t x_prime, y_prime, x_prime2, y_prime2;

		x_prime = fusiongene_iter->second.a;

		y_prime = fusiongene_iter->second.b;

		x_prime2 = fusiongene_iter->second.c;

		y_prime2 = fusiongene_iter->second.d;

		for (fusion_iter = fusiongene_iter->second.fusions.begin(); fusion_iter != fusiongene_iter->second.fusions.end(); ++fusion_iter)
		{
			size_t newstart, newend;

			newstart = convert_point_coordinate((*fusion_iter)->start, fusiongene_iter->second.gene1_start_ext, fusiongene_iter->second.gene1_end_ext, x_prime, y_prime, fusiongene_iter->second.strand1);

			newend = convert_point_coordinate((*fusion_iter)->end, fusiongene_iter->second.gene2_start_ext, fusiongene_iter->second.gene2_end_ext, x_prime2, y_prime2, fusiongene_iter->second.strand2);

			ofs_gene <<(*fusion_iter)->fusionline << endl;

			ofs_gene << fusiongene_iter->second.fusiongenename <<'\t' << newstart << '\t' << newend<< '\t' << "++"<<endl;
		}

		ofs_gene << endl;
	}
}

int main(int argc, char** argv)
{
	if (argc < 3)
	{
		cout << "fusionjunc extend_length sam_file output_sam_file gene_file output_gene_file chromdir outout_junc_seq_file" <<endl;
		exit(0);
	}

	char* fusionjunc = argv[1];

	int length = atoi(argv[2]);

	char* sam_file = argv[3];

	char* output_sam_file = argv[4];

	char* gene_file = argv[5];

	char* output_gene_file = argv[6];

	char* chromdir =argv[7];

	char* outout_junc_seq_file = argv[8];

	char* fusion_gene_profile = argv[9];

	vector<Fusion> fusions;

	cout << "read_fusion"<<endl;

	read_fusion(fusionjunc, fusions);

	vector<Exon> exons;

	cout << "read_annotation"<<endl;

	read_annotation(gene_file, exons, output_gene_file);

	cout << "exons:"<< exons.size()<<endl;

	map<string, Gene> genes;

	cout << "exon2gene"<<endl;

	map<string, vector<Gene* > > sorted_genes;

	size_t max_gene_len = exon2gene(exons, genes, length, output_gene_file, sorted_genes);

	map<string, FusionGene> fusiongenes;

	vector<pair<string, size_t> > syn_chrom_name_length;

	cout << "locate_corresponding_gene"<<endl;

	locate_corresponding_gene(fusions, genes, fusiongenes, syn_chrom_name_length, sorted_genes, length);

	cout << "fusiongene2chromseq"<<endl;

	fusiongene2chromseq(fusiongenes, chromdir, outout_junc_seq_file);

	cout << "gene2fusiongene"<<endl;

	gene2fusiongene(fusiongenes, output_gene_file);

	cout << "fusiongene2profile"<<endl;

	fusiongene2profile(fusiongenes, fusion_gene_profile);

	map<string, vector<GeneGenomeLoc> > genes_genome_locs;

	cout << "gene2genome_coordinate"<<endl;

	gene2genome_coordinate(fusiongenes, genes_genome_locs);

	cout << "extract_and_convert_sam"<<endl;

	extract_and_convert_sam(genes_genome_locs, sam_file, output_sam_file, max_gene_len, syn_chrom_name_length);
}