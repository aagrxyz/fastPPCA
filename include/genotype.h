#ifndef GENOTYPE_H
#define GENOTYPE_H
#include <bits/stdc++.h>

class genotype {
	std::vector< std::vector <bool> > msb;
	std::vector< std::vector <bool> > lsb;
	std::vector<int> columnsum;
	public:
		
		int Nsnp,Nindv,Nsegments_hori,segment_size_hori,segment_size_ver,Nsegments_ver;
		std::vector< std::vector<int> > p;
		std::vector< std::vector<int> > q;
		
		void read_genotype_mailman (std::string filename);
		void read_genotype_naive (std::string filename);
		float get_geno(int snpindex,int indvindex);
		std::vector<float> get_geno_row(int snpindex);
		std::vector<float> get_geno_col(int indvindex);
		float get_col_mean(int snpindex);

		float get_col_std(int snpindex);

		bool is_observed(int snpindex,int indvindex);
		

};

#endif