#ifndef GENOTYPE_H
#define GENOTYPE_H
#include <bits/stdc++.h>
#include "storage.h"

class genotype {
	std::vector< std::vector <bool> > msb;
	std::vector< std::vector <bool> > lsb;
	std::vector<int> columnsum;
	public:
		
		int Nsnp,Nindv,Nsegments_hori,segment_size_hori,segment_size_ver,Nsegments_ver;
		int Nbits_hori,Nbits_ver;
		int Nelements_hori,Nelements_ver;
		std::vector< std::vector<int> > p;
		std::vector< std::vector<int> > q;

		std::vector< std::vector<unsigned> > p_eff;
		std::vector< std::vector<unsigned> > q_eff;

		std::vector< std::vector<int> > not_O_j;
		std::vector< std::vector<int> > not_O_i;
		
		void read_genotype_mailman_missing (std::string filename);
		void read_genotype_mailman (std::string filename);
		void read_genotype_eff (std::string filename);		
		void read_genotype_naive (std::string filename);

		float get_geno(int snpindex,int indvindex,bool var_normalize);
		std::vector<float> get_geno_row(int snpindex);
		std::vector<float> get_geno_col(int indvindex);


		float get_col_mean(int snpindex);
		float get_col_sum(int snpindex);		
		float get_col_std(int snpindex);

		bool is_observed(int snpindex,int indvindex);
		

};

#endif