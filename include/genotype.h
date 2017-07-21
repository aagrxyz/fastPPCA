#ifndef GENOTYPE_H
#define GENOTYPE_H
#include <bits/stdc++.h>
#include "storage.h"

class genotype {
	std::vector< std::vector <bool> > msb;
	std::vector< std::vector <bool> > lsb;
	std::vector<int> columnsum;
	std::vector<double> columnmeans;
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
		
		void init_means(bool is_missing);
		void read_genotype_mailman (std::string filename,bool allow_missing);
		void read_genotype_eff (std::string filename,bool allow_missing);		
		void read_genotype_naive (std::string filename,bool allow_missing);


		double get_geno(int snpindex,int indvindex,bool var_normalize);


		double get_col_mean(int snpindex);
		void update_col_mean(int snpindex,double value);
		double get_col_sum(int snpindex);		
		double get_col_std(int snpindex);		

};

#endif