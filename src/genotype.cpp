#include <bits/stdc++.h>
#include "genotype.h"

using namespace std;

void genotype::read_genotype_naive (std::string filename){
	FILE* fp;
	fp= fopen(filename.c_str(),"r");
	int j=0;
	int i=0;
	char ch;
	vector <bool> m;
	vector <bool> l;
	int sum=0;
	int snp_file,indv_file;
	int rd = fscanf(fp,"%d %d\n",&snp_file,&indv_file);
    do{
		int rd = fscanf(fp,"%c",&ch);
		if(ch=='\n'){
			i++;
			msb.push_back(m);
			lsb.push_back(l);
			m.clear();
			l.clear();
			columnsum.push_back(sum);
			sum=0;
			j=0;
		}
		else{
			int val = int(ch-'0');
			if(val==0){
				l.push_back(false);
				m.push_back(false);
			}
			else if(val==1){
				sum+=1;
				l.push_back(true);
				m.push_back(false);
			}
			else if(val==2){
				sum+=2;
				l.push_back(false);
				m.push_back(true);
			}
			else{
				l.push_back(true);
				m.push_back(true);
			}
			j++;
		}
	}while(!feof(fp));
	Nsnp = msb.size()-1;
	Nindv = msb[0].size();
	assert(snp_file==Nsnp);
	assert(indv_file==Nindv);
}

void genotype::read_genotype_mailman (std::string filename){
	FILE* fp;
	fp= fopen(filename.c_str(),"r");
	int j=0;
	int i=0;
	char ch;
	// Calculating the sizes and other stuff for genotype matrix
	int rd = fscanf(fp,"%d %d\n",&Nsnp,&Nindv);
	segment_size_hori = ceil(log(Nindv)/log(3));
	segment_size_ver = ceil(log(Nsnp)/log(3));
	Nsegments_hori = ceil(Nsnp*1.0/(segment_size_hori*1.0));
	Nsegments_ver = ceil(Nindv*1.0/(segment_size_ver*1.0));
	p.resize(Nsegments_hori,std::vector<int>(Nindv));
	q.resize(Nsegments_ver,std::vector<int>(Nsnp));
	int sum=0;
    do{
		int rd = fscanf(fp,"%c",&ch);
		if(ch=='\n'){
			i++;
			columnsum.push_back(sum);
			sum=0;
			j=0;
		}
		else{
			int val = int(ch-'0');
			int horiz_seg_no = i/segment_size_hori ;
			int ver_seg_no = j/segment_size_ver ;
			if(val==0){
				p[horiz_seg_no][j] = 3 * p[horiz_seg_no][j] ;
				q[ver_seg_no][i] = 3 * q[ver_seg_no][i];
			}
			else if(val==1){
				sum+=1;

				p[horiz_seg_no][j] = 3 * p[horiz_seg_no][j]  + 1;
				q[ver_seg_no][i] = 3*q[ver_seg_no][i] + 1;
			}
			else if(val==2){
				sum+=2;
				p[horiz_seg_no][j] = 3 * p[horiz_seg_no][j] + 2 ;
				q[ver_seg_no][i] = 3*q[ver_seg_no][i] + 2;
			}
			else{
				cout<<"Invalid entry in Genotype Matrix"<<endl;
				exit(-1);
			}
			j++;
		}
	}while(!feof(fp));
	i--;
}

float genotype::get_geno(int snpindex,int indvindex){
	float m = msb[snpindex][indvindex];
	float l = lsb[snpindex][indvindex];
	if ((m*2+l)==3.0)
		return 0.0;
	else
		return ( (m*2.0+l) - get_col_mean(snpindex));
}

vector<float> genotype::get_geno_row(int snpindex){
	vector<float> v;
	float mean = get_col_mean(snpindex);
	for(int i=0;i<msb[snpindex].size();i++){
		float m = msb[snpindex][i];
		float l = lsb[snpindex][i];
		v.push_back((m*2.0+l)-mean);
	}
	return v;
}

vector<float> genotype::get_geno_col(int indvindex){
	vector<float> v;
	for(int i=0;i<Nsnp;i++){
		float m = msb[i][indvindex];
		float l = lsb[i][indvindex];
		v.push_back((m*2.0+l)-get_col_mean(i));
	}
	return v;
}

float genotype::get_col_mean(int snpindex){
	float temp = columnsum[snpindex]*1.0 / Nindv ;
	return temp;
}

float genotype::get_col_std(int snpindex){
	float p_i = get_col_mean(snpindex);
	float temp = sqrt(p_i*(1-(0.5*p_i))) ; 
	return temp;
}


bool genotype::is_observed(int snpindex,int indvindex){
	if(msb[snpindex][indvindex] && lsb[snpindex][indvindex])
		return false;
	else 
		return true;
}

