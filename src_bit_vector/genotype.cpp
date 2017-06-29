#include <bits/stdc++.h>
#include "genotype.h"

using namespace std;

void genotype::read_genotype_naive (char *filename)
{
	FILE* fp;
	fp= fopen(filename,"r");
	int j=0;
	int i=0;
	char ch;
	vector <bool> m;
	vector <bool> l;
	int sum=0;
	int temp1,temp2;
	int rd = fscanf(fp,"%d %d\n",&temp1,&temp2);
    do
    {
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
	assert(temp1==Nsnp);
	assert(temp2==Nindv);
}

void genotype::read_genotype_mailman (char *filename)
{
	FILE* fp;
	fp= fopen(filename,"r");
	int j=0;
	int i=0;
	char ch;
	// Calculating the sizes and other stuff for genotype matrix
	int rd = fscanf(fp,"%d %d\n",&Nsnp,&Nindv);
	segment_size_hori = ceil(log2(Nindv));
	segment_size_ver = ceil(log2(Nsnp));
	Nsegments_hori = ceil(Nsnp*1.0/(segment_size_hori*1.0));
	Nsegments_ver = ceil(Nindv*1.0/(segment_size_ver*1.0));
	p_msb.resize(Nsegments_hori,std::vector< std::vector<bool> >(Nindv));
	p_lsb.resize(Nsegments_hori,std::vector< std::vector<bool> >(Nindv));
	q_msb.resize(Nsegments_ver,std::vector< std::vector<bool> >(Nsnp));
	q_lsb.resize(Nsegments_ver,std::vector< std::vector<bool> >(Nsnp));
	int sum=0;

    do
    {
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
				p_msb[horiz_seg_no][j].push_back(false);
				p_lsb[horiz_seg_no][j].push_back(false); 
				q_msb[ver_seg_no][i].push_back(false); 				
				q_lsb[ver_seg_no][i].push_back(false); 		
			}
			else if(val==1){
				sum+=1;
				p_msb[horiz_seg_no][j].push_back(false);
				p_lsb[horiz_seg_no][j].push_back(true);
				q_msb[ver_seg_no][i].push_back(false); 
				q_lsb[ver_seg_no][i].push_back(true);
			}
			else if(val==2){
				sum+=2;
				p_msb[horiz_seg_no][j].push_back(true);
				p_lsb[horiz_seg_no][j].push_back(false);
				q_msb[ver_seg_no][i].push_back(true);
				q_lsb[ver_seg_no][i].push_back(false);
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

float genotype::get_geno(int snpindex,int indvindex)
{
	float m = msb[snpindex][indvindex];
	float l = lsb[snpindex][indvindex];
	if ((m*2+l)==3.0)
		return 0.0;
	else
		return ( (m*2.0+l) - get_col_mean(snpindex));
}

vector<float> genotype::get_geno_row(int snpindex)
{
	vector<float> v;
	float mean = get_col_mean(snpindex);
	for(int i=0;i<msb[snpindex].size();i++)
	{
		float m = msb[snpindex][i];
		float l = lsb[snpindex][i];
		v.push_back((m*2.0+l)-mean);
	}
	return v;
}

vector<float> genotype::get_geno_col(int indvindex)
{
	vector<float> v;
	for(int i=0;i<Nsnp;i++)
	{
		float m = msb[i][indvindex];
		float l = lsb[i][indvindex];
		v.push_back((m*2.0+l)-get_col_mean(i));
	}
	return v;
}

float genotype::get_col_mean(int snpindex)
{
	float temp = columnsum[snpindex]*1.0 / Nindv ;
	return temp;
}

bool genotype::is_observed(int snpindex,int indvindex)
{
	if(msb[snpindex][indvindex] && lsb[snpindex][indvindex])
		return false;
	else 
		return true;
}