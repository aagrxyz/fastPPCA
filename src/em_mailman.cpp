#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include "genotype.h"
#include "mailman.h"
#include "helper.h"

using namespace Eigen;
using namespace std;
clock_t total_begin = clock();

int MAX_ITER;
genotype g;
int k,p,n;

MatrixXf c; //(p,k)
MatrixXf x; //(k,n)
MatrixXf v; //(p,k)

options command_line_opts;

bool debug = false;
bool check_accuracy = false;

MatrixXf get_evec(MatrixXf &c)
{
	JacobiSVD<MatrixXf> svd(c, ComputeThinU | ComputeThinV);
	MatrixXf c_orth(k,p);
	MatrixXf data(k,n);
	c_orth = (svd.matrixU()).transpose();
	for(int n_iter=0;n_iter<n;n_iter++)
	{
		for(int k_iter=0;k_iter<k;k_iter++)
		{
			float res=0;
			for(int p_iter=0;p_iter<p;p_iter++)
				res+= c_orth(k_iter,p_iter)*(g.get_geno(p_iter,n_iter));
			data(k_iter,n_iter)=res;
		}
	}
	MatrixXf means(k,1);
	for(int i=0;i<k;i++)
	{
		float sum=0.0;
		for(int j=0;j<n;j++)
			sum+=data(i,j);
		means(i,0)=sum/(n*1.0);
	}
	data = data - (means*(MatrixXf::Constant(1,n,1)));
	MatrixXf cov(k,k);
	cov = data*(data.transpose())*(1.0/(n));
	JacobiSVD<MatrixXf> svd_cov(cov, ComputeThinU | ComputeThinV);
	MatrixXf to_return(p,k);
	to_return =(c_orth.transpose())*svd_cov.matrixU() ;
	return to_return;
}


float get_accuracy(MatrixXf &u)
{
	
	MatrixXf temp(k,k);
	temp = (u.transpose()) * v ;
	float accuracy = 0.0;
	for(int j=0;j<k;j++)
	{
		float sum=0.0;
		for(int i=0;i<k;i++)
			sum += temp(i,j)*temp(i,j);
		accuracy += sqrt(sum);
	}
	return (accuracy/k);
}

MatrixXf get_reference_evec()
{

	MatrixXf y_m(p,n);
	for(int i=0;i<p;i++)
	{
		for(int j=0;j<n;j++)
			y_m(i,j) = g.get_geno(i,j);
	}
	MatrixXf cov(p,p);
	printf("Calculating covariance\n");
	cov = y_m*(y_m.transpose())*(1.0/n);
	printf("Calculating SVD\n");
	JacobiSVD<MatrixXf> svd_cov(cov, ComputeThinU | ComputeThinV);
	MatrixXf to_return(p,k);
	MatrixXf U(p,k);
	U = svd_cov.matrixU();
	for(int i=0;i<p;i++)
	{
		for(int j=0;j<k;j++)
			to_return(i,j) = U(i,j);
	}

	return to_return;
}
int main(int argc, char const *argv[])
{
	clock_t io_begin = clock();
	
	parse_args(argc,argv);
	g.read_genotype_mailman((char*)command_line_opts.GENOTYPE_FILE_PATH);	
	MAX_ITER =  command_line_opts.max_iterations ; 
	k = command_line_opts.num_of_evec ;
	debug = command_line_opts.debugmode ;
	check_accuracy = command_line_opts.getaccuracy;
	p = g.Nsnp;
	n = g.Nindv;
	srand((unsigned int) time(0));
	
	clock_t io_end = clock();

	c = MatrixXf::Random(p,k);
	
	ofstream c_file;
	if(debug){
		c_file.open("cvals_orig.txt");
		c_file<<c<<endl;
		c_file.close();
		printf("Read Matrix\n");
	}

	// TODO: CHECK ACCURACY with Mailman
	
	check_accuracy = false;

	if(check_accuracy){
		v = get_reference_evec(); 
	 	printf("Computed Reference Eigen Vectors\n");	
	}
	double *yint_m = new double[(int)pow(2,g.segment_size_hori)];
	double *yint_e = new double[(int)pow(2,g.segment_size_ver)];

	cout<<"Running on Dataset of "<<g.Nsnp<<" SNPs and "<<g.Nindv<<" Individuals"<<endl;
	if(check_accuracy)
		cout<<endl<<"Iterations vs accuracy"<<endl;
	
	clock_t it_begin = clock();
	for(int i=0;i<MAX_ITER;i++)
	{
		if(debug)
			printf("Iteration %d -- E Step\n",i);

		MatrixXf c_temp(p,k);
		c_temp = (((c.transpose()*c).inverse())*(c.transpose())).transpose();

		for(int k_iter=0;k_iter<k;k_iter++)
		{
			double *y_msb = new double[g.segment_size_ver];
			double *y_lsb = new double[g.segment_size_ver];
			int seg_iter;
			MatrixXf c_temp_col(p,1);
			c_temp_col = c_temp.col(k_iter);

			for(seg_iter=0; seg_iter < g.Nsegments_ver; seg_iter++)
			{
				if(seg_iter==g.Nsegments_ver-1)
				{
					if(g.Nindv%g.segment_size_ver!=0)
					{
						double *y_msb_final = new double[g.Nindv%g.segment_size_ver];
						double *y_lsb_final = new double[g.Nindv%g.segment_size_ver];
						mailman::fastmultiply(g.Nindv%g.segment_size_ver,g.Nsnp,g.q_msb[seg_iter],c_temp_col,yint_e,y_msb_final);
						mailman::fastmultiply(g.Nindv%g.segment_size_ver,g.Nsnp,g.q_lsb[seg_iter],c_temp_col,yint_e,y_lsb_final);
						for(int n_iter=seg_iter*g.segment_size_ver ; n_iter<seg_iter*g.segment_size_ver + g.Nindv%g.segment_size_ver  && n_iter<g.Nindv ; n_iter++)
							x(k_iter,n_iter) = 2*y_msb_final[n_iter-(seg_iter*g.segment_size_ver)] + y_lsb_final[n_iter-(seg_iter*g.segment_size_ver)];
						break;
					}

				}

				mailman::fastmultiply(g.segment_size_ver,g.Nsnp,g.q_msb[seg_iter],c_temp_col,yint_e,y_msb);
				mailman::fastmultiply(g.segment_size_ver,g.Nsnp,g.q_lsb[seg_iter],c_temp_col,yint_e,y_lsb);
				int n_base = seg_iter*g.segment_size_ver; 
				for(int n_iter=n_base; (n_iter<n_base+g.segment_size_ver) && (n_iter<g.Nindv) ; n_iter++ ) 
					x(k_iter,n_iter) = 2*y_msb[n_iter-n_base] + y_lsb[n_iter-n_base];

			}
		}
		double *sums_elements = new double[k];
    	memset (sums_elements, 0, k * sizeof(int));

		for(int k_iter=0;k_iter<k;k_iter++)
		{
			double sum_to_calc=0.0;
			for(int p_iter=0;p_iter<p;p_iter++)
				sum_to_calc += g.get_col_mean(p_iter)*c_temp(p_iter,k_iter);
			sums_elements[k_iter] = sum_to_calc;
		}

		MatrixXf to_subtract_e(n,k);
		for(int k_iter=0;k_iter<k;k_iter++)
		{
			for(int n_iter=0;n_iter<n;n_iter++)
				x(k_iter,n_iter) = x(k_iter,n_iter) - sums_elements[k_iter];
		}


		if(debug)
			printf("Iteration %d -- M Step\n",i);

		MatrixXf x_temp(n,k);
		x_temp = (x.transpose()) * ((x*(x.transpose())).inverse());

		double *sum_x_temp = new double[k];
		for(int k_iter=0;k_iter<k;k_iter++)
		{
			double *y_msb = new double[g.segment_size_hori];
			double *y_lsb = new double[g.segment_size_hori];
			int seg_iter;
			MatrixXf x_temp_col(n,1);
			x_temp_col = x_temp.col(k_iter);
			sum_x_temp[k_iter]=x_temp_col.sum();

			for(seg_iter=0;seg_iter<g.Nsegments_hori;seg_iter++)
			{
				if(seg_iter==g.Nsegments_hori-1)
				{
					if(g.Nsnp%g.segment_size_hori!=0)
					{
						double *y_msb_final = new double[g.Nsnp%g.segment_size_hori];
						double *y_lsb_final = new double[g.Nsnp%g.segment_size_hori];
						mailman::fastmultiply(g.Nsnp%g.segment_size_hori,g.Nindv,g.p_msb[seg_iter],x_temp_col,yint_m,y_msb_final);
						mailman::fastmultiply(g.Nsnp%g.segment_size_hori,g.Nindv,g.p_lsb[seg_iter],x_temp_col,yint_m,y_lsb_final);
						for(int p_iter=seg_iter*g.segment_size_hori;p_iter<seg_iter*g.segment_size_hori + g.Nsnp%g.segment_size_hori  && p_iter<g.Nsnp;p_iter++)
							c(p_iter,k_iter) = 2*y_msb_final[p_iter-(seg_iter*g.segment_size_hori)] + y_lsb_final[p_iter-(seg_iter*g.segment_size_hori)];
						break;
					}
				}

				mailman::fastmultiply(g.segment_size_hori,g.Nindv,g.p_msb[seg_iter],x_temp_col,yint_m,y_msb);
				mailman::fastmultiply(g.segment_size_hori,g.Nindv,g.p_lsb[seg_iter],x_temp_col,yint_m,y_lsb);
				int p_base = seg_iter*g.segment_size_hori; 
				for(int p_iter=p_base; (p_iter<p_base+g.segment_size_hori) && (p_iter<g.Nsnp) ; p_iter++ ) 
					c(p_iter,k_iter) = 2*y_msb[p_iter-p_base] + y_lsb[p_iter-p_base];
			}			
		}

		for(int p_iter=0;p_iter<p;p_iter++)
		{
			for(int k_iter=0;k_iter<k;k_iter++)
				c(p_iter,k_iter) = c(p_iter,k_iter) - (g.get_col_mean(p_iter)*sum_x_temp[k_iter]);
		}
		if(check_accuracy){
			MatrixXf eigenvectors(p,k);
			eigenvectors = get_evec(c);
			cout<<"Iteration "<<i<<" "<<get_accuracy(eigenvectors)<<endl;	
		}
	}
	clock_t it_end = clock();
	c_file.open("cvals_end.txt");
	c_file<<c<<endl;
	c_file.close();
	ofstream x_file;
	x_file.open("xvals.txt");
	x_file<<x<<endl;
	x_file.close();
	clock_t total_end = clock();
	double io_time = double(io_end - io_begin) / CLOCKS_PER_SEC;
	double avg_it_time = double(it_end - it_begin) / (MAX_ITER * 1.0 * CLOCKS_PER_SEC);
	double total_time = double(total_end - total_begin) / CLOCKS_PER_SEC;
	cout<<"IO Time:  "<< io_time << "\nAVG Iteration Time:  "<<avg_it_time<<"\nTotal runtime:   "<<total_time<<endl;
	return 0;
}
