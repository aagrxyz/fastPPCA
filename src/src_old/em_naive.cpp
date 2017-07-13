#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/QR>
#include "genotype.h"
#include "helper.h"

using namespace Eigen;
using namespace std;
clock_t total_begin = clock();

int MAX_ITER;
genotype g;
int k,p,n;
int k_orig;

MatrixXd c; //(p,k)
MatrixXd x; //(k,n)
MatrixXd v; //(p,k)

options command_line_opts;

bool debug = false;
bool check_accuracy = false;
bool var_normalize=false;
int accelerated_em=0;
double convergence_limit;


void multiply_y_pre(MatrixXd &op, int Ncol_op ,MatrixXd &res){

	for(int p_iter=0;p_iter<p;p_iter++)
	{
		for(int k_iter=0;k_iter<Ncol_op;k_iter++)
		{
			float temp=0;
			for(int n_iter=0;n_iter<n;n_iter++)
				temp+= g.get_geno(p_iter,n_iter,var_normalize)*op(n_iter,k_iter);
			res(p_iter,k_iter)=temp;
		}
	}
}

void multiply_y_post(MatrixXd &op, int Nrow_op ,MatrixXd &res){

	for(int n_iter=0;n_iter<n;n_iter++)
	{
		for(int k_iter=0;k_iter<Nrow_op;k_iter++)
		{
			float temp=0;
			for(int p_iter=0;p_iter<p;p_iter++)
				temp+= op(k_iter,p_iter)*(g.get_geno(p_iter,n_iter,var_normalize));
			res(k_iter,n_iter)=temp;
		}
	}
}



MatrixXd get_evec(MatrixXd &c){
	JacobiSVD<MatrixXd> svd(c, ComputeThinU | ComputeThinV);
	MatrixXd c_orth(k,p);
	MatrixXd data(k,n);
	c_orth = (svd.matrixU()).transpose();
	multiply_y_post(c_orth,k,data);

	MatrixXd means(k,1);
	for(int i=0;i<k;i++){
		float sum=0.0;
		for(int j=0;j<n;j++)
			sum+=data(i,j);
		means(i,0)=sum/(n*1.0);
	}
	data = data - (means*(MatrixXd::Constant(1,n,1)));
	MatrixXd cov(k,k);
	cov = data*(data.transpose())*(1.0/(n));
	JacobiSVD<MatrixXd> svd_cov(cov, ComputeThinU | ComputeThinV);
	MatrixXd to_return(p,k);
	to_return =(c_orth.transpose())*svd_cov.matrixU() ;
	return to_return;
}


float get_accuracy(MatrixXd &u){
	
	MatrixXd temp(k,k);
	temp = (u.transpose()) * v ;
	float accuracy = 0.0;
	for(int j=0;j<k;j++){
		float sum=0.0;
		for(int i=0;i<k;i++)
			sum += temp(i,j)*temp(i,j);
		accuracy += sqrt(sum);
	}
	return (accuracy/k);
}

MatrixXd get_reference_evec()
{

	MatrixXd y_m(p,n);
	for(int i=0;i<p;i++){
		for(int j=0;j<n;j++)
			y_m(i,j) = g.get_geno(i,j,var_normalize);
	}
	MatrixXd cov(p,p);
	if(debug)
		printf("Calculating covariance\n");
	cov = y_m*(y_m.transpose())*(1.0/n);
	if(debug)
		printf("Calculating SVD\n");
	JacobiSVD<MatrixXd> svd_cov(cov, ComputeThinU | ComputeThinV);
	MatrixXd to_return(p,k);
	MatrixXd U(p,k);
	U = svd_cov.matrixU();
	for(int i=0;i<p;i++){
		for(int j=0;j<k;j++)
			to_return(i,j) = U(i,j);
	}

	return to_return;
}

pair<double,double> get_error_norm(MatrixXd &c)
{
	HouseholderQR<MatrixXd> qr(c);
	MatrixXd Q;
	Q = qr.householderQ() * MatrixXd::Identity(p,k);
	MatrixXd q_t(k,p);
	q_t = Q.transpose();
	MatrixXd b(k,n);
	multiply_y_post(q_t,k,b);
	JacobiSVD<MatrixXd> b_svd(b, ComputeThinU | ComputeThinV);
	
	MatrixXd u_l,d_l,v_l; 
	u_l = Q * b_svd.matrixU();
	v_l = b_svd.matrixV();
	d_l = MatrixXd::Zero(k,k);
	for(int kk=0;kk<k; kk++)
		d_l(kk,kk) = (b_svd.singularValues())(kk);
	
	MatrixXd u_k,v_k,d_k;
	u_k = u_l.leftCols(k_orig);
	v_k = v_l.leftCols(k_orig);
	d_k = MatrixXd::Zero(k_orig,k_orig);
	for(int kk =0 ; kk < k_orig ; kk++)
		d_k(kk,kk)  =(b_svd.singularValues())(kk);

	MatrixXd b_l,b_k;
	b_l = u_l * d_l * (v_l.transpose());
	b_k = u_k * d_k * (v_k.transpose());

	MatrixXd e_l(p,n);
	MatrixXd e_k(p,n);
	for(int p_iter=0;p_iter<p;p_iter++)
	{
		for(int n_iter=0;n_iter<n;n_iter++){
			e_l(p_iter,n_iter) = g.get_geno(p_iter,n_iter,var_normalize) - b_l(p_iter,n_iter);
			e_k(p_iter,n_iter) = g.get_geno(p_iter,n_iter,var_normalize) - b_k(p_iter,n_iter);
		}
	}

	double ek_norm = e_k.norm();
	double el_norm = e_l.norm();
	return make_pair(ek_norm,el_norm);

}



MatrixXd run_EM(MatrixXd &c_orig){
	MatrixXd c_temp(k,p);
	MatrixXd c_new(p,k);
	c_temp = ( (c_orig.transpose()*c_orig).inverse() ) * (c_orig.transpose());

	MatrixXd x_fn(k,n);
	multiply_y_post(c_temp,k,x_fn);

	MatrixXd x_temp(n,k);
	x_temp = (x_fn.transpose()) * ((x_fn*(x_fn.transpose())).inverse());
	multiply_y_pre(x_temp,k,c_new);
	return c_new;
}


int main(int argc, char const *argv[])
{

	clock_t io_begin = clock();
	pair<double,double> prev_error = make_pair(0.0,0.0);
	double prevnll=0.0;
	
	parse_args(argc,argv);
	g.read_genotype_naive(command_line_opts.GENOTYPE_FILE_PATH);	
	MAX_ITER =  command_line_opts.max_iterations ; 
	k_orig = command_line_opts.num_of_evec ;
	debug = command_line_opts.debugmode ;
	check_accuracy = command_line_opts.getaccuracy;
	k = k_orig + command_line_opts.l;
	p = g.Nsnp;
	n = g.Nindv;
	c.resize(p,k);
	x.resize(k,n);
	v.resize(p,k);
	convergence_limit = command_line_opts.convergence_limit;
	srand((unsigned int) time(0));
	
	clock_t io_end = clock();

	c = MatrixXd::Random(p,k);
	
	ofstream c_file;
	if(debug){
		c_file.open((string(command_line_opts.OUTPUT_PATH)+string("cvals_orig_naive.txt")).c_str());
		c_file<<c<<endl;
		c_file.close();
		printf("Read Matrix\n");
	}
	if(check_accuracy){
		// v = get_reference_evec(); 
	 	printf("Computed Reference Eigen Vectors\n");	
	}
	cout<<"Running on Dataset of "<<g.Nsnp<<" SNPs and "<<g.Nindv<<" Individuals"<<endl;
	if(check_accuracy)
		cout<<endl<<"Iterations vs accuracy"<<endl;
	clock_t it_begin = clock();
	for(int i=0;i<MAX_ITER;i++)
	{
		MatrixXd c1,c2,cint,r,v;
		double a,nll;
		if(accelerated_em!=0){
			c1 = run_EM(c);
			c2 = run_EM(c1);
			r = c1-c;
			v = (c2-c1) - r;
			a = -1.0 * r.norm() / (v.norm()) ;
			if(accelerated_em==1){
				if(a>-1){
					a=-1;
					cint=c2;
				}
				else {
					cint = c - 2*a*r + a*a*v;
					nll = get_error_norm(cint).second;
					if(i>0){
						while(nll>prevnll && a<-1){
							a = 0.5 * (a-1);
							cint = c - 2*a*r +(a*a*v);
							nll = get_error_norm(cint).second;
						}
					}
				}
				c = cint;
			}
			else if(accelerated_em==2){
				cint = c - 2*a*r + a*a*v;
				c = run_EM(cint);				
			}
		}
		else
			c = run_EM(c);
		
		pair<double,double> e = get_error_norm(c);
		prevnll = e.second;
		if(check_accuracy)
			cout<<"Iteration "<<i+1<<"  "<<prev_error.first - e.first<<"  "<<prev_error.second - e.second<<endl;
		prev_error = e;
		
	}
	clock_t it_end = clock();
	c_file.open((string(command_line_opts.OUTPUT_PATH)+string("cvals_end_naive.txt")).c_str());
	c_file<<c<<endl;
	c_file.close();
	ofstream x_file;
	x_file.open((string(command_line_opts.OUTPUT_PATH) + string("xvals_naive.txt")).c_str());
	x_file<<x<<endl;
	x_file.close();
	clock_t total_end = clock();
	double io_time = double(io_end - io_begin) / CLOCKS_PER_SEC;
	double avg_it_time = double(it_end - it_begin) / (MAX_ITER * 1.0 * CLOCKS_PER_SEC);
	double total_time = double(total_end - total_begin) / CLOCKS_PER_SEC;
	cout<<"IO Time:  "<< io_time << "\nAVG Iteration Time:  "<<avg_it_time<<"\nTotal runtime:   "<<total_time<<endl;
	return 0;
}
