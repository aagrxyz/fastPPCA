/** 
 All of this code is written by Aman Agrawal 
 (Indian Institute of Technology, Delhi)
*/

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

/* Redundant Functions 
MatrixXd get_evec(MatrixXd &c)
{
	JacobiSVD<MatrixXd> svd(c, ComputeThinU | ComputeThinV);
	MatrixXd c_orth(k,p);
	MatrixXd data(k,n);
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
	MatrixXd means(k,1);
	for(int i=0;i<k;i++)
	{
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

float get_accuracy(MatrixXd &u)
{
	
	MatrixXd temp(k,k);
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

MatrixXd get_reference_evec()
{

	MatrixXd y_m(p,n);
	for(int i=0;i<p;i++)
	{
		for(int j=0;j<n;j++)
			y_m(i,j) = g.get_geno(i,j);
	}
	MatrixXd cov(p,p);
	printf("Calculating covariance\n");
	cov = y_m*(y_m.transpose())*(1.0/n);
	printf("Calculating SVD\n");
	JacobiSVD<MatrixXd> svd_cov(cov, ComputeThinU | ComputeThinV);
	MatrixXd to_return(p,k);
	MatrixXd U(p,k);
	U = svd_cov.matrixU();
	for(int i=0;i<p;i++)
	{
		for(int j=0;j<k;j++)
			to_return(i,j) = U(i,j);
	}

	return to_return;
}

*/


void multiply_y_pre(MatrixXd &op, int Ncol_op ,MatrixXd &res){
	double *sum_op = new double[Ncol_op];
	double *yint_m = new double[(int)pow(3,g.segment_size_hori)];

	for(int k_iter=0;k_iter<Ncol_op;k_iter++){
		double *y = new double[g.segment_size_hori];
		int seg_iter;
		MatrixXd op_col(n,1);
		op_col = op.col(k_iter);
		sum_op[k_iter]=op_col.sum();

		for(seg_iter=0;seg_iter<g.Nsegments_hori;seg_iter++){
			if(seg_iter==g.Nsegments_hori-1){
				if(g.Nsnp%g.segment_size_hori!=0){
					double *y_final = new double[g.Nsnp%g.segment_size_hori];
					mailman::fastmultiply2(g.Nsnp%g.segment_size_hori,g.Nindv,g.p[seg_iter],op_col,yint_m,y_final);
					for(int p_iter=seg_iter*g.segment_size_hori;p_iter<seg_iter*g.segment_size_hori + g.Nsnp%g.segment_size_hori  && p_iter<g.Nsnp;p_iter++)
						res(p_iter,k_iter) = y_final[p_iter-(seg_iter*g.segment_size_hori)];
					break;
				}
			}
			mailman::fastmultiply2(g.segment_size_hori,g.Nindv,g.p[seg_iter],op_col,yint_m,y);
			int p_base = seg_iter*g.segment_size_hori; 
			for(int p_iter=p_base; (p_iter<p_base+g.segment_size_hori) && (p_iter<g.Nsnp) ; p_iter++ ) 
				res(p_iter,k_iter) = y[p_iter-p_base];
		}			
	}

	for(int p_iter=0;p_iter<p;p_iter++){
		for(int k_iter=0;k_iter<Ncol_op;k_iter++){
			res(p_iter,k_iter) = res(p_iter,k_iter) - (g.get_col_mean(p_iter)*sum_op[k_iter]);
			if(var_normalize)
				res(p_iter,k_iter) = res(p_iter,k_iter)/(g.get_col_std(p_iter));
		}
	}	
}

void multiply_y_post(MatrixXd &op_orig, int Nrows_op, MatrixXd &res){
	double *yint_e = new double[(int)pow(3,g.segment_size_ver)];

	MatrixXd op;
	op = op_orig.transpose();

	if(var_normalize){
		for(int p_iter=0;p_iter<p;p_iter++){
			for(int k_iter=0;k_iter<Nrows_op;k_iter++)
				op(p_iter,k_iter) = op(p_iter,k_iter) / (g.get_col_std(p_iter));
		}
	}

	int Ncol_op = Nrows_op;
	for(int k_iter=0;k_iter<Ncol_op;k_iter++){
		double *y = new double[g.segment_size_ver];
		int seg_iter;
		MatrixXd op_col(p,1);
		op_col = op.col(k_iter);

		for(seg_iter=0; seg_iter < g.Nsegments_ver; seg_iter++){
			if(seg_iter==g.Nsegments_ver-1){
				if(g.Nindv%g.segment_size_ver!=0){
					double *y_final = new double[g.Nindv%g.segment_size_ver];
					mailman::fastmultiply2(g.Nindv%g.segment_size_ver,g.Nsnp,g.q[seg_iter],op_col,yint_e,y_final);
					for(int n_iter=seg_iter*g.segment_size_ver ; n_iter<seg_iter*g.segment_size_ver + g.Nindv%g.segment_size_ver  && n_iter<g.Nindv ; n_iter++)
						res(k_iter,n_iter) = y_final[n_iter-(seg_iter*g.segment_size_ver)];
					break;
				}

			}
			mailman::fastmultiply2(g.segment_size_ver,g.Nsnp,g.q[seg_iter],op_col,yint_e,y);
			int n_base = seg_iter*g.segment_size_ver; 
			for(int n_iter=n_base; (n_iter<n_base+g.segment_size_ver) && (n_iter<g.Nindv) ; n_iter++ ) 
				res(k_iter,n_iter) = y[n_iter-n_base];

		}
	}
	double *sums_elements = new double[Ncol_op];
	memset (sums_elements, 0, Nrows_op * sizeof(int));

	for(int k_iter=0;k_iter<Ncol_op;k_iter++){
		double sum_to_calc=0.0;
		for(int p_iter=0;p_iter<p;p_iter++)
			sum_to_calc += g.get_col_mean(p_iter)*op(p_iter,k_iter);
		sums_elements[k_iter] = sum_to_calc;
	}

	for(int k_iter=0;k_iter<Ncol_op;k_iter++){
		for(int n_iter=0;n_iter<n;n_iter++)
			res(k_iter,n_iter) = res(k_iter,n_iter) - sums_elements[k_iter];
	}
}

pair<double,double> get_error_norm(MatrixXd &c){
	HouseholderQR<MatrixXd> qr(c);
	MatrixXd Q;
	Q = qr.householderQ() * MatrixXd::Identity(p,k);
	MatrixXd q_t(k,p);
	q_t = Q.transpose();
	MatrixXd b(k,n);
	multiply_y_post(q_t,k,b);
	JacobiSVD<MatrixXd> b_svd(b, ComputeThinU | ComputeThinV);
	MatrixXd u_l,d_l,v_l; 
	u_l = b_svd.matrixU();
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

	double temp_k=0.0;
	double temp_l=0.0;
	for(int k_iter=0;k_iter<k;k_iter++){
		for(int n_iter=0;n_iter<n;n_iter++){
			temp_k += b_k(k_iter,n_iter)*b(k_iter,n_iter);
			temp_l += b_l(k_iter,n_iter)*b(k_iter,n_iter);
		}
	}

	double b_knorm = b_k.norm();
	double b_lnorm = b_l.norm();
	double norm_k = (b_knorm*b_knorm) - (2*temp_k);
	double norm_l = (b_lnorm*b_lnorm) - (2*temp_l);	
	return make_pair(norm_k,norm_l);
}

void print_vals(){

	HouseholderQR<MatrixXd> qr(c);
	MatrixXd Q;
	Q = qr.householderQ() * MatrixXd::Identity(p,k);
	MatrixXd q_t(k,p);
	q_t = Q.transpose();
	MatrixXd b(k,n);
	multiply_y_post(q_t,k,b);
	JacobiSVD<MatrixXd> b_svd(b, ComputeThinU | ComputeThinV);
	MatrixXd u_l; 
	u_l = b_svd.matrixU();
	MatrixXd v_l;
	v_l = b_svd.matrixV();
	MatrixXd u_k;
	MatrixXd v_k,d_k;
	u_k = u_l.leftCols(k_orig);
	v_k = v_l.leftCols(k_orig);

	ofstream evec_file;
	evec_file.open((string(command_line_opts.OUTPUT_PATH)+string("evecs_mailman.txt")).c_str());
	evec_file<< Q*u_k << endl;
	evec_file.close();
	ofstream eval_file;
	eval_file.open((string(command_line_opts.OUTPUT_PATH)+string("evals_mailman.txt")).c_str());
	for(int kk =0 ; kk < k_orig ; kk++)
		eval_file << (b_svd.singularValues())(kk)<<endl;
	eval_file.close();

	d_k = MatrixXd::Zero(k_orig,k_orig);
	for(int kk =0 ; kk < k_orig ; kk++)
		d_k(kk,kk)  =(b_svd.singularValues())(kk);

	MatrixXd x_k;
	x_k = d_k * (v_k.transpose());
	if(debug){
		ofstream c_file;
		c_file.open((string(command_line_opts.OUTPUT_PATH)+string("cvals_mailman.txt")).c_str());
		c_file<<c<<endl;
		c_file.close();
		ofstream x_file;
		x_file.open((string(command_line_opts.OUTPUT_PATH) + string("xvals_mailman.txt")).c_str());
		x_file<<x_k<<endl;
		x_file.close();
	}

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

int main(int argc, char const *argv[]){
	
	clock_t io_begin = clock();
	pair<double,double> prev_error = make_pair(0.0,0.0);
	double prevnll=0.0;

	parse_args(argc,argv);

	g.read_genotype_mailman(command_line_opts.GENOTYPE_FILE_PATH);	
	MAX_ITER =  command_line_opts.max_iterations ; 
	k_orig = command_line_opts.num_of_evec ;
	debug = command_line_opts.debugmode ;
	check_accuracy = command_line_opts.getaccuracy;
	var_normalize = command_line_opts.var_normalize;
	accelerated_em = command_line_opts.accelerated_em;
	k = k_orig + command_line_opts.l;
	p = g.Nsnp;
	n = g.Nindv;
	convergence_limit = command_line_opts.convergence_limit;
	srand((unsigned int) time(0));
	c.resize(p,k);
	x.resize(k,n);
	v.resize(p,k);

	clock_t io_end = clock();

	//TODO: Initialization of c

	c = MatrixXd::Random(p,k);
	
	ofstream c_file;
	if(debug){
		c_file.open((string(command_line_opts.OUTPUT_PATH)+string("cvals_orig_mailman.txt")).c_str());
		c_file<<c<<endl;
		c_file.close();
		printf("Read Matrix\n");
	}

	cout<<"Running on Dataset of "<<g.Nsnp<<" SNPs and "<<g.Nindv<<" Individuals"<<endl;
	
	clock_t it_begin = clock();
	for(int i=0;i<MAX_ITER;i++){

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

	print_vals();
		
	clock_t total_end = clock();
	double io_time = double(io_end - io_begin) / CLOCKS_PER_SEC;
	double avg_it_time = double(it_end - it_begin) / (MAX_ITER * 1.0 * CLOCKS_PER_SEC);
	double total_time = double(total_end - total_begin) / CLOCKS_PER_SEC;
	cout<<"IO Time:  "<< io_time << "\nAVG Iteration Time:  "<<avg_it_time<<"\nTotal runtime:   "<<total_time<<endl;
	return 0;
}