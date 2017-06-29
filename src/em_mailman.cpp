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
double convergence_limit;

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


void multiply_y_pre(MatrixXd &op, int Ncol_op ,MatrixXd &res)
{
	double *sum_op = new double[Ncol_op];
	double *yint_m = new double[(int)pow(2,g.segment_size_hori)];
	for(int k_iter=0;k_iter<Ncol_op;k_iter++)
	{
		double *y_msb = new double[g.segment_size_hori];
		double *y_lsb = new double[g.segment_size_hori];
		int seg_iter;
		MatrixXd op_col(n,1);
		op_col = op.col(k_iter);
		sum_op[k_iter]=op_col.sum();

		for(seg_iter=0;seg_iter<g.Nsegments_hori;seg_iter++)
		{
			if(seg_iter==g.Nsegments_hori-1)
			{
				if(g.Nsnp%g.segment_size_hori!=0)
				{
					double *y_msb_final = new double[g.Nsnp%g.segment_size_hori];
					double *y_lsb_final = new double[g.Nsnp%g.segment_size_hori];
					mailman::fastmultiply(g.Nsnp%g.segment_size_hori,g.Nindv,g.p_msb[seg_iter],op_col,yint_m,y_msb_final);
					mailman::fastmultiply(g.Nsnp%g.segment_size_hori,g.Nindv,g.p_lsb[seg_iter],op_col,yint_m,y_lsb_final);
					for(int p_iter=seg_iter*g.segment_size_hori;p_iter<seg_iter*g.segment_size_hori + g.Nsnp%g.segment_size_hori  && p_iter<g.Nsnp;p_iter++)
						res(p_iter,k_iter) = 2*y_msb_final[p_iter-(seg_iter*g.segment_size_hori)] + y_lsb_final[p_iter-(seg_iter*g.segment_size_hori)];
					break;
				}
			}
			mailman::fastmultiply(g.segment_size_hori,g.Nindv,g.p_msb[seg_iter],op_col,yint_m,y_msb);
			mailman::fastmultiply(g.segment_size_hori,g.Nindv,g.p_lsb[seg_iter],op_col,yint_m,y_lsb);
			int p_base = seg_iter*g.segment_size_hori; 
			for(int p_iter=p_base; (p_iter<p_base+g.segment_size_hori) && (p_iter<g.Nsnp) ; p_iter++ ) 
				c(p_iter,k_iter) = 2*y_msb[p_iter-p_base] + y_lsb[p_iter-p_base];
		}			
	}

	for(int p_iter=0;p_iter<p;p_iter++)
	{
		for(int k_iter=0;k_iter<Ncol_op;k_iter++)
			res(p_iter,k_iter) = res(p_iter,k_iter) - (g.get_col_mean(p_iter)*sum_op[k_iter]);
	}	
}


void multiply_y_post(MatrixXd &op_orig, int Nrows_op, MatrixXd &res)
{
	double *yint_e = new double[(int)pow(2,g.segment_size_ver)];
	MatrixXd op;
	op = op_orig.transpose();
	int Ncol_op = Nrows_op;
	for(int k_iter=0;k_iter<Ncol_op;k_iter++)
	{
		double *y_msb = new double[g.segment_size_ver];
		double *y_lsb = new double[g.segment_size_ver];
		int seg_iter;
		MatrixXd op_col(p,1);
		op_col = op.col(k_iter);

		for(seg_iter=0; seg_iter < g.Nsegments_ver; seg_iter++)
		{
			if(seg_iter==g.Nsegments_ver-1)
			{
				if(g.Nindv%g.segment_size_ver!=0)
				{
					double *y_msb_final = new double[g.Nindv%g.segment_size_ver];
					double *y_lsb_final = new double[g.Nindv%g.segment_size_ver];
					mailman::fastmultiply(g.Nindv%g.segment_size_ver,g.Nsnp,g.q_msb[seg_iter],op_col,yint_e,y_msb_final);
					mailman::fastmultiply(g.Nindv%g.segment_size_ver,g.Nsnp,g.q_lsb[seg_iter],op_col,yint_e,y_lsb_final);
					for(int n_iter=seg_iter*g.segment_size_ver ; n_iter<seg_iter*g.segment_size_ver + g.Nindv%g.segment_size_ver  && n_iter<g.Nindv ; n_iter++)
						res(k_iter,n_iter) = 2*y_msb_final[n_iter-(seg_iter*g.segment_size_ver)] + y_lsb_final[n_iter-(seg_iter*g.segment_size_ver)];
					break;
				}

			}
			mailman::fastmultiply(g.segment_size_ver,g.Nsnp,g.q_msb[seg_iter],op_col,yint_e,y_msb);
			mailman::fastmultiply(g.segment_size_ver,g.Nsnp,g.q_lsb[seg_iter],op_col,yint_e,y_lsb);
			int n_base = seg_iter*g.segment_size_ver; 
			for(int n_iter=n_base; (n_iter<n_base+g.segment_size_ver) && (n_iter<g.Nindv) ; n_iter++ ) 
				res(k_iter,n_iter) = 2*y_msb[n_iter-n_base] + y_lsb[n_iter-n_base];

		}
	}
	double *sums_elements = new double[Ncol_op];
	memset (sums_elements, 0, Nrows_op * sizeof(int));

	for(int k_iter=0;k_iter<Ncol_op;k_iter++)
	{
		double sum_to_calc=0.0;
		for(int p_iter=0;p_iter<p;p_iter++)
			sum_to_calc += g.get_col_mean(p_iter)*op(p_iter,k_iter);
		sums_elements[k_iter] = sum_to_calc;
	}

	for(int k_iter=0;k_iter<Ncol_op;k_iter++)
	{
		for(int n_iter=0;n_iter<n;n_iter++)
			res(k_iter,n_iter) = res(k_iter,n_iter) - sums_elements[k_iter];
	}
}

double get_error_norm()
{
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
	d_k = MatrixXd::Zero(k_orig,k_orig);
	for(int kk =0 ; kk < k_orig ; kk++)
		d_k(kk,kk)  =(b_svd.singularValues())(kk);

	MatrixXd b_k;
	b_k = u_k * d_k * (v_k.transpose());
	double sum_norm_temp=0.0;
	for(int k_iter=0;k_iter<k;k_iter++)
	{
		for(int n_iter=0;n_iter<n;n_iter++)
			sum_norm_temp+= b_k(k_iter,n_iter)*b(k_iter,n_iter);
	}
	double b_knorm = b_k.norm();
	double ef_norm = (b_knorm*b_knorm) - (2*sum_norm_temp);
	return ef_norm;
}

void print_vals()
{

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
	if(debug){
		ofstream c_file;
		c_file.open((string(command_line_opts.OUTPUT_PATH)+string("cvals_mailman.txt")).c_str());
		c_file<<c<<endl;
		c_file.close();
		ofstream x_file;
		x_file.open((string(command_line_opts.OUTPUT_PATH) + string("xvals_mailman.txt")).c_str());
		x_file<<x<<endl;
		x_file.close();
	}

}


int main(int argc, char const *argv[])
{
	clock_t io_begin = clock();

	double prev_error = 0.0;
	
	parse_args(argc,argv);
	g.read_genotype_mailman(command_line_opts.GENOTYPE_FILE_PATH);	
	MAX_ITER =  command_line_opts.max_iterations ; 
	k_orig = command_line_opts.num_of_evec ;
	debug = command_line_opts.debugmode ;
	check_accuracy = command_line_opts.getaccuracy;
	k = k_orig + command_line_opts.l;
	p = g.Nsnp;
	n = g.Nindv;
	convergence_limit = command_line_opts.convergence_limit;
	srand((unsigned int) time(0));
	c.resize(p,k);
	x.resize(k,n);
	v.resize(p,k);

	clock_t io_end = clock();

	c = MatrixXd::Random(p,k);
	
	ofstream c_file;
	if(debug){
		c_file.open((string(command_line_opts.OUTPUT_PATH)+string("cvals_orig_mailman.txt")).c_str());
		c_file<<c<<endl;
		c_file.close();
		printf("Read Matrix\n");
	}

	cout<<"Running on Dataset of "<<g.Nsnp<<" SNPs and "<<g.Nindv<<" Individuals"<<endl;

	if(check_accuracy)
		cout<<endl<<"Iterations vs accuracy"<<endl;
	
	clock_t it_begin = clock();
	for(int i=0;i<MAX_ITER;i++)
	{
		if(debug)
			printf("Iteration %d -- E Step\n",i);

		MatrixXd c_temp(p,k);
		c_temp = ( (c.transpose()*c).inverse() ) * (c.transpose());

		multiply_y_post(c_temp,k,x);
		
		if(debug)
			printf("Iteration %d -- M Step\n",i);

		MatrixXd x_temp(n,k);
		x_temp = (x.transpose()) * ((x*(x.transpose())).inverse());

		multiply_y_pre(x_temp,k,c);

		if(check_accuracy){

			double e = get_error_norm();
			cout<<"Iteration "<<i+1<<" "<<e<<"  "<<abs(e-prev_error)<<endl;
			if( abs((e- prev_error))< convergence_limit)
				break;
			else
				prev_error = e;
		}
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
