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
#include "time.h"
#include <Eigen/QR>
#include "storage.h"

#ifdef SSE_SUPPORT
	#define fastmultiply fastmultiply_sse
	#define fastmultiply_pre fastmultiply_pre_sse
#else
	#define fastmultiply fastmultiply_normal
	#define fastmultiply_pre fastmultiply_pre_normal
#endif

using namespace Eigen;
using namespace std;

// Storing in RowMajor Form
typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdr;

//Intermediate Variables
int blocksize;
double *partialsums;
double *sum_op;
double *yint_e;
double *yint_m;
double **y_e;
double **y_m;


struct timespec t0;
struct timespec
elapsed (){
  struct timespec ts;
  clock_gettime (CLOCK_REALTIME, &ts);
  if (ts.tv_nsec < t0.tv_nsec)
    {
      ts.tv_nsec = 1000000000 + ts.tv_nsec - t0.tv_nsec;
      ts.tv_sec--;
    }
  ts.tv_sec -= t0.tv_sec;
  return (ts);
}

int timelog (const char* message){
  struct timespec ts = elapsed ();
  return (printf ("[%06ld.%09ld] %s\n", ts.tv_sec, ts.tv_nsec, message));
}

void * malloc_double_align(size_t n, unsigned int a /*alignment*/, double *& output){
    void *adres=NULL;
    void *adres2=NULL;
    adres=malloc(n*sizeof(double)+a);
    size_t adr=(size_t)adres;
    size_t adr2=adr+a-(adr&(a-1u)); 	// a valid address for a alignment
    adres2=(void * ) adr2;
    output=(double *)adres2;
    return adres;                		// pointer to be used in free()
}


clock_t total_begin = clock();

int MAX_ITER;
genotype g;
int k,p,n;
int k_orig;

MatrixXdr c; //(p,k)
MatrixXdr x; //(k,n)
MatrixXdr v; //(p,k)
MatrixXdr means; //(p,1)

options command_line_opts;

bool debug = false;
bool check_accuracy = false;
bool var_normalize=false;
int accelerated_em=0;
double convergence_limit;
bool memory_efficient = false;
bool missing=false;

bool fast_mode = true;

bool text_version = false;
void print_timenl () {
	clock_t c = clock();
	double t = double(c) / CLOCKS_PER_SEC;
	cout << "Time = " << t << endl ;	
}

void print_time () {
	clock_t c = clock();
	double t = double(c) / CLOCKS_PER_SEC;
	cout << "Time = " << t  << " : ";	
}

void multiply_y_pre_fast(MatrixXdr &op, int Ncol_op ,MatrixXdr &res,bool subtract_means){
	
	for(int k_iter=0;k_iter<Ncol_op;k_iter++){
		sum_op[k_iter]=op.col(k_iter).sum();
	}

	if(debug){
		print_time (); 
		cout <<"Starting mailman on premultiply"<<endl;
		cout << "Nops = " << Ncol_op << "\t" <<g.Nsegments_hori << endl;
		cout << "Segment size = " << g.segment_size_hori << endl;
		cout << "Matrix size = " <<g.segment_size_hori<<"\t" <<g.Nindv << endl;
		cout << "op = " <<  op.rows () << "\t" << op.cols () << endl;
		cout << g.segment_size_ver << "\t" << g.Nsegments_ver << endl;
	}


	for(int seg_iter=0;seg_iter<g.Nsegments_hori-1;seg_iter++){
		//TODO: Memory Effecient SSE FastMultipy
		// if(memory_efficient)
		// 	mailman::fastmultiply_memory_eff(g.segment_size_hori,g.Nindv,Ncol_op,g.p_eff[seg_iter],op,yint_m,partialsums,y_m,g.Nbits_hori);
		// else
		
		mailman::fastmultiply(g.segment_size_hori,g.Nindv,Ncol_op,g.p[seg_iter],op,yint_m,partialsums,y_m);
		
		int p_base = seg_iter*g.segment_size_hori; 
		for(int p_iter=p_base; (p_iter<p_base+g.segment_size_hori) && (p_iter<g.Nsnp) ; p_iter++ ){
			for(int k_iter=0;k_iter<Ncol_op;k_iter++) 
				res(p_iter,k_iter) = y_m[p_iter-p_base][k_iter];
		}
	}

	int last_seg_size = (g.Nsnp%g.segment_size_hori !=0 ) ? g.Nsnp%g.segment_size_hori : g.segment_size_hori;
	//TODO: Memory Effecient SSE FastMultipy	
	// if(memory_efficient)
	// 	mailman::fastmultiply_memory_eff(last_seg_size,g.Nindv,Ncol_op,g.p_eff[g.Nsegments_hori-1],op,yint_m,partialsums,y_m,g.Nbits_hori);
	// else
	mailman::fastmultiply(last_seg_size,g.Nindv,Ncol_op,g.p[g.Nsegments_hori-1],op,yint_m,partialsums,y_m);	
	
	if(debug){
		print_time (); 
		cout <<"Ending mailman on premultiply"<<endl;
	}

	int p_base = (g.Nsegments_hori-1)*g.segment_size_hori;
	for(int p_iter=p_base; (p_iter<p_base+g.segment_size_hori) && (p_iter<g.Nsnp) ; p_iter++){
		for(int k_iter=0;k_iter<Ncol_op;k_iter++) 
			res(p_iter,k_iter) = y_m[p_iter-p_base][k_iter];
	}


	if(!subtract_means)
		return;
	for(int p_iter=0;p_iter<p;p_iter++){
		for(int k_iter=0;k_iter<Ncol_op;k_iter++){
			res(p_iter,k_iter) = res(p_iter,k_iter) - (g.get_col_mean(p_iter)*sum_op[k_iter]);
			if(var_normalize)
				res(p_iter,k_iter) = res(p_iter,k_iter)/(g.get_col_std(p_iter));
		}
	}	
}

void multiply_y_post_fast(MatrixXdr &op_orig, int Nrows_op, MatrixXdr &res,bool subtract_means){
	MatrixXdr op;
	op = op_orig.transpose();

	if(var_normalize){
		for(int p_iter=0;p_iter<p;p_iter++){
			for(int k_iter=0;k_iter<Nrows_op;k_iter++)
				op(p_iter,k_iter) = op(p_iter,k_iter) / (g.get_col_std(p_iter));
		}
	}

	if(debug){
		print_time (); cout <<"Starting mailman on postmultiply"<<endl;
		cout << "Nops = " << Nrows_op << "\t" <<g.Nsegments_ver << endl;
		cout << "Segment size = " << g.segment_size_ver << endl;
		cout << "Matrix size = " <<g.segment_size_ver <<"\t" <<g.Nsnp << endl;
	}

	int Ncol_op = Nrows_op;


	int seg_iter;
	for(seg_iter=0;seg_iter<g.Nsegments_hori-1;seg_iter++){
		mailman::fastmultiply_pre(g.segment_size_hori,g.Nindv,Ncol_op, seg_iter * g.segment_size_hori, g.p[seg_iter],op,yint_e,partialsums,y_e);
	}
	int last_seg_size = (g.Nsnp%g.segment_size_hori !=0 ) ? g.Nsnp%g.segment_size_hori : g.segment_size_hori;
	mailman::fastmultiply_pre(last_seg_size,g.Nindv,Ncol_op, seg_iter * g.segment_size_hori, g.p[seg_iter],op,yint_e,partialsums,y_e);


	if(debug){
		print_time (); 
		cout <<"Ending mailman on postmultiply"<<endl;
	}

	for(int n_iter=0; n_iter<n; n_iter++)  {
		for(int k_iter=0;k_iter<Ncol_op;k_iter++) {
			res(k_iter,n_iter) = y_e[n_iter][k_iter];
			y_e[n_iter][k_iter] = 0;
		}
	}
	
	if(!subtract_means)
		return;
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

void multiply_y_pre_naive(MatrixXdr &op, int Ncol_op ,MatrixXdr &res){
	for(int p_iter=0;p_iter<p;p_iter++){
		for(int k_iter=0;k_iter<Ncol_op;k_iter++){
			double temp=0;
			for(int n_iter=0;n_iter<n;n_iter++)
				temp+= g.get_geno(p_iter,n_iter,var_normalize)*op(n_iter,k_iter);
			res(p_iter,k_iter)=temp;
		}
	}
}

void multiply_y_post_naive(MatrixXdr &op, int Nrows_op ,MatrixXdr &res){
	for(int n_iter=0;n_iter<n;n_iter++){
		for(int k_iter=0;k_iter<Nrows_op;k_iter++){
			double temp=0;
			for(int p_iter=0;p_iter<p;p_iter++)
				temp+= op(k_iter,p_iter)*(g.get_geno(p_iter,n_iter,var_normalize));
			res(k_iter,n_iter)=temp;
		}
	}
}

void multiply_y_post(MatrixXdr &op, int Nrows_op ,MatrixXdr &res,bool subtract_means){
    if(fast_mode)
        multiply_y_post_fast(op,Nrows_op,res,subtract_means);
    else
        multiply_y_post_naive(op,Nrows_op,res);
}

void multiply_y_pre(MatrixXdr &op, int Ncol_op ,MatrixXdr &res,bool subtract_means){
    if(fast_mode)
        multiply_y_pre_fast(op,Ncol_op,res,subtract_means);
    else
        multiply_y_pre_naive(op,Ncol_op,res);
}

pair<double,double> get_error_norm(MatrixXdr &c){
	HouseholderQR<MatrixXdr> qr(c);
	MatrixXdr Q;
	Q = qr.householderQ() * MatrixXdr::Identity(p,k);
	MatrixXdr q_t(k,p);
	q_t = Q.transpose();
	MatrixXdr b(k,n);
	multiply_y_post(q_t,k,b,true);
	JacobiSVD<MatrixXdr> b_svd(b, ComputeThinU | ComputeThinV);
	MatrixXdr u_l,d_l,v_l; 
	if(fast_mode)
        u_l = b_svd.matrixU();
    else
        u_l = Q * b_svd.matrixU();
	v_l = b_svd.matrixV();
	d_l = MatrixXdr::Zero(k,k);
	for(int kk=0;kk<k; kk++)
		d_l(kk,kk) = (b_svd.singularValues())(kk);
	
	MatrixXdr u_k,v_k,d_k;
	u_k = u_l.leftCols(k_orig);
	v_k = v_l.leftCols(k_orig);
	d_k = MatrixXdr::Zero(k_orig,k_orig);
	for(int kk =0 ; kk < k_orig ; kk++)
		d_k(kk,kk)  =(b_svd.singularValues())(kk);

	MatrixXdr b_l,b_k;
    b_l = u_l * d_l * (v_l.transpose());
    b_k = u_k * d_k * (v_k.transpose());

    if(fast_mode){
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
    else{
        MatrixXdr e_l(p,n);
        MatrixXdr e_k(p,n);
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
}

MatrixXdr run_EM_not_missing(MatrixXdr &c_orig){
	if(debug){
		print_time ();
		cout << "Enter: run_EM_not_missing" << endl;
	}
	MatrixXdr c_temp(k,p);
	MatrixXdr c_new(p,k);
	c_temp = ( (c_orig.transpose()*c_orig).inverse() ) * (c_orig.transpose());
	if(debug)
		print_timenl ();
	MatrixXdr x_fn(k,n);
	multiply_y_post(c_temp,k,x_fn,true);
	if(debug)
		print_timenl ();
	MatrixXdr x_temp(n,k);
	x_temp = (x_fn.transpose()) * ((x_fn*(x_fn.transpose())).inverse());
	multiply_y_pre(x_temp,k,c_new,true);
	if(debug){
		print_time ();
		cout << "Exiting: run_EM_not_missing" << endl;
	}
	return c_new;
}

MatrixXdr run_EM_missing(MatrixXdr &c_orig){
	
	MatrixXdr c_new(p,k);

	MatrixXdr mu(k,n);
	
	// E step
	MatrixXdr c_temp(k,k);
	c_temp = c_orig.transpose() * c_orig;

	MatrixXdr T(k,n);
	MatrixXdr c_fn;
	c_fn = c_orig.transpose();
	multiply_y_post(c_fn,k,T,false);

	MatrixXdr M_temp(k,1);
	for(int k_iter=0;k_iter<k;k_iter++){
		double sum=0.0;
		for(int p_iter=0;p_iter<p;p_iter++){
			sum+= c_orig(p_iter,k_iter)*g.get_col_mean(p_iter);
		}
		M_temp(k_iter,0) = sum;
	}

	for(int j=0;j<n;j++){
		MatrixXdr D(k,k);
		MatrixXdr M_to_remove(k,1);
		D = MatrixXdr::Zero(k,k);
		M_to_remove = MatrixXdr::Zero(k,1);
		for(int i=0;i<g.not_O_j[j].size();i++){
			int idx = g.not_O_j[j][i];
			D = D + (c_orig.row(idx).transpose() * c_orig.row(idx));
			M_to_remove = M_to_remove + (c_orig.row(idx).transpose()*g.get_col_mean(idx));
		}
		mu.col(j) = (c_temp-D).inverse() * ( T.col(j) - M_temp + M_to_remove);
	}

	// M step

	MatrixXdr mu_temp(k,k);
	mu_temp = mu * mu.transpose();
	MatrixXdr T1(p,k);
	MatrixXdr mu_fn;
	mu_fn = mu.transpose();
	multiply_y_pre(mu_fn,k,T1,false);
	MatrixXdr mu_sum(k,1);
	mu_sum = MatrixXdr::Zero(k,1);
	for(int j=0;j<n;j++)
		mu_sum += mu.col(j);
	
	for(int i=0;i<p;i++){
		MatrixXdr D(k,k);
		MatrixXdr mu_to_remove(k,1);
		D = MatrixXdr::Zero(k,k);
		mu_to_remove = MatrixXdr::Zero(k,1);
		for(int j=0;j<g.not_O_i[i].size();j++){
			int idx = g.not_O_i[i][j];
			D = D + (mu.col(idx) * mu.col(idx).transpose());
			mu_to_remove = mu_to_remove + (mu.col(idx));
		}
		c_new.row(i) = (((mu_temp-D).inverse()) * (T1.row(i).transpose() -  ( g.get_col_mean(i) * (mu_sum-mu_to_remove)))).transpose();
		double mean;
		mean = g.get_col_sum(i);
		mean = mean -  (c_orig.row(i)*(mu_sum-mu_to_remove))(0,0);
		mean = mean * 1.0 / (n-g.not_O_i[i].size());
		g.update_col_mean(i,mean);
	}
	return c_new;
}

MatrixXdr run_EM(MatrixXdr &c_orig){
	
	if(missing)
		return run_EM_missing(c_orig);
	else
		return run_EM_not_missing(c_orig);
}

void print_vals(){

	HouseholderQR<MatrixXdr> qr(c);
	MatrixXdr Q;
	Q = qr.householderQ() * MatrixXdr::Identity(p,k);
	MatrixXdr q_t(k,p);
	q_t = Q.transpose();
	MatrixXdr b(k,n);
	multiply_y_post(q_t,k,b,true);
	JacobiSVD<MatrixXdr> b_svd(b, ComputeThinU | ComputeThinV);
	MatrixXdr u_l; 
	u_l = b_svd.matrixU();
	MatrixXdr v_l;
	v_l = b_svd.matrixV();
	MatrixXdr u_k;
	MatrixXdr v_k,d_k;
	u_k = u_l.leftCols(k_orig);
	v_k = v_l.leftCols(k_orig);

	ofstream evec_file;
	evec_file.open((string(command_line_opts.OUTPUT_PATH)+string("evecs.txt")).c_str());
	evec_file<< std::setprecision(15) << Q*u_k << endl;
	evec_file.close();
	ofstream eval_file;
	eval_file.open((string(command_line_opts.OUTPUT_PATH)+string("evals.txt")).c_str());
	for(int kk =0 ; kk < k_orig ; kk++)
		eval_file << std::setprecision(15)<< (b_svd.singularValues())(kk)<<endl;
	eval_file.close();

	ofstream proj_file;
	proj_file.open((string(command_line_opts.OUTPUT_PATH) + string("projections.txt")).c_str());
	proj_file << std::setprecision(15)<< v_k<<endl;
	proj_file.close();
	if(debug){
		ofstream c_file;
		c_file.open((string(command_line_opts.OUTPUT_PATH)+string("cvals.txt")).c_str());
		c_file<<c<<endl;
		c_file.close();
		
		d_k = MatrixXdr::Zero(k_orig,k_orig);
		for(int kk =0 ; kk < k_orig ; kk++)
		d_k(kk,kk)  =(b_svd.singularValues())(kk);
		MatrixXdr x_k;
		x_k = d_k * (v_k.transpose());
		ofstream x_file;
		x_file.open((string(command_line_opts.OUTPUT_PATH) + string("xvals.txt")).c_str());
		x_file<<x_k.transpose()<<endl;
		x_file.close();
	}
}

int main(int argc, char const *argv[]){

	clock_t io_begin = clock();
    clock_gettime (CLOCK_REALTIME, &t0);

	pair<double,double> prev_error = make_pair(0.0,0.0);
	double prevnll=0.0;

	parse_args(argc,argv);

	//TODO: Memory Effecient Version
	// memory_efficient = command_line_opts.memory_efficient;
	text_version = command_line_opts.text_version;
    fast_mode = command_line_opts.fast_mode;
	missing = command_line_opts.missing;

	
	if(fast_mode){
		if(text_version)
			g.read_genotype_mailman(command_line_opts.GENOTYPE_FILE_PATH,missing);
		else
			g.read_plink(command_line_opts.GENOTYPE_FILE_PATH,missing);	
	}
	else
		g.read_genotype_naive(command_line_opts.GENOTYPE_FILE_PATH,missing);	

    MAX_ITER =  command_line_opts.max_iterations ; 
	k_orig = command_line_opts.num_of_evec ;
	debug = command_line_opts.debugmode ;
	check_accuracy = command_line_opts.getaccuracy;
	var_normalize = command_line_opts.var_normalize;
	accelerated_em = command_line_opts.accelerated_em;
	k = k_orig + command_line_opts.l;
	k = (int)ceil(k/10.0)*10;
	command_line_opts.l = k - k_orig;
	p = g.Nsnp;
	n = g.Nindv;
	convergence_limit = command_line_opts.convergence_limit;
	bool toStop=false;
	if(convergence_limit!=-1)
		toStop=true;
	srand((unsigned int) time(0));
	c.resize(p,k);
	x.resize(k,n);
	v.resize(p,k);
	means.resize(p,1);

	clock_t io_end = clock();

	//TODO: Initialization of c
	c = MatrixXdr::Random(p,k);


	// Initial intermediate data structures
	blocksize = k;
	int hsegsize = g.segment_size_hori; 	// = log_3(n)
	int hsize = pow(3,hsegsize);		 
	int vsegsize = g.segment_size_ver; 		// = log_3(p)
	int vsize = pow(3,vsegsize);		 

	partialsums = new double [blocksize];
	sum_op = new double[blocksize];

	yint_e = new double [hsize*blocksize];
	yint_m = new double [hsize*blocksize];
	memset (yint_m, 0, hsize*blocksize * sizeof(double));
	memset (yint_e, 0, hsize*blocksize * sizeof(double));

	y_e  = new double*[g.Nindv];
	for (int i = 0 ; i < g.Nindv ; i++) {
		y_e[i] = new double[blocksize];
		memset (y_e[i], 0, blocksize * sizeof(double));
	}

	y_m = new double*[hsegsize];
	for (int i = 0 ; i < hsegsize ; i++)
		y_m[i] = new double[blocksize];

	for(int i=0;i<p;i++)
		means(i,0) = g.get_col_mean(i);
		

	ofstream c_file;
	if(debug){
		c_file.open((string(command_line_opts.OUTPUT_PATH)+string("cvals_orig.txt")).c_str());
		c_file<<c<<endl;
		c_file.close();
		printf("Read Matrix\n");
	}

	cout<<"Running on Dataset of "<<g.Nsnp<<" SNPs and "<<g.Nindv<<" Individuals"<<endl;

	
	clock_t it_begin = clock();
	for(int i=0;i<MAX_ITER;i++){

		MatrixXdr c1,c2,cint,r,v;
		double a,nll;
		if(debug){
			print_time (); 
			cout << "*********** Begin epoch " << i << "***********" << endl;
		}
		if(accelerated_em!=0){
			if(debug){
				print_time();
				cout << "Before EM" << endl;
			}
			c1 = run_EM(c);
			c2 = run_EM(c1);
			if(debug){
				print_time(); 
				cout << "After EM but before acceleration" << endl;
			}
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
				c = cint;
				// c = run_EM(cint);				
			}
		}
		else{
			c = run_EM(c);
		}
		
		if ( accelerated_em == 1 || check_accuracy || toStop) {
			pair<double,double> e = get_error_norm(c);
			prevnll = e.second;
			if(check_accuracy) 
				cout<<"Iteration "<<i+1<<"  "<<e.first<<"  "<<e.second<<endl;
			if(abs(e.first-prev_error.first)<=convergence_limit){
				cout<<"Breaking after "<<i+1<<" iterations"<<endl;
				break;
			}
			prev_error = e;
		}
		if(debug){
			print_time (); 
			cout << "*********** End epoch " << i << "***********" << endl;
		}

	}
	clock_t it_end = clock();

    print_vals();
		
	clock_t total_end = clock();
	double io_time = double(io_end - io_begin) / CLOCKS_PER_SEC;
	double avg_it_time = double(it_end - it_begin) / (MAX_ITER * 1.0 * CLOCKS_PER_SEC);
	double total_time = double(total_end - total_begin) / CLOCKS_PER_SEC;
	cout<<"IO Time:  "<< io_time << "\nAVG Iteration Time:  "<<avg_it_time<<"\nTotal runtime:   "<<total_time<<endl;


	delete[] partialsums;
	delete[] sum_op;
	delete[] yint_e; 
	delete[] yint_m;

	for (int i  = 0 ; i < hsegsize; i++)
		delete[] y_m [i]; 
	delete[] y_m;

	for (int i  = 0 ; i < g.Nindv; i++)
		delete[] y_e[i]; 
	delete[] y_e;

	return 0;
}
