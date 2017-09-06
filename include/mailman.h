#ifndef MAILMAN_H

#define MAILMAN_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <bits/stdc++.h>
#include "storage.h"

namespace mailman {

	void fastmultiply ( int m, int n, std::vector<int> &p, Eigen::MatrixXd &x, double *yint, double *y) {

		memset (yint, 0, pow(2,m) * sizeof(yint));

		for (int i = 0 ; i < n; i++)
			yint[p[i]] += x(i,0);

		for (int j  = 0 ;  j < m ; j++)  {
			int d = pow(2.,m-j-1);
			double c = 0 ; 
			for (int i = 0 ; i < d; i++) { 
		        double z = yint[i+d];
				yint[i] = yint[i] + z;
				c += z;
			}
			y[j] = c;
		}
	}



	void fastmultiply2 (int m, int n , std::vector<int> &p, Eigen::MatrixXd &x, double *yint, double *y, int col){

		memset (yint, 0, pow(3,m) * sizeof(yint));

		for (int i = 0 ; i < n; i++)
			yint[p[i]] += x(i,col);

		int d = pow(3,m);

		for (int j  = 0 ;  j < m ; j++)  {
			// int d = pow(3.,m-j-1);
			d = d/3;
			double c = 0 ; 
			for (int i = 0 ; i < d; i++) { 
		        double z1 = yint[i+d];
		        double z2 = yint[i+(2*d)];
				yint[i] = yint[i] + z1 + z2;
				c += (z1 + 2*z2);
			}
			y[j] = c;
		}

	}



	void fastmultiply3 (int m, int n , std::vector<unsigned> &p, Eigen::MatrixXd &x, double *yint, double *y,int Nbits,int col = 0){

		memset (yint, 0, pow(3,m) * sizeof(yint));

		for (int i = 0 ; i < n; i++){
			int temp = extract_from_arr(i,Nbits,p);
			yint[temp] += x(i,col);
		}

		for (int j  = 0 ;  j < m ; j++)  {
			int d = pow(3.,m-j-1);
			double c = 0 ; 
			for (int i = 0 ; i < d; i++) { 
		        double z1 = yint[i+d];
		        double z2 = yint[i+(2*d)];
				yint[i] = yint[i] + z1 + z2;
				c += (z1 + 2*z2);
			}
			y[j] = c;
		}

	}


/*
	void fastmultiply2 (int m, int n , std::vector<int> &p, Eigen::MatrixXd &x, double *yint1, double *y, int col = 0){

		int m0  = pow(3*1.,m);
//		Eigen::VectorXd yint (m0);
		Eigen::VectorXd yint =Eigen::VectorXd::Zero(m0);

		for (int i = 0 ; i < n; i++)
			yint(p[i]) += x(i,col);

		for (int j  = 0 ;  j < m ; j++)  {
			int d = pow(3.,m-j-1);
			yint.head(d) = yint.segment (d,d)+yint.segment(2*d,d);
			double c  = yint.segment (d,d).sum()+ 2* yint.segment(2*d,d).sum();
			y[j] = c;
		}

	}*/

/*
	void fastmultiply2 (int m, int n , int k, std::vector<int> &p, Eigen::MatrixXd &x, double **yint, double **y){

//		memset (yint, 0, pow(3,m) * sizeof(yint)  * k);
		double *c = new double[k];

		for (int i = 0 ; i < n; i++) { 
			for (int col = 0 ; col <  k ; col++) {
				yint[p[i]][col] += x(i,col);
			}
		}

		for (int j  = 0 ;  j < m ; j++)  {
			int d = pow(3.,m-j-1);
			for (int col = 0 ; col < k; col++)
				c[col] = 0 ;
			for (int i = 0 ; i < d; i++) { 
				for (int col = 0; col < k; col++) {
		        double z1 = yint[i+d][col];
		        double z2 = yint[i+(2*d)][col];
				yint[i][col] = yint[i][col] + z1 + z2;
				c[col] += (z1 + 2*z2);
				}
			}
			for (int col = 0 ; col < k; col++)
				y[j][col] = c[col];
		}

	}*/





}


#endif
