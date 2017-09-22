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
		int d = pow(3,m);
		

		for (int j  = 0 ;  j < m ; j++)  {
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

	// Efficient versions working on full matrix
	void fastmultiply2 (int m, int n , int k, std::vector<int> &p, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &x, double *yint, double *c, double **y){
		
		for (int i = 0 ; i < n; i++)  {
			int l = p[i]  ;
			for (int j = 0 ; j < k ; j ++)
				yint[l*k + j] += x(i,j);

		}

		int d = pow(3,m);
		for (int j  = 0 ;  j < m ; j++)  {
			d =d /3;
			for (int l = 0; l < k ; l++)
				c [l] = 0 ; 
			for (int i = 0 ; i < d; i++) { 
				for (int l = 0; l < k ; l++){
					double z1 = yint[l + (i + d)*k];
					double z2 = yint[l + (i +2*d)*k];
					yint[l+(i+d)*k] = 0;
					yint[l+(i+2*d)*k] = 0 ;
					yint[l+i*k] = yint[l+i*k] + z1 + z2;
					c[l] += (z1 + 2*z2);
				}
			}
			for (int l = 0; l < k ; l++)
				y[j][l] = c[l];
		}
		for (int l = 0; l < k ; l++)
			yint[l] = 0;
	}

	void fastmultiply3 (int m, int n , int k, std::vector<unsigned> &p, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &x, double *yint, double *c, double **y,int Nbits){
		
		for (int i = 0 ; i < n; i++)  {
			int l = extract_from_arr(i,Nbits,p);
			for (int j = 0 ; j < k ; j ++)
				yint[l*k + j] += x(i,j);

		}

		int d = pow(3,m);
		for (int j  = 0 ;  j < m ; j++)  {
			d =d /3;
			for (int l = 0; l < k ; l++)
				c [l] = 0 ; 
			for (int i = 0 ; i < d; i++) { 
				for (int l = 0; l < k ; l++){
					double z1 = yint[l + (i + d)*k];
					double z2 = yint[l + (i +2*d)*k];
					yint[l+(i+d)*k] = 0;
					yint[l+(i+2*d)*k] = 0 ;
					yint[l+i*k] = yint[l+i*k] + z1 + z2;
					c[l] += (z1 + 2*z2);
				}
			}
			for (int l = 0; l < k ; l++)
				y[j][l] = c[l];
		}
		for (int l = 0; l < k ; l++)
			yint[l] = 0;
	}


}


#endif
