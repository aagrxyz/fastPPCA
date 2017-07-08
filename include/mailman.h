#ifndef MAILMAN_H

#define MAILMAN_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <bits/stdc++.h>

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

	void fastmultiply2 (int m, int n , std::vector<int> &p, Eigen::MatrixXd &x, double *yint, double *y){

		memset (yint, 0, pow(3,m) * sizeof(yint));

		for (int i = 0 ; i < n; i++)
			yint[p[i]] += x(i,0);

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

}


#endif
