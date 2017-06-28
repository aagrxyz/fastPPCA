#ifndef MAILMAN_H

#define MAILMAN_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <bits/stdc++.h>

namespace mailman {

	void fastmultiply ( int m, int n, std::vector<int> &p, Eigen::MatrixXf &x, double *yint, double *y) {

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

}


#endif
