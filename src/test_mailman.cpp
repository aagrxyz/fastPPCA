#include <bits/stdc++.h>
#include "mailman.h"
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/LU"
#include "Eigen/SVD"

using namespace std;

int main()
{

    vector<int> p;
    int n = 20;
    int m = ceil(log(n)/log(3));
    p.resize(n);
    for(int i=0;i<n;i++)
        p[i] = (5*i+1)%((int)pow(3,m));
    

    Eigen::MatrixXd x(n,1);
    Eigen::MatrixXd A(m,n);
    A = Eigen::MatrixXd::Zero(m,n);

    for(int i=0;i<n;i++)
        x(i,0) = i*0.33*i*i;
    
    double *y_int = new double[(int)pow(3,m)];

    double *y = new double[m];

    mailman::fastmultiply2(m,n,p,x,y_int,y);

    for(int i=0;i<n;i++){
        int res = p[i];
        int j=m-1;
        while(res>0 && j>=0){
            A(j,i) = res%3;
            res = res/3;
            j--;
        }
    }


    cout<<"A is "<<endl;
    cout<<A<<endl;


    cout<<"p is"<<endl;
    for(int i=0;i<n;i++)
        cout<<p[i]<<endl;
    
    

    cout<<"x is"<<endl;

    cout<<x<<endl;

    cout<<"y is"<<endl;
    for(int i=0;i<m;i++)
        cout<<y[i]<<endl;
    
    return 0;
}