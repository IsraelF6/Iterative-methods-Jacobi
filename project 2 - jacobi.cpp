#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;
 

void jacobi(int n, double* d, double* l,double* u, double* b, double* x0)
{
	const double eps = 1.0e-5;
	double* AX=new double[n];
	double* x_new=new double[n];
	double* Minus_Tj=new double[n];
	double err=0; 
	double* residual=new double[n];
	int N=100;
	
	
	for (int k=1; k<=N; k++)
	{
		for (int i=1; i<n; i++)
		{
			AX[i]=0;						//where A = d - l - u
			AX[i] -= (l[i]*x0[i-1]);
		}
				
		for (int i=0; i<n-1; i++)
		{
			AX[0]=0;
			AX[i] -=  (u[i]*x0[i+1]);
		}
	
		for (int i=0; i<n; i++)
			Minus_Tj[i] = AX[i]/d[i];		//Tj = (l+u)/d

		for (int i=0; i<n; i++)
			AX[i] += (d[i]*x0[i]);		//AX = (d-l-u)x

		for (int i=0; i<n; i++)
			x_new[i] = ((b[i] / d[i]) - Minus_Tj[i]);			// new x vector

		for (int i=0; i<n; i++)
		{
			residual[i] = b[i] - AX[i];				// norm
			err += residual[i]*residual[i];		
		}
		err = sqrt(err);
			cout << "err = " << err << endl;

		if (err < eps)			
		{
			cout << " x = ( ";
			for (int i=0; i<n; i++)
				cout << x_new[i] << " ";
			cout << ")" << endl << "After " << k << " iterations " << endl;
			break;		
		}
		
		else
		{
			for (int i=0; i<n; i++)
				x0[i] = x_new[i];
		}
	
		if (k==N)
		{
			cout << "the procedure was successful," << endl;	//error bound not met
			cout << "x = ( ";
				for (int i=0; i<n; i++)
					cout << x_new[i] << " ";
				cout << ")" << endl << "After " << k << " iterations. " << endl;
		}
	}
}
	

int main() 
{
	int n = 5;
	
	double* d=new double[n];
	double* u=new double[n];
	double* l=new double[n];
	double* XO=new double[n];
	double* b=new double[n];
	
	//matrix d
	for (int i=0; i<n; i++)
		d[i] = 4;

	//matrix u
	for (int i=0; i<n-1; i++)
		u[i] = -1;

	//matrix l
	for (int i=1; i<n; i++)
		l[i] = -1;	

	// matrix b
	for (int i=0; i<n; i++)
	{
		if (2*i<n)
		b[i] = 1;
		else
		b[i]=0;
	}
	
	//set x0={0...0}
	for (int i=0; i<n; i++)			
				XO[i] = 0;

	jacobi(n, d, l, u , b, XO); 

	return 0;
}