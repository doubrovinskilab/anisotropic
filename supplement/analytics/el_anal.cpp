/*
Computing the deformation field of an elastic sheet harboring 
a rectangular contractile domain subjected to isotropic uniform
stress.

Computation done by summing analytically derived Fourier-expansion
of the solution.
*/

#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

const double Pi = 3.1415926535897932385;

//-----------------------------------------------------------------------------------------------------------
// Fourier coefficients of the expansion - derived using computer algebra (Maple)
inline double phiBx(int n, int m, double x0, double y0, double Lx, double Ly, double phi0);
inline double phiAy(int n, int m, double x0, double y0, double Lx, double Ly, double phi0);

double Bx(int n, int m, double x0, double y0, double Lx, double Ly, double phi0, double eta, double alpha, double beta, double t);
double Ay(int n, int m, double x0, double y0, double Lx, double Ly, double phi0, double eta, double alpha, double beta, double t);
//-----------------------------------------------------------------------------------------------------------

int main()
{
	//=======================================================================================================
	//
	const double Lx = 50.0; // Embedding domain length
	const double Ly = 50.0; // Embedding domain width

	// contractile domain dimensions: Lx-x0*2 by Ly-y0*2
	const double x0 = 20.0; // Distance between the "left-most" boundary of the embedding domain and that of the contractile domain
	const double y0 = 24.0; // Distance between the "bottom" boundary of the embedding domain and that of the contractile domain

	const double phi0 = 1.50/2; // Magnitude of active stress

	const double eta = 1.0; // Viscous drag on the sheet

	const double E = 1.0; // Young's modulus
	double sigma = 0.20; // Poisson's ratio

	double alpha = E/2/(1.0 + sigma);
	double beta =  E/2/(1.0 - sigma);
	//=======================================================================================================


	//=======================================================================================================
	//
	int i=0, j=0, n=0, m=0;
	int ii = 0;

	double x = 0;
	double y = 0;

	const int NumTPoints = 100/2;

	double* ux_ = new double[NumTPoints]; // x-component of the displacement field (at a particular position) at different time-points
	double* uy_ = new double[NumTPoints]; // y-component of the displacement field (at a particular position) at different time-points

	ofstream file1("outt.m");
	//=======================================================================================================

	for(i=NumTPoints-1; i<NumTPoints; i++) // Loop time-points (only one point looped since asymptotic value is being calsulated)
	{
		cout<<ii<<"\t"<<i<<endl;

		double t_ = i*100.0 ;

		// Left-most point on the contractile domain along the "horizontal" mid-cross-section
		x = x0;
		y = Ly/2;

		ux_[i] = 0;
		uy_[i] = 0;

		// sum 40000x40000 lowest Fourier-modes
		for(n=1; n<40000; n++)
		{
			for(m=0; m<40000; m++)
			{
				ux_[i] += Bx(n,m,x0,y0,Lx,Ly,phi0,eta,alpha,beta,t_)*sin(Pi*n*x/Lx)*cos(Pi*m*y/Ly);
			}
		}

		// Bottom-most point on the contractile domain along the "vertical" mid-cross-section
		x = Lx/2;
		y = y0;

		// sum 40000x40000 lowest Fourier-modes
		for(n=0; n<40000; n++)
		{
			for(m=1; m<40000; m++)
			{
				uy_[i] += Ay(n,m,x0,y0,Lx,Ly,phi0,eta,alpha,beta,t_)*cos(Pi*n*x/Lx)*sin(Pi*m*y/Ly);
			}
		}
	}

	// Output the dimensions of the contractile domain in the deformed state
	file1<<"ux=[ux; ";
	for(i=NumTPoints-1; i<NumTPoints; i++)
	{
		file1<<Lx-2*(x0+ux_[i])<<" ";
	}
	file1<<"];"<<endl;

	file1<<"uy=[uy; ";
	for(i=NumTPoints-1; i<NumTPoints; i++)
	{
		file1<<Ly-2*(y0+uy_[i])<<" ";
	}
	file1<<"];"<<endl;

	return 1;
}

inline double phiBx(int n, int m, double x0, double y0, double Lx, double Ly, double phi0)
{
	if( m!=0 )
	{
		return -4.0*(sin(0.3141592653589793E1*n*x0/Lx)-sin(0.3141592653589793E1*n*(Lx-x0)/Lx))*sin(0.3141592653589793E1*m/Ly*y0)*(pow(-1.0,1.0*m)+1.0)/0.3141592653589793E1/m*phi0/Lx;

	}
	else
	{
		return 2.0*(sin(0.3141592653589793E1*n*x0/Lx)-sin(0.3141592653589793E1*n*(Lx-x0)/Lx))*(Ly-2.0*y0)*phi0/Lx/Ly;

	}
}

inline double phiAy(int n, int m, double x0, double y0, double Lx, double Ly, double phi0)
{
	if(n!=0)
	{
		return -4.0*(sin(0.3141592653589793E1*m/Ly*y0)-sin((Ly-y0)*0.3141592653589793E1*m/Ly))*sin(0.3141592653589793E1*n*x0/Lx)*(pow(-1.0,1.0*n)+1.0)/0.3141592653589793E1/n*phi0/Ly;


	}
	else
	{
		return 2.0*(sin(0.3141592653589793E1*m/Ly*y0)-sin((Ly-y0)*0.3141592653589793E1*m/Ly))*(Lx-2.0*x0)*phi0/Lx/Ly;
	}
}

double Bx(int n, int m, double x0, double y0, double Lx, double Ly, double phi0, double eta, double alpha, double beta, double t)
{
	double phiA = phiAy(n,m,x0,y0,Lx,Ly,phi0);
	double phiB = phiBx(n,m,x0,y0,Lx,Ly,phi0);

return -Ly*Ly*Lx*Lx*(n*n*Ly*Ly*exp(-0.3141592653589793E1*0.3141592653589793E1*eta*(m*m*Lx*Lx+n*n*Ly*Ly)*(beta+alpha)*t/(Lx*Lx)/(Ly*Ly))*
alpha*phiB-phiB*Ly*Ly*alpha*n*n+Ly*Lx*beta*n*m*phiA-n*beta*Ly*exp(-
0.3141592653589793E1*0.3141592653589793E1*eta*alpha*(m*m*Lx*Lx+n*n*Ly*Ly)/(Lx*
Lx)/(Ly*Ly)*t)*m*Lx*phiA-n*Ly*exp(-0.3141592653589793E1*0.3141592653589793E1*
eta*alpha*(m*m*Lx*Lx+n*n*Ly*Ly)/(Lx*Lx)/(Ly*Ly)*t)*m*Lx*phiA*alpha+n*Ly*exp(-
0.3141592653589793E1*0.3141592653589793E1*eta*(m*m*Lx*Lx+n*n*Ly*Ly)*(beta+alpha
)*t/(Lx*Lx)/(Ly*Ly))*alpha*Lx*m*phiA+beta*exp(-0.3141592653589793E1*
0.3141592653589793E1*eta*alpha*(m*m*Lx*Lx+n*n*Ly*Ly)/(Lx*Lx)/(Ly*Ly)*t)*m*m*Lx*
Lx*phiB-phiB*m*m*Lx*Lx*beta-phiB*m*m*Lx*Lx*alpha+exp(-0.3141592653589793E1*
0.3141592653589793E1*eta*alpha*(m*m*Lx*Lx+n*n*Ly*Ly)/(Lx*Lx)/(Ly*Ly)*t)*m*m*Lx*
Lx*phiB*alpha)/alpha/(0.3141592653589793E1*0.3141592653589793E1)/(beta+alpha)/
pow(m*m*Lx*Lx+n*n*Ly*Ly,2.0);

}

double Ay(int n, int m, double x0, double y0, double Lx, double Ly, double phi0, double eta, double alpha, double beta, double t)
{
	double phiA = phiAy(n,m,x0,y0,Lx,Ly,phi0);
	double phiB = phiBx(n,m,x0,y0,Lx,Ly,phi0);

	return -Lx*Lx*Ly*Ly*(-beta*n*n*Ly*Ly*phiA+n*n*Ly*Ly*beta*exp(-0.3141592653589793E1*0.3141592653589793E1*eta*alpha*(m*m*Lx*Lx+n*n*Ly*Ly)/(Lx*
Lx)/(Ly*Ly)*t)*phiA+n*n*Ly*Ly*exp(-0.3141592653589793E1*0.3141592653589793E1*
eta*alpha*(m*m*Lx*Lx+n*n*Ly*Ly)/(Lx*Lx)/(Ly*Ly)*t)*phiA*alpha-alpha*n*n*Ly*Ly*
phiA+Lx*Ly*beta*n*m*phiB-n*Ly*beta*exp(-0.3141592653589793E1*
0.3141592653589793E1*eta*alpha*(m*m*Lx*Lx+n*n*Ly*Ly)/(Lx*Lx)/(Ly*Ly)*t)*m*phiB*
Lx+n*Ly*Lx*m*exp(-0.3141592653589793E1*0.3141592653589793E1*eta*(m*m*Lx*Lx+n*n*
Ly*Ly)*(beta+alpha)*t/(Lx*Lx)/(Ly*Ly))*alpha*phiB-n*Ly*exp(-
0.3141592653589793E1*0.3141592653589793E1*eta*alpha*(m*m*Lx*Lx+n*n*Ly*Ly)/(Lx*
Lx)/(Ly*Ly)*t)*m*phiB*Lx*alpha-phiA*Lx*Lx*alpha*m*m+Lx*Lx*m*m*exp(-
0.3141592653589793E1*0.3141592653589793E1*eta*(m*m*Lx*Lx+n*n*Ly*Ly)*(beta+alpha
)*t/(Lx*Lx)/(Ly*Ly))*alpha*phiA)/alpha/(0.3141592653589793E1*
0.3141592653589793E1)/(beta+alpha)/pow(m*m*Lx*Lx+n*n*Ly*Ly,2.0);
}
