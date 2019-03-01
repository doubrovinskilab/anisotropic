/*
Finite difference code for simulating anisotropic contraction

Output written to ./outpt and ./outt.m 
Output files in ./outpt/ store deformation components at given time points
./outt.m stores the maximal x- and y-displacement at a given time-point
*/

#include <iostream>
#include <fstream>

#include <string.h>
#include <math.h>

using namespace std;

char* itoa(int val, int base) ;

int main()
{
	//=============================================================================================================
	// Parameters
	double E = 1.0 ; // Young's modulus
	const double eta = 1.0 ; // Viscous drag

	double sigma = 0.20 ; // Poisson's ratio

	const double alpha = E/2/(1+sigma) ;
	const double beta = E/2/(1-sigma) ;

	const double Lx = 50.0 ; // Embedding domain length
	const double Ly = 50.0 ; // Embedding domain width
	//=============================================================================================================


	//=============================================================================================================
	// Auxiliaries
	const int NumXPoints = 800; // Number of discretization points
	const int NumYPoints = 800;

	const double dt = 0.000010 *5*5 ; // Euler time step

	int i=0, j=0, k=0, m=0 ;

	double dx = Lx/(NumXPoints-1) ;
	double dy = Ly/(NumYPoints-1) ;

	double** ux = new double*[NumXPoints]; // Displacement field
	double** uy = new double*[NumXPoints];

	double** fx = new double*[NumXPoints]; // Active force
	double** fy = new double*[NumXPoints];

	double** dux = new double*[NumXPoints]; // Time-derivative of the displacement field
	double** duy = new double*[NumXPoints];

	for(i=0; i<NumXPoints; i++)
	{
		ux[i] = new double[NumYPoints] ;
		uy[i] = new double[NumYPoints] ;

		dux[i] = new double[NumYPoints] ;
		duy[i] = new double[NumYPoints] ;

		fx[i] = new double[NumYPoints] ;
		fy[i] = new double[NumYPoints] ;
	}

	ofstream* fl = new ofstream("./outt.m");
	(*fl)<<"Ux=[];"<<endl;
	(*fl)<<"Uy=[];"<<endl;
	fl->close();
	delete fl;
	//=============================================================================================================

	//=============================================================================================================
	//
	const double L_x = 10.0 ; // Contractile domain length
	const double L_y = 2.0 ; // Contractile domain width

	const double xi = .30 /4 ; // The width of the "delta-peak" force at the contractile domain boundary

	const double phi0 = 7.50/2 ; // Magnitude of the force peak at the boundary of the contractile domain

	const double x0 = Lx/2-L_x/2 ;
	const double y0 = Ly/2-L_y/2 ;

	double x = 0, y = 0;

	double r1 = 0, r2 = 0 ;

	for(i=0; i<NumXPoints; i++) // Active force at the boundary of the contractile domain
	{
		for(j=0; j<NumYPoints; j++)
		{
			fx[i][j] = 0;
			fy[i][j] = 0;

			x = i*dx ;
			y = j*dy ;

			r1 = x - x0 ;
			r2 = x - (Lx-x0) ;

			if( (y>y0) && (y<(Ly-y0)) )
			{
				fx[i][j] +=  phi0 * exp( -r1*r1/(xi*xi) ) ;
				fx[i][j] += -phi0 * exp( -r2*r2/(xi*xi) ) ;
			}

			r1 = y - y0 ;
			r2 = y - (Ly-y0) ;

			if( (x>x0) && (x<(Lx-x0)) )
			{
				fy[i][j] +=  phi0 * exp( -r1*r1/(xi*xi) ) ;
				fy[i][j] += -phi0 * exp( -r2*r2/(xi*xi) ) ;
			}
		}
	}

	for(i=0; i<NumXPoints; i++) // Initial displacement field set to zero
	{
		for(j=0; j<NumYPoints; j++)
		{
			ux[i][j] = 0;
			uy[i][j] = 0;
		}
	}
	//=============================================================================================================

	for(i=0; ;i++) // Time-evolution
	{
		if( !(i/10000.0-i/10000) )
		{
			cout<<"Iteration # "<<i<<endl;
		}

		for(j=0+1; j<NumXPoints-1; j++)
		{
			for(k=0+1; k<NumYPoints-1; k++)
			{
				int PrvX = j-1;
				int NxtX = j+1;

				int PrvY = k-1;
				int NxtY = k+1;

				dux[j][k] = 0;
				duy[j][k] = 0;

				//==================================================================================================================
				// Finite difference discretization of the governing equations (second-order centered differences)
				dux[j][k] += alpha * (ux[NxtX][k]-2*ux[j][k]+ux[PrvX][k])/(dx*dx) ;
				dux[j][k] += alpha * (ux[j][NxtY]-2*ux[j][k]+ux[j][PrvY])/(dy*dy) ;

				dux[j][k] += beta * (ux[NxtX][k]-2*ux[j][k]+ux[PrvX][k])/(dx*dx) ;
				dux[j][k] += beta * ( (uy[NxtX][NxtY]-uy[NxtX][PrvY])/(2*dy) - (uy[PrvX][NxtY]-uy[PrvX][PrvY])/(2*dy) )/(2*dx) ;


				duy[j][k] += alpha * (uy[NxtX][k]-2*uy[j][k]+uy[PrvX][k])/(dx*dx) ;
				duy[j][k] += alpha * (uy[j][NxtY]-2*uy[j][k]+uy[j][PrvY])/(dy*dy) ;

				duy[j][k] += beta * ( (ux[NxtX][NxtY]-ux[NxtX][PrvY])/(2*dy) - (ux[PrvX][NxtY]-ux[PrvX][PrvY])/(2*dy) )/(2*dx) ;
				duy[j][k] += beta * (uy[j][NxtY]-2*uy[j][k]+uy[j][PrvY])/(dy*dy) ;

				dux[j][k] += fx[j][k] ;
				duy[j][k] += fy[j][k] ;
				//==================================================================================================================
			}
		}

		for(j=0+1; j<NumXPoints-1; j++) // Euler forward step
		{
			for(k=0+1; k<NumYPoints-1; k++)
			{
				ux[j][k] += dt* dux[j][k] ;
				uy[j][k] += dt* duy[j][k] ;
			}
		}

		int qq = 1000 ;
		if( !(i/double(qq)-i/qq) ) // Output the maximal values of the x- and y-components of the displacement (to judge convergence)
		{
			fl = new ofstream("./outt.m",ios::app);

			double Ux = 0;
			double Uy = 0;

			for(j=0; j<NumXPoints; j++)
			{
				for(k=0; k<NumYPoints; k++)
				{
					if( ux[j][k]>Ux )
					{
						Ux = ux[j][k];
					}

					if( uy[j][k]>Uy )
					{
						Uy = uy[j][k];
					}
				}
			}

			(*fl)<<"Ux=[Ux "<<Ux<<"];"<<endl;
			(*fl)<<"Uy=[Uy "<<Uy<<"];"<<endl;

			fl->close();
			delete fl;
		}

		int q = 5000 ;
		if( !(i/double(q)-i/q) ) // Output the fields every q time-steps
		{
			char str1[80]="./outpt/out";
			char str3[80]=".m";
	
			strncat(str1,itoa(i/q,10),10);
			strncat(str1,str3,10);

			ofstream * file = new ofstream(str1);

			// file->precision(20);

			(*file)<<"xs=[];"<<endl;
			for(j=0; j<NumXPoints; j++)
			{
				(*file)<<"xs=[xs; ";
				for(k=0; k<NumYPoints; k++)
				{
					(*file)<<j*dx<<" ";
				}
				(*file)<<"];"<<endl;
			}

			(*file)<<"ys=[];"<<endl;
			for(j=0; j<NumXPoints; j++)
			{
				(*file)<<"ys=[ys; ";
				for(k=0; k<NumYPoints; k++)
				{
					(*file)<<k*dy<<" ";
				}
				(*file)<<"];"<<endl;
			}

			(*file)<<"ux=[];"<<endl;
			for(j=0; j<NumXPoints; j++)
			{
				(*file)<<"ux=[ux; ";
				for(k=0; k<NumYPoints; k++)
				{
					(*file)<<ux[j][k]<<" ";
				}
				(*file)<<"];"<<endl;
			}

			(*file)<<"uy=[];"<<endl;
			for(j=0; j<NumXPoints; j++)
			{
				(*file)<<"uy=[uy; ";
				for(k=0; k<NumYPoints; k++)
				{
					(*file)<<uy[j][k]<<" ";
				}
				(*file)<<"];"<<endl;
			}

			file->close();
			delete file;
		}
	}

	return 1;
}

char* itoa(int val, int base) {
     static char buf[32] = {0};
     int i = 30;
     	if(val==0)
     	{
     	buf[i]='0';
	return &buf[i];
	}
     for (; val && i; --i, val /= base)
         buf[i] = "0123456789abcdef"[val % base];
     return &buf[i+1];
}
