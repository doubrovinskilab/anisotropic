/*
Simulation of an ablation experiment
The mesh is stored in x.dat y.dat tri_1.dat tri_2.dat tri_3.dat,
Rows in x.dat and y.dat are vertex coordinates
Rows in tri_1.dat tri_2.dat tri_3.dat are indices of vertices making up a triangle

After prescribed time, vertices within a circular "hole" are effectively removed
to simulate a laser ablation

Output is written into ./output folder
*/

#include <iostream>
#include <fstream>
using namespace std;
#include <math.h>
#include <stdlib.h>
#include <string.h>

char* itoa(int val, int base);


int main()
{
	//=============================================================================================
	// Parameters

	const double mu_0 = 3.5; // tension within the contractile domain
	const double E = 50.0 ; // spring constant (determining the Young's modulus)

	const double dt = 0.00001/2 ; // time step

	const double Tabl = 5.0; // time of "ablation"
	//=============================================================================================

	//=============================================================================================
	// Auxiliaries

	double L0 = 0;

	double* x;
	double* y;

	double* fx;
	double* fy;

	int* tri_1;
	int* tri_2;
	int* tri_3;

	int i=0, j=0, k=0, m=0;

	char str1[100];
	char str2[100];
	char str3[100];

	ifstream* rd1;
	ifstream* rd2;
	ifstream* rd3;

	int NumPoints = 0;

	//---------------------------------------------------------------------------------------------
	// Read files with the triangluated mesh
	/* 
	After files are read, (x[tri_1[i]] y[tri_1[i]]) (x[tri_2[i]] y[tri_2[i]]) (x[tri_3[i]] y[tri_3[i]])
	are the coordinates of the vertices of the i:th triagle (traversed counterclockwise)
	*/
	rd1 = new ifstream("./meshfiles/x.dat");
	while( rd1->getline(str1,100) )
	{
		NumPoints++;
	}
	rd1->close();
	delete rd1;

	cout<<NumPoints<<endl;

	x = new double[NumPoints];
	y = new double[NumPoints];

	fx = new double[NumPoints];
	fy = new double[NumPoints];

	rd1 = new ifstream("./meshfiles/x.dat");
	rd2 = new ifstream("./meshfiles/y.dat");

	for(i=0; i<NumPoints; i++)
	{
		rd1->getline(str1,100);
		rd2->getline(str2,100);

		x[i] = atof(str1);
		y[i] = atof(str2);
	}
	rd1->close();
	delete rd1;

	rd2->close();
	delete rd2;

	int NumTriangles = 0;
	rd1 = new ifstream("./meshfiles/tri_1.dat");
	while( rd1->getline(str1,100) )
	{
		NumTriangles++;
	}
	rd1->close();
	delete rd1;

	cout<<NumTriangles<<endl;

	double** L0s_ = new double*[NumPoints];

	for(i=0; i<NumPoints; i++)
	{
		L0s_[i] = new double[NumPoints];
	}

	tri_1 = new int[NumTriangles];
	tri_2 = new int[NumTriangles];
	tri_3 = new int[NumTriangles];

	rd1 = new ifstream("./meshfiles/tri_1.dat");
	rd2 = new ifstream("./meshfiles/tri_2.dat");
	rd3 = new ifstream("./meshfiles/tri_3.dat");

	for(i=0; i<NumTriangles; i++)
	{
		rd1->getline(str1,100);
		rd2->getline(str2,100);
		rd3->getline(str3,100);

		tri_1[i] = int( atof(str1) )-1;
		tri_2[i] = int( atof(str2) )-1;
		tri_3[i] = int( atof(str3) )-1;
	}
	//---------------------------------------------------------------------------------------------

	double vx_hat1=0, vy_hat1=0;

	double vx1=0, vy1=0, vz1=0;
	double vx2=0, vy2=0, vz2=0;

	double lambda = 0;

	//---------------------------------------------------------------------------------------------
	// Store initial (rest) lengths of all edges in L0s_ 
	for(i=0; i<NumTriangles; i++)
	{
		//-----------------------------------------------------------------------------------------
		//
		int idx_1 = tri_2[i];
		int idx_2 = tri_1[i];

		vx1 = x[tri_2[i]]-x[tri_1[i]];
		vy1 = y[tri_2[i]]-y[tri_1[i]];

		L0s_[idx_1][idx_2] = sqrt(vx1*vx1+vy1*vy1) ;
		//-----------------------------------------------------------------------------------------

		//-----------------------------------------------------------------------------------------
		//
		idx_1 = tri_3[i];
		idx_2 = tri_2[i];

		vx1 = x[tri_3[i]]-x[tri_2[i]];
		vy1 = y[tri_3[i]]-y[tri_2[i]];

		L0s_[idx_1][idx_2] = sqrt(vx1*vx1+vy1*vy1) ;
		//-----------------------------------------------------------------------------------------

		//-----------------------------------------------------------------------------------------
		//
		idx_1 = tri_1[i];
		idx_2 = tri_3[i];

		vx1 = x[tri_1[i]]-x[tri_3[i]];
		vy1 = y[tri_1[i]]-y[tri_3[i]];

		L0s_[idx_1][idx_2] = sqrt(vx1*vx1+vy1*vy1) ;
		//-----------------------------------------------------------------------------------------
	}
	//---------------------------------------------------------------------------------------------

	//---------------------------------------------------------------------------------------------
	// Specify the contractile domain
	/*
	WhichTeiangles[i] set to 1 iff all vertices of i:th 
	triangle fall within the conractile domain
	*/
	double* WhichTriangles = new double[NumTriangles];

	for(i=0; i<NumTriangles; i++)
	{
		WhichTriangles[i] = 0;

		if( !( (x[tri_1[i]]>-0.605) && (x[tri_1[i]]<0.605) ) ) continue;
		if( !( (x[tri_2[i]]>-0.605) && (x[tri_2[i]]<0.605) ) ) continue;
		if( !( (x[tri_3[i]]>-0.605) && (x[tri_3[i]]<0.605) ) ) continue;

		if( !( (y[tri_1[i]]>-0.3) && (y[tri_1[i]]<0.3) ) ) continue;
		if( !( (y[tri_2[i]]>-0.3) && (y[tri_2[i]]<0.3) ) ) continue;
		if( !( (y[tri_3[i]]>-0.3) && (y[tri_3[i]]<0.3) ) ) continue;

		WhichTriangles[i] = 1;
	}

	/*
	EdgesUnderTension[i][j]=1 iff i:th and j:th vertices joint by an edge
	EdgeUnderTension_1[i] EdgeUnderTension_2[i] are vertices adjacent to i:th edge
	NumEdgesUnderTension the number of edges
	*/

	int** EdgesUnderTension = new int*[NumPoints];
	for(i=0; i<NumPoints; i++)
	{
		EdgesUnderTension[i] = new int[NumPoints];
	}

	for(i=0; i<NumPoints; i++)
	{
		for(j=0; j<NumPoints; j++)
		{
			EdgesUnderTension[i][j] = 0;
		}
	}

	int NumEdgesUnderTension = 0;
	int* EdgeUnderTension_1 = new int[NumTriangles*3];
	int* EdgeUnderTension_2 = new int[NumTriangles*3];

	for(i=0; i<NumTriangles; i++)
	{
		if( WhichTriangles[i] )
		{
			int vrtx_1 = tri_1[i];
			int vrtx_2 = tri_2[i];
			int vrtx_3 = tri_3[i];

			if( ! EdgesUnderTension[vrtx_1][vrtx_2] )
			{
				EdgesUnderTension[vrtx_1][vrtx_2] = 1;
				EdgesUnderTension[vrtx_2][vrtx_1] = 1;

				EdgeUnderTension_1[NumEdgesUnderTension] = vrtx_1;
				EdgeUnderTension_2[NumEdgesUnderTension] = vrtx_2;

				NumEdgesUnderTension++;
			}

			if( ! EdgesUnderTension[vrtx_2][vrtx_3] )
			{
				EdgesUnderTension[vrtx_2][vrtx_3] = 1;
				EdgesUnderTension[vrtx_3][vrtx_2] = 1;

				EdgeUnderTension_1[NumEdgesUnderTension] = vrtx_2;
				EdgeUnderTension_2[NumEdgesUnderTension] = vrtx_3;

				NumEdgesUnderTension++;
			}

			if( ! EdgesUnderTension[vrtx_3][vrtx_1] )
			{
				EdgesUnderTension[vrtx_3][vrtx_1] = 1;
				EdgesUnderTension[vrtx_1][vrtx_3] = 1;

				EdgeUnderTension_1[NumEdgesUnderTension] = vrtx_1;
				EdgeUnderTension_2[NumEdgesUnderTension] = vrtx_3;

				NumEdgesUnderTension++;
			}
		}
	}
	//---------------------------------------------------------------------------------------------

	//---------------------------------------------------------------------------------------------
	// The hole created by the ablation
	// IsInHole[i]=1 iff i:th vertex was removed by ablation, initialized to zeros
	int NumPointsInHole = 0;
	int* PointsInHole = new int[NumPoints];
	int* IsInHole = new int[NumPoints];

	for(i=0; i<NumPoints; i++)
	{
		IsInHole[i] = 0;
	}
	//---------------------------------------------------------------------------------------------

	//=============================================================================================

	//=============================================================================================
	//=============================================================================================

	for(i=0; ;i++)
	{
		if( !(i/1000.0-i/1000) ) {
			cout<<"Iteration # "<<i<<endl;
		}

		#include "outpt.h" // produce output

		for(j=0; j<NumPoints; j++) // initialize forces at the vertices to zero
		{
			fx[j]=0;
			fy[j]=0;
		}

		//-----------------------------------------------------------------------------------------
		// After ablation at time Tabl, remove the ablated points within the created "hole"
		if( i*dt>Tabl && (!NumPointsInHole) )
		{
			for(j=0; j<NumPoints; j++)
			{
				double r = sqrt( x[j]*x[j] + y[j]*y[j] );

				if( r<0.1/4 )
				{
					PointsInHole[NumPointsInHole] = j;
					IsInHole[j] = 1;

					NumPointsInHole++;
				}
			}
		}
		//-----------------------------------------------------------------------------------------

		//-----------------------------------------------------------------------------------------
		// Elasticity
		for(j=0; j<NumTriangles; j++)
		{
			if( IsInHole[tri_1[j]] || IsInHole[tri_2[j]] || IsInHole[tri_3[j]] ) continue; // disregard points within the "hole" after ablation

			//-------------------------------------------------------------------------------------
			L0 = L0s_[tri_2[j]][tri_1[j]];

			vx1 = x[tri_2[j]]-x[tri_1[j]];
			vy1 = y[tri_2[j]]-y[tri_1[j]];

			lambda = sqrt( vx1*vx1+vy1*vy1 );

			vx_hat1 = vx1/lambda;
			vy_hat1 = vy1/lambda;

			fx[tri_1[j]] += -E*(L0-lambda) *vx_hat1 ;
			fy[tri_1[j]] += -E*(L0-lambda) *vy_hat1 ;

			fx[tri_2[j]] +=  E*(L0-lambda) *vx_hat1 ;
			fy[tri_2[j]] +=  E*(L0-lambda) *vy_hat1 ;
			//-------------------------------------------------------------------------------------

			//-------------------------------------------------------------------------------------
			L0 = L0s_[tri_3[j]][tri_2[j]];

			vx1 = x[tri_3[j]]-x[tri_2[j]];
			vy1 = y[tri_3[j]]-y[tri_2[j]];

			lambda = sqrt( vx1*vx1+vy1*vy1 );

			vx_hat1 = vx1/lambda;
			vy_hat1 = vy1/lambda;

			fx[tri_2[j]] += -E*(L0-lambda)*vx_hat1 ;
			fy[tri_2[j]] += -E*(L0-lambda)*vy_hat1 ;

			fx[tri_3[j]] +=  E*(L0-lambda)*vx_hat1 ;
			fy[tri_3[j]] +=  E*(L0-lambda)*vy_hat1 ;
			//-------------------------------------------------------------------------------------

			//-------------------------------------------------------------------------------------
			L0 = L0s_[tri_1[j]][tri_3[j]];

			vx1 = x[tri_1[j]]-x[tri_3[j]];
			vy1 = y[tri_1[j]]-y[tri_3[j]];

			lambda = sqrt( vx1*vx1+vy1*vy1 );

			vx_hat1 = vx1/lambda;
			vy_hat1 = vy1/lambda;

			fx[tri_3[j]] += -E*(L0-lambda)*vx_hat1 ;
			fy[tri_3[j]] += -E*(L0-lambda)*vy_hat1 ;

			fx[tri_1[j]] +=  E*(L0-lambda)*vx_hat1 ;
			fy[tri_1[j]] +=  E*(L0-lambda)*vy_hat1 ;
			//-------------------------------------------------------------------------------------
		}
		//-----------------------------------------------------------------------------------------

		//-----------------------------------------------------------------------------------------
		// Active stress
		for(j=0; j<NumEdgesUnderTension; j++)
		{
			int vrtx_1 = EdgeUnderTension_1[j];
			int vrtx_2 = EdgeUnderTension_2[j];

			if( IsInHole[vrtx_1] || IsInHole[vrtx_2] ) continue; // disregard points within the "hole" after the ablation

			double x1 = x[vrtx_1];
			double y1 = y[vrtx_1];

			double x2 = x[vrtx_2];
			double y2 = y[vrtx_2];

			double r2 = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1);
			double r = sqrt(r2);

			vx1 = (x2-x1)/r;
			vy1 = (y2-y1)/r;

			fx[vrtx_1] +=  mu_0 * vx1;
			fy[vrtx_1] +=  mu_0 * vy1;

			fx[vrtx_2] += -mu_0 * vx1;
			fy[vrtx_2] += -mu_0 * vy1;
		}
		//-----------------------------------------------------------------------------------------

		//-----------------------------------------------------------------------------------------
		// Euler step
		// Points at the boundaries not updated
		for(j=0; j<NumPoints; j++)
		{
			if( x[j]< -(1.0-0.0010) ) continue;
			if( x[j]>  (1.0-0.0010) ) continue;

			if( y[j]< -(1.0-0.0010) ) continue;
			if( y[j]>  (1.0-0.0010) ) continue;

			x[j] += fx[j]*dt;
			y[j] += fy[j]*dt;
		}
		//-----------------------------------------------------------------------------------------
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
