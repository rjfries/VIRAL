#include <mpi.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sstream>
#include <string>
#include <cstring>
#include <float.h>
#include <vector>
#include <sys/types.h>
#include <unistd.h>


#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h> 
#include <gsl/gsl_bspline.h> 
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>


#ifdef __unix__  
#include <fpu_control.h>
#include <fenv.h> 
#endif


using namespace std;

#include "maindef.h"
#include "thermo.h"
#include "alloc.h"
#include "init.h"
#include "utils.h"
#include "pack.h"
#include "boundary.h"
#include "debug.h"
#include "Source.h"
#include "write.h"
#include "root.h"
#ifdef SHAS
#include "hydroshasta.h"
#endif

#ifdef ZAL
#include "hydrozalesak.h"
#endif

#ifdef KT
#include "hydroKT.h"
#endif

void FirstOrder(GRID HydroGrid, double ts, double tau);
void SecondOrder(GRID HydroGrid, double ts, double tau);

int main(int argc, char* argv[])
{
#ifdef __unix__   
	feenableexcept( FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW); 
#endif
		
	
	MPI_Init (&argc, &argv); 
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &size); 
	int i,j,k;
	
	ts = TS;  //set time step	
	tau = TAUSTART;

	int l=1;	
	
	clock_t start =clock();	
	double printFreq = 0.5; //print out every printFreq fm/c
	double tauPrint = TAUSTART;
	
	double sourceFreq = 0.050; //print out source term every sourceFreq fm/c
	double sourcePrint = TAUSTART;
	
	double derFreq = 1; //print out source term every sourceFreq fm/c
	double derPrint = TAUSTART;

	init(tau,ts);

	WriteResultsXYCom(tau,HydroGrid);	
	//WriteSourceXY(tau,HydroGrid);
	
		
	while(1)
	{
	//	FirstOrder(HydroGrid, ts, tau);
		 SecondOrder(HydroGrid, ts, tau);
		
		tau += ts;

		CheckRoot( HydroGrid ,  tau);
		
		if(rank==root)
			cout<<NP*sizeof(cell)*XCM*YCM*ZCM/(1024*1024) <<" mega bytes"<<endl;
			
		if(rank==root)
		{
			cout<<"Time Step is "<<ts <<" @ "; 
			cout<<std::fixed<<std::setprecision(13)<<"TAU -->"<<tau;
			cout<<std::fixed<<std::setprecision(5)<<" This is tau step no. "<<l;
			clock_t now = clock();
			cout<<std::fixed<<std::setw(5)<<std::setprecision(5)<<" at time "<< (((double)(now-start))/CLOCKS_PER_SEC)<< " seconds"<<endl<<endl<<endl<<endl; 	
			fflush(stdout);
		}	

		if( fabs( (tau-tauPrint) - printFreq) < 1e-6)
		{	tauPrint += printFreq;		WriteResultsXYCom(tau,HydroGrid);}
		
		if( fabs( (tau-sourcePrint) - sourceFreq) < 1e-6 )
		{	sourcePrint += sourceFreq;				}//WriteSourceXY(tau-ts,HydroGrid);}

		if( fabs( (tau-ts-derPrint) - derFreq) < 1e-6 )
			derPrint += derFreq;	
						
		l++;
		
		if(tau>25.3)
			break;
	}
	

	if(!rank)
		cout<<endl<<"*************************************** The End ******************************************************"<<endl;
		
	MPI_Finalize();
 }
 
 
 

void FirstOrder(GRID HydroGrid, double ts, double tau)
{
	
	int a =0;
	CheckPhysics(HydroGrid, 0);
	
	CalcSource(HydroGrid, tau, ts);		//for entire grid
	CalcCentreFlux(HydroGrid, tau,ts);
	hydroExplicit(HydroGrid, tau, ts); 			//gets update for Tmunu evrywhere excluding boundary region
	UpdatePrimaryVariablesAtEndOfTimeStep( HydroGrid,ts); //Updates Tmunu everywhere excluding boundary region
	RootSearchForEnVelUsingDerivatives(HydroGrid, tau+ts );		//Finds En,P,V's everywhere excluding boundary region
	
	DebugMSG(HydroGrid);
		
	pack(HydroGrid);  //Exchanges En&Vel  P,Tmunu at the cell interfaces
	boundary(HydroGrid);  // Outflowing boundary condition for all the 9 quantities (4+4+1)	

}



void SOBackUpPrimaryVariables(GRID HydroGrid);
void SORestorePrimaryVariables(GRID HydroGrid);

void SecondOrder(GRID HydroGrid, double ts, double tau)
{
	
	CheckPhysics(HydroGrid, 0);
				
	SOBackUpPrimaryVariables(HydroGrid);
		
	CalcSource(HydroGrid, tau, ts/2);	
	CalcCentreFlux(HydroGrid, tau,ts);
	hydroExplicit(HydroGrid, tau, ts/2); 		
	UpdatePrimaryVariablesAtEndOfTimeStep( HydroGrid,ts/2);		//update primary vars, en,v,P to tau+ts/2		
	RootSearchForEnVelUsingDerivatives(HydroGrid, tau+ts );
	 	
	/**************************************************************/
	
	pack(HydroGrid);
	boundary(HydroGrid);	
	
	CalcSource(HydroGrid, tau+ts/2, ts/2);	
	CalcCentreFlux(HydroGrid, tau+ts/2, ts/2);
	
	/**************************************************************/
	
	
	SORestorePrimaryVariables(HydroGrid);	
	
	hydroExplicit(HydroGrid, tau, ts);	
	UpdatePrimaryVariablesAtEndOfTimeStep( HydroGrid,ts);
	RootSearchForEnVelUsingDerivatives(HydroGrid, tau+ts );
	
	
	/**************************************************************/
	
	
	pack(HydroGrid);
	boundary(HydroGrid);	

	DebugMSG(HydroGrid);
}



void SOBackUpPrimaryVariables(GRID HydroGrid)
{
	int i,j,k;
	
	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		HydroGrid[i][j][k].BackUp[0] = HydroGrid[i][j][k].T00;
		HydroGrid[i][j][k].BackUp[1] = HydroGrid[i][j][k].T10;
		HydroGrid[i][j][k].BackUp[2] = HydroGrid[i][j][k].T20;
		HydroGrid[i][j][k].BackUp[3] = HydroGrid[i][j][k].T30;
		
		HydroGrid[i][j][k].BackUp[4] = HydroGrid[i][j][k].En;
		HydroGrid[i][j][k].BackUp[5] = HydroGrid[i][j][k].Vx;
		HydroGrid[i][j][k].BackUp[6] = HydroGrid[i][j][k].Vy;
		HydroGrid[i][j][k].BackUp[7] = HydroGrid[i][j][k].Ve;
		
		HydroGrid[i][j][k].BackUp[8] = HydroGrid[i][j][k].P;
	}
}



void SORestorePrimaryVariables(GRID HydroGrid)
{
	int i,j,k;
	
	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		HydroGrid[i][j][k].T00 = HydroGrid[i][j][k].BackUp[0];
		HydroGrid[i][j][k].T10 = HydroGrid[i][j][k].BackUp[1];
		HydroGrid[i][j][k].T20 = HydroGrid[i][j][k].BackUp[2];
		HydroGrid[i][j][k].T30 = HydroGrid[i][j][k].BackUp[3];
		
		HydroGrid[i][j][k].En = HydroGrid[i][j][k].BackUp[4];
		HydroGrid[i][j][k].Vx = HydroGrid[i][j][k].BackUp[5];
		HydroGrid[i][j][k].Vy = HydroGrid[i][j][k].BackUp[6];
		HydroGrid[i][j][k].Ve = HydroGrid[i][j][k].BackUp[7];
		
		HydroGrid[i][j][k].P = HydroGrid[i][j][k].BackUp[8];
	}
}
