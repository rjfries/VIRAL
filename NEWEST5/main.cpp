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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>


#ifdef __unix__  
#include <fpu_control.h>
#include <fenv.h> 
#endif


using namespace std;

#include "maindef.h"
#include "thermo.h"
#include "alloc.h"
#include "init.h"

#ifdef GINIT
#include "ginit.h"
#endif

#include "utils.h"
#include "pack.h"
#include "boundary.h"
#include "debug.h"
#include "Source.h"
#include "write.h"
#include "multiroot.h"
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
	//~ feenableexcept( FE_INVALID|FE_DIVBYZERO); 
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
	double printFreq = PFREQ; //print out every printFreq fm/c
	double tauPrint = TAUSTART;
	
	double sourceFreq = 0.050; //print out source term every sourceFreq fm/c
	double sourcePrint = TAUSTART;
	
	double derFreq = 1; //print out source term every sourceFreq fm/c
	double derPrint = TAUSTART;

	init(tau,ts);

	WriteResultsXY(tau,HydroGrid);	
	//~ WriteResultsXYCom(tau,HydroGrid);	
	//~ WriteTempXYCom(tau,HydroGrid);	
	//~ WriteSourceXY(tau,HydroGrid);
	
#ifdef BJORKEN
		ofstream myfile;
		myfile.open ("temp.txt");
#endif		

    
	while(1)
	{
		 
		//~ FirstOrder(HydroGrid, ts, tau);
		SecondOrder(HydroGrid, ts, tau);
		
#ifdef BJORKEN
		myfile<<tau<< "  "<< fmtoMev(FT(HydroGrid[XCM/2][YCM/2][ZCM/2].En , HydroGrid[XCM/2][YCM/2][ZCM/2].r))/1000<<endl;
#endif
		 
		
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
		{	tauPrint += printFreq;	WriteResultsXY(tau,HydroGrid);}//WriteTempXYCom(tau,HydroGrid) ;}//WriteResultsXYCom(tau,HydroGrid);}//	WriteTempXYCom(tau,HydroGrid) ;}//
		
		if( fabs( (tau-sourcePrint) - sourceFreq) < 1e-6 )
		{	sourcePrint += sourceFreq;				}//WriteSourceXY(tau-ts,HydroGrid);}

		if( fabs( (tau-ts-derPrint) - derFreq) < 1e-6 )
			derPrint += derFreq;	
						
		l++;
		
		
#ifdef GUBSER
		if(tau>11)
			break;
#endif
		if(tau>25.3)
			break;
	}
	

	if(!rank)
		cout<<endl<<"*************************************** The End ******************************************************"<<endl;
		
	MPI_Finalize();
 }
 
 
 

void FirstOrder(GRID HydroGrid, double ts, double tau)
{
	
	CheckPhysics(HydroGrid, 0);
	
	CalcSource(HydroGrid, tau, ts);		//for entire grid	

#ifdef KT
	CalcCentreFlux(HydroGrid, tau);
#endif


	hydroExplicit(HydroGrid, tau, ts); 			//gets update for PV,pi and PI everywhere excluding boundary region
	UpdatePrimaryVariablesAtEndOfTimeStep( HydroGrid, tau, ts); //update PV's, pi and PI to "tau+ts"	from "tau"
	MultiRootSearchForEnVelUsingDerivatives(HydroGrid, tau+ts );	//Finds En,P,V's, 4vel, everywhere excluding boundary region
	
	DebugMSG(HydroGrid);
	pack(HydroGrid);  //Exchanges En&Vel, pi and PI and updates P,4vel,Tmunu at the cell interfaces
	boundary(HydroGrid);  // Outflowing boundary condition for everything
}


/*
 * 
 * 
 * SECOND ORDER
 * 
 * 
 */
void SOBackUpPrimaryVariables(GRID HydroGrid);
void SORestorePrimaryVariables(GRID HydroGrid);


void SecondOrder(GRID HydroGrid, double ts, double tau)
{
	
	CheckPhysics(HydroGrid, 0);
	
				
	SOBackUpPrimaryVariables(HydroGrid);
	
	CalcSource(HydroGrid, tau, ts);	  //calculating source after a previous "ts" timestep at time "tau". "ts" needed for time derivative of 4 velocity	
#ifdef KT		
	CalcCentreFlux(HydroGrid, tau);  
#endif
	
	hydroExplicit(HydroGrid, tau, ts/2); 	
	
	
	UpdatePrimaryVariablesAtEndOfTimeStep( HydroGrid, tau, ts/2);		//update primary vars, pi and PI to "tau+ts/2"	from "tau"
	MultiRootSearchForEnVelUsingDerivatives(HydroGrid , tau + ts/2);   // find en,vel at tau+ts/2and update 4 vel and Pressure. Before updating 4 vel save the current one (at tau) to prevu[4]
    CheckPhysics(HydroGrid, 1);
	/**************************************************************/
	
	pack(HydroGrid );    // update ghost cells
	boundary(HydroGrid);	//implement boundary condition	
	CalcSource(HydroGrid, tau+ts/2, ts/2);	//calculating source after a previous "ts/2" timestep now at time "tau+ts/2"
#ifdef KT		
	CalcCentreFlux(HydroGrid, tau+ts/2);  // calculate centrafluxes at the middle of time step. for second order accuracy in time integration
#endif	
	/**************************************************************/
	
	
	SORestorePrimaryVariables(HydroGrid);	 // now revert back the PV, pi , PI , en, vel and 4 vel to tau
	hydroExplicit(HydroGrid, tau, ts);	// now advance full time step using the middle of time step values of source terms and centra fluxes
	UpdatePrimaryVariablesAtEndOfTimeStep( HydroGrid,tau, ts);  // //update primary vars, pi and PI to "tau+ts" from "tau"
	MultiRootSearchForEnVelUsingDerivatives(HydroGrid , tau + ts); // find en,vel at tau+ts and update 4 vel and Pressure to "tau+ts". Before updating 4 vel save the current one (again at tau) to prevu[4].
	CheckPhysics(HydroGrid, 2);
	
	/**************************************************************/

	
	pack(HydroGrid);
	boundary(HydroGrid);	

	DebugMSG(HydroGrid);
}



void SOBackUpPrimaryVariables(GRID HydroGrid)
{
	int i,j,k,l;
	
	
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
		
		HydroGrid[i][j][k].BackUp[8] = HydroGrid[i][j][k].u[0];
		HydroGrid[i][j][k].BackUp[9] = HydroGrid[i][j][k].u[1];
		HydroGrid[i][j][k].BackUp[10] = HydroGrid[i][j][k].u[2];
		HydroGrid[i][j][k].BackUp[11] = HydroGrid[i][j][k].u[3];
		
		HydroGrid[i][j][k].BackUp[12] = HydroGrid[i][j][k].prevu[0];
		HydroGrid[i][j][k].BackUp[13] = HydroGrid[i][j][k].prevu[1];
		HydroGrid[i][j][k].BackUp[14] = HydroGrid[i][j][k].prevu[2];
		HydroGrid[i][j][k].BackUp[15] = HydroGrid[i][j][k].prevu[3];
		
		HydroGrid[i][j][k].BackUp[4*VARN] = HydroGrid[i][j][k].P;
		
		for(l=0;l<Npi;l++)
			HydroGrid[i][j][k].BackUp[4*VARN+1+l] = HydroGrid[i][j][k].pi[l];
			
		HydroGrid[i][j][k].BackUp[4*VARN+1+Npi] = HydroGrid[i][j][k].PI;
	}
}



void SORestorePrimaryVariables(GRID HydroGrid)
{
	int i,j,k,l;
	
	
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
		
		HydroGrid[i][j][k].u[0] = HydroGrid[i][j][k].BackUp[8];
		HydroGrid[i][j][k].u[1] = HydroGrid[i][j][k].BackUp[9];
		HydroGrid[i][j][k].u[2] = HydroGrid[i][j][k].BackUp[10];
		HydroGrid[i][j][k].u[3] = HydroGrid[i][j][k].BackUp[11];
		
		HydroGrid[i][j][k].prevu[0] = HydroGrid[i][j][k].BackUp[12];
		HydroGrid[i][j][k].prevu[1] = HydroGrid[i][j][k].BackUp[13];
		HydroGrid[i][j][k].prevu[2] = HydroGrid[i][j][k].BackUp[14];
		HydroGrid[i][j][k].prevu[3] = HydroGrid[i][j][k].BackUp[15];
		
		HydroGrid[i][j][k].P = HydroGrid[i][j][k].BackUp[4*VARN];
		
		for(l=0;l<Npi;l++)
			HydroGrid[i][j][k].pi[l] = HydroGrid[i][j][k].BackUp[4*VARN+1+l];
			
		HydroGrid[i][j][k].PI  = HydroGrid[i][j][k].BackUp[4*VARN+1+Npi];
	}
}
