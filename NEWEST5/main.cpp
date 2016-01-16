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


#ifdef SHAS
	#define CON
#endif
#define GEVFM 0.1973 
#define PIE 3.141592653589793


using namespace std; 
#include "s95p.h" 
#include "maindef.h"


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



#include "time.h"

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
	
#ifdef LBI
	WriteResultsXY(tau,HydroGrid);	
#else
	WriteResults(tau,HydroGrid);	 
#endif	
	
#if defined BJORKEN || BULKTEST
		ofstream myfile;
		myfile.open ("temp.txt");
#endif		

    
	while(1)
	{
	
#ifdef RK1	 
		FirstOrder(HydroGrid, tau, ts);
#endif		
#ifdef RK2	 
		TVDRK2(HydroGrid, tau, ts);
#endif		
#ifdef RK3	 
		TVDRK3(HydroGrid, tau, ts);
#endif		

#ifdef BJORKEN
		myfile<<tau<< "  "<< fmtoMev(FT(HydroGrid[XCM/2][YCM/2][ZCMA/2].En))/1000<<endl;
#endif
		 
		
#ifdef BULKTEST
		myfile<<tau<< "  "<<  HydroGrid[XCM/2][YCM/2][ZCMA/2].PI<<endl;
#endif
		 
		
		tau += ts;

		CheckRoot( HydroGrid ,  tau);
		
		//~ if(rank==root)
		//~ cout<<NP*sizeof(cell)*XCM*YCM*ZCMA/(1024*1024) <<" mega bytes"<<endl;
			
		double tmaxMev = 1000*MaxTempGev(HydroGrid) ;
		
		if(rank==root)
		{
			cout<<"Time Step is "<<ts <<" @ "; 
			cout<<std::fixed<<std::setprecision(4)<<" TempMax is "<< tmaxMev<<" MeV at TAU -->"<<tau;
			cout<<std::fixed<<std::setprecision(5)<<" This is tau step no. "<<l;
			clock_t now = clock();
			cout<<std::fixed<<std::setw(5)<<std::setprecision(5)<<" at time "<< (((double)(now-start))/CLOCKS_PER_SEC)<< " seconds"<<endl<<endl<<endl<<endl; 	
			fflush(stdout);
		}	

#ifdef LBI
		if( fabs( (tau-tauPrint) - printFreq) < 1e-6)
		{	tauPrint += printFreq;	WriteResultsXY(tau,HydroGrid);}//WriteTempXYCom(tau,HydroGrid) ;}//WriteResultsXYCom(tau,HydroGrid);}//	WriteTempXYCom(tau,HydroGrid) ;}//
		if( tau<(0.25+1E-6)   &&  fabs(100*tau - int(100*tau+1E-6))<1E-6 )
		{	WriteResultsXY(tau,HydroGrid);}
#else
		if( fabs( (tau-tauPrint) - printFreq) < 1e-6)// || fabs( (tau-(TAUSTART+TS))) < 1e-6 )
		{	tauPrint += printFreq;	 WriteResults(tau,HydroGrid);}
		
	
#endif
		l++;
		
		
#ifdef GUBSER
		if(tau>9)
			break;
#endif

#ifdef GINIT
		if(tmaxMev<120)
			break;
#endif

#if defined BJORKEN || BULKTEST
		if(tau>10)
			break;
#endif

		if(tau>25.3)
			break;
	}
	
#if defined BJORKEN || BULKTEST 
		myfile.close();
#endif		

	
	FreeFromHeap();
	if(!rank)
		cout<<endl<<"***********************************  Memory Freed from Heap    ***********************************"<<endl;
	
	MPI_Barrier(mpi_grid);
	
	if(!rank)
		cout<<endl<<"*************************************** The End ******************************************************"<<endl;
	MPI_Finalize();
 }
