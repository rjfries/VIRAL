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


#ifdef SHAS
	#define CON
#endif

using namespace std;

#include "maindef.h"

#ifdef S95P
#include "s95p.h"
#endif
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
		myfile<<tau<< "  "<< fmtoMev(FT(HydroGrid[XCM/2][YCM/2][ZCM/2].En , HydroGrid[XCM/2][YCM/2][ZCM/2].r))/1000<<endl;
#endif
		 
		
		tau += ts;

		CheckRoot( HydroGrid ,  tau);
		
		//~ if(rank==root)
			//~ cout<<NP*sizeof(cell)*XCM*YCM*ZCM/(1024*1024) <<" mega bytes"<<endl;
			
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

		if( fabs( (tau-tauPrint) - printFreq) < 1e-6)
		{	tauPrint += printFreq;	WriteResultsXY(tau,HydroGrid);}//WriteTempXYCom(tau,HydroGrid) ;}//WriteResultsXYCom(tau,HydroGrid);}//	WriteTempXYCom(tau,HydroGrid) ;}//
		
		if( fabs( (tau-sourcePrint) - sourceFreq) < 1e-6 )
		{	sourcePrint += sourceFreq;				}//WriteSourceXY(tau-ts,HydroGrid);}

		if( fabs( (tau-ts-derPrint) - derFreq) < 1e-6 )
			derPrint += derFreq;	
						
		l++;
		
		
#ifdef GUBSER
		if(tau>9)
			break;
#endif

#ifdef GINIT
		if(tmaxMev<120)
			break;
#endif

		if(tau>25.3)
			break;
	}
	

	if(!rank)
		cout<<endl<<"*************************************** The End ******************************************************"<<endl;
		
	MPI_Finalize();
 }
