#define PIE 3.141592653589793


/*  
 *    x massless quark degrees of freedom
 * 	  gq = 2·2·3 = 12 for quarks
 *    gg = 16 for gluons
 * 	  
 * 	DOF = x*gq*(7/8) + (2*8)
 */

#define DOF 42.25  //x=2.5


#define FACTOR (3.0*DOF*PIE*PIE/90.0)

#ifdef IDEAL
#define SCALE_VIS 1e-6
#define SCALE_TPI 1e-6
#define SCALE_VIS_BULK 1e-6
#define SCALE_TPI_BULK 1e-6
#endif


inline double EOS(double en , double r)
{
	double ret;
	
	ret = (en/3.0);//*WOODSAXON(r, 0.5, 18); 
	
#ifdef S95P
	ret = s95p_p(en);
#endif

	return(ret);
}

inline double DPDE(double en , double r)
{
	double ret;
	
	ret = (1.0/3.0);//*WOODSAXON(r, 0.5, 18);

#ifdef S95P
	ret = s95p_a(en);
#endif
		
	return(ret);
}

inline double FEnFromTemp(double temp)
{
	double ret;
	
	ret = FACTOR*temp*temp*temp*temp;

 
	return ret;
}

inline double FT(double en , double r)
{
	double ret;
	
	ret = pow(en/FACTOR,0.25);//*WOODSAXON(r, 0.5, 18);

#ifdef S95P
	ret = s95p_T(en);
#endif		
		
	return(ret);
}



inline double FS( double en, double Pr, double T)
{
	double ret;
	
	ret = ((en+Pr)*pow(T,-1));// * WOODSAXON(r, width, loc);
	
#ifdef S95P
	ret = s95p_s(en);
#endif		
		
	return(ret);
}




//SHEAR VISCOSITY


//shear viscosity terms
#ifndef IDEAL
#define SCALE_VIS 1
#define SCALE_TPI 1
#define SCALE_VIS_BULK 1
#define SCALE_TPI_BULK 1
#endif

inline double Feta(  double s, double en, double r)
{
	double ret;	
	double loc = 6;
	double width = 0.2;
	double scale;
	
	scale = SCALE_VIS;
	
	ret = (scale*(s/(4.0*PIE)));//*WOODSAXON(r, width, loc);
	
	
	#ifdef BJORKEN
		ret  = 0.2*s;
	#endif
	
	#ifdef GUBSER
		ret  = 0.2*s;
	#endif
	
	return(ret);
	
}

inline double Ftaupi( double eta , double  p, double en, double r)
{
	double ret;	
	double loc = 16;
	double width = 1;
	
	double scale = SCALE_TPI;
			
	ret = (1.5*scale*eta/p);// * WOODSAXON(r, width, loc);
	
	#ifdef BJORKEN
	ret  = 0.01;
	#endif
	
	
#ifdef GUBSER
	double c=5;
	ret  = c*eta/(en+p);
#endif
	
	
	return(ret);
}





// BULK VISCOSTY


inline double FZeta( double s, double en, double r)
{
	double ret;
	
	ret = (SCALE_VIS_BULK*(s/(4.0*PIE)) );
		
	return(ret);
}

inline double FtauPI( double zeta , double  p, double en, double r)
{
	double ret;
	
	ret = (1.5*SCALE_TPI_BULK* zeta/p );
		
	return(ret);
}
