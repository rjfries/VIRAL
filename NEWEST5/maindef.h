

/*    x massless quark degrees of freedom
 * 	  gq = 2·2·3 = 12 for quarks
 *    gg = 16 for gluons
 * 	  
 * 	DOF = x*gq*(7/8) + (2*8)
 */
#define DOF 42.25  //x=2.5
#define FACTOR (3.0*DOF*PIE*PIE/90.0) 


#define OFF 0.5  //(0.5 for correct typecasting)
#define XCM ((int)((XL/XS)+2*BORDER+OFF))
#define YCM ((int)((YL/YS)+2*BORDER+OFF))

#define XCMA ((int)(XL/XS+OFF))
#define YCMA ((int)(YL/YS+OFF))
  


#if defined LBI 
	#define ZCM 1
	#define ZCMA 1
#else 
	#define ZCM ((int)(2*(ZL/ZS)+OFF+1+2*BORDER)) 
	#define ZCMA ((int)(2*(ZL/ZS)+OFF+1)) 
#endif

double ts;
double tau;
 



#if !defined(GINIT) && !defined(BJORKEN)&& !defined(GUBSER) && !defined(BULKTEST)  && !defined(FLUCT) 
	//~ #define NSINIT
	//~ #define ZEROINIT   
	//~ #define LBI

#ifdef SHAS
	#define CON
#endif

	#define NPX  1
	#define NPY  1 
	#define NP (NPX*NPY*1)

	#define TAUSTART 0.60 
	#define TS 0.002
	#define XL 4
	#define XS 0.1
	#define YL 4
	#define YS 0.1
	
#ifdef LBI
	#define ZL 0
	#define ZS 0.1
#else
	#define ZL 5
	#define ZS 0.1
#endif	
	
	#define PFREQ 0.01
	#define FREQ ((int)1)
	#define FREQZ ((int)1)
	
	inline double EOS(double en )                                       {return (en/3.0                        );}	
	inline double DPDE(double en )                                      {return (1.0/3.0                       );}	 
	inline double FEnFromTemp(double temp)                              {return (FACTOR*pow(temp,4)            );}
	inline double FT(double en)                                         {return (pow(en/FACTOR,0.25)           );}
	inline double FS( double en, double Pr, double T)                   {return ((en+Pr)*pow(T,-1)             );}

	#define SCALE_VIS 1
	#define SCALE_TPI 1	
	inline double Feta( double s, double en)                            {return (SCALE_VIS*(s/(4.0*PIE))       );}
	inline double Ftaupi( double eta , double  p, double en)			{return (1.5*SCALE_TPI*(eta/p)         );}

#ifdef BULK
	#define SCALE_BULK_VIS 1
	#define SCALE_BULK_TPI 1	
	inline double FZeta( double s, double en)                           {return (SCALE_BULK_VIS*(s/(4.0*PIE))  );}
	inline double FtauPI(double zeta , double  p, double en)            {return (1.5*SCALE_TPI_BULK*zeta/p     );}
#endif

	
#endif







#if defined FLUCT 
	#define LBI

#ifdef SHAS
	#define CON
#endif

	#define NPX  4
	#define NPY  4 
	#define NP (NPX*NPY*1)

	#define TAUSTART 0.60 
	#define TS 0.002
	#define XL 5
	#define XS 0.1
	#define YL 5
	#define YS 0.1
	 
#ifdef LBI
	#define ZL 0
	#define ZS 0.1
#else
	#define ZL 5
	#define ZS 0.1
#endif	
	
	#define PFREQ 0.01
	#define FREQ ((int)1)
	#define FREQZ ((int)1)
	
	inline double EOS(double en )                                       {return (en/3.0                        );}	
	inline double DPDE(double en )                                      {return (1.0/3.0                       );}	 
	inline double FEnFromTemp(double temp)                              {return (FACTOR*pow(temp,4)            );}
	inline double FT(double en)                                         {return (pow(en/FACTOR,0.25)           );}
	inline double FS( double en, double Pr, double T)                   {return ((en+Pr)*pow(T,-1)             );}

	#define SCALE_VIS 0.01
	#define SCALE_TPI 1	
	inline double Feta( double s, double en)                            {return (SCALE_VIS*(s/(4.0*PIE))       );}
	inline double Ftaupi( double eta , double  p, double en)			{return (1.5*SCALE_TPI*(eta/p)         );}

#ifdef BULK
	#define SCALE_BULK_VIS 1
	#define SCALE_BULK_TPI 1	
	inline double FZeta( double s, double en)                           {return (SCALE_BULK_VIS*(s/(4.0*PIE))  );}
	inline double FtauPI(double zeta , double  p, double en)            {return (1.5*SCALE_TPI_BULK*zeta/p     );}
#endif

	
#endif











#if defined BULKTEST
	#define BULK  
#ifdef SHAS
	#define CON
#endif
	
	#define NPX  1
	#define NPY  1
	#define NP (NPX*NPY*1)

	#define TAUSTART 0.60 
	#define TS 0.01
	#define XL 1
	#define XS 0.1
	#define YL 1
	#define YS 0.1 
	#define ZL 2
	#define ZS 0.1
	 

	#define PFREQ 0.1	
	#define FREQ ((int)2)
	#define FREQZ ((int)1) 
	
	inline double EOS(double en )                                       {return (en/3.0                        );}	
	inline double DPDE(double en )                                      {return (1.0/3.0                       );}	 
	inline double FEnFromTemp(double temp)                              {return (FACTOR*pow(temp,4)            );}
	inline double FT(double en)                                         {return (pow(en/FACTOR,0.25)           );}
	inline double FS( double en, double Pr, double T)                   {return ((en+Pr)*pow(T,-1)             );}

	#define SCALE_VIS 1
	#define SCALE_TPI 1	
	inline double Feta( double s, double en)                            {return (0.2*s                         );}
	inline double Ftaupi( double eta , double  p, double en)			{return (0.01                          );}

#ifdef BULK
	#define SCALE_BULK_VIS 1
	#define SCALE_BULK_TPI 1	
	inline double FZeta( double s, double en)                           {return (1                             );}
	inline double FtauPI(double zeta , double  p, double en)            {return (0.1						   );}
#endif


#endif



#if defined(BJORKEN)

#ifdef SHAS
	#define CON
#endif

	#define NPX  1
	#define NPY  1 
	#define NP (NPX*NPY*1)

	#define TAUSTART 0.60 
	#define TS 0.01
	#define XL 1
	#define XS 0.1
	#define YL 1
	#define YS 0.1
	#define ZL 2
	#define ZS 0.1
	  
	#define PFREQ 0.2
	#define FREQ ((int)1)
	#define FREQZ ((int)1)	
	
	inline double EOS(double en )                                       {return (en/3.0                        );}	
	inline double DPDE(double en )                                      {return (1.0/3.0                       );}	 
	inline double FEnFromTemp(double temp)                              {return (FACTOR*pow(temp,4)            );}
	inline double FT(double en)                                         {return (pow(en/FACTOR,0.25)           );}
	inline double FS( double en, double Pr, double T)                   {return ((en+Pr)*pow(T,-1)             );}

	#define SCALE_VIS 1
	#define SCALE_TPI 1	
	inline double Feta( double s, double en)                            {return (0.2*s                         );}
	inline double Ftaupi( double eta , double  p, double en)			{return (0.01                          );}

#ifdef BULK
	#define SCALE_BULK_VIS 1
	#define SCALE_BULK_TPI 1	
	inline double FZeta( double s, double en)                           {return (SCALE_BULK_VIS*(s/(4.0*PIE))  );}
	inline double FtauPI(double zeta , double  p, double en)            {return (1.5*SCALE_TPI_BULK*zeta/p     );}
#endif

#endif



#if defined(GINIT)  
	//~ #define BULK
	#define S95P
	//~ #define VORT 

#ifdef SHAS
		#define CON
#endif

	#define PEDESTAL 0 	
	
	#define NOS  6 //per side	for ginit.h
	#define NPX  1
	#define NPY  1
	#define NP (NPX*NPY*1)

	#define TAUSTART 0.1
	#define TS 0.002
	#define XL 4
	#define XS 0.1
	#define YL 4
	#define YS 0.1
	
#ifdef LBI
	#define ZL 0
	#define ZS 0.1
#else
	#define ZL 5
	#define ZS 0.1
#endif	

	#define PFREQ 0.01
	#define FREQ ((int)1)
	#define FREQZ ((int)1)
	
	inline double EOS(double en )                                       {return (s95p_p(en)                    );}	
	inline double DPDE(double en )                                      {return (s95p_a(en)                    );}	 
	inline double FEnFromTemp(double temp)                              {return (FACTOR*pow(temp,4)            );}
	inline double FT(double en)                                         {return (s95p_T(en)                    );}
	inline double FS( double en, double Pr, double T)                   {return (s95p_s(en)                    );}

	#define SCALE_VIS 1
	#define SCALE_TPI 1	
	inline double Feta( double s, double en)                            {return (SCALE_VIS*(s/(4.0*PIE))       );}
	inline double Ftaupi( double eta , double  p, double en)			{return (1.5*SCALE_TPI*(eta/p)         );}


#ifdef BULK
	#define SCALE_BULK_VIS 0.01
	#define SCALE_BULK_TPI 1	
	inline double FZeta( double s, double en)                           {return (SCALE_BULK_VIS*(s/(4.0*PIE))  );}
	inline double FtauPI(double zeta , double  p, double en)            {return (1.5*SCALE_BULK_TPI*zeta/p     );}
#endif


#endif



#if defined(GUBSER) 
	#define LBI
	
	
	#define NPX  6
	#define NPY  6
	#define NP (NPX*NPY*1)

	#define TAUSTART 1.00 
	#define TS 0.005
	#define XL 5
	#define XS 0.05
	#define YL 5
	#define YS 0.05 
	#define ZL 0
	#define ZS 0.1 


	#define PFREQ 0.1		
	#define FREQ ((int)2)
	#define FREQZ ((int)1) 
	
	inline double EOS(double en )                                       {return (en/3.0                        );}	
	inline double DPDE(double en )                                      {return (1.0/3.0                       );}	 
	inline double FEnFromTemp(double temp)                              {return (FACTOR*pow(temp,4)            );}
	inline double FT(double en)                                         {return (pow(en/FACTOR,0.25)           );}
	inline double FS( double en, double Pr, double T)                   {return ((en+Pr)*pow(T,-1)             );}

	#define SCALE_VIS 1
	#define SCALE_TPI 1	
	inline double Feta( double s, double en)                            {return (0.2*s                         );}
	inline double Ftaupi( double eta , double  p, double en)			{return (5*eta/(en+p) 			       );}


#ifdef BULK
	#define SCALE_BULK_VIS 0.01
	#define SCALE_BULK_TPI 1	
	inline double FZeta( double s, double en)                           {return (SCALE_BULK_VIS*(s/(4.0*PIE))  );}
	inline double FtauPI(double zeta , double  p, double en)            {return (1.5*SCALE_TPI_BULK*zeta/p     );}
#endif

#endif 
	
	
	









	
#define DIM 3
#define WENOP 2
#define WENOEPS 1E-6




#ifdef SHAS
#define BORDER 3
#endif

#ifdef ZAL
#define BORDER 4
#endif

#ifdef KT
#define GMINV 1.1
#define BORDER 3
#endif



#define VARN 4
#define Npi 5
#define NPI 1

#define PACKVAR    (4*VARN+Npi+NPI+1)//23
#define BACKUPVAR (4*VARN+Npi+NPI+1) //23
#define SVAR (VARN+Npi+NPI)  //15
#define TEMPVAR 15


#if defined CON
	#define EVAR 2
#else
	#define EVAR 1
#endif

typedef struct
{
	bool RELEVANT;
	double X,Y,eta,r;

//for debugging
	double temp[TEMPVAR];
	
	double T00,T10,T20,T30;
	double En,P,Vx,Vy,Ve,Temp;
	double u[VARN];
	double prevu[VARN], du[VARN][VARN];	//4 velocity and its derivative

//for viscosity	
	double pi[Npi], PI; 
	double nspi[Npi], nsPI;	
	
//for source terms	
	double Source [SVAR];

//for backup variable	
	double BackUp[BACKUPVAR];


//Common For all hydro schemes
	double Var[SVAR],  Result[SVAR],  PartialResult[SVAR];
	
//For time integration
	double Var0[SVAR], Var1[SVAR],  Var2[SVAR],  Var3[SVAR];
	double L0[SVAR],  L1[SVAR],  L2[SVAR];
	
	
#if defined SHAS || defined ZAL
	double UTD,Ubar,A,Ac;
	double NVx[SVAR],NVy[SVAR],NVz[SVAR];
#endif
	

#ifdef VORT
	double Vort[Npi];
#endif

#if defined KT
	double Fx[SVAR], Fy[SVAR], Fz[SVAR];  //centered fluxes
	
	double FxLX[SVAR], FxRX[SVAR];
	double FyLY[SVAR], FyRY[SVAR];
	double FzLZ[SVAR], FzRZ[SVAR];  //centered fluxes
	
	double VarLX[SVAR], VarRX[SVAR];
	double VarLY[SVAR], VarRY[SVAR];
	double VarLZ[SVAR], VarRZ[SVAR];
	
 
	double Ax[EVAR],Ay[EVAR],Az[EVAR];//all or 2 of the center eigen values
	
	double AxLX[EVAR], AxRX[EVAR]; //reconstruct all of the center eigen values
	double AyLY[EVAR], AyRY[EVAR]; //reconstruct all of the center eigen values
	double AzLZ[EVAR], AzRZ[EVAR]; //reconstruct all of the center eigen values
 
  
	double AxLXMAX, AxRXMAX;
	double AyLYMAX, AyRYMAX;
	double AzLZMAX, AzRZMAX;
	
	
	double fluxT[SVAR];	
#endif

}cell;

double BMax;

typedef cell (*GRID) [YCM][ZCM];
GRID HydroGrid;


typedef double (*PCKX) [BORDER][YCMA][ZCMA];
typedef double (*PCKY) [BORDER][XCMA][ZCMA];
typedef double (*PCKZ) [BORDER][XCMA][YCMA];
typedef double (*PCKXY) [BORDER][BORDER][ZCMA];

 

#define WOODSAXON(r , width, loc)   (1.0/(1 + exp( (r - loc)/width) ) )
#define HeaviSideTheta(num)   ( (num>=0)?1:0)


void CalcNS(GRID HydroGrid, double tau, double ts);


inline double fmtoMev(double temp)
{
	double ret;
	ret = temp*GEVFM; //converts from 1/fm to Gev	
	return 1000*ret;	
}


#define DECLTmu0      double T00 = HydroGrid[i][j][k].T00;double T10 = HydroGrid[i][j][k].T10;double T20 = HydroGrid[i][j][k].T20;double T30 = HydroGrid[i][j][k].T30
#define DECLePPIa     double e = HydroGrid[i][j][k].En;double P = HydroGrid[i][j][k].P;double PI = HydroGrid[i][j][k].PI;double a = DPDE(HydroGrid[i][j][k].P )
#define DECLp5u4      double p1 = HydroGrid[i][j][k].pi[0];double p2 = HydroGrid[i][j][k].pi[1];double p3 = HydroGrid[i][j][k].pi[2];double p4 = HydroGrid[i][j][k].pi[3];double p5 = HydroGrid[i][j][k].pi[4];double u0 = HydroGrid[i][j][k].u[0];double u1 = HydroGrid[i][j][k].u[1];double u2 = HydroGrid[i][j][k].u[2];double u3 = HydroGrid[i][j][k].u[3]
#define DECLp5        double p1 = HydroGrid[i][j][k].pi[0];double p2 = HydroGrid[i][j][k].pi[1];double p3 = HydroGrid[i][j][k].pi[2];double p4 = HydroGrid[i][j][k].pi[3];double p5 = HydroGrid[i][j][k].pi[4]
#define DECLu4        double u0 = HydroGrid[i][j][k].u[0]; double u1 = HydroGrid[i][j][k].u[1]; double u2 = HydroGrid[i][j][k].u[2];double u3 = HydroGrid[i][j][k].u[3]
#define DECLcoord     double X = HydroGrid[i][j][k].X;     double Y = HydroGrid[i][j][k].Y;     double eta = HydroGrid[i][j][k].eta;double r = HydroGrid[i][j][k].r
#define A1 (-((2*p3*u1*u2 + p1*pow(u1,2) + p2*pow(u2,2) + pow(tau,2)*(2*(p4*u1 + p5*u2)*u3 - (p1 + p2)*pow(u3,2)))*pow(-pow(u0,2) + pow(tau,2)*pow(u3,2),-1)))      
#define A2 ((p1*u1 + p3*u2 + p4*u3*pow(tau,2))*pow(u0,-1))             
#define A3 ((p3*u1 + p2*u2 + p5*u3*pow(tau,2))*pow(u0,-1))             
#define A4 (pow(u0,-1)*(2*p3*u1*u2*u3 + (p4*u1 + p5*u2 - (p1 + p2)*u3)*pow(u0,2) + p1*u3*pow(u1,2) + p2*u3*pow(u2,2) + p4*u1*pow(tau,2)*pow(u3,2) + p5*u2*pow(tau,2)*pow(u3,2))*pow(pow(u0,2) - pow(tau,2)*pow(u3,2),-1))     
#define A5 (-(pow(tau,-2)*(2*p3*u1*u2 + 2*(p4*u1 + p5*u2)*u3*pow(tau,2) - (p1 + p2)*pow(u0,2) + p1*pow(u1,2) + p2*pow(u2,2))*pow(-pow(u0,2) + pow(tau,2)*pow(u3,2),-1)))     
