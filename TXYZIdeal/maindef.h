

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



#if  !defined(FLUCT) 
#ifdef SHAS
	#define CON
#endif

	#define NPX  1
	#define NPY  1
	#define NP (NPX*NPY*1)

	#define TAUSTART 0 
	#define TS 0.004
	#define XL 0.4
	#define XS 0.1
	#define YL 0.4
	#define YS 0.1
	
	#define ZL 4
	#define ZS 0.02

	#define PFREQ 0.2
	#define FREQ ((int)1)
	#define FREQZ ((int)1)
	
	inline double EOS(double en )                                       {return (en/3.0                        );}	
	inline double DPDE(double en )                                      {return (1.0/3.0                       );}	 
	inline double FEnFromTemp(double temp)                              {return (FACTOR*pow(temp,4)            );}
	inline double FT(double en)                                         {return (pow(en/FACTOR,0.25)           );}
	inline double FS( double en, double Pr, double T)                   {return ((en+Pr)*pow(T,-1)             );}
#endif




#if defined(FLUCT)
	#define LBI
#ifdef SHAS
	#define CON
#endif
	#define NPX  4
	#define NPY  4 
	#define NP (NPX*NPY*1)

	#define TAUSTART 0.0 
	#define TS 0.01
	#define XL 4
	#define XS 0.1
	#define YL 4
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
#endif

#define DIM 3
#define WENOP 2
#define WENOEPS 1E-6

#ifdef SHAS
#define BORDER 3
#endif

#ifdef KT
#define GMINV 1.1
#define BORDER 3
#endif



#define VARN 4 

#define PACKVAR    (4*VARN + 1) 
#define SVAR (VARN )  
 
#define EVAR 1
typedef struct
{
	bool RELEVANT;
	double X,Y,eta,r; 
	double T00,T10,T20,T30;
	double En,P,Vx,Vy,Ve,Temp;
	double u[VARN];
	double prevu[VARN], du[VARN][VARN];	//4 velocity and its derivative
 
//for source terms	
	double Source [SVAR];
 

//Common For all hydro schemes
	double Var[SVAR],  Result[SVAR],  PartialResult[SVAR];
	
//For time integration
	double Var0[SVAR], Var1[SVAR],  Var2[SVAR],  Var3[SVAR];
	double L0[SVAR],  L1[SVAR],  L2[SVAR];
	
	
#if defined SHAS  
	double UTD,Ubar,A,Ac;
	double NVx[SVAR],NVy[SVAR],NVz[SVAR];
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


inline double fmtoMev(double temp)
{
	double ret;
	ret = temp*GEVFM; //converts from 1/fm to Gev	
	return 1000*ret;	
}


#define DECLTmu0      double T00 = HydroGrid[i][j][k].T00;double T10 = HydroGrid[i][j][k].T10;double T20 = HydroGrid[i][j][k].T20;double T30 = HydroGrid[i][j][k].T30
#define DECLePa     double e = HydroGrid[i][j][k].En;double P = HydroGrid[i][j][k].P; double a = DPDE(HydroGrid[i][j][k].P )
#define DECLu4        double u0 = HydroGrid[i][j][k].u[0]; double u1 = HydroGrid[i][j][k].u[1]; double u2 = HydroGrid[i][j][k].u[2];double u3 = HydroGrid[i][j][k].u[3]
#define DECLcoord     double X = HydroGrid[i][j][k].X;     double Y = HydroGrid[i][j][k].Y;     double eta = HydroGrid[i][j][k].eta;double r = HydroGrid[i][j][k].r
