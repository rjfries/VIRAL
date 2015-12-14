double ts;
double tau;

#define GEVFM 0.1973 

#if !defined(GINIT) && !defined(BJORKEN)&& !defined(GUBSER) 
//~ #define NSINIT
//~ #define ZEROINIT   
	#define NPX  4
	#define NPY  4
	#define NPZ  1
	#define NP (NPX*NPY*NPZ)

	#define TAUSTART 0.60 
	#define TS 0.01
	#define XL 6
	#define XS 0.1
	#define YL 6
	#define YS 0.1
	#define ZL 2
	#define ZS 0.1
	
	#define PFREQ 0.1
	#define FREQ ((int)2)
	#define FREQZ ((int)1)
#endif



#if defined(BJORKEN) 	
//~ #define NSINIT
//~ #define ZEROINIT

	#define NPX  1
	#define NPY  1
	#define NPZ  1
	#define NP (NPX*NPY*NPZ)

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
#endif



#if defined(GINIT) 
	#define LBI
	#define BULK
	#define S95P
	#define VORT
	#define FIX
	
	
	//~ #define NSINIT
	//~ #define ZEROINIT
	
	#define PEDESTAL 0 	
	
	#define NOS  6 //per side	for ginit.h
	#define NPX  10
	#define NPY  10
	#define NPZ  1
	#define NP (NPX*NPY*NPZ)

	#define TAUSTART 0.1
	#define TS 0.002
	#define XL 4
	#define XS 0.1
	#define YL 4
	#define YS 0.1
	#define ZL 0
	#define ZS 0.1
	
	#define PFREQ 0.1
	#define FREQ ((int)2)
	#define FREQZ ((int)1)
#endif



#if defined(GUBSER) 
//~ #define NSINIT
//~ #define ZEROINIT
	#define LBI
	
	
	#define NPX  4
	#define NPY  4
	#define NPZ  1
	#define NP (NPX*NPY*NPZ)

	#define TAUSTART 1.00 
	#define TS 0.005
	#define XL 5
	#define XS 0.05
	#define YL 5
	#define YS 0.05
	#define ZL 0
	#define ZS 0.05
	
	#define PFREQ 0.1		
	#define FREQ ((int)2)
	#define FREQZ ((int)1)
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

#define OFF 0.5  //(0.5 for correct typecasting)



#define XCM ((int)((XL/XS)+2*BORDER+OFF))
#define YCM ((int)((YL/YS)+2*BORDER+OFF))

#define XCMA ((int)(XL/XS+OFF))
#define YCMA ((int)(YL/YS+OFF))


#if !defined LBI
#define ZCM ((int)(2*(ZL/ZS)+2*BORDER+OFF+1)) 
#define ZCMA ((int)(2*(ZL/ZS)+OFF+1)) 
#endif

#if defined LBI
#define ZCM  1
#define ZCMA 1
#endif

#define VARN 4
#define Npi 10
#define NPI 1

#define PACKVAR     (4*VARN+Npi+NPI+1)//28
#define BACKUPVAR   (4*VARN+Npi+NPI+1) //28
#define SVAR (VARN+Npi+NPI)  //15
#define TEMPVAR 15


#if defined CON
	#define EVAR 2
#else
	//~ #define EVAR SVAR
	#define EVAR 1
#endif

typedef struct
{
	bool RELEVANT;
	double X,Y,eta,r;	//TODO Later

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

//#define WOODSAXON(r , width, loc)   (1.0)



inline double ratio(double num, double denom)
{
	if(fabs(num)<1e-14)
	{
		num=0;
		denom=1;
	}
	else if(num > 1e-14 &&	fabs(denom)<1e-14)
	{
		num =1e14;
		denom=1;
	}
	else if(num < -1e-14 &&  fabs(denom)<1e-14)
	{
		num = -1e14;
		denom=1;
	}

	return(num/denom);

}


void CalcNS(GRID HydroGrid, double tau, double ts);

inline double fmtoMev(double temp)
{
	double ret;
	
	ret = temp*GEVFM; //converts from 1/fm to Gev
	
	return 1000*ret;	
}

#define DECLTmu0      double T00 = HydroGrid[i][j][k].T00;  double T10 = HydroGrid[i][j][k].T10;double T20 = HydroGrid[i][j][k].T20;double T30 = HydroGrid[i][j][k].T30
#define DECLePPIa     double e = HydroGrid[i][j][k].En;  double P = HydroGrid[i][j][k].P;double PI = HydroGrid[i][j][k].PI;double a = DPDE(HydroGrid[i][j][k].P,HydroGrid[i][j][k].r)
#define DECLp10u4     double p1 = HydroGrid[i][j][k].pi[0];double p2 = HydroGrid[i][j][k].pi[1];double p3 = HydroGrid[i][j][k].pi[2];double p4 = HydroGrid[i][j][k].pi[3];double p5 = HydroGrid[i][j][k].pi[4];double p6 = HydroGrid[i][j][k].pi[5];double p7 = HydroGrid[i][j][k].pi[6];double p8 = HydroGrid[i][j][k].pi[7];double p9 = HydroGrid[i][j][k].pi[8];double p10 = HydroGrid[i][j][k].pi[9];double u0 = HydroGrid[i][j][k].u[0];double u1 = HydroGrid[i][j][k].u[1];double u2 = HydroGrid[i][j][k].u[2];double u3 = HydroGrid[i][j][k].u[3]
#define DECLp10       double p1 = HydroGrid[i][j][k].pi[0];double p2 = HydroGrid[i][j][k].pi[1];double p3 = HydroGrid[i][j][k].pi[2];double p4 = HydroGrid[i][j][k].pi[3];double p5 = HydroGrid[i][j][k].pi[4];double p6 = HydroGrid[i][j][k].pi[5];double p7 = HydroGrid[i][j][k].pi[6];double p8 = HydroGrid[i][j][k].pi[7];double p9 = HydroGrid[i][j][k].pi[8];double p10 = HydroGrid[i][j][k].pi[9]
#define DECLu4        double u0 = HydroGrid[i][j][k].u[0]; double u1 = HydroGrid[i][j][k].u[1]; double u2 = HydroGrid[i][j][k].u[2];double u3 = HydroGrid[i][j][k].u[3]
#define DECLcoord     double X = HydroGrid[i][j][k].X;     double Y = HydroGrid[i][j][k].Y;     double eta = HydroGrid[i][j][k].eta;double r = HydroGrid[i][j][k].r
