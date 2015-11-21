double ts;
double tau;
#define GEVFM 0.1973 

#if  !defined(BJORKEN)&& !defined(GUBSER)&& !defined(STAT)
	#define NPX  4
	#define NPY  4

	#define NPZ  1
	#define NP (NPX*NPY*NPZ)
	#define TAUSTART 0.10 
	#define XL 6
	#define XS 0.1

	#define YL 6
	#define YS 0.1

	#define ZL 0.1
	#define ZS 0.1
#endif

#ifdef STAT		

	#define NPX  2
	#define NPY  2
	#define NPZ  1

	#define NP (NPX*NPY*NPZ)
	#define TAUSTART 0.60 
	#define XL 1
	#define XS 0.1

	#define YL 1
	#define YS 0.1

	#define ZL 3
	#define ZS 0.1
#endif


#ifdef BJORKEN		

	#define NPX  1
	#define NPY  1
	#define NPZ  1

	#define NP (NPX*NPY*NPZ)
	#define TAUSTART 0.60 
	#define XL 1
	#define XS 0.1

	#define YL 1
	#define YS 0.1

	#define ZL 2
	#define ZS 0.1
#endif




#ifdef GUBSER

	#define NPX  4
	#define NPY  4

	#define NPZ  1
	#define NP (NPX*NPY*NPZ)

	#define TAUSTART 1.00 
	#define XL 6
	#define XS 0.1

	#define YL 6
	#define YS 0.1

	#define ZL 0.1
	#define ZS 0.1
#endif 
	
	
#define DIM 3
#define WENOP 2
#define WENOEPS 1E-6


#define FREQ ((int)1)
#define FREQZ ((int)1)



#define TS 0.01
#define ETASTART 0

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
#define ZCM ((int)(2*(ZL/ZS)+2*BORDER+OFF+1)) 





#define XCMA ((int)(XL/XS+OFF))
#define YCMA ((int)(YL/YS+OFF))
#define ZCMA ((int)(2*(ZL/ZS)+OFF+1)) 


#define VARN 4

#define PACKVAR (2*VARN+1)//9
#define BACKUPVAR (2*VARN+1)//9
typedef struct
{
	bool RELEVANT;
	
	double X,Y,eta,r;	
	double En,Vx,Vy,Ve;
	double NVx[4], NVy[4], NVe[4];
	double Fx[4], Fy[4], Fe[4];
	double P;
	double T00,T10,T20,T30;
	double Source [VARN];
	double BufA[VARN],Result[VARN],PartialResult[VARN];
	double Buf,Temp,Fun;
	double BackUp[BACKUPVAR];
	
#if defined SHAS || defined ZAL
	double UTD,Ubar,A,Ac;
#endif

#if defined KT
	double  fluxL, fluxR, fluxT;
	double velL, velR;
	double varL, varR;
	double u0L, u0R;
	double left, right,VAR;
#endif

}cell;

typedef cell (*GRID) [YCM][ZCM];
GRID HydroGrid;


typedef double (*ARRXY) [YCM];
typedef double (*PCKX) [BORDER][YCMA][ZCMA];
typedef double (*PCKY) [BORDER][XCMA][ZCMA];
typedef double (*PCKZ) [BORDER][XCMA][YCMA];


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
