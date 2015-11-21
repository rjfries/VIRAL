void hydroExplicit(GRID HydroGrid, double tau, double taustep);
void fvX(GRID HydroGrid, double tau);
void fvY(GRID HydroGrid, double tau);
void fvZ(GRID HydroGrid, double tau);


#define WENOP 2
#define WENOEPS 1E-6


int method;
int type;



void CalcCentreFlux(GRID HydroGrid, double tau)
{
	int i,j,k;
	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		DECLp10u4;
		DECLePPIa;

#if defined CON			
		HydroGrid[i][j][k].Fx[0] =  tau*(       p2 + u1*u0*(e+P-PI)       );
		HydroGrid[i][j][k].Fx[1] =  tau*( HydroGrid[i][j][k].T10*u1/u0    );
		HydroGrid[i][j][k].Fx[2] =  tau*( HydroGrid[i][j][k].T20*u1/u0    );
		HydroGrid[i][j][k].Fx[3] =  tau*( HydroGrid[i][j][k].T30*u1/u0    );
#else
   		HydroGrid[i][j][k].Fx[0] =  tau*(  p2 + u1*u0*(e+P-PI)            );
		HydroGrid[i][j][k].Fx[1] =  tau*(  p5 + (P-PI) + u1*u1*(e+P-PI)   );
		HydroGrid[i][j][k].Fx[2] =  tau*(  p6 + u1*u2*(e+P-PI)            );
		HydroGrid[i][j][k].Fx[3] =  tau*(  p7 + u1*u3*(e+P-PI)            );
#endif

        HydroGrid[i][j][k].Fx[4]=  (p1*u1)/u0;
        HydroGrid[i][j][k].Fx[5]=  (p2*u1)/u0;
        HydroGrid[i][j][k].Fx[6]=  (p3*u1)/u0;
        HydroGrid[i][j][k].Fx[7]=  (p4*u1)/u0;
        HydroGrid[i][j][k].Fx[8]=  (p5*u1)/u0;
        HydroGrid[i][j][k].Fx[9]=  (p6*u1)/u0;
        HydroGrid[i][j][k].Fx[10]=  (p7*u1)/u0;
        HydroGrid[i][j][k].Fx[11]=  (p8*u1)/u0;
        HydroGrid[i][j][k].Fx[12]=  (p9*u1)/u0;
        HydroGrid[i][j][k].Fx[13]=  (p10*u1)/u0;
        HydroGrid[i][j][k].Fx[14]=  (PI*u1)/u0;

#if  defined CON
		HydroGrid[i][j][k].Fy[0] =  tau*(    p3 + u2*u0*(e+P-PI)            );
		HydroGrid[i][j][k].Fy[1] =  tau*(  HydroGrid[i][j][k].T10*u2/u0      );
		HydroGrid[i][j][k].Fy[2] =  tau*(  HydroGrid[i][j][k].T20*u2/u0      );
		HydroGrid[i][j][k].Fy[3] =  tau*(  HydroGrid[i][j][k].T30*u2/u0      );
#else

		HydroGrid[i][j][k].Fy[0] =  tau*(  p3 + u2*u0*(e+P-PI)            );
		HydroGrid[i][j][k].Fy[1] =  tau*(  p6 + u2*u1*(e+P-PI)            );
		HydroGrid[i][j][k].Fy[2] =  tau*(  p8 + (P-PI) + u2*u2*(e+P-PI)   );
		HydroGrid[i][j][k].Fy[3] =  tau*(  p9 + u2*u3*(e+P-PI)            );
#endif		
		HydroGrid[i][j][k].Fy[4] =  (p1*u2)/u0;
        HydroGrid[i][j][k].Fy[5] =  (p2*u2)/u0;
        HydroGrid[i][j][k].Fy[6] =  (p3*u2)/u0;
        HydroGrid[i][j][k].Fy[7] =  (p4*u2)/u0;
        HydroGrid[i][j][k].Fy[8] =  (p5*u2)/u0;
        HydroGrid[i][j][k].Fy[9] =  (p6*u2)/u0;
        HydroGrid[i][j][k].Fy[10]=  (p7*u2)/u0;
        HydroGrid[i][j][k].Fy[11]=  (p8*u2)/u0;
        HydroGrid[i][j][k].Fy[12]=  (p9*u2)/u0;
        HydroGrid[i][j][k].Fy[13]=  (p10*u2)/u0;
        HydroGrid[i][j][k].Fy[14]=  (PI*u2)/u0;
	

#if defined CON			
		HydroGrid[i][j][k].Fz[0] =  tau*(  p4 + u3*u0*(e+P-PI)                );
		HydroGrid[i][j][k].Fz[1] =  tau*(  HydroGrid[i][j][k].T10*u3/u0       );
		HydroGrid[i][j][k].Fz[2] =  tau*(  HydroGrid[i][j][k].T20*u3/u0       );
		HydroGrid[i][j][k].Fz[3] =  tau*(  HydroGrid[i][j][k].T30*u3/u0       );	
#else	
		HydroGrid[i][j][k].Fz[0] =  tau*(  p4 + u3*u0*(e+P-PI)                     );
		HydroGrid[i][j][k].Fz[1] =  tau*(  p7 + u3*u1*(e+P-PI)                     );
		HydroGrid[i][j][k].Fz[2] =  tau*(  p9 + u3*u2*(e+P-PI)                     );
		HydroGrid[i][j][k].Fz[3] =  tau*(  p10 + ((P-PI)/(tau*tau)) + u3*u3*(e+P-PI)     );	
#endif
	
		HydroGrid[i][j][k].Fz[4] =  (p1*u3)/u0;
        HydroGrid[i][j][k].Fz[5] =  (p2*u3)/u0;
        HydroGrid[i][j][k].Fz[6] =  (p3*u3)/u0;
        HydroGrid[i][j][k].Fz[7] =  (p4*u3)/u0;
        HydroGrid[i][j][k].Fz[8] =  (p5*u3)/u0;
        HydroGrid[i][j][k].Fz[9] =  (p6*u3)/u0;
        HydroGrid[i][j][k].Fz[10]=  (p7*u3)/u0;
        HydroGrid[i][j][k].Fz[11]=  (p8*u3)/u0;
        HydroGrid[i][j][k].Fz[12]=  (p9*u3)/u0;
        HydroGrid[i][j][k].Fz[13]=  (p10*u3)/u0;
        HydroGrid[i][j][k].Fz[14]=  (PI*u3)/u0;		
	}	
}


#define GMMX(VAR )\
void GMMX_##VAR(GRID HydroGrid , int NVAR  ){\
	int i,j,k,l;\
	for(l=0;l<NVAR;l++)\
	for(j=jl;j<jr;j++)\
	for(k=kl;k<kr;k++)\
	for(i=1;i<XCM-1;i++){\
		double qm1,qc,qp1;\
		qm1=HydroGrid[i-1][j][k].VAR[l];\
		qc=HydroGrid[i][j][k].VAR[l];\
		qp1=HydroGrid[i+1][j][k].VAR[l];\
		double slope1 = genminmod(qm1,qc,qp1,XS, GMINV);\
		double num = (qc  - qm1 );\
		double denom = (qp1 - qc );\
		double r = ratio(num,denom);\
		slope1 *= phi(r,-1);\
		HydroGrid[i][j][k].VAR##LX[l]   = qc + slope1*(XS/2);\
		HydroGrid[i-1][j][k].VAR##RX[l] = qc - slope1*(XS/2);\
	}\
}

#define GMMY(VAR )\
void GMMY_##VAR(GRID HydroGrid , int NVAR    ){\
	int i,j,k,l;\
	for(l=0;l<NVAR;l++)\
	for(i=il;i<ir;i++)\
	for(k=kl;k<kr;k++)\
	for(j=1;j<YCM-1;j++){\
		double qm1,qc,qp1,qp2;\
		qm1=HydroGrid[i][j-1][k].VAR[l];\
		qc=HydroGrid[i][j][k].VAR[l];\
		qp1=HydroGrid[i][j+1][k].VAR[l];\
		double slope1 = genminmod(qm1,qc,qp1,YS,GMINV);\
		double num = (qc  - qm1 );\
		double denom = (qp1 - qc );\
		double r = ratio(num,denom);\
		slope1 *= phi(r,-1);\
		HydroGrid[i][j][k].VAR##LY[l]   = qc + slope1*(YS/2);\
		HydroGrid[i][j-1][k].VAR##RY[l] = qc - slope1*(YS/2);\
	}\
}


#define GMMZ(VAR )\
void GMMZ_##VAR(GRID HydroGrid , int NVAR  ){\
	int i,j,k,l;\
	for(l=0;l<NVAR;l++)\
	for(i=il;i<ir;i++)\
	for(j=jl;j<jr;j++)\
	for(k=1;k<ZCM-1;k++){\
		double qm1,qc,qp1,qp2;\
		qm1=HydroGrid[i][j][k-1].VAR[l];\
		qc=HydroGrid[i][j][k].VAR[l];\
		qp1=HydroGrid[i][j][k+1].VAR[l];\
		double slope1 = genminmod(qm1,qc,qp1,ZS,GMINV);\
		double num = (qc  - qm1 );\
		double denom = (qp1 - qc );\
		double r = ratio(num,denom);\
		slope1 *= phi(r,-1);\	
		HydroGrid[i][j][k].VAR##LZ[l]   = qc + slope1*(ZS/2);\
		HydroGrid[i][j][k-1].VAR##RZ[l] = qc - slope1*(ZS/2);\
	}\
}







#define VLRX(VAR )\
void VLRX_##VAR( GRID HydroGrid , int NVAR   ){\
	double eps=1e-12,Sl, Sr,S;\
	int i,j,k,l,c;\
	for(l=0;l<NVAR;l++)\
	for(j=jl;j<jr;j++)\
	for(k=kl;k<kr;k++)\
	for(i=1;i<XCM-1;i++){\
		double qm1,qc,qp1;\
		qm1=HydroGrid[i-1][j][k].VAR[l];\
		qc=HydroGrid[i][j][k].VAR[l];\
		qp1=HydroGrid[i+1][j][k].VAR[l];\
		Sl = (qc-qm1)/XS; Sr = (qp1-qc)/XS;\
		S = (SGN(Sl)+SGN(Sr))* ((fabs(Sl)* fabs(Sr))/ (fabs(Sl)+fabs(Sr)+eps));\
		HydroGrid[i][j][k].VAR##LX[l] =  qc + S*(XS/2);\
		HydroGrid[i-1][j][k].VAR##RX[l] = qc - S*(XS/2);\			
	}\
}


#define VLRY(VAR )\
void VLRY_##VAR( GRID HydroGrid , int NVAR   ){\
	double eps=1e-12,Sl, Sr,S;\
	int i,j,k,l,c;\
	for(l=0;l<NVAR;l++)\	
	for(i=il;i<ir;i++)\
	for(k=kl;k<kr;k++)\
	for(j=1;j<YCM-1;j++){\
		double qm1,qc,qp1;\
		qm1=HydroGrid[i][j-1][k].VAR[l];\
		qc=HydroGrid[i][j][k].VAR[l];\
		qp1=HydroGrid[i][j+1][k].VAR[l];\
		Sl = (qc-qm1)/YS; Sr = (qp1-qc)/YS;\
		S = (SGN(Sl)+SGN(Sr))* ((fabs(Sl)* fabs(Sr))/ (fabs(Sl)+fabs(Sr)+eps));\
		HydroGrid[i][j][k].VAR##LY[l] =  qc + S*(YS/2);\
		HydroGrid[i][j-1][k].VAR##RY[l]= qc -  S*(YS/2);\
	}\
}

#define VLRZ(VAR )\
void VLRZ_##VAR( GRID HydroGrid  , int NVAR  ){\
	double eps=1e-12,Sl, Sr,S;\
	int i,j,k,l,c;\
	for(l=0;l<NVAR;l++)\	
	for(i=il;i<ir;i++)\
	for(j=jl;j<jr;j++)\
	for(k=1;k<ZCM-1;k++){\
		double qm1,qc,qp1;\
		qm1=HydroGrid[i][j][k-1].VAR[l];\
		qc=HydroGrid[i][j][k].VAR[l];\
		qp1=HydroGrid[i][j][k+1].VAR[l];\
		Sl = (qc-qm1)/ZS;Sr = (qp1-qc)/ZS;\
		S = (SGN(Sl)+SGN(Sr))* ((fabs(Sl)* fabs(Sr))/ (fabs(Sl)+fabs(Sr)+eps));\
		HydroGrid[i][j][k] .VAR##LZ[l]=  qc + S*(ZS/2);\
		HydroGrid[i][j][k-1].VAR##RZ[l]= qc -  S*(ZS/2);\
	}\
}

#define WENOX(VAR )\
void WENOX_##VAR(GRID HydroGrid , int NVAR  ){\
	double w[3], q[3], d[3], alpha[3];\
	double wt[3], qt[3], dt[3], alphat[3], beta[4];\
	double p,eps=WENOEPS,sum;\
	int i,j,k,l,c;\	
	d[0]= dt[2]=0.3;  d[1]= dt[1]=0.6;  d[2]= dt[0]=0.1;\
	for(l=0;l<NVAR;l++)\
	for(j=jl;j<jr;j++)\
	for(k=kl;k<kr;k++)\
	for(i=1;i<XCM-1;i++){\
		double qm2,qm1,qc,qp1,qp2,left=0,right=0;\			
		qc=HydroGrid[i][j][k].VAR[l]; qm1=HydroGrid[i-1][j][k].VAR[l]; qp1=HydroGrid[i+1][j][k].VAR[l];\
		if(i==1){qm2=qm1; qp2=HydroGrid[i+2][j][k].VAR[l];}\
		else if(i==XCM-2){qm2=HydroGrid[i-2][j][k].VAR[l]; qp2=qp1;}\
		else {qm2=HydroGrid[i-2][j][k].VAR[l]; qp2=HydroGrid[i+2][j][k].VAR[l];}\
		HydroGrid[i][j][k].VAR##LX[l] = genWENOL(qm2,qm1,qc,qp1,qp2);\
		HydroGrid[i-1][j][k].VAR##RX[l] = genWENOR(qm2,qm1,qc,qp1,qp2);\
	}\
} 

#define WENOY(VAR )\
void WENOY_##VAR(GRID HydroGrid , int NVAR  ){\
	double w[3], q[3], d[3], alpha[3];\
	double wt[3], qt[3], dt[3], alphat[3], beta[4];\
	double p,eps=WENOEPS,sum;\
	int i,j,k,l,c;\	
	d[0]= dt[2]=0.3;  d[1]= dt[1]=0.6;  d[2]= dt[0]=0.1;\
	for(l=0;l<NVAR;l++)\
	for(i=il;i<ir;i++)\
	for(k=kl;k<kr;k++)\
	for(j=1;j<YCM-1;j++){\
		double qm2,qm1,qc,qp1,qp2,left=0,right=0;\			
		qc=HydroGrid[i][j][k].VAR[l]; qm1=HydroGrid[i][j-1][k].VAR[l]; qp1=HydroGrid[i][j+1][k].VAR[l];\
		if(j==1){qm2=qm1; qp2=HydroGrid[i][j+2][k].VAR[l];}\
		else if(j==YCM-2){qm2=HydroGrid[i][j-2][k].VAR[l]; qp2=qp1;}\
		else {qm2=HydroGrid[i][j-2][k].VAR[l]; qp2=HydroGrid[i][j+2][k].VAR[l];}\
		HydroGrid[i][j][k].VAR##LY[l] = genWENOL(qm2,qm1,qc,qp1,qp2);\
		HydroGrid[i][j-1][k].VAR##RY[l] = genWENOR(qm2,qm1,qc,qp1,qp2);\
	}\
} 
		
		
#define WENOZ(VAR )\
void WENOZ_##VAR(GRID HydroGrid , int NVAR  ){\
	double w[3], q[3], d[3], alpha[3];\
	double wt[3], qt[3], dt[3], alphat[3], beta[4];\
	double p,eps=WENOEPS,sum;\
	int i,j,k,l,c;\	
	d[0]= dt[2]=0.3;  d[1]= dt[1]=0.6;  d[2]= dt[0]=0.1;\
	for(l=0;l<NVAR;l++)\
	for(i=il;i<ir;i++)\
	for(j=jl;j<jr;j++)\
	for(k=1;k<ZCM-1;k++){\
		double qm2,qm1,qc,qp1,qp2,left=0,right=0;\			
		qc=HydroGrid[i][j][k].VAR[l]; qm1=HydroGrid[i][j][k-1].VAR[l]; qp1=HydroGrid[i][j][k+1].VAR[l];\
		if(k==1){qm2=qm1; qp2=HydroGrid[i][j][k+2].VAR[l];}\
		else if(k==ZCM-2){qm2=HydroGrid[i][j][k-2].VAR[l]; qp2=qp1;}\
		else {qm2=HydroGrid[i][j][k-2].VAR[l]; qp2=HydroGrid[i][j][k+2].VAR[l];}\
		HydroGrid[i][j][k].VAR##LZ[l] = genWENOL(qm2,qm1,qc,qp1,qp2);\
		HydroGrid[i][j][k-1].VAR##RZ[l] = genWENOR(qm2,qm1,qc,qp1,qp2);\
	}\
} 
		
WENOX(Fx );
WENOY(Fy );
WENOZ(Fz );
WENOX(Var);
WENOY(Var);
WENOZ(Var);
WENOX(Ax );
WENOY(Ay );
WENOZ(Az ); 

				
VLRX(Fx );
VLRY(Fy );
VLRZ(Fz );
VLRX(Var);
VLRY(Var);
VLRZ(Var);
VLRX(Ax );
VLRY(Ay );
VLRZ(Az ); 	
			

			
				
GMMX(Fx );
GMMY(Fy );
GMMZ(Fz );
GMMX(Var);
GMMY(Var);
GMMZ(Var); 
GMMX(Ax );
GMMY(Ay );
GMMZ(Az ); 	
		

		
void FindEigenValuesX(GRID HydroGrid, double tau)
{

	int i,j,k;
	
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		DECLePPIa;
		DECLp10u4;
		DECLTmu0;
		

#if defined CON
		HydroGrid[i][j][k].Ax[0] = u1/u0;
		HydroGrid[i][j][k].Ax[1] = (u1/u0) 
					+  (  
							(   -u0*(2*p2*u0+P*u1)   + a*(e*u0*u1+2*p2*(-1+u0*u0-2*u1*u1) ))
							/
							(     ( (e+P-PI)*(u0*u0)*(-u0*u0 + a*(-1+u0*u0) )  )           )
					   ) 
					; 
		
#else
		HydroGrid[i][j][k].Ax[0] = u1/u0;
		continue;
		
		double buf[SVAR] = {0,0,0,0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,PI};
		double denom = tau*(e+P-PI)*u0*u0*(-u0*u0+a*(-1+u0*u0));
		double col[SVAR] = { 
			(1+a)*u1*u0, 
			-u0*u0 - a + a*u0*u0 - 2*a*u1*u1,
			-2*a*u1*u2,
			-2*a*u1*u3*tau*tau,
			-(1+a)*u0*u1*tau,
			(u0*u0 + a - a*u0*u0 + 2*a*u1*u1)*tau,
			2*a*u1*u2*tau,
			2*a*u1*u3*tau*tau*tau,
			0,0,0,0,0,0,
			-u0*u1*tau
		};
		gsl_matrix *matrix;
	  
		matrix = gsl_matrix_calloc(SVAR, SVAR);
	 

		for(int m = 4; m < SVAR; m++)
		for(int n = 0; n < SVAR; n++)
			gsl_matrix_set(matrix, m, n,  (buf[m]*col[n])/denom);
		
		
		for(int m = 4; m < SVAR; m++)
			gsl_matrix_set( 
							matrix, m, m,  
							gsl_matrix_get(matrix,m,m) + (u1/u0)
						);
						
		
		gsl_matrix_set(matrix,0,1,1);
		denom = u0*u0*u0+a*u0-a*u0*u0*u0;
		
		int l=1;
		gsl_matrix_set(matrix,l,0,(  a*u0*(-1 + 2*pow(u0,2)) - (1 + a)*u0*pow(u1,2)     )/denom);
		gsl_matrix_set(matrix,l,1,(  2*u1*(pow(u0,2) + a*(1 - 2*pow(u0,2) + pow(u1,2)))     )/denom);
		gsl_matrix_set(matrix,l,2,(  2*a*u2*(-pow(u0,2) + pow(u1,2))    )/denom);
		gsl_matrix_set(matrix,l,3,(  2*a*u3*pow(tau,2)*(-pow(u0,2) + pow(u1,2))     )/denom);
		gsl_matrix_set(matrix,l,4,(  tau*u0*(pow(u1,2) + a*(1 - 2*pow(u0,2) + pow(u1,2)))     )/denom);
		gsl_matrix_set(matrix,l,5,(  -2*tau*u1*(pow(u0,2) + a*(1 - 2*pow(u0,2) + pow(u1,2)))     )/denom);
		gsl_matrix_set(matrix,l,6,(  2*a*tau*(u0 - u1)*(u0 + u1)*u2     )/denom);
		gsl_matrix_set(matrix,l,7,(  2*a*(u0 - u1)*(u0 + u1)*u3*pow(tau,3)     )/denom);
		gsl_matrix_set(matrix,l,8,       tau);
		gsl_matrix_set(matrix,l,14,( tau*u0*(-pow(u0,2) + pow(u1,2))   )/denom);


		l=2;
		gsl_matrix_set(matrix,l,0,(  -((1 + a)*u0*u1*u2)     )/denom);
		gsl_matrix_set(matrix,l,1,( u2*(a + pow(u0,2) - a*pow(u0,2) + 2*a*pow(u1,2))      )/denom);
		gsl_matrix_set(matrix,l,2,( u1*(a + pow(u0,2) - a*pow(u0,2) + 2*a*pow(u2,2))      )/denom);
		gsl_matrix_set(matrix,l,3,( 2*a*u1*u2*u3*pow(tau,2)      )/denom);
		gsl_matrix_set(matrix,l,4,( (1 + a)*tau*u0*u1*u2      )/denom);
		gsl_matrix_set(matrix,l,5,( tau*u2*(-pow(u0,2) + a*(-1 + pow(u0,2) - 2*pow(u1,2)))      )/denom);
		gsl_matrix_set(matrix,l,6,( tau*u1*(-pow(u0,2) + a*(-1 + pow(u0,2) - 2*pow(u2,2)))      )/denom);
		gsl_matrix_set(matrix,l,7,( -2*a*u1*u2*u3*pow(tau,3)      )/denom);
		gsl_matrix_set(matrix,l,9,tau);
		gsl_matrix_set(matrix,l,14,(  tau*u0*u1*u2   )/denom);
		
		
		l=3;
		gsl_matrix_set(matrix,l,0,(  -((1 + a)*u0*u1*u3)     )/denom);
		gsl_matrix_set(matrix,l,1,(  u3*(a + pow(u0,2) - a*pow(u0,2) + 2*a*pow(u1,2))     )/denom);
		gsl_matrix_set(matrix,l,2,( 2*a*u1*u2*u3      )/denom);
		gsl_matrix_set(matrix,l,3,(  u1*(a + pow(u0,2) - a*pow(u0,2) + 2*a*pow(tau,2)*pow(u3,2))     )/denom);
		gsl_matrix_set(matrix,l,4,( (1 + a)*tau*u0*u1*u3      )/denom);
		gsl_matrix_set(matrix,l,5,( tau*u3*(-pow(u0,2) + a*(-1 + pow(u0,2) - 2*pow(u1,2)))      )/denom);
		gsl_matrix_set(matrix,l,6,( -2*a*tau*u1*u2*u3      )/denom);
		gsl_matrix_set(matrix,l,7,(  -(tau*u1*(a + pow(u0,2) - a*pow(u0,2) + 2*a*pow(tau,2)*pow(u3,2)))     )/denom);
		gsl_matrix_set(matrix,l,10,tau);
		gsl_matrix_set(matrix,l,14,(   u0*u1*u2*tau     )/denom);
		
		
		
	
		gsl_eigen_nonsymm_workspace	*workspace = gsl_eigen_nonsymm_alloc(SVAR);	
		gsl_vector_complex  *eig = gsl_vector_complex_alloc(SVAR);	
		gsl_eigen_nonsymm(matrix, eig, workspace);	 
		gsl_complex c_element;
		for(int m = 0; m < SVAR; m++)
		{
			c_element = gsl_vector_complex_get(eig, m);				
			HydroGrid[i][j][k].Ax[m] = fabs(GSL_REAL(c_element));
			//~ cout<<m<<std::scientific<<" -->"<<GSL_REAL(c_element)<<"+ i * "<<GSL_IMAG(c_element)<<endl;	
			
			if(fabs(GSL_IMAG(c_element) >1e-6))
				HydroGrid[i][j][k].Ax[m] = u1/u0;
				
			HydroGrid[i][j][k].temp[m] = HydroGrid[i][j][k].Ax[m];
			
		}
		gsl_eigen_nonsymm_free(workspace);
		gsl_matrix_free(matrix);
		gsl_vector_complex_free(eig);
#endif

	}
}
		
void FindEigenValuesY(GRID HydroGrid, double tau)
{
	int i,j,k;
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		DECLePPIa;
		DECLp10u4;
		DECLTmu0;
		

#if defined CON
		HydroGrid[i][j][k].Ay[0]=u2/u0;
		HydroGrid[i][j][k].Ay[1]=(u2/u0) 
					+  (  (-u0*(2*p3*u0+P*u2)   + a*(e*u0*u2+2*p3*(-1+u0*u0-2*u2*u2) ))
							/
						(  ( (e+P-PI)*(u0*u0)*(-u0*u0 + a*(-1+u0*u0) )  ) )
						) 
					; 
#else
		HydroGrid[i][j][k].Ay[0] = u2/u0;
		continue;
		
		double buf[SVAR] = {0,0,0,0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,PI};
		double denom = tau*(e+P-PI)*u0*u0*(-u0*u0+a*(-1+u0*u0));
		double col[SVAR] = { 
			(1+a)*u2*u0, 
			-2*a*u1*u2,
			-u0*u0 - a + a*u0*u0 - 2*a*u2*u2,
			-2*a*u2*u3*tau*tau,
			-(1+a)*u0*u2*tau,
			2*a*u1*u2*tau,
			(u0*u0 + a - a*u0*u0 + 2*a*u2*u2)*tau,
			2*a*u2*u3*tau*tau*tau,
			0,0,0,0,0,0,
			-u0*u2*tau
		};

		gsl_matrix *matrix = gsl_matrix_calloc(SVAR, SVAR);


		for(int m = 4; m < SVAR; m++)
		for(int n = 0; n < SVAR; n++)
			gsl_matrix_set(matrix, m, n,  (buf[m]*col[n])/denom);
		
		
		for(int m = 4; m < SVAR; m++)
			gsl_matrix_set( 
							matrix, m, m,  
							gsl_matrix_get(matrix,m,m) + (u2/u0)
						);
		
		gsl_matrix_set(matrix,0,2,1);
		denom = u0*u0*u0+a*u0-a*u0*u0*u0;
	
		int l=1;
		gsl_matrix_set(matrix,l,0,(-((1 + a)*u0*u1*u2)        )/denom);
		gsl_matrix_set(matrix,l,1,( u2*(a + pow(u0,2) - a*pow(u0,2) + 2*a*pow(u1,2))        )/denom);
		gsl_matrix_set(matrix,l,2,( u1*(a + pow(u0,2) - a*pow(u0,2) + 2*a*pow(u2,2))      )/denom);
		gsl_matrix_set(matrix,l,3,(  2*a*u1*u2*u3*pow(tau,2)      )/denom);
		gsl_matrix_set(matrix,l,4,( (1 + a)*tau*u0*u1*u2    )/denom);
		gsl_matrix_set(matrix,l,5,( tau*u2*(-pow(u0,2) + a*(-1 + pow(u0,2) - 2*pow(u1,2)))  )/denom);
		gsl_matrix_set(matrix,l,6,( tau*u1*(-pow(u0,2) + a*(-1 + pow(u0,2) - 2*pow(u2,2)))   )/denom);
		gsl_matrix_set(matrix,l,7,( -2*a*u1*u2*u3*pow(tau,3)      )/denom);
		gsl_matrix_set(matrix,l,9,       tau);
		gsl_matrix_set(matrix,l,14,( tau*u0*u1*u2   )/denom);


		l=2;
		gsl_matrix_set(matrix,l,0,(a*u0*(-1 + 2*pow(u0,2)) - (1 + a)*u0*pow(u2,2)       )/denom);
		gsl_matrix_set(matrix,l,1,( -2*a*u1*(u0 - u2)*(u0 + u2)      )/denom);
		gsl_matrix_set(matrix,l,2,( 2*u2*(pow(u0,2) + a*(1 - 2*pow(u0,2) + pow(u2,2)))      )/denom);
		gsl_matrix_set(matrix,l,3,( -2*a*(u0 - u2)*(u0 + u2)*u3*pow(tau,2)      )/denom);
		gsl_matrix_set(matrix,l,4,(tau*u0*(pow(u2,2) + a*(1 - 2*pow(u0,2) + pow(u2,2)))       )/denom);
		gsl_matrix_set(matrix,l,5,( 2*a*tau*u1*(u0 - u2)*(u0 + u2)      )/denom);
		gsl_matrix_set(matrix,l,6,( -2*tau*u2*(pow(u0,2) + a*(1 - 2*pow(u0,2) + pow(u2,2)))      )/denom);
		gsl_matrix_set(matrix,l,7,( 2*a*(u0 - u2)*(u0 + u2)*u3*pow(tau,3)      )/denom);
		gsl_matrix_set(matrix,l,11,             tau);
		gsl_matrix_set(matrix,l,14,(tau*u0*(-pow(u0,2) + pow(u2,2))       )/denom);
		
		
		l=3;
		gsl_matrix_set(matrix,l,0,( -((1 + a)*u0*u2*u3)      )/denom);
		gsl_matrix_set(matrix,l,1,( 2*a*u1*u2*u3      )/denom);
		gsl_matrix_set(matrix,l,2,(u3*(a + pow(u0,2) - a*pow(u0,2) + 2*a*pow(u2,2))       )/denom);
		gsl_matrix_set(matrix,l,3,( u2*(a + pow(u0,2) - a*pow(u0,2) + 2*a*pow(tau,2)*pow(u3,2))      )/denom);
		gsl_matrix_set(matrix,l,4,(  (1 + a)*tau*u0*u2*u3     )/denom);
		gsl_matrix_set(matrix,l,5,(  -2*a*tau*u1*u2*u3     )/denom);
		gsl_matrix_set(matrix,l,6,(tau*u3*(-pow(u0,2) + a*(-1 + pow(u0,2) - 2*pow(u2,2)))       )/denom);
		gsl_matrix_set(matrix,l,7,(-(tau*u2*(a + pow(u0,2) - a*pow(u0,2) + 2*a*pow(tau,2)*pow(u3,2)))       )/denom);
		gsl_matrix_set(matrix,l,12,           tau);
		gsl_matrix_set(matrix,l,14,( tau*u0*u2*u3     )/denom);



		gsl_eigen_nonsymm_workspace	*workspace = gsl_eigen_nonsymm_alloc(SVAR);	
		gsl_vector_complex  *eig = gsl_vector_complex_alloc(SVAR);	
		gsl_eigen_nonsymm(matrix, eig, workspace);	 
		gsl_complex c_element;
		for(int m = 0; m < SVAR; m++)
		{
			c_element = gsl_vector_complex_get(eig, m);
			//~ cout<<m<<std::scientific<<" -->"<<GSL_REAL(c_element)<<"+ i* "<<GSL_IMAG(c_element)<<endl;			
			if(fabs(GSL_IMAG(c_element) >1e-6))
				HydroGrid[i][j][k].Ay[m] = u2/u0;
				
			HydroGrid[i][j][k].temp[m] = HydroGrid[i][j][k].Ay[m];

		}
		gsl_eigen_nonsymm_free(workspace);
		gsl_matrix_free(matrix);
		gsl_vector_complex_free(eig);
#endif


	}
}

			
void FindEigenValuesZ(GRID HydroGrid, double tau)
{

	int i,j,k;
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		DECLePPIa;
		DECLp10u4;
		DECLTmu0;
		
		
#if defined CON
		HydroGrid[i][j][k].Az[0]=u3/u0;
		HydroGrid[i][j][k].Az[1]= (u3/u0) +
							( 
							  ( (-P + p1 + a*(e+p1)-(1+a)*p1*PI)*u3   ) 
							/ ( (e+P-PI)*(u0)*(-u0*u0+a*(-1+u0*u0) )  ) 
							) ;		
#else
		HydroGrid[i][j][k].Az[0] = u3/u0;
		continue;
		

		double buf[SVAR] = {0,0,0,0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,PI};
		double denom = tau*(e+P-PI)*u0*u0*(-u0*u0+a*(-1+u0*u0));
		double col[SVAR] = { 
			(1+a)*u3*u0, 
			-2*a*u1*u3,
			-2*a*u2*u3,
			-u0*u0 - a + a*u0*u0 - 2*a*u3*u3*tau*tau,
			-(1+a)*u0*u3*tau,
			2*a*u1*u3*tau,
			2*a*u2*u3*tau,
			( u0*u0 + a - a*u0*u0 + 2*a*u3*u3*tau*tau)*tau,
			0,0,0,0,0,0,
			-u0*u3*tau
		};
		
		gsl_matrix *matrix = gsl_matrix_calloc(SVAR, SVAR);
	 

		for(int m = 4; m < SVAR; m++)
		for(int n = 0; n < SVAR; n++)
			gsl_matrix_set(matrix, m, n,  (buf[m]*col[n])/denom);
		
		
		for(int m = 4; m < SVAR; m++)
			gsl_matrix_set( 
							matrix, m, m,  
							gsl_matrix_get(matrix,m,m) + (u3/u0)
						);
						
						
		gsl_matrix_set(matrix,0,3,1);
		denom = u0*u0*u0+a*u0-a*u0*u0*u0;
	
		int l=1;
		gsl_matrix_set(matrix,l,0,(   -((1 + a)*u0*u1*u3)     )/denom);
		gsl_matrix_set(matrix,l,1,(   u3*(a + pow(u0,2) - a*pow(u0,2) + 2*a*pow(u1,2))     )/denom);
		gsl_matrix_set(matrix,l,2,(   2*a*u1*u2*u3     )/denom);
		gsl_matrix_set(matrix,l,3,(  u1*(a + pow(u0,2) - a*pow(u0,2) + 2*a*pow(tau,2)*pow(u3,2))      )/denom);
		gsl_matrix_set(matrix,l,4,(   (1 + a)*tau*u0*u1*u3     )/denom);
		gsl_matrix_set(matrix,l,5,(   tau*u3*(-pow(u0,2) + a*(-1 + pow(u0,2) - 2*pow(u1,2)))     )/denom);
		gsl_matrix_set(matrix,l,6,(  -2*a*tau*u1*u2*u3      )/denom);
		gsl_matrix_set(matrix,l,7,( -(tau*u1*(a + pow(u0,2) - a*pow(u0,2) + 2*a*pow(tau,2)*pow(u3,2)))      )/denom);
		gsl_matrix_set(matrix,l,10,     tau   );
		gsl_matrix_set(matrix,l,14,(  tau*u0*u1*u3     )/denom);


		l=2;
		gsl_matrix_set(matrix,l,0,(   -((1 + a)*u0*u2*u3)    )/denom);
		gsl_matrix_set(matrix,l,1,(   2*a*u1*u2*u3    )/denom);
		gsl_matrix_set(matrix,l,2,(   u3*(a + pow(u0,2) - a*pow(u0,2) + 2*a*pow(u2,2))    )/denom);
		gsl_matrix_set(matrix,l,3,(   u2*(a + pow(u0,2) - a*pow(u0,2) + 2*a*pow(tau,2)*pow(u3,2))    )/denom);
		gsl_matrix_set(matrix,l,4,(   (1 + a)*tau*u0*u2*u3    )/denom);
		gsl_matrix_set(matrix,l,5,(   -2*a*tau*u1*u2*u3    )/denom);
		gsl_matrix_set(matrix,l,6,(   tau*u3*(-pow(u0,2) + a*(-1 + pow(u0,2) - 2*pow(u2,2)))    )/denom);
		gsl_matrix_set(matrix,l,7,(   -(tau*u2*(a + pow(u0,2) - a*pow(u0,2) + 2*a*pow(tau,2)*pow(u3,2)))    )/denom);
		gsl_matrix_set(matrix,l,12,     tau);
		gsl_matrix_set(matrix,l,14,(  tau*u0*u2*u3     )/denom);
		
		
		l=3;
		gsl_matrix_set(matrix,l,0,(a*u0*pow(tau,-2)*(-1 + 2*pow(u0,2)) + u0*(-2 + a*(-3 + 2*pow(u0,2)))*pow(u3,2)      )/denom);
		gsl_matrix_set(matrix,l,1,( 2*a*u1*(-(pow(tau,-2)*pow(u0,2)) - (-2 + pow(u0,2))*pow(u3,2))      )/denom);
		gsl_matrix_set(matrix,l,2,( 2*a*u2*(-(pow(tau,-2)*pow(u0,2)) - (-2 + pow(u0,2))*pow(u3,2))      )/denom);
		gsl_matrix_set(matrix,l,3,(  2*u3*(a + pow(u0,2) - 2*a*pow(u0,2) - a*pow(tau,2)*(-2 + pow(u0,2))*pow(u3,2))     )/denom);
		gsl_matrix_set(matrix,l,4,( pow(tau,-1)*(a*(u0 - 2*pow(u0,3)) + u0*pow(tau,2)*(2 + a*(3 - 2*pow(u0,2)))*pow(u3,2))      )/denom);
		gsl_matrix_set(matrix,l,5,( 2*a*u1*pow(tau,-1)*(pow(u0,2) + pow(tau,2)*(-2 + pow(u0,2))*pow(u3,2))      )/denom);
		gsl_matrix_set(matrix,l,6,( 2*a*u2*pow(tau,-1)*(pow(u0,2) + pow(tau,2)*(-2 + pow(u0,2))*pow(u3,2))      )/denom);
		gsl_matrix_set(matrix,l,7,(  2*tau*u3*(-pow(u0,2) + a*(-1 + 2*pow(u0,2) + pow(tau,2)*(-2 + pow(u0,2))*pow(u3,2)))     )/denom);
		gsl_matrix_set(matrix,l,13,     tau);
		gsl_matrix_set(matrix,l,14,( -(pow(tau,-1)*pow(u0,3)) - tau*u0*(-2 + pow(u0,2))*pow(u3,2)     )/denom);
		

		gsl_eigen_nonsymm_workspace	*workspace = gsl_eigen_nonsymm_alloc(SVAR);	
		gsl_vector_complex  *eig = gsl_vector_complex_alloc(SVAR);	
		gsl_eigen_nonsymm(matrix, eig, workspace);	 
		gsl_complex c_element;
		
		for(int m = 0; m < SVAR; m++)
		{
			c_element = gsl_vector_complex_get(eig, m);
			
			if(fabs(GSL_IMAG(c_element) >1e-6))
				HydroGrid[i][j][k].Az[m] = u3/u0;
				
			HydroGrid[i][j][k].temp[m] = HydroGrid[i][j][k].Az[m];		
		}
		gsl_eigen_nonsymm_free(workspace);
		gsl_matrix_free(matrix);
		gsl_vector_complex_free(eig);
#endif

	}
}



#define GMMTYPE 0	
#define VLRTYPE  1	
#define WENOTYPE 2
	
void Reconstruct(GRID HydroGrid)
{
	//~ int type = GMMTYPE;
	//~ int type = VLRTYPE;
	int type = WENOTYPE;

	switch(type)
	{
		case GMMTYPE:
		
			GMMX_Fx(HydroGrid , SVAR);
			GMMY_Fy(HydroGrid , SVAR);
			GMMZ_Fz(HydroGrid , SVAR);
			GMMX_Var(HydroGrid , SVAR);
			GMMY_Var(HydroGrid , SVAR);
			GMMZ_Var(HydroGrid , SVAR);
			
			 
			GMMX_Ax(HydroGrid , EVAR);
			GMMY_Ay(HydroGrid , EVAR);
			GMMZ_Az(HydroGrid , EVAR);
		 
			break;
			
		case VLRTYPE:
		
			VLRX_Fx(HydroGrid , SVAR);
			VLRY_Fy(HydroGrid , SVAR);
			VLRZ_Fz(HydroGrid , SVAR);
			VLRX_Var(HydroGrid , SVAR);
			VLRY_Var(HydroGrid , SVAR);
			VLRZ_Var(HydroGrid , SVAR);
			
			
			VLRX_Ax( HydroGrid , EVAR );
			VLRY_Ay( HydroGrid , EVAR );
			VLRZ_Az( HydroGrid , EVAR );
			
			
			break;
			
		case WENOTYPE:
		
			WENOX_Fx(HydroGrid , SVAR);
			WENOY_Fy(HydroGrid , SVAR);
			WENOZ_Fz(HydroGrid , SVAR);
			WENOX_Var(HydroGrid , SVAR);
			WENOY_Var(HydroGrid , SVAR);
			WENOZ_Var(HydroGrid , SVAR);
			
			WENOX_Ax( HydroGrid , EVAR );
			WENOY_Ay( HydroGrid , EVAR );
			WENOZ_Az( HydroGrid , EVAR );
						
			break;
	}
	
}
		
void FindMaxEigenValueXYZ(GRID HydroGrid)
{
	int i,j,k,l;
	
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		double maxlx=0,maxrx=0;
		double maxly=0,maxry=0;
		double maxlz=0,maxrz=0;
		
		for(l=0;l<EVAR;l++)
		{
			if(fabs(HydroGrid[i][j][k].AxLX[l]) > maxlx)
				maxlx = fabs(HydroGrid[i][j][k].AxLX[l]);
			
			if(fabs(HydroGrid[i][j][k].AxRX[l]) > maxrx)
				maxrx = fabs(HydroGrid[i][j][k].AxRX[l]);
				
			
			if(fabs(HydroGrid[i][j][k].AyLY[l]) > maxly)
				maxly = fabs(HydroGrid[i][j][k].AyLY[l]);
			
			if(fabs(HydroGrid[i][j][k].AyRY[l]) > maxry)
				maxry = fabs(HydroGrid[i][j][k].AyRY[l]);
				
				
			if(fabs(HydroGrid[i][j][k].AzLZ[l]) > maxly)
				maxly = fabs(HydroGrid[i][j][k].AzLZ[l]);
			
			if(fabs(HydroGrid[i][j][k].AzRZ[l]) > maxry)
				maxry = fabs(HydroGrid[i][j][k].AzRZ[l]);
		}
				
		HydroGrid[i][j][k].AxLXMAX = maxlx;
		HydroGrid[i][j][k].AxRXMAX = maxrx;
		HydroGrid[i][j][k].AyLYMAX = maxly;
		HydroGrid[i][j][k].AyRYMAX = maxry; 
		HydroGrid[i][j][k].AzLZMAX = maxlz;
		HydroGrid[i][j][k].AzRZMAX = maxrz;
	}
}
			

void CFL(GRID HydroGrid, double tau, double taustep)
{
	int i,j,k,l;
	
	double maxflowgrid=0;
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		double xflow = MAX(   HydroGrid[i][j][k].AxLXMAX  ,   HydroGrid[i][j][k].AxRXMAX  );
		double yflow = MAX(   HydroGrid[i][j][k].AyLYMAX  ,   HydroGrid[i][j][k].AyRYMAX  );
		
		double ttt= sqrt(xflow*xflow+yflow*yflow);
		
		if(ttt>maxflowgrid)
			maxflowgrid = ttt;
	}
	
	double timesug = XS/(4*maxflowgrid);
	double mintime;
	
    MPI_Reduce(&timesug, &mintime, 1, MPI_DOUBLE,MPI_MIN, 0, MPI_COMM_WORLD);
    
	if(!rank && mintime<TS)
		cout<<std::scientific<<"The time step based on CFL should be "<<mintime<<endl;
}
			



void CopyPrimaryVariablesToVar(GRID HydroGrid, double tau)
{

	for( int i=0; i<XCM ; i++)
	for( int j=0; j<YCM ; j++)
	for( int k=0; k<ZCM ; k++)
	{
		HydroGrid[i][j][k].Var[0]= tau*HydroGrid[i][j][k].T00;
		HydroGrid[i][j][k].Var[1]= tau*HydroGrid[i][j][k].T10;
		HydroGrid[i][j][k].Var[2]= tau*HydroGrid[i][j][k].T20;
		HydroGrid[i][j][k].Var[3]= tau*HydroGrid[i][j][k].T30;
		
		for(int l=0;l<Npi;l++)
			HydroGrid[i][j][k].Var[VARN+l]=  HydroGrid[i][j][k].pi[l];
		
		HydroGrid[i][j][k].Var[VARN+Npi]=  HydroGrid[i][j][k].PI;
	}
}


void ClearResultVariable(GRID HydroGrid)
{
	for( int i=0; i<XCM ; i++)
	for( int j=0; j<YCM ; j++)
	for( int k=0; k<ZCM ; k++)
	for( int l=0; l<SVAR; l++)
		HydroGrid[i][j][k].Result[l] = 0;
}

void AddPartialResultToFinalResult(GRID HydroGrid)
{
	int i,j,k;
	int l;
	
	for(l=0;l<SVAR;l++)
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
		HydroGrid[i][j][k].Result[l] += HydroGrid[i][j][k].PartialResult[l];
}


void hydroExplicit(GRID HydroGrid, double tau, double taustep)
{
	ClearResultVariable( HydroGrid);
	CopyPrimaryVariablesToVar( HydroGrid , tau);
	
	FindEigenValuesX(HydroGrid, tau);
	FindEigenValuesY(HydroGrid, tau);

#if !defined LBI
	FindEigenValuesZ(HydroGrid, tau);
#endif


	Reconstruct(HydroGrid);
	

	FindMaxEigenValueXYZ(HydroGrid);
	
	
	fvX( HydroGrid, tau );
	AddPartialResultToFinalResult( HydroGrid);
	
	fvY( HydroGrid, tau );
	AddPartialResultToFinalResult( HydroGrid);
	

#if !defined LBI
	fvZ( HydroGrid, tau);
	AddPartialResultToFinalResult( HydroGrid);	
#endif
	
	CFL(HydroGrid, tau ,  taustep);
}




//The actual hydro code


void fvX(GRID HydroGrid, double tau)
{

	int i,j,k,l;
 
	 	
	for( l=0; l<SVAR; l++)
	for( j=jl; j<jr ; j++)
	for( k=kl; k<kr ; k++) 
	for( i=il-1; i<ir ; i++)
		HydroGrid[i][j][k].fluxT[l] =   0.5*( HydroGrid[i][j][k].FxLX[l] + HydroGrid[i][j][k].FxRX[l] )
				- 0.5 *(  MAX(  HydroGrid[i][j][k].AxLXMAX , HydroGrid[i][j][k].AxRXMAX  ) ) * ( HydroGrid[i][j][k].VarRX[l] - HydroGrid[i][j][k].VarLX[l] ); 
		
	
	for( l=0; l<SVAR; l++)
	for( j=jl;j<jr;j++)
	for( k=kl;k<kr;k++)
	for( i=il;i<ir;i++)
		HydroGrid[i][j][k].PartialResult[l] =  -(1.0/XS)*(HydroGrid[i][j][k].fluxT[l] - HydroGrid[i-1][j][k].fluxT[l]);

}

void fvY(GRID HydroGrid, double tau)
{
	int i,j,k,l;	
 
	 	
	for( l=0; l<SVAR; l++)
	for( i=il; i< ir;  i++)
	for( k=kl; k< kr;  k++)
	for( j=jl-1 ; j<jr; j++)
		HydroGrid[i][j][k].fluxT[l] =   0.5*( HydroGrid[i][j][k].FyLY[l] + HydroGrid[i][j][k].FyRY[l] )
				- 0.5 *(  MAX( fabs( HydroGrid[i][j][k].AyLYMAX  ), fabs(  HydroGrid[i][j][k].AyRYMAX ) ) ) * ( HydroGrid[i][j][k].VarRY[l] - HydroGrid[i][j][k].VarLY[l] ); 
		
	
	for( l=0; l<SVAR; l++)
	for( i=il; i< ir;  i++)
	for( k=kl; k< kr;  k++)
	for( j=jl; j<jr; j++)
		HydroGrid[i][j][k].PartialResult[l] =  -(1.0/YS)*(HydroGrid[i][j][k].fluxT[l] - HydroGrid[i][j-1][k].fluxT[l]);

}

void fvZ(GRID HydroGrid,  double tau)
{
	int i,j,k,l;	
 
	 	
	for( l=0; l<SVAR; l++)
	for( i=il; i< ir; i++)
	for( j=jl; j< jr; j++)
	for( k=kl-1; k< kr; k++) 
		HydroGrid[i][j][k].fluxT[l] =   0.5*( HydroGrid[i][j][k].FzLZ[l] + HydroGrid[i][j][k].FzRZ[l] )
				- 0.5 *(  MAX( fabs( HydroGrid[i][j][k].AzLZMAX   ), fabs( HydroGrid[i][j][k].AzRZMAX  ) ) ) * ( HydroGrid[i][j][k].VarRZ[l] - HydroGrid[i][j][k].VarLZ[l] ); 
		
	
	for( l=0; l<SVAR; l++)
	for( i=il; i< ir; i++)
	for( j=jl; j< jr; j++)
	for( k=kl; k< kr; k++)
		HydroGrid[i][j][k].PartialResult[l] =  -(1.0/ZS)*(HydroGrid[i][j][k].fluxT[l] - HydroGrid[i][j][k-1].fluxT[l]);
}
