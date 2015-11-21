void hydroExplicit(GRID HydroGrid, double tau, double taustep);
void fvX(GRID HydroGrid, double taustep, double tau);
void fvY(GRID HydroGrid,  double taustep, double tau);
void fvZ(GRID HydroGrid,  double taustep, double tau);


#define WENOP 2
#define WENOEPS 1E-6


int method;
int type;




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
		

		
void FindEigenValuesX(GRID HydroGrid, double tau, double taustep)
{

	int i,j,k;
	
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		DECLePPIa;
		DECLp5u4;
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
	 

#endif

	}
}
		
void FindEigenValuesY(GRID HydroGrid, double tau, double taustep)
{
	int i,j,k;
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		DECLePPIa;
		DECLp5u4;
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

#endif


	}
}

			
void FindEigenValuesZ(GRID HydroGrid, double tau, double taustep)
{

	int i,j,k;
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		DECLePPIa;
		DECLp5u4;
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
#endif

	}
}



#define GMMTYPE 0	
#define VLRTYPE  1	
#define WENOTYPE 2
	
void Reconstruct(GRID HydroGrid, double tau, double taustep)
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
		
void FindMaxEigenValueXYZ(GRID HydroGrid, double tau, double taustep)
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
	
	cout<<std::scientific<<"The time step based on CFL should be "<< XS/(4*maxflowgrid)<<endl;
}
			




void hydroExplicit(GRID HydroGrid, double tau, double taustep)
{
	ClearResultVariable( HydroGrid);
	CopyPrimaryVariablesToVar( HydroGrid , tau);
	
	FindEigenValuesX(HydroGrid, tau,  taustep);
	FindEigenValuesY(HydroGrid, tau,  taustep);

#if !defined LBI
	FindEigenValuesZ(HydroGrid, tau,  taustep);
#endif


	Reconstruct(HydroGrid, tau ,  taustep);
	

	FindMaxEigenValueXYZ(HydroGrid, tau ,  taustep);
	
	
	fvX( HydroGrid,  tau ,  taustep);
	AddPartialResultToFinalResult( HydroGrid);
	
	fvY( HydroGrid, tau ,  taustep);
	AddPartialResultToFinalResult( HydroGrid);
	

#if !defined LBI
	fvZ( HydroGrid,  tau ,  taustep);
	AddPartialResultToFinalResult( HydroGrid);
	
#endif
	
	//~ CFL(HydroGrid, tau ,  taustep);
}




//The actual hydro code


void fvX(GRID HydroGrid, double tau, double taustep)
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
		HydroGrid[i][j][k].PartialResult[l] =  -(taustep/XS)*(HydroGrid[i][j][k].fluxT[l] - HydroGrid[i-1][j][k].fluxT[l]);

}

void fvY(GRID HydroGrid, double tau, double taustep)
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
		HydroGrid[i][j][k].PartialResult[l] =  -(taustep/YS)*(HydroGrid[i][j][k].fluxT[l] - HydroGrid[i][j-1][k].fluxT[l]);

}

void fvZ(GRID HydroGrid,  double tau, double taustep)
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
		HydroGrid[i][j][k].PartialResult[l] =  -(taustep/ZS)*(HydroGrid[i][j][k].fluxT[l] - HydroGrid[i][j][k-1].fluxT[l]);
}
