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
		DECLu4;
		DECLePa;
 
		HydroGrid[i][j][k].Fx[0] =  (  u1*u0*(e+P)            );
		HydroGrid[i][j][k].Fx[1] =  (  P + u1*u1*(e+P)        );
		HydroGrid[i][j][k].Fx[2] =  (  u1*u2*(e+P)            );
		HydroGrid[i][j][k].Fx[3] =  (  u1*u3*(e+P)            );  

 
		HydroGrid[i][j][k].Fy[0] =  (  u2*u0*(e+P)            );
		HydroGrid[i][j][k].Fy[1] =  (  u2*u1*(e+P)            );
		HydroGrid[i][j][k].Fy[2] =  (  P + u2*u2*(e+P)        );
		HydroGrid[i][j][k].Fy[3] =  (  u2*u3*(e+P)            );  

 
  		HydroGrid[i][j][k].Fz[0] =  (  u3*u0*(e+P)            );
		HydroGrid[i][j][k].Fz[1] =  (  u3*u1*(e+P)            );
		HydroGrid[i][j][k].Fz[2] =  (  u3*u2*(e+P)            );
		HydroGrid[i][j][k].Fz[3] =  (  P + u3*u3*(e+P)        );	  
	}	
}




#define WENOX(VAR )\
void WENOX_##VAR(GRID HydroGrid , int NVAR  ){\  
	int i,j,k,l;\	 
	for(l=0;l<NVAR;l++)\
	for(j=jl;j<jr;j++)\
	for(k=kl;k<kr;k++)\
	for(i=1;i<XCM;i++){\
		double qm2,qm1,qc,qp1,qp2;\			
		qc=HydroGrid[i][j][k].VAR[l];\
		if(i>1 && i<XCM-2){qm1=HydroGrid[i-1][j][k].VAR[l];qm2=HydroGrid[i-2][j][k].VAR[l];qp1=HydroGrid[i+1][j][k].VAR[l];qp2=HydroGrid[i+2][j][k].VAR[l];}\
		else if(i==0){qm2=qm1=qc;qp1=HydroGrid[i+1][j][k].VAR[l];qp2=HydroGrid[i+2][j][k].VAR[l];}\
		else if(i==1){qm2=qm1=HydroGrid[i-1][j][k].VAR[l]; qp1=HydroGrid[i+1][j][k].VAR[l];qp2=HydroGrid[i+2][j][k].VAR[l];}\
		else if(i==XCM-2){qm2=HydroGrid[i-2][j][k].VAR[l]; qm1=HydroGrid[i-1][j][k].VAR[l];qp1=qp2=HydroGrid[i+1][j][k].VAR[l];}\
		else if(i==XCM-1){qm2=HydroGrid[i-2][j][k].VAR[l]; qm1=HydroGrid[i-1][j][k].VAR[l];qp1=qp2=qc;}\
		HydroGrid[i][j][k].VAR##LX[l] = genWENOL(qm2,qm1,qc,qp1,qp2);HydroGrid[i-1][j][k].VAR##RX[l] = genWENOR(qm2,qm1,qc,qp1,qp2);\
	}\
} 

#define WENOY(VAR )\
void WENOY_##VAR(GRID HydroGrid , int NVAR  ){\  
	int i,j,k,l;\	 
	for(l=0;l<NVAR;l++)\
	for(i=il;i<ir;i++)\
	for(k=kl;k<kr;k++)\
	for(j=1;j<YCM;j++){\
		double qm2,qm1,qc,qp1,qp2;\			
		qc=HydroGrid[i][j][k].VAR[l];\
		if(j>1 && j<YCM-2){qm1=HydroGrid[i][j-1][k].VAR[l];qm2=HydroGrid[i][j-2][k].VAR[l];qp1=HydroGrid[i][j+1][k].VAR[l];qp2=HydroGrid[i][j+2][k].VAR[l];}\
		else if(j==0){qm2=qm1=qc;qp1=HydroGrid[i][j+1][k].VAR[l];qp2=HydroGrid[i][j+2][k].VAR[l];}\
		else if(j==1){qm2=qm1=HydroGrid[i][j-1][k].VAR[l]; qp1=HydroGrid[i][j+1][k].VAR[l];qp2=HydroGrid[i][j+2][k].VAR[l];}\
		else if(j==YCM-2){qm2=HydroGrid[i][j-2][k].VAR[l]; qm1=HydroGrid[i][j-1][k].VAR[l];qp1=qp2=HydroGrid[i][j+1][k].VAR[l];}\
		else if(j==YCM-1){qm2=HydroGrid[i][j-2][k].VAR[l]; qm1=HydroGrid[i][j-1][k].VAR[l];qp1=qp2=qc;}\
		HydroGrid[i][j][k].VAR##LY[l] = genWENOL(qm2,qm1,qc,qp1,qp2);	HydroGrid[i][j-1][k].VAR##RY[l] = genWENOR(qm2,qm1,qc,qp1,qp2);\
	}\
} 
		
		
#if !defined LBI

#define WENOZ(VAR )\
void WENOZ_##VAR(GRID HydroGrid , int NVAR  ){\
	int i,j,k,l;\	
	for(l=0;l<NVAR;l++)\
	for(i=il;i<ir;i++)\
	for(j=jl;j<jr;j++)\
	for(k=1;k<ZCM;k++){\
		double qm2,qm1,qc,qp1,qp2;\			
		qc=HydroGrid[i][j][k].VAR[l];\
		if(k>1 && k<ZCM-2){qm1=HydroGrid[i][j][k-1].VAR[l];qm2=HydroGrid[i][j][k-2].VAR[l];qp1=HydroGrid[i][j][k+1].VAR[l];qp2=HydroGrid[i][j][k+2].VAR[l];}\
		else if(k==0){qm2=qm1=qc;qp1=HydroGrid[i][j][k+1].VAR[l];qp2=HydroGrid[i][j][k+2].VAR[l];}\
		else if(k==1){qm2=qm1=HydroGrid[i][j][k-1].VAR[l]; qp1=HydroGrid[i][j][k+1].VAR[l];qp2=HydroGrid[i][j][k+2].VAR[l];}\
		else if(k==ZCM-2){qm2=HydroGrid[i][j][k-2].VAR[l]; qm1=HydroGrid[i][j][k-1].VAR[l];qp1=qp2=HydroGrid[i][j][k+1].VAR[l];}\
		else if(k==ZCM-1){qm2=HydroGrid[i][j][k-2].VAR[l]; qm1=HydroGrid[i][j][k-1].VAR[l];qp1=qp2=qc;}\ 
		HydroGrid[i][j][k].VAR##LZ[l] = genWENOL(qm2,qm1,qc,qp1,qp2);	HydroGrid[i][j][k-1].VAR##RZ[l] = genWENOR(qm2,qm1,qc,qp1,qp2);\
	}\
} 
	
WENOZ(Fz );
WENOZ(Var);
WENOZ(Az ); 
#endif
		
		
WENOX(Fx );
WENOY(Fy );
WENOX(Var);
WENOY(Var);
WENOX(Ax );
WENOY(Ay );
  	
		

		
void FindEigenValuesX(GRID HydroGrid, double tau)
{

	int i,j,k;
	
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		DECLu4; 
		HydroGrid[i][j][k].Ax[0] = u1/u0;  
	}
}
		
void FindEigenValuesY(GRID HydroGrid, double tau)
{
	int i,j,k;
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
	{ 
		DECLu4; 
		HydroGrid[i][j][k].Ay[0] = u2/u0;  
	}
}

#if !defined LBI
void FindEigenValuesZ(GRID HydroGrid, double tau)
{

	int i,j,k;
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		DECLu4;
		HydroGrid[i][j][k].Az[0] = u3/u0; 

	}
}
#endif
 
	
void Reconstruct(GRID HydroGrid)
{ 
	WENOX_Fx(HydroGrid , SVAR);
	WENOX_Var(HydroGrid , SVAR);	
	WENOX_Ax( HydroGrid , EVAR );
	
	WENOY_Fy(HydroGrid , SVAR);
	WENOY_Var(HydroGrid , SVAR);
	WENOY_Ay( HydroGrid , EVAR );
	
#if !defined LBI	
	WENOZ_Fz(HydroGrid , SVAR);
	WENOZ_Var(HydroGrid , SVAR);
	WENOZ_Az( HydroGrid , EVAR );
#endif

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
		
#if !defined LBI	
		double maxlz=0,maxrz=0;
#endif				
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
				
#if !defined LBI	
			if(fabs(HydroGrid[i][j][k].AzLZ[l]) > maxlz)
				maxlz = fabs(HydroGrid[i][j][k].AzLZ[l]); 
				
			if(fabs(HydroGrid[i][j][k].AzRZ[l]) > maxrz)
				maxrz = fabs(HydroGrid[i][j][k].AzRZ[l]);
#endif		

		}
				
		HydroGrid[i][j][k].AxLXMAX = maxlx;
		HydroGrid[i][j][k].AxRXMAX = maxrx;
		HydroGrid[i][j][k].AyLYMAX = maxly;
		HydroGrid[i][j][k].AyRYMAX = maxry; 
#if !defined LBI	
		HydroGrid[i][j][k].AzLZMAX = maxlz;
		HydroGrid[i][j][k].AzRZMAX = maxrz;
#endif		

	}
}
		 


void CopyPrimaryVariablesToVar(GRID HydroGrid, double tau)
{

	for( int i=0; i<XCM ; i++)
	for( int j=0; j<YCM ; j++)
	for( int k=0; k<ZCM ; k++)
	{
		HydroGrid[i][j][k].Var[0]= HydroGrid[i][j][k].T00;
		HydroGrid[i][j][k].Var[1]= HydroGrid[i][j][k].T10;
		HydroGrid[i][j][k].Var[2]= HydroGrid[i][j][k].T20;
		HydroGrid[i][j][k].Var[3]= HydroGrid[i][j][k].T30;
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
	
	//~ CFL(HydroGrid, tau ,  taustep);
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
	for( i=il; i<ir;  i++)
	for( k=kl; k<kr; k++)
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
				- 0.5*(MAX(fabs( HydroGrid[i][j][k].AzLZMAX), fabs( HydroGrid[i][j][k].AzRZMAX)))*(HydroGrid[i][j][k].VarRZ[l]-HydroGrid[i][j][k].VarLZ[l]); 
		
	
	for( l=0; l<SVAR; l++)
	for( i=il; i< ir; i++)
	for( j=jl; j< jr; j++)
	for( k=kl; k< kr; k++)
		HydroGrid[i][j][k].PartialResult[l] =  -(1.0/ZS)*(HydroGrid[i][j][k].fluxT[l] - HydroGrid[i][j][k-1].fluxT[l]);
		
}
