void hydroExplicit(GRID HydroGrid, double tau, double taustep);
void fvX(GRID HydroGrid, double taustep, double tau);
void fvY(GRID HydroGrid,  double taustep, double tau);
void fvZ(GRID HydroGrid,  double taustep, double tau);


#define WENOP 2
#define WENOEPS 1E-6


int method;
int type;

void ClearLR(GRID HydroGrid){
	
	int i,j,k;
	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		HydroGrid[i][j][k].left = 0; 
		HydroGrid[i][j][k].right = 0; 
		HydroGrid[i][j][k].VAR = 0; 
	}
}

void hydroExplicit(GRID HydroGrid, double tau, double taustep)
{
	
	ClearResultVariable( HydroGrid);
	CopyPrimaryVariablesToBufA( HydroGrid);

	fvX( HydroGrid, taustep,  tau);	
	AddPartialResultToFinalResult( HydroGrid);

	fvY( HydroGrid, taustep,  tau);		
	AddPartialResultToFinalResult( HydroGrid);

	fvZ( HydroGrid, taustep,  tau);
	AddPartialResultToFinalResult( HydroGrid);
}


inline void ClearFluxVelVar(GRID HydroGrid)
{
	int  i, j,k;
	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		HydroGrid[i][j][k].fluxL=0;
		HydroGrid[i][j][k].fluxR=0;
		HydroGrid[i][j][k].fluxT=0;
		HydroGrid[i][j][k].velL=0;
		HydroGrid[i][j][k].velR=0;
		HydroGrid[i][j][k].varL=0;
		HydroGrid[i][j][k].varR=0;
		HydroGrid[i][j][k].PartialResult[0]=0;
		HydroGrid[i][j][k].PartialResult[1]=0;
		HydroGrid[i][j][k].PartialResult[2]=0;
		HydroGrid[i][j][k].PartialResult[3]=0;
	}
}





/*
  **
  **
  **Reconstruction procedures
  **
  **
  */
  
  
void  WenoShuX(GRID HydroGrid)
{
	double w[3], q[3], d[3], alpha[3];
	double wt[3], qt[3], dt[3], alphat[3];

	double   beta[4];

	double p,eps=WENOEPS,sum;
	int i,j,k,c;
	
	
	d[0]= dt[2]=0.3;
	d[1]= dt[1]=0.6 ;
	d[2]= dt[0]=0.1 ;
	
/*	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)*/
	
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		for(i=2;i<XCM-2;i++)
		{
			double qm2,qm1,qc,qp1,qp2;
			
			qm2=HydroGrid[i-2][j][k].VAR;
			qm1=HydroGrid[i-1][j][k].VAR;
			qc=HydroGrid[i][j][k].VAR;
			qp1=HydroGrid[i+1][j][k].VAR;
			qp2=HydroGrid[i+2][j][k].VAR;

			beta[0]=(13./12.)*pow(qc - 2*qp1 + qp2 , 2) + 0.25*pow(3*qc - 4*qp1 + qp2 , 2);
			beta[1]=(13./12.)*pow(qm1 - 2*qc + qp1 , 2) + 0.25*pow(qm1 - qp1 , 2);
			beta[2]=(13./12.)*pow(qm2 - 2*qm1 + qc , 2) + 0.25*pow(qm2 - 4*qm1 + 3*qc , 2);


			
			q[0]= (2*qc  + 5*qp1 - 1*qp2)/6.0;
			q[1]= (-1*qm1 + 5*qc + 2*qp1)/6.0;
			q[2]= (2*qm2 - 7*qm1 + 11*qc)/6.0;

			for(c=0;c<3;c++)
				alpha[c]= d[c]/pow(eps + beta[c],WENOP);

			sum=0;
			for(c=0;c<3;c++)
				sum += alpha[c];

			for(c=0;c<3;c++)
				w[c]= alpha[c]/sum;
			
			HydroGrid[i][j][k].left=0;
			for(c=0;c<3;c++)
				HydroGrid[i][j][k].left += w[c]*q[c];


			
			qt[0]= (11*qc - 7*qp1 + 2*qp2)/6.0;
			qt[1]= (2*qm1 + 5*qc - 1*qp1)/6.0;
			qt[2]= (-1*qm2 + 5*qm1 + 2*qc)/6.0;

			for(c=0;c<3;c++)
				alphat[c]= dt[c]/pow(eps + beta[c],WENOP);
			
			sum=0;
			for(c=0;c<3;c++)
				sum += alphat[c];
			
			for(c=0;c<3;c++)
				wt[c]= alphat[c]/sum;
			
			HydroGrid[i-1][j][k].right=0;
			for(c=0;c<3;c++)
				HydroGrid[i-1][j][k].right += wt[c]*qt[c];
			
		}
	}
}

void  WenoShuY( GRID HydroGrid )
{
	double w[3], q[3], d[3], alpha[3];
	double wt[3], qt[3], dt[3], alphat[3];

	double   beta[4];

	double p,eps=WENOEPS,sum;
	int i,j,k,c;
	
	
	d[0]= dt[2]=0.3;
	d[1]= dt[1]=0.6;
	d[2]= dt[0]=0.1;
	
	
/*	for(i=il;i<ir;i++)
	for(k=kl;k<kr;k++)*/
	for(i=0;i<XCM;i++)
	for(k=0;k<ZCM;k++)
	{
		for(j=2;j<YCM-2;j++)
		{
			double qm2,qm1,qc,qp1,qp2;
			
			qm2=HydroGrid[i][j-2][k].VAR;
			qm1=HydroGrid[i][j-1][k].VAR;
			qc =HydroGrid[i][j][k].VAR;
			qp1=HydroGrid[i][j+1][k].VAR;
			qp2=HydroGrid[i][j+2][k].VAR;

			beta[0]=(13./12.)*pow(qc - 2*qp1 + qp2 , 2) + 0.25*pow(3*qc - 4*qp1 + qp2 , 2);
			beta[1]=(13./12.)*pow(qm1 - 2*qc + qp1 , 2) + 0.25*pow(qm1 - qp1 , 2);
			beta[2]=(13./12.)*pow(qm2 - 2*qm1 + qc , 2) + 0.25*pow(qm2 - 4*qm1 + 3*qc , 2);


			
			q[0]= (2*qc + 5*qp1 - 1*qp2)/6.0;
			q[1]= (-1*qm1 + 5*qc + 2*qp1)/6.0;
			q[2]= (2*qm2 - 7*qm1 + 11*qc)/6.0;

			for(c=0;c<3;c++)
				alpha[c]= d[c]/pow(eps + beta[c],WENOP);

			sum=0;
			for(c=0;c<3;c++)
				sum += alpha[c];

			for(c=0;c<3;c++)
				w[c]= alpha[c]/sum;
			
			HydroGrid[i][j][k].left=0;
			for(c=0;c<3;c++)
				HydroGrid[i][j][k].left+= w[c]*q[c];


			
			qt[0]= (11*qc   - 7 *qp1 + 2*qp2)/6.0;
			qt[1]= (2*qm1 + 5*qc - 1*qp1)/6.0;
			qt[2]= (-1*qm2 + 5*qm1 + 2*qc)/6.0;

			for(c=0;c<3;c++)
				alphat[c]= dt[c]/pow(eps + beta[c],WENOP);
			
			sum=0;
			for(c=0;c<3;c++)
				sum += alphat[c];
			
			for(c=0;c<3;c++)
				wt[c]= alphat[c]/sum;
			
			HydroGrid[i][j-1][k].right=0;
			for(c=0;c<3;c++)
				HydroGrid[i][j-1][k].right += wt[c]*qt[c];
			
		}
	}
}


void  WenoShuZ(GRID HydroGrid )
{
	double w[3], q[3], d[3], alpha[3];
	double wt[3], qt[3], dt[3], alphat[3];

	double   beta[4];

	double p,eps=WENOEPS,sum;
	int i,j,k,c;
/*
	ClearLR(HydroGrid);
	return;
	*/
	
	d[0]= dt[2]=0.3;
	d[1]= dt[1]=0.6 ;
	d[2]= dt[0]=0.1 ;
	/*
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)*/
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	{

		for(k=2;k<ZCM-2;k++)
		{
			double qm2,qm1,qc,qp1,qp2;
			
			qm2=HydroGrid[i][j][k-2].VAR;
			qm1=HydroGrid[i][j][k-1].VAR;
			qc=HydroGrid[i][j][k].VAR;
			qp1=HydroGrid[i][j][k+1].VAR;
			qp2=HydroGrid[i][j][k+2].VAR;

			beta[0]=(13./12.)*pow(qc - 2*qp1 + qp2 , 2) + 0.25*pow(3*qc - 4*qp1 + qp2 , 2);
			beta[1]=(13./12.)*pow(qm1 - 2*qc + qp1 , 2) + 0.25*pow(qm1 - qp1 , 2);
			beta[2]=(13./12.)*pow(qm2 - 2*qm1 + qc , 2) + 0.25*pow(qm2 - 4*qm1 + 3*qc , 2);


			
			q[0]= (2*qc + 5*qp1 - 1*qp2)/6.0;
			q[1]= (-1*qm1 + 5*qc + 2*qp1)/6.0;
			q[2]= (2*qm2 - 7*qm1 + 11*qc)/6.0;

			for(c=0;c<3;c++)
				alpha[c]= d[c]/pow(eps + beta[c],WENOP);

			sum=0;
			for(c=0;c<3;c++)
				sum += alpha[c];

			for(c=0;c<3;c++)
				w[c]= alpha[c]/sum;
			
			HydroGrid[i][j][k].left=0;
			for(c=0;c<3;c++)
				HydroGrid[i][j][k].left+= w[c]*q[c];


			
			qt[0]= (11*qc   - 7 *qp1 + 2*qp2)/6.0;
			qt[1]= (2*qm1 + 5*qc - 1*qp1)/6.0;
			qt[2]= (-1*qm2 + 5*qm1 + 2*qc)/6.0;

			for(c=0;c<3;c++)
				alphat[c]= dt[c]/pow(eps + beta[c],WENOP);
			
			sum=0;
			for(c=0;c<3;c++)
				sum += alphat[c];
			
			for(c=0;c<3;c++)
				wt[c]= alphat[c]/sum;
			
			HydroGrid[i][j][k-1].right=0;
			for(c=0;c<3;c++)
				HydroGrid[i][j][k-1].right += wt[c]*qt[c];
			
		}
	}
}



inline void  ReconstructX(GRID HydroGrid, int mode, int var )
{

	int i,j,k;
	
	switch(mode)
	{
		case 0:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
					HydroGrid[i][j][k].VAR = HydroGrid[i][j][k].Vx; 
					
			break;
			
		case 1:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
					HydroGrid[i][j][k].VAR = HydroGrid[i][j][k].BufA[var];
					
			break;		
			
		case 2:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
					HydroGrid[i][j][k].VAR = HydroGrid[i][j][k].Fx[var];
					
			break;
						
		case 3:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
					HydroGrid[i][j][k].VAR = 1.0/sqrt(1- HydroGrid[i][j][k].Vx*HydroGrid[i][j][k].Vx- HydroGrid[i][j][k].Vy*HydroGrid[i][j][k].Vy- tau*tau* HydroGrid[i][j][k].Ve*HydroGrid[i][j][k].Ve );
					
			break;
	}

	WenoShuX( HydroGrid); 
	
	switch(mode)
	{
		case 0:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
				{
					HydroGrid[i][j][k].velL  = HydroGrid[i][j][k].left; 
					HydroGrid[i][j][k].velR  = HydroGrid[i][j][k].right; 
				}
					
			break;
			
		case 1:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
				{
					HydroGrid[i][j][k].varL= HydroGrid[i][j][k].left; 
					HydroGrid[i][j][k].varR = HydroGrid[i][j][k].right; 
				}
			break;	
				
		case 2:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
				{
					HydroGrid[i][j][k].fluxL= HydroGrid[i][j][k].left; 
					HydroGrid[i][j][k].fluxR = HydroGrid[i][j][k].right; 
				}
			break;		
				
		case 3:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
				{
					HydroGrid[i][j][k].u0L = HydroGrid[i][j][k].left; 
					HydroGrid[i][j][k].u0R = HydroGrid[i][j][k].right; 
				}
			break;
	}
	ClearLR(HydroGrid);
}


inline void  ReconstructY(GRID HydroGrid, int mode, int var   )
{

	int i,j,k;
	
	switch(mode)
	{
		case 0:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
					HydroGrid[i][j][k].VAR = HydroGrid[i][j][k].Vy ; 
					
			break;
			
		case 1:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
					HydroGrid[i][j][k].VAR = HydroGrid[i][j][k].BufA[var] ;
					
			break;
			
						
		case 2:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
					HydroGrid[i][j][k].VAR = HydroGrid[i][j][k].Fy[var] ;
			
			break;
						
		case 3:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
					HydroGrid[i][j][k].VAR = 1.0/sqrt(1- HydroGrid[i][j][k].Vx*HydroGrid[i][j][k].Vx- HydroGrid[i][j][k].Vy*HydroGrid[i][j][k].Vy- tau*tau*HydroGrid[i][j][k].Ve*HydroGrid[i][j][k].Ve );
					
			break;
	
	}

	WenoShuY( HydroGrid);		
	
	switch(mode)
	{
		case 0:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
				{
					HydroGrid[i][j][k].velL = HydroGrid[i][j][k].left ; 
					HydroGrid[i][j][k].velR = HydroGrid[i][j][k].right ; 
				}
					
			break;
			
		case 1:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
				{
					HydroGrid[i][j][k].varL = HydroGrid[i][j][k].left; 
					HydroGrid[i][j][k].varR = HydroGrid[i][j][k].right; 
				}
			break;		
				
		case 2:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
				{
					HydroGrid[i][j][k].fluxL = HydroGrid[i][j][k].left; 
					HydroGrid[i][j][k].fluxR = HydroGrid[i][j][k].right; 
				}
			break;	
				
		case 3:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
				{
					HydroGrid[i][j][k].u0L = HydroGrid[i][j][k].left; 
					HydroGrid[i][j][k].u0R = HydroGrid[i][j][k].right; 
				}
			break;
	}
	ClearLR(HydroGrid);
	
}

inline void  ReconstructZ(GRID HydroGrid, int mode, int var   )
{

	int i,j,k;
	
	switch(mode)
	{
		case 0:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
					HydroGrid[i][j][k].VAR  = HydroGrid[i][j][k].Ve ; 
					
			break;
			
		case 1:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
					HydroGrid[i][j][k].VAR  = HydroGrid[i][j][k].BufA[var] ;					
			break;			
			
		case 2:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
					HydroGrid[i][j][k].VAR  = HydroGrid[i][j][k].Fe[var] ;
			break;
										
		case 3:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
					HydroGrid[i][j][k].VAR = 1.0/sqrt(1- HydroGrid[i][j][k].Vx*HydroGrid[i][j][k].Vx- HydroGrid[i][j][k].Vy*HydroGrid[i][j][k].Vy-tau*tau* HydroGrid[i][j][k].Ve*HydroGrid[i][j][k].Ve );
			break;
	}

	WenoShuZ( HydroGrid);
	
	switch(mode)
	{
		case 0:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
				{
					HydroGrid[i][j][k].velL = HydroGrid[i][j][k].left; 
					HydroGrid[i][j][k].velR = HydroGrid[i][j][k].right; 
				}
					
			break;
			
		case 1:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
				{
					HydroGrid[i][j][k].varL = HydroGrid[i][j][k].left; 
					HydroGrid[i][j][k].varR = HydroGrid[i][j][k].right; 
				}
			break;
				
		case 2:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
				{
					HydroGrid[i][j][k].fluxL = HydroGrid[i][j][k].left; 
					HydroGrid[i][j][k].fluxR = HydroGrid[i][j][k].right; 
				}
			break;	
				
		case 3:
				for(i=0;i<XCM;i++)
				for(j=0;j<YCM;j++)
				for(k=0;k<ZCM;k++)
				{
					HydroGrid[i][j][k].u0L = HydroGrid[i][j][k].left; 
					HydroGrid[i][j][k].u0R = HydroGrid[i][j][k].right; 
				}
			break;
	}
	ClearLR(HydroGrid);
}



//The actual hydro code


void fvX(GRID HydroGrid, double taustep, double tau)
{

	int i,j,k,l;

	ClearFluxVelVar(HydroGrid);
	
	/*
	 * 
	 * Reconstruct  Mode
	 *    0 -> for vel
	 *    1 -> for var  
	 *    2 -> for centerd flux
	 *    3 -> for u0
	 *    Var --> which variable from BufA
	 */ 
	
	ReconstructX(HydroGrid,0,0);	
	ReconstructX(HydroGrid,3,0);

	for(l=0; l< VARN;l++)
	{
		
		ReconstructX(HydroGrid,1,l);	
		ReconstructX(HydroGrid,2,l);	
		
		//Slope limiter on the variable	
		
		//~ for( j=jl; j<jr ; j++)
		//~ for( k=kl; k<kr ; k++) 
		for( j=0; j<YCM ; j++)
		for( k=0; k<ZCM ; k++)
		for( i=2; i<XCM-2 ; i++)
		{
			 
			 
			double f = 2;
			HydroGrid[i][j][k].fluxT = 0.5*( HydroGrid[i][j][k].fluxL + HydroGrid[i][j][k].fluxR )
				- 0.5 *(MAX(fabs(f*HydroGrid[i][j][k].velL ),fabs(f*HydroGrid[i][j][k].velR ))) * (HydroGrid[i][j][k].varR - HydroGrid[i][j][k].varL ); 
								
			//~ double a = 1.0/3.0;
			//~ 
			//~ double u0L = HydroGrid[i][j][k].u0L;
			//~ double ukL = u0L*HydroGrid[i][j][k].velL;
			//~ double u0R = HydroGrid[i][j][k].u0R;
			//~ double ukR = u0R*HydroGrid[i][j][k].velR;
			//~ 
			//~ double lambdaL = ( fabs( (1-a)*(u0L*ukL)) + sqrt( a*( u0L*u0L - ukL*ukL - a*(u0L*u0L-ukL*ukL-1) ) ) ) / (u0L*u0L - a*(u0L*u0L-1) )   ;
			//~ double lambdaR = ( fabs( (1-a)*(u0R*ukR)) + sqrt( a*( u0R*u0R - ukR*ukR - a*(u0R*u0R-ukR*ukR-1) ) ) ) / (u0R*u0R - a*(u0R*u0R-1) )   ;
			//~ 
			//~ HydroGrid[i][j][k].fluxT = 0.5*( HydroGrid[i][j][k].fluxL + HydroGrid[i][j][k].fluxR )
				//~ - 0.5 * ( MAX(fabs(lambdaL),fabs(lambdaR) ) ) * ( HydroGrid[i][j][k].varR - HydroGrid[i][j][k].varL );
			//~ if( 0.5*(HydroGrid[i][j][k].Vx + HydroGrid[i+1][j][k].Vx) > 0)
				//~ HydroGrid[i][j][k].fluxT =  HydroGrid[i][j][k].fluxL;
			//~ else
				//~ HydroGrid[i][j][k].fluxT =  HydroGrid[i][j][k].fluxR; 
		}
		
			
		for( j=0; j<YCM   ; j++)
		for( k=0; k<ZCM   ; k++)
		for( i=il; i<ir   ; i++)
		{
			double diff = HydroGrid[i][j][k].fluxT - HydroGrid[i-1][j][k].fluxT;			
		
			HydroGrid[i][j][k].PartialResult[l] =  -(taustep/XS)*(diff);
		}		
	}
}

void fvY(GRID HydroGrid, double taustep, double tau)
{
	int i,j,k,l;	
	
	
	ClearFluxVelVar(HydroGrid);
	
	/*
	 * 
	 * Reconstruct  Mode
	 *    0 -> for vel
	 *    1 -> for var  
	 *    2 -> for centerd flux
	 *    3 -> for u0
	 * Var --> which variable from BufA
	 * 
	 */ 
	
	
	ReconstructY(HydroGrid,0,0);	
	ReconstructY(HydroGrid,3,0);			
	
	
	for(l=0; l<VARN; l++)
	{
		
		ReconstructY(HydroGrid,1,l);	
		ReconstructY(HydroGrid,2,l);
/*
		for( i=il;   i< ir;  i++)
		for( k=kl;   k< kr;  k++)*/
	
		
	
		for( i=0;  i<XCM ; i++)
		for( k=0;  k<ZCM ; k++)
		for( j=2 ; j<YCM-2; j++)
		{
			
			
			
			double f = 2;
			HydroGrid[i][j][k].fluxT = 0.5*( HydroGrid[i][j][k].fluxL + HydroGrid[i][j][k].fluxR ) 
				- 0.5 *(MAX(fabs(f*HydroGrid[i][j][k].velL ),fabs(f*HydroGrid[i][j][k].velR ))) * (HydroGrid[i][j][k].varR - HydroGrid[i][j][k].varL ); 		
			
				
			//~ if(i==7  && l==1)
				//~ cout<<std::scientific<<HydroGrid[i][j][k].fluxT<<endl;
					
			/*	
			
			double a = DPDE(HydroGrid[i][j][k].En, HydroGrid[i][j][k].r);
			
			double u0L = HydroGrid[i][j][k].u0L;
			double ukL = u0L*HydroGrid[i][j][k].velL;
			double u0R = HydroGrid[i][j][k].u0R;
			double ukR = u0R*HydroGrid[i][j][k].velR;
			
			double lambdaL = ( fabs( (1-a)*(u0L*ukL)) + sqrt( a*( u0L*u0L - ukL*ukL - a*(u0L*u0L-ukL*ukL-1) ) ) ) / (u0L*u0L - a*(u0L*u0L-1) )   ;
			double lambdaR = ( fabs( (1-a)*(u0R*ukR)) + sqrt( a*( u0R*u0R - ukR*ukR - a*(u0R*u0R-ukR*ukR-1) ) ) ) / (u0R*u0R - a*(u0R*u0R-1) )   ;
			
			HydroGrid[i][j][k].fluxT = 0.5*( HydroGrid[i][j][k].fluxL + HydroGrid[i][j][k].fluxR )
				- 0.5 * ( MAX(fabs(lambdaL),fabs(lambdaR) ) ) * ( HydroGrid[i][j][k].varR - HydroGrid[i][j][k].varL );
			
			if( 0.5*(HydroGrid[i][j][k].Vy + HydroGrid[i][j+1][k].Vy) > 0)
				HydroGrid[i][j][k].fluxT =  HydroGrid[i][j][k].fluxL;
			else
				HydroGrid[i][j][k].fluxT =  HydroGrid[i][j][k].fluxR;
			*/	
		}


	


		for( i=0;  i<XCM ; i++)
		for( k=0;  k<ZCM ; k++)
		for( j=1 ; j<YCM ; j++)
			HydroGrid[i][j][k].PartialResult[l] =  -(taustep/YS)*(HydroGrid[i][j][k].fluxT  - HydroGrid[i][j-1][k].fluxT ) ;
			
		
		
	}
}

		
void fvZ(GRID HydroGrid, double taustep, double tau)
{
	int i,j,k,l;	


	ClearFluxVelVar(HydroGrid);
	
	/*
	 * 
	 * Reconstruct  Mode
	 *    0 -> for vel
	 *    1 -> for var  
	 *    2 -> for centerd flux
	 *    3 -> for u0
	 * Var --> which variable from BufA
	 * 
	 */ 
	
	
	ReconstructZ(HydroGrid,0,0);
	ReconstructZ(HydroGrid,3,0);
	
	for(l=0; l< VARN;l++)
	{

		ReconstructZ(HydroGrid,1,l);	
		ReconstructZ(HydroGrid,2,l);	

		
/*
		for( i=il; i< ir; i++)
		for( j=jl; j< jr; j++)
		for( k=2; k< ZCM-2; k++)*/
		for( i=0; i< XCM; i++)
		for( j=0; j< YCM; j++)
		for( k=1; k< ZCM-1; k++)

		{
			double f=2;
			HydroGrid[i][j][k].fluxT  = 0.5*( HydroGrid[i][j][k].fluxL + HydroGrid[i][j][k].fluxR )
				- 0.5 *(MAX(fabs(f*HydroGrid[i][j][k].velL ),fabs(f*HydroGrid[i][j][k].velR ))) * (HydroGrid[i][j][k].varR - HydroGrid[i][j][k].varL);
		
		/*	double a = 1.0/3.0;
			
			double u0L = HydroGrid[i][j][k].u0L;
			double ukL = u0L*HydroGrid[i][j][k].velL;
			double u0R = HydroGrid[i][j][k].u0R;
			double ukR = u0R*HydroGrid[i][j][k].velR;
			
			double lambdaL = (1.0/tau)*( fabs( (1-a)*(u0L*ukL)) + sqrt( a*( u0L*u0L - ukL*ukL - a*(u0L*u0L-ukL*ukL-1) ) ) ) / (u0L*u0L - a*(u0L*u0L-1) )   ;
			double lambdaR = (1.0/tau)*( fabs( (1-a)*(u0R*ukR)) + sqrt( a*( u0R*u0R - ukR*ukR - a*(u0R*u0R-ukR*ukR-1) ) ) ) / (u0R*u0R - a*(u0R*u0R-1) )   ;
			
			
			HydroGrid[i][j][k].fluxT = 0.5*( HydroGrid[i][j][k].fluxL + HydroGrid[i][j][k].fluxR )
				- 0.5 * ( MAX(fabs(lambdaL),fabs(lambdaR) ) ) * ( HydroGrid[i][j][k].varR - HydroGrid[i][j][k].varL );
				
				*/							
			/*if( 0.5*(HydroGrid[i][j][k].Ve + HydroGrid[i][j][k+1].Ve) > 0)
				HydroGrid[i][j][k].fluxT =  HydroGrid[i][j][k].fluxL;
			else
				HydroGrid[i][j][k].fluxT =  HydroGrid[i][j][k].fluxR;
			*/
		}


/*
		for( i=il; i< ir; i++)
		for( j=jl; j< jr; j++)
		for( k=kl; k< kr; k++)*/
	
		for( i=0; i< XCM; i++)
		for( j=0; j< YCM; j++)
		for( k=kl; k< kr; k++)
		{
			HydroGrid[i][j][k].PartialResult[l]  =  -(taustep/ZS)*(HydroGrid[i][j][k].fluxT  - HydroGrid[i][j][k-1].fluxT);
			
			//~ if( fabs(HydroGrid[i][j][k].PartialResult[l]) > 1e-10 && k==k0)
				//~ cout<<std::scientific<<"jehrvfbkcnlskdncs   "<<l<< "     "<<k<< "     "<<fabs(HydroGrid[i][j][k].PartialResult[l])<<endl;
		
		}
	}
}




