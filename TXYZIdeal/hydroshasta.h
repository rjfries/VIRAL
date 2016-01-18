void hydroExplicit(GRID HydroGrid, double tau, double taustep);
void shastaX(GRID HydroGrid, double taustep);
void shastaY(GRID HydroGrid,  double taustep);
void shastaZ(GRID HydroGrid,  double taustep);



void CalcNumVel(GRID HydroGrid, double tau)
{
	int i,j,k,l;
	for(l=0;l<SVAR;l++)
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		DECLu4;
		DECLePa;
		
		if(l==0)
		{
			HydroGrid[i][j][k].NVx[l] = (   u1*u0*(e+P) )/HydroGrid[i][j][k].T00;
			HydroGrid[i][j][k].NVy[l] = (   u2*u0*(e+P) )/HydroGrid[i][j][k].T00;
			HydroGrid[i][j][k].NVz[l] = (   u3*u0*(e+P) )/HydroGrid[i][j][k].T00;		
		}
		else
		{
			HydroGrid[i][j][k].NVx[l] = u1/u0; 
			HydroGrid[i][j][k].NVy[l] = u2/u0;
			HydroGrid[i][j][k].NVz[l] = u3/u0;		
		}
	}
}


void CopyPrimaryVariablesToVar(GRID HydroGrid, double tau)
{
	int i,j,k,l;

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
	int i,j,k,l;
	
	for( int l=0; l<SVAR; l++)
	for( int i=0; i<XCM ; i++)
	for( int j=0; j<YCM ; j++)
	for( int k=0; k<ZCM ; k++)
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

#define SHASTA 1

void hydroExplicit(GRID HydroGrid, double tau, double taustep)
{
	
	ClearResultVariable( HydroGrid);
	CopyPrimaryVariablesToVar(HydroGrid,tau);

	shastaX( HydroGrid, taustep);	
	AddPartialResultToFinalResult(HydroGrid);

	shastaY( HydroGrid, taustep);
	AddPartialResultToFinalResult( HydroGrid);
 
#if !defined LBI
	shastaZ( HydroGrid, taustep);
	AddPartialResultToFinalResult( HydroGrid);
#endif

}


void shastaX(GRID HydroGrid, double taustep)
{

//counters
	int i,j,k,l;

	//
	//Variables used during calculation of lower order solution
	//
	double Qplus ,Qminus;
	double lambdax;
	
	lambdax = taustep/XS;
	
	//
	//Variables for doing anitdiffusion and prelimiting
	//
	double signa,temp,term1,term3,min;
	double factor = SHASTA;

	for(l=0; l< SVAR;l++)
	{
		for( j=jl; j<jr ; j++)  
		for( k=kl; k<kr ; k++)
		{
			for( i=1; i<XCM-1;i++)
			{
				double Vx,Vxp1,Vxm1;
				
				Vx = HydroGrid[i][j][k].NVx[l];
				Vxp1 = HydroGrid[i+1][j][k].NVx[l];
				Vxm1 = HydroGrid[i-1][j][k].NVx[l];
								
				Qplus  = (0.5 - Vx*lambdax ) /  ( 1.0 + (Vxp1- Vx)*lambdax);
 				Qminus = (0.5 + Vx*lambdax ) /  ( 1.0 - (Vxm1 - Vx)*lambdax);
				
				HydroGrid[i][j][k].UTD = 0.5 *( (Qplus*Qplus*(HydroGrid[i+1][j][k].Var[l]- HydroGrid[i][j][k].Var[l])) 
					+ (Qminus*Qminus*(HydroGrid[i-1][j][k].Var[l]-HydroGrid[i][j][k].Var[l] )) )
					+ (Qplus+Qminus)*HydroGrid[i][j][k].Var[l];
			}

			for(i=2;i<XCM-3;i++)
			{
				signa =1;
			 	//factor=SHASTAV(HydroGrid[i][j][k].En);
				temp = factor*0.125*(HydroGrid[i+1][j][k].UTD - HydroGrid[i][j][k].UTD);

				if(temp<0)
					signa=-1;

				term1 =  signa*(HydroGrid[i+2][j][k].UTD - HydroGrid[i+1][j][k].UTD);		
				term3 =  signa*(HydroGrid[i][j][k].UTD - HydroGrid[i-1][j][k].UTD);
				min = MIN(term1,fabs(temp),term3);

				HydroGrid[i][j][k].A=signa*MAX(0,min);
			}
			
			for(i=il;i<ir;i++)
				HydroGrid[i][j][k].PartialResult[l] = ( -HydroGrid[i][j][k].Var[l] + HydroGrid[i][j][k].UTD - HydroGrid[i][j][k].A + HydroGrid[i-1][j][k].A )/taustep;
		}
	}
}

void shastaY(GRID HydroGrid, double taustep)
{
 
	int i,j,k,l;

	//
	//Variables used during calculation of lower order solution
	//
	double Qplus ,Qminus;
	double lambday;
	
	lambday = taustep/YS;
	
	//
	//Variables for doing anitdiffusion and prelimiting
	//
	double signa,temp,term1,term3,min;
	double factor = SHASTA;

	for(l=0; l< SVAR;l++)
	{		
		for( i=il; i< ir; i++)
		for( k=kl; k< kr; k++)
		{
			for( j=1; j<YCM-1;j++)
			{
				double Vy,Vyp1,Vym1;				
				 
				Vy = HydroGrid[i][j][k].NVy[l];
				Vyp1 = HydroGrid[i][j+1][k].NVy[l];
				Vym1 = HydroGrid[i][j-1][k].NVy[l];
				 
				
				Qplus  = (0.5 - Vy*lambday) / ( 1.0 + (Vyp1 - Vy)*lambday);
				Qminus = (0.5 + Vy*lambday) / ( 1.0 - (Vym1 - Vy)*lambday);
				
				HydroGrid[i][j][k].UTD  = 0.5 *( (Qplus*Qplus*(HydroGrid[i][j+1][k].Var[l]-HydroGrid[i][j][k].Var[l]))
					+ (Qminus*Qminus*(HydroGrid[i][j-1][k].Var[l]-HydroGrid[i][j][k].Var[l]) ) ) 
					+ (Qplus+Qminus)*HydroGrid[i][j][k].Var[l];
			}

			for(j=2;j<YCM-3;j++)
			{
				signa =1;
			 	//factor=SHASTAV(HydroGrid[i][j][k].En);
				temp = factor*0.125*(HydroGrid[i][j+1][k].UTD - HydroGrid[i][j][k].UTD);
		
				if(temp<0)
					signa=-1;
		
				term1 =  signa*(HydroGrid[i][j+2][k].UTD - HydroGrid[i][j+1][k].UTD);		
				term3 =  signa*(HydroGrid[i][j][k].UTD - HydroGrid[i][j-1][k].UTD);
				min = MIN(term1,fabs(temp),term3);
			
				HydroGrid[i][j][k].A=signa*MAX(0,min);
			}
		
			for(j=jl;j<jr;j++)
				HydroGrid[i][j][k].PartialResult[l] = ( -HydroGrid[i][j][k].Var[l] + HydroGrid[i][j][k].UTD - HydroGrid[i][j][k].A + HydroGrid[i][j-1][k].A )/taustep;
		}
	}
}

void shastaZ(GRID HydroGrid, double taustep)
{

	int i,j,k,l;


	
	double Qplus ,Qminus;
	double lambdaz;
	lambdaz = taustep/ZS;

	
	double signa,temp,term1,term3,min;
	double factor = SHASTA;
	
	for(l=0; l< SVAR;l++)
	{

		for( i=il; i< ir; i++)
		for( j=jl; j< jr; j++)
		{
			for( k=1; k<ZCM-1;k++)
			{
				double Vz = HydroGrid[i][j][k].NVz[l];
				double Vzp1 = HydroGrid[i][j][k+1].NVz[l];
				double Vzm1 = HydroGrid[i][j][k-1].NVz[l];
				
				
				Qplus  = (0.5 - Vz*lambdaz) / ( 1.0 + (Vzp1 - Vz)*lambdaz);
				Qminus = (0.5 + Vz*lambdaz) / ( 1.0 - (Vzm1 - Vz)*lambdaz);
				
				HydroGrid[i][j][k].UTD	= 0.5 *( (Qplus*Qplus*(HydroGrid[i][j][k+1].Var[l]-HydroGrid[i][j][k].Var[l]))
						+ (Qminus*Qminus*(HydroGrid[i][j][k-1].Var[l]-HydroGrid[i][j][k].Var[l]) ) ) 
						+ (Qplus+Qminus)*HydroGrid[i][j][k].Var[l];

			}
			
			for(k=2;k<ZCM-3;k++)
			{
				signa =1;
				//factor=SHASTAV(HydroGrid[i][j][k].En);
				temp = factor*0.125*(HydroGrid[i][j][k+1].UTD - HydroGrid[i][j][k].UTD);
		
				if(temp<0)
					signa=-1;
		
				term1 =  signa*(HydroGrid[i][j][k+2].UTD - HydroGrid[i][j][k+1].UTD);		
				term3 =  signa*(HydroGrid[i][j][k].UTD - HydroGrid[i][j][k-1].UTD);
				min = MIN(term1,fabs(temp),term3);
		
				HydroGrid[i][j][k].A = signa*MAX(0,min);
			}
		
			for(k=kl;k<kr;k++)
				HydroGrid[i][j][k].PartialResult[l] = ( -HydroGrid[i][j][k].Var[l] + HydroGrid[i][j][k].UTD - HydroGrid[i][j][k].A + HydroGrid[i][j][k-1].A )/taustep;
		}
	}
}


