void hydroExplicit(GRID HydroGrid, double tau, double taustep);
void shastaX(GRID HydroGrid, double taustep);
void shastaY(GRID HydroGrid,  double taustep);
void shastaZ(GRID HydroGrid,  double taustep);


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
				double Vx = HydroGrid[i][j][k].Vx;
				double Vxp1 = HydroGrid[i+1][j][k].Vx;
				double Vxm1 = HydroGrid[i-1][j][k].Vx;
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
				HydroGrid[i][j][k].PartialResult[l] = ( -HydroGrid[i][j][k].Var[l] + HydroGrid[i][j][k].UTD 
											- HydroGrid[i][j][k].A + HydroGrid[i-1][j][k].A );
		}
	}
}

void shastaY(GRID HydroGrid, double taustep)
{

//counters
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

	for(l=0; l< VARN;l++)
	{		
		for( i=il; i< ir; i++)
		for( k=kl; k< kr; k++)
		{
			for( j=1; j<YCM-1;j++)
			{
				double Vy = HydroGrid[i][j][k].Vy;
				double Vyp1 = HydroGrid[i][j+1][k].Vy;
				double Vym1 = HydroGrid[i][j-1][k].Vy;
				
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
				HydroGrid[i][j][k].PartialResult[l] = ( -HydroGrid[i][j][k].Var[l] + HydroGrid[i][j][k].UTD 
											- HydroGrid[i][j][k].A + HydroGrid[i][j-1][k].A );
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
	
	for(l=0; l< VARN;l++)
	{

		for( i=il; i< ir; i++)
		for( j=jl; j< jr; j++)
		{
			for( k=1; k<ZCM-1;k++)
			{
				double Ve = HydroGrid[i][j][k].Ve;
				double Vep1 = HydroGrid[i][j][k+1].Ve;
				double Vem1 = HydroGrid[i][j][k-1].Ve;
				
				
				Qplus  = (0.5 - Ve*lambdaz) / ( 1.0 + (Vep1 - Ve)*lambdaz);
				Qminus = (0.5 + Ve*lambdaz) / ( 1.0 - (Vem1 - Ve)*lambdaz);
				
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
				HydroGrid[i][j][k].PartialResult[l] = ( -HydroGrid[i][j][k].Var[l] + HydroGrid[i][j][k].UTD 
											- HydroGrid[i][j][k].A + HydroGrid[i][j][k-1].A );
		}
	}
}


