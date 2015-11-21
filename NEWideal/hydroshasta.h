void hydroExplicit(GRID HydroGrid, double tau, double taustep);
void shastaX(GRID HydroGrid, double taustep);
void shastaY(GRID HydroGrid,  double taustep);
void shastaZ(GRID HydroGrid,  double taustep);


#define SHASTA 1
 
void hydroExplicit(GRID HydroGrid, double tau, double taustep)
{
	
	ClearResultVariable( HydroGrid);
	CopyPrimaryVariablesToBufA( HydroGrid);

	shastaX( HydroGrid, taustep);	
	AddPartialResultToFinalResult( HydroGrid);

	shastaY( HydroGrid, taustep);
	AddPartialResultToFinalResult( HydroGrid);
 
	shastaZ( HydroGrid, taustep);
	AddPartialResultToFinalResult( HydroGrid);
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

	for(l=0; l< VARN;l++)
	{
		for( j=jl; j<jr ; j++)
		for( k=kl; k<kr ; k++)
		{
			for( i=1; i<XCM-1;i++)
			{
				Qplus  = (0.5 - HydroGrid[i][j][k].NVx[l]*lambdax ) /  ( 1.0 + (HydroGrid[i+1][j][k].NVx[l] - HydroGrid[i][j][k].NVx[l])*lambdax);
 				Qminus = (0.5 + HydroGrid[i][j][k].NVx[l]*lambdax ) /  ( 1.0 - (HydroGrid[i-1][j][k].NVx[l] - HydroGrid[i][j][k].NVx[l])*lambdax);
				
				HydroGrid[i][j][k].UTD = 0.5 *( (Qplus*Qplus*(HydroGrid[i+1][j][k].BufA[l]- HydroGrid[i][j][k].BufA[l])) 
					+ (Qminus*Qminus*(HydroGrid[i-1][j][k].BufA[l]-HydroGrid[i][j][k].BufA[l] )) )
					+ (Qplus+Qminus)*HydroGrid[i][j][k].BufA[l];
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
				HydroGrid[i][j][k].PartialResult[l] = ( -HydroGrid[i][j][k].BufA[l] + HydroGrid[i][j][k].UTD 
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
				Qplus  = (0.5 - HydroGrid[i][j][k].NVy[l]*lambday) / ( 1.0 + (HydroGrid[i][j+1][k].NVy[l] - HydroGrid[i][j][k].NVy[l])*lambday);
				Qminus = (0.5 + HydroGrid[i][j][k].NVy[l]*lambday) / ( 1.0 - (HydroGrid[i][j-1][k].NVy[l] - HydroGrid[i][j][k].NVy[l])*lambday);
				
				HydroGrid[i][j][k].UTD  = 0.5 *( (Qplus*Qplus*(HydroGrid[i][j+1][k].BufA[l]-HydroGrid[i][j][k].BufA[l]))
					+ (Qminus*Qminus*(HydroGrid[i][j-1][k].BufA[l]-HydroGrid[i][j][k].BufA[l]) ) ) 
					+ (Qplus+Qminus)*HydroGrid[i][j][k].BufA[l];
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
				HydroGrid[i][j][k].PartialResult[l] = ( -HydroGrid[i][j][k].BufA[l] + HydroGrid[i][j][k].UTD 
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
				Qplus  = (0.5 - HydroGrid[i][j][k].NVe[l]*lambdaz) / ( 1.0 + (HydroGrid[i][j][k+1].NVe[l] - HydroGrid[i][j][k].NVe[l])*lambdaz);
				Qminus = (0.5 + HydroGrid[i][j][k].NVe[l]*lambdaz) / ( 1.0 - (HydroGrid[i][j][k-1].NVe[l] - HydroGrid[i][j][k].NVe[l])*lambdaz);
				
				HydroGrid[i][j][k].UTD	= 0.5 *( (Qplus*Qplus*(HydroGrid[i][j][k+1].BufA[l]-HydroGrid[i][j][k].BufA[l]))
						+ (Qminus*Qminus*(HydroGrid[i][j][k-1].BufA[l]-HydroGrid[i][j][k].BufA[l]) ) ) 
						+ (Qplus+Qminus)*HydroGrid[i][j][k].BufA[l];

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
				HydroGrid[i][j][k].PartialResult[l] = ( -HydroGrid[i][j][k].BufA[l] + HydroGrid[i][j][k].UTD 
											- HydroGrid[i][j][k].A + HydroGrid[i][j][k-1].A );
		}
	}
}


