void StoreInHeap();
void AllocatePack();

PCKX leftX ,rightX;
PCKY leftY ,rightY;
PCKZ leftZ ,rightZ;

void AllocatePack()
{

	int i,j,k,l;

	leftX = new double[PACKVAR][BORDER][YCMA][ZCMA];  //certified contigous .. first element at [0][0][0] and last at [max][max][max]
	rightX = new double[PACKVAR][BORDER][YCMA][ZCMA];
	leftY = new double[PACKVAR][BORDER][XCMA][ZCMA];
	rightY = new double[PACKVAR][BORDER][XCMA][ZCMA];
	leftZ = new double[PACKVAR][BORDER][XCMA][YCMA];
	rightZ = new double[PACKVAR][BORDER][XCMA][YCMA];

	for(l=0;l<PACKVAR;l++)
	for(i=0;i<BORDER;i++)
	for(j=0;j<YCMA;j++)
	for(k=0;k<ZCMA;k++)
	{
		(leftX)[l][i][j][k] = 0; 
		(rightX)[l][i][j][k] = 0; 
	}

	for(l=0;l<PACKVAR;l++)
	for(i=0;i<BORDER;i++)
	for(j=0;j<XCMA;j++)
	for(k=0;k<ZCMA;k++)
	{
		(leftY)[l][i][j][k] = 0; 
		(rightY)[l][i][j][k] = 0; 
	}


	
	for(l=0;l<PACKVAR;l++)
	for(i=0;i<BORDER;i++)
	for(j=0;j<XCMA;j++)
	for(k=0;k<YCMA;k++)
	{
		(leftZ)[l][i][j][k] = 0; 
		(rightZ)[l][i][j][k] = 0; 
	}		

}



void StoreInHeap()
{

	HydroGrid = new cell  [XCM][YCM][ZCM];

	for(int i=0;i<XCM;i++)
	for(int j=0;j<YCM;j++)
	for(int k=0;k<ZCM;k++)
	{
		HydroGrid[i][j][k].En=0;
		HydroGrid[i][j][k].Vx=0;
		HydroGrid[i][j][k].Vy=0;
		HydroGrid[i][j][k].Ve=0;
		HydroGrid[i][j][k].P=0;

		
		HydroGrid[i][j][k].T00=0;
		HydroGrid[i][j][k].T10=0;
		HydroGrid[i][j][k].T20=0;
		HydroGrid[i][j][k].T30=0;

		for(int l=0;l<VARN;l++)
			HydroGrid[i][j][k].Source[l]=0;




		
		for(int l=0;l<VARN;l++)
		{
			HydroGrid[i][j][k].BufA[l]=0;
			HydroGrid[i][j][k].Result[l]=0;
			HydroGrid[i][j][k].PartialResult[l]=0;
		}


#if defined SHAS || defined ZAL
	HydroGrid[i][j][k].UTD=0;
	HydroGrid[i][j][k].Ubar=0;
	HydroGrid[i][j][k].A=0;
	HydroGrid[i][j][k].Ac=0;
#endif

#if defined KT
	HydroGrid[i][j][k].fluxL=0;
	HydroGrid[i][j][k].fluxR=0;
	
	HydroGrid[i][j][k].fluxT=0;

	HydroGrid[i][j][k].velL=0;
	HydroGrid[i][j][k].velR=0;
	
	HydroGrid[i][j][k].varL=0;
	HydroGrid[i][j][k].varR=0;
	
	HydroGrid[i][j][k].left=0;
	HydroGrid[i][j][k].right=0;
	HydroGrid[i][j][k].VAR=0;
#endif

	}
		
	
}


