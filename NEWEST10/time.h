
void CalcL0(GRID HydroGrid )
{
	int i,j,k,l;


	for(l=0;l<SVAR;l++)
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
		HydroGrid[i][j][k].L0[l] =  (HydroGrid[i][j][k].Result[l] + HydroGrid[i][j][k].Source[l] ); 
}


void CalcL1(GRID HydroGrid )
{
	int i,j,k,l;


	for(l=0;l<SVAR;l++)
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
		HydroGrid[i][j][k].L1[l] =  (HydroGrid[i][j][k].Result[l] + HydroGrid[i][j][k].Source[l] ); 
}

void CalcL2(GRID HydroGrid )
{
	int i,j,k,l;


	for(l=0;l<SVAR;l++)
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
		HydroGrid[i][j][k].L2[l] =  (HydroGrid[i][j][k].Result[l] + HydroGrid[i][j][k].Source[l] ); 
}

void CalcVar0(GRID HydroGrid )
{
	int i,j,k,l;


	for(l=0;l<SVAR;l++)
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
		HydroGrid[i][j][k].Var0[l] =  HydroGrid[i][j][k].Var[l]; 
}

void CalcVar1(GRID HydroGrid,  double ts)
{
	int i,j,k,l;


	for(l=0;l<SVAR;l++)
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
		HydroGrid[i][j][k].Var1[l] =  HydroGrid[i][j][k].Var0[l]+ts*(HydroGrid[i][j][k].L0[l]); 
}

void CalcVar2RK2(GRID HydroGrid,  double ts)
{
	int i,j,k,l;


	for(l=0;l<SVAR;l++)
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
		HydroGrid[i][j][k].Var2[l] =  0.5*HydroGrid[i][j][k].Var0[l]+  0.5*HydroGrid[i][j][k].Var1[l]+ 0.5*ts*HydroGrid[i][j][k].L1[l]; 
}

void CalcVar2RK3(GRID HydroGrid,  double ts)
{
	int i,j,k,l;


	for(l=0;l<SVAR;l++)
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
		HydroGrid[i][j][k].Var2[l] =  0.75*HydroGrid[i][j][k].Var0[l]+  0.25*HydroGrid[i][j][k].Var1[l]+ 0.25*ts*HydroGrid[i][j][k].L1[l]; 
}


void CalcVar3(GRID HydroGrid,  double ts)
{
	int i,j,k,l;


	for(l=0;l<SVAR;l++)
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
		HydroGrid[i][j][k].Var3[l] =  (1.0/3.0)*HydroGrid[i][j][k].Var0[l]+  (2.0/3.0)*HydroGrid[i][j][k].Var2[l]+ (2.0/3.0)*ts*HydroGrid[i][j][k].L2[l]; 
}

void UpdatePrimaryVariablesFromVar1( GRID HydroGrid, double tau )
{
	int i,j,k,l;


	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
	{
		HydroGrid[i][j][k].T00 =  (HydroGrid[i][j][k].Var1[0])/(tau);								
		HydroGrid[i][j][k].T10 =  (HydroGrid[i][j][k].Var1[1])/(tau);								
		HydroGrid[i][j][k].T20 =  (HydroGrid[i][j][k].Var1[2])/(tau);								
		HydroGrid[i][j][k].T30 =  (HydroGrid[i][j][k].Var1[3])/(tau);	

		for(l=0;l<Npi;l++)
			HydroGrid[i][j][k].pi[l] =  ( HydroGrid[i][j][k].Var1[VARN+l]);
			
		HydroGrid[i][j][k].PI = (HydroGrid[i][j][k].Var1[VARN+Npi]);
	}
}

void UpdatePrimaryVariablesFromVar2( GRID HydroGrid, double tau )
{
	int i,j,k,l;


	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
	{
		HydroGrid[i][j][k].T00 =  (HydroGrid[i][j][k].Var2[0])/(tau);								
		HydroGrid[i][j][k].T10 =  (HydroGrid[i][j][k].Var2[1])/(tau);								
		HydroGrid[i][j][k].T20 =  (HydroGrid[i][j][k].Var2[2])/(tau);								
		HydroGrid[i][j][k].T30 =  (HydroGrid[i][j][k].Var2[3])/(tau);	

		for(l=0;l<Npi;l++)
			HydroGrid[i][j][k].pi[l] =  ( HydroGrid[i][j][k].Var2[VARN+l]);
			
		HydroGrid[i][j][k].PI = (HydroGrid[i][j][k].Var2[VARN+Npi]);
	}
}
 

void UpdatePrimaryVariablesFromVar3( GRID HydroGrid, double tau )
{
	int i,j,k,l;


	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
	{
		HydroGrid[i][j][k].T00 =  (HydroGrid[i][j][k].Var3[0])/(tau );								
		HydroGrid[i][j][k].T10 =  (HydroGrid[i][j][k].Var3[1])/(tau );								
		HydroGrid[i][j][k].T20 =  (HydroGrid[i][j][k].Var3[2])/(tau );								
		HydroGrid[i][j][k].T30 =  (HydroGrid[i][j][k].Var3[3])/(tau );	

		for(l=0;l<Npi;l++)
			HydroGrid[i][j][k].pi[l] =  ( HydroGrid[i][j][k].Var3[VARN+l]);
			
		HydroGrid[i][j][k].PI = (HydroGrid[i][j][k].Var3[VARN+Npi]);
	}
}


		
void UpdatePrevU( GRID HydroGrid)
{
	int i,j,k,l;


	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
	{
		HydroGrid[i][j][k].prevu[0] = HydroGrid[i][j][k].u[0];
		HydroGrid[i][j][k].prevu[1] = HydroGrid[i][j][k].u[1];
		HydroGrid[i][j][k].prevu[2] = HydroGrid[i][j][k].u[2];
		HydroGrid[i][j][k].prevu[3] = HydroGrid[i][j][k].u[3];
	}
}


void FirstOrder(GRID HydroGrid, double tau, double ts)
{
	
	CheckPhysics(HydroGrid, 0);		
	
	
	
	CalcSource(HydroGrid, tau, ts); 	
#ifdef KT
	CalcCentreFlux(HydroGrid, tau);
#endif
#ifdef SHAS
	CalcNumVel(HydroGrid, tau);
#endif
	hydroExplicit(HydroGrid, tau, ts); 			//gets update for PV,pi and PI everywhere excluding boundary region

	CalcL0(HydroGrid);
	CalcVar0( HydroGrid);
	CalcVar1( HydroGrid,ts);	
	UpdatePrimaryVariablesFromVar1(HydroGrid,  tau+ts); //update PV's, pi and PI to "tau+ts"	from "tau"	
	UpdatePrevU(HydroGrid);
	RootSearchForEnVelUsingDerivatives(HydroGrid, tau+ts );	//Finds En,P,V's, 4vel, 3vel everywhere excluding boundary region
	
		
	pack(HydroGrid);  //Exchanges En&Vel, pi and PI and updates P,4vel,Tmunu at the cell interfaces
	boundary(HydroGrid);  // Outflowing boundary condition for everything
	
	
	DebugMSG(HydroGrid);
}


void TVDRK2(GRID HydroGrid, double tau, double ts)
{
	
	CheckPhysics(HydroGrid, 0);				
	
	
	
	CalcSource(HydroGrid, tau, ts);	  	
#ifdef KT		
	CalcCentreFlux(HydroGrid, tau);  
#endif
#ifdef SHAS
	CalcNumVel(HydroGrid, tau);
#endif	
	hydroExplicit(HydroGrid, tau, ts/2); //ts/2 is not used in KT and in SHAS it is used 	
	CalcL0(HydroGrid);
	CalcVar0( HydroGrid);	
	CalcVar1( HydroGrid,ts);	
	
	UpdatePrimaryVariablesFromVar1(HydroGrid,  tau+ts);  
	UpdatePrevU(HydroGrid);
	RootSearchForEnVelUsingDerivatives(HydroGrid , tau + ts);    
	
	
	pack(HydroGrid );     
	boundary(HydroGrid);	 	
	CalcSource(HydroGrid, tau+ts, ts);	
#ifdef KT		
	CalcCentreFlux(HydroGrid, tau+ts);  
#endif	
#ifdef SHAS
	CalcNumVel(HydroGrid, tau+ts);
#endif
	hydroExplicit(HydroGrid, tau+ts, ts/2);
	CalcL1(HydroGrid);
	CalcVar2RK2( HydroGrid,ts);	
	
	
	UpdatePrimaryVariablesFromVar2(HydroGrid,  tau+ts);  
    RootSearchForEnVelUsingDerivatives(HydroGrid , tau+ts); 
    	

		
	pack(HydroGrid);
	boundary(HydroGrid);	
	DebugMSG(HydroGrid);
}




void TVDRK3(GRID HydroGrid, double tau, double ts)
{
	CheckPhysics(HydroGrid, 0);		
	
			
	
	CalcSource(HydroGrid, tau, ts);	  	
#ifdef KT		
	CalcCentreFlux(HydroGrid, tau);  
#endif
#ifdef SHAS
	CalcNumVel(HydroGrid, tau);
#endif	
	hydroExplicit(HydroGrid, tau, ts/3.0); //ts/2 is not used in KT and in SHAS it is used 	
	CalcL0(HydroGrid);
	CalcVar0( HydroGrid);	
	CalcVar1( HydroGrid,ts);	
	
	UpdatePrimaryVariablesFromVar1(HydroGrid,  tau+ts);  
	UpdatePrevU(HydroGrid);
	RootSearchForEnVelUsingDerivatives(HydroGrid , tau + ts);    
	
	
	
	
	
	pack(HydroGrid );     
	boundary(HydroGrid);	 	
	CalcSource(HydroGrid, tau+ts, ts);	
#ifdef KT		
	CalcCentreFlux(HydroGrid, tau+ts);  
#endif	
#ifdef SHAS
	CalcNumVel(HydroGrid, tau+ts);
#endif
	hydroExplicit(HydroGrid, tau+ts, ts/3.0);
	CalcL1(HydroGrid);
	CalcVar2RK3( HydroGrid,ts);	
		
	UpdatePrimaryVariablesFromVar2(HydroGrid,  tau+ts/2.0);  
    RootSearchForEnVelUsingDerivatives(HydroGrid , tau + ts/2.0); 
    	




	
	
	pack(HydroGrid );     
	boundary(HydroGrid);	 	
	CalcSource(HydroGrid, tau+ts/2.0, ts/2.0);	
#ifdef KT		
	CalcCentreFlux(HydroGrid, tau+ts/2.0);  
#endif	
#ifdef SHAS
	CalcNumVel(HydroGrid, tau+ts/2.0);
#endif
	hydroExplicit(HydroGrid, tau+ts/2.0, ts/3.0);
	CalcL2(HydroGrid);
	CalcVar3( HydroGrid,ts);	
		
	UpdatePrimaryVariablesFromVar3(HydroGrid,  tau+ts);  
    RootSearchForEnVelUsingDerivatives(HydroGrid , tau+ts); 
    	
		
	pack(HydroGrid);
	boundary(HydroGrid);	
	DebugMSG(HydroGrid);
}


