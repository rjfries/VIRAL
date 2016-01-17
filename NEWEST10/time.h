
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


void Rescalepi(GRID HydroGrid, double tau)
{
	
	int i,j,k,l;
 
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		DECLp10u4;
		DECLePPIa;
		double T[4][4];
		double pi[4][4] = {{p1,p2,p3,p4},{p2,p5,p6,p7},{p3,p6,p8,p9},{p4,p7,p9,p10}};
		T[0][0] = -P + ((e + P)*u0*u0);
		T[0][1] = ((e + P)*u0*u1);
		T[0][2] = ((e + P)*u0*u2);
		T[0][3] = ((e + P)*u0*u3);
		
		T[1][0] = T[0][1];
		T[1][1] = P + ((e + P)*u1*u1);
		T[1][2] = ((e + P)*u1*u2);
		T[1][3] = ((e + P)*u1*u3);
		
		T[2][0] = T[0][2];;
		T[2][1] = ((e + P)*u2*u1);
		T[2][2] = P + ((e + P)*u2*u2);
		T[2][3] = ((e + P)*u2*u3); 
	
		T[3][0] = T[0][3]; 
		T[3][1] = T[1][3]; 
		T[3][2] = T[2][3]; 
		T[3][3] = ((e + P)*u3*u3 + P/(tau*tau) ); 
	
		double maxpi=0,maxT=0;
		
        for (int m = 0; m < 4; m++)
        for (int n= 0; n < 4; n++)
			if(fabs(pi[m][n])>maxpi) maxpi=fabs(pi[m][n]);
			

        for (int m = 0; m < 4; m++)
        for (int n= 0; n < 4; n++)
			if(fabs(T[m][n])>maxT) maxT=fabs(T[m][n]);
			
			maxT=fmax( maxT , P + (e+P)*(u1*u1+u2*u2+u3*u3) );

		if (maxT / maxpi < 1.0) 
		{ 
			HydroGrid[i][j][k].pi[0] = 0.1*p1*maxT/maxpi;
			HydroGrid[i][j][k].pi[1] = 0.1*p2*maxT/maxpi;
			HydroGrid[i][j][k].pi[2] = 0.1*p3*maxT/maxpi;
			HydroGrid[i][j][k].pi[3] = 0.1*p4*maxT/maxpi;
			HydroGrid[i][j][k].pi[4] = 0.1*p5*maxT/maxpi;
			HydroGrid[i][j][k].pi[5] = 0.1*p6*maxT/maxpi;
			HydroGrid[i][j][k].pi[6] = 0.1*p7*maxT/maxpi;
			HydroGrid[i][j][k].pi[7] = 0.1*p8*maxT/maxpi;
			HydroGrid[i][j][k].pi[8] = 0.1*p9*maxT/maxpi;
			HydroGrid[i][j][k].pi[9] = 0.1*p10*maxT/maxpi;
		}
		 
		if (fabs(PI) > P) 
		{
			if (PI != 0.) HydroGrid[i][j][k].PI=0.1 * PI/fabs(PI)*P;
		}
	 }

}








void FixOne(GRID HydroGrid, double tau)
{
	
	int i,j,k,l;
  
	if(BMax>0.98)
	{
		for(i=il;i<ir;i++)
		for(j=jl;j<jr;j++)
		for(k=kl;k<kr;k++)
		{			 
			HydroGrid[i][j][k].Ve =  0; 
		
			for(l=0;l<Npi;l++)
				HydroGrid[i][j][k].pi[l]=0;
			
			HydroGrid[i][j][k].PI=0;
			
		 }
	 }

}


void FixTwo(GRID HydroGrid, double tau)
{
	
	int i,j,k,l;
  
	if(BMax>0.98)
	{
		for(i=il;i<ir;i++)
		for(j=jl;j<jr;j++)
		for(k=kl;k<kr;k++)
		{			 
			HydroGrid[i][j][k].Ve =  0;  
			
		 }
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
	
	
#ifdef FIX
	//~ Rescalepi(HydroGrid,tau);	
	//~ FixOne(HydroGrid,tau);	
	FixTwo(HydroGrid,tau);	
#endif

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


#ifdef FIX
	//~ Rescalepi(HydroGrid,tau);	
	//~ FixOne(HydroGrid,tau);	
	FixTwo(HydroGrid,tau);	
#endif
	
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

#ifdef FIX
	//~ Rescalepi(HydroGrid,tau);	
	//~ FixOne(HydroGrid,tau);	
	FixTwo(HydroGrid,tau);	
#endif

	DebugMSG(HydroGrid);
}


