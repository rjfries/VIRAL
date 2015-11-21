
void CalcL0(GRID HydroGrid )
{
	int i,j,k,l;


	for(l=0;l<SVAR;l++)
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
		HydroGrid[i][j][k].L0[l] =  (HydroGrid[i][j][k].Result[l] + HydroGrid[i][j][k].Source[l] ); 
}

void CalcVar1(GRID HydroGrid,  double ts)
{
	int i,j,k,l;


	for(l=0;l<SVAR;l++)
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
		HydroGrid[i][j][k].Var1[l] =  HydroGrid[i][j][k].Var[l]+ts*(HydroGrid[i][j][k].L0[l]); 
}


void FirstOrderUpdatePrimaryVariables( GRID HydroGrid, double tau,double ts)
{
	int i,j,k,l;


	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
	{
		HydroGrid[i][j][k].T00 =  (HydroGrid[i][j][k].Var1[0])/(tau+ts);								
		HydroGrid[i][j][k].T10 =  (HydroGrid[i][j][k].Var1[1])/(tau+ts);								
		HydroGrid[i][j][k].T20 =  (HydroGrid[i][j][k].Var1[2])/(tau+ts);								
		HydroGrid[i][j][k].T30 =  (HydroGrid[i][j][k].Var1[3])/(tau+ts);	

		for(l=0;l<Npi;l++)
			HydroGrid[i][j][k].pi[l] =  ( HydroGrid[i][j][k].Var1[VARN+l]);
			
		HydroGrid[i][j][k].PI = (HydroGrid[i][j][k].Var1[VARN+Npi]);
	}
}


void FirstOrder(GRID HydroGrid, double tau, double ts)
{
	
	CheckPhysics(HydroGrid, 0);	
	
	
	
	
	
	
	
	CalcSource(HydroGrid, tau, ts);		//for entire grid	
#ifdef KT
	CalcCentreFlux(HydroGrid, tau);
#endif
#ifdef SHAS
	CalcNumVel(HydroGrid, tau);
#endif

	hydroExplicit(HydroGrid, tau, ts); 			//gets update for PV,pi and PI everywhere excluding boundary region

	CalcL0(HydroGrid);
	CalcVar1( HydroGrid,ts);
	
	FirstOrderUpdatePrimaryVariables(HydroGrid,  tau, ts); //update PV's, pi and PI to "tau+ts"	from "tau"	
	RootSearchForEnVelUsingDerivatives(HydroGrid, tau+ts );	//Finds En,P,V's, 4vel, everywhere excluding boundary region
	
	DebugMSG(HydroGrid);
	pack(HydroGrid);  //Exchanges En&Vel, pi and PI and updates P,4vel,Tmunu at the cell interfaces
	boundary(HydroGrid);  // Outflowing boundary condition for everything
}



void SOBackUpPrimaryVariables(GRID HydroGrid);
void SORestorePrimaryVariables(GRID HydroGrid);


void TVDRK2(GRID HydroGrid, double tau, double ts)
{
	
	CheckPhysics(HydroGrid, 0);				
	SOBackUpPrimaryVariables(HydroGrid);	
	
	
	
	CalcSource(HydroGrid, tau, ts);	  	
#ifdef KT		
	CalcCentreFlux(HydroGrid, tau);  
#endif
#ifdef SHAS
	CalcNumVel(HydroGrid, tau);
#endif
	
	hydroExplicit(HydroGrid, tau, ts/2); 	
	
	
	//~ UpdatePrimaryVariablesAtEndOfTimeStep( HydroGrid, tau, ts/2);		//update primary vars, pi and PI to "tau+ts/2"	from "tau"
	RootSearchForEnVelUsingDerivatives(HydroGrid , tau + ts/2);   // find en,vel at tau+ts/2and update 4 vel and Pressure. Before updating 4 vel save the current one (at tau) to prevu[4]
    CheckPhysics(HydroGrid, 1);
	/**************************************************************/
	
	pack(HydroGrid );    // update ghost cells
	boundary(HydroGrid);	//implement boundary condition	
	CalcSource(HydroGrid, tau+ts/2, ts/2);	//calculating source after a previous "ts/2" timestep now at time "tau+ts/2"
#ifdef KT		
	CalcCentreFlux(HydroGrid, tau+ts/2);  // calculate centrafluxes at the middle of time step. for second order accuracy in time integration
#endif	
#ifdef SHAS
	CalcNumVel(HydroGrid, tau);
#endif

	/**************************************************************/
	
	
	SORestorePrimaryVariables(HydroGrid);	 // now revert back the PV, pi , PI , en, vel and 4 vel to tau
	hydroExplicit(HydroGrid, tau, ts);	// now advance full time step using the middle of time step values of source terms and centra fluxes
	//~ UpdatePrimaryVariablesAtEndOfTimeStep( HydroGrid,tau, ts);  // //update primary vars, pi and PI to "tau+ts" from "tau"
    RootSearchForEnVelUsingDerivatives(HydroGrid , tau + ts); // find en,vel at tau+ts and update 4 vel and Pressure to "tau+ts". Before updating 4 vel save the current one (again at tau) to prevu[4].
    CheckPhysics(HydroGrid, 2);
	
	/**************************************************************/

	
	pack(HydroGrid);
	boundary(HydroGrid);	

	DebugMSG(HydroGrid);
}



void SOBackUpPrimaryVariables(GRID HydroGrid)
{
	int i,j,k,l;
	
	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		HydroGrid[i][j][k].BackUp[0] = HydroGrid[i][j][k].T00;
		HydroGrid[i][j][k].BackUp[1] = HydroGrid[i][j][k].T10;
		HydroGrid[i][j][k].BackUp[2] = HydroGrid[i][j][k].T20;
		HydroGrid[i][j][k].BackUp[3] = HydroGrid[i][j][k].T30;
		
		HydroGrid[i][j][k].BackUp[4] = HydroGrid[i][j][k].En;
		HydroGrid[i][j][k].BackUp[5] = HydroGrid[i][j][k].Vx;
		HydroGrid[i][j][k].BackUp[6] = HydroGrid[i][j][k].Vy;
		HydroGrid[i][j][k].BackUp[7] = HydroGrid[i][j][k].Ve;
		
		HydroGrid[i][j][k].BackUp[8] = HydroGrid[i][j][k].u[0];
		HydroGrid[i][j][k].BackUp[9] = HydroGrid[i][j][k].u[1];
		HydroGrid[i][j][k].BackUp[10] = HydroGrid[i][j][k].u[2];
		HydroGrid[i][j][k].BackUp[11] = HydroGrid[i][j][k].u[3];
		
		HydroGrid[i][j][k].BackUp[12] = HydroGrid[i][j][k].prevu[0];
		HydroGrid[i][j][k].BackUp[13] = HydroGrid[i][j][k].prevu[1];
		HydroGrid[i][j][k].BackUp[14] = HydroGrid[i][j][k].prevu[2];
		HydroGrid[i][j][k].BackUp[15] = HydroGrid[i][j][k].prevu[3];
		
		HydroGrid[i][j][k].BackUp[4*VARN] = HydroGrid[i][j][k].P;
		
		for(l=0;l<Npi;l++)
			HydroGrid[i][j][k].BackUp[4*VARN+1+l] = HydroGrid[i][j][k].pi[l];
			
		HydroGrid[i][j][k].BackUp[4*VARN+1+Npi] = HydroGrid[i][j][k].PI;
	}
}



void SORestorePrimaryVariables(GRID HydroGrid)
{
	int i,j,k,l;
	
	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		HydroGrid[i][j][k].T00 = HydroGrid[i][j][k].BackUp[0];
		HydroGrid[i][j][k].T10 = HydroGrid[i][j][k].BackUp[1];
		HydroGrid[i][j][k].T20 = HydroGrid[i][j][k].BackUp[2];
		HydroGrid[i][j][k].T30 = HydroGrid[i][j][k].BackUp[3];
		
		HydroGrid[i][j][k].En = HydroGrid[i][j][k].BackUp[4];
		HydroGrid[i][j][k].Vx = HydroGrid[i][j][k].BackUp[5];
		HydroGrid[i][j][k].Vy = HydroGrid[i][j][k].BackUp[6];
		HydroGrid[i][j][k].Ve = HydroGrid[i][j][k].BackUp[7];
		
		HydroGrid[i][j][k].u[0] = HydroGrid[i][j][k].BackUp[8];
		HydroGrid[i][j][k].u[1] = HydroGrid[i][j][k].BackUp[9];
		HydroGrid[i][j][k].u[2] = HydroGrid[i][j][k].BackUp[10];
		HydroGrid[i][j][k].u[3] = HydroGrid[i][j][k].BackUp[11];
		
		HydroGrid[i][j][k].prevu[0] = HydroGrid[i][j][k].BackUp[12];
		HydroGrid[i][j][k].prevu[1] = HydroGrid[i][j][k].BackUp[13];
		HydroGrid[i][j][k].prevu[2] = HydroGrid[i][j][k].BackUp[14];
		HydroGrid[i][j][k].prevu[3] = HydroGrid[i][j][k].BackUp[15];
		
		HydroGrid[i][j][k].P = HydroGrid[i][j][k].BackUp[4*VARN];
		
		for(l=0;l<Npi;l++)
			HydroGrid[i][j][k].pi[l] = HydroGrid[i][j][k].BackUp[4*VARN+1+l];
			
		HydroGrid[i][j][k].PI  = HydroGrid[i][j][k].BackUp[4*VARN+1+Npi];
	}
}
/*
 * 
 * 
 * Fourth ORDER
 * 
 * 
 */
/*
void rk4adder1(ARR phi[VARN], double factor, ARR k0[VARN], ARR result[VARN])
{
	int l, i, j,k;
	
	for(l=0;l<VARN;l++)
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
		result[l][i][j][k] = phi[l][i][j][k] + factor*k0[l][i][j][k];

}

void rk4adder2( ARR phi[VARN],ARR k1[VARN],ARR k2[VARN],ARR k3[VARN],ARR k4[VARN], ARR result[VARN])
{
	int l, i, j,k;
	double f1 = 1.0/6.0;
	double f2 = 2.0/6.0;
	
	
	for(l=0;l<VARN;l++)
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
		result[l][i][j][k] = phi[l][i][j][k] + f1*(k1[l][i][j][k]+k4[l][i][j][k]) +  f2*(k2[l][i][j][k]+ k3[l][i][j][k]);

}


void ruku4()
{



	CalcSource(HydroGrid, tau, ts);		//for entire grid	
	CalcCentreFlux(HydroGrid, tau);
	
	hydroExplicit(HydroGrid, tau, ts); 			//gets update for PV,pi and PI everywhere excluding boundary region
	UpdatePrimaryVariablesAtEndOfTimeStep( HydroGrid, tau, ts); //update PV's, pi and PI to "tau+ts"	from "tau"
	RootSearchForEnVelUsingDerivatives(HydroGrid, tau+ts );	//Finds En,P,V's, 4vel, everywhere excluding boundary region
	pack(HydroGrid);  //Exchanges En&Vel, pi and PI and updates P,4vel,Tmunu at the cell interfaces
	boundary(HydroGrid);  // Outflowing boundary condition for everything
	
	
	
	//RK4 STEP1	
	pack(VARS); 	
	boundary(VARS); 		
	Calc_Source(VARS,tau,Source);
	DebugCFL(VARS,Source);	// for debug messages
	Cal_NumVel(VARS,NumVel); 
	hydroExplicit(VARS,Result,NumVel,ts );
	AddListOfVariables( Source, Result, Source);
	FactorTimesListOfVariables(ts, Source, k1);



	//RK4  STEP2
	rk4adder1(VARS,0.5,k1, k0);
	pack(k0);		
	boundary(k0);		
	Calc_Source(k0,tau+ts/2,Source);
	DebugCFL(k0,Source);  // for debug messages
	Cal_NumVel(k0,NumVel); 
	hydroExplicit(k0,Result,NumVel,ts );
	AddListOfVariables( Source, Result, Source);
	FactorTimesListOfVariables(ts, Source, k2);
	

	//RK4  STEP3
	rk4adder1(VARS,0.5,k2, k0);
	pack(k0);		
	boundary(k0);		
	Calc_Source(k0,tau+ts/2,Source);
	DebugCFL(k0,Source);  // for debug messages
	Cal_NumVel(k0,NumVel); 
	hydroExplicit(k0,Result,NumVel,ts );
	AddListOfVariables( Source, Result, Source);
	FactorTimesListOfVariables(ts, Source, k3);


	//RK4  STEP4
	rk4adder1(VARS,1,k3, k0);
	pack(k0);		
	boundary(k0);		
	Calc_Source(k0,tau+ts,Source);
	DebugCFL(k0,Source);  // for debug messages
	Cal_NumVel(k0,NumVel); 
	hydroExplicit(k0,Result,NumVel,ts );
	AddListOfVariables( Source, Result, Source);
	FactorTimesListOfVariables(ts, Source, k4);


	rk4adder2(VARS,k1,k2,k3,k4,VARS);	

}
*/
