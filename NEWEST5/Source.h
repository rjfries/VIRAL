void CopyPrimaryVariablesToVar(GRID HydroGrid, double tau)
{

	for( int i=0; i<XCM ; i++)
	for( int j=0; j<YCM ; j++)
	for( int k=0; k<ZCM ; k++)
	{
		HydroGrid[i][j][k].Var[0]= tau*HydroGrid[i][j][k].T00;
		HydroGrid[i][j][k].Var[1]= tau*HydroGrid[i][j][k].T10;
		HydroGrid[i][j][k].Var[2]= tau*HydroGrid[i][j][k].T20;
		HydroGrid[i][j][k].Var[3]= tau*HydroGrid[i][j][k].T30;
		
		for(int l=0;l<Npi;l++)
			HydroGrid[i][j][k].Var[VARN+l]=  HydroGrid[i][j][k].pi[l];
			
			
		
		
		HydroGrid[i][j][k].Var[VARN+Npi]=  HydroGrid[i][j][k].PI;
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

void UpdatePrimaryVariablesAtEndOfTimeStep(GRID HydroGrid, double tau, double ts)
{
	int i,j,k,l;

/*********************************/		
	double KST=1;
	double KSpi=1;
	double KSPI=1;
/*********************************/


	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		
		HydroGrid[i][j][k].T00 =  (tau*HydroGrid[i][j][k].T00 
								+ HydroGrid[i][j][k].Result[0]
								+ KST*HydroGrid[i][j][k].Source[0]*ts)/(tau+ts);
								
		HydroGrid[i][j][k].T10 =  (tau*HydroGrid[i][j][k].T10 
								+ HydroGrid[i][j][k].Result[1])/(tau+ts);
								
		HydroGrid[i][j][k].T20 =  (tau*HydroGrid[i][j][k].T20 
								+ HydroGrid[i][j][k].Result[2])/(tau+ts);
								
		HydroGrid[i][j][k].T30 =  (tau*HydroGrid[i][j][k].T30 
								+ HydroGrid[i][j][k].Result[3]
								+ KST*HydroGrid[i][j][k].Source[3]*ts)/(tau+ts);	

		for(l=0;l<Npi;l++)
			HydroGrid[i][j][k].pi[l] +=  ( HydroGrid[i][j][k].Result[VARN+l] + KSpi*HydroGrid[i][j][k].Source[VARN+l]*ts);
			
		
		/* Fixes for pi's */
		//~ DECLp10u4;
		//~ DECLePPIa;
		//~ double rhomax=1;
		//~ double zeta0=0.1;
		//~ 
		//~ double piupCdown = p1*p1 - 2*p2*p2 -2*p3*p3+ p5*p5 + 2*p6*p6 + p8*p8 - 2*p4*p4*tau*tau  + 2*p7*p7*tau*tau + 2*p9*p9*tau*tau + p10*p10*tau*tau*tau*tau+1e-6;
		//~ double trace = p1 - p5 - p8 - p10*tau*tau;
		//~ double row0  = p1*u0 - p2*u1 - p3*u2 - p4*u3*tau*tau;
		//~ double row1  = p2*u0 - p5*u1 - p6*u2 - p7*u3*tau*tau;
		//~ double row2  = p3*u0 - p6*u1 - p8*u2 - p9*u3*tau*tau;
		//~ double row3  = p4*u0 - p7*u1 - p9*u2 - p10*u3*tau*tau;
		//~ 
		//~ double r[6] = {
						//~ sqrt(piupCdown)/(rhomax*sqrt(e*e+3*P*P) ),
						//~ trace /(rhomax*zeta0*sqrt(piupCdown) ),
						//~ row0 /(rhomax*zeta0*sqrt(piupCdown) ),
						//~ row1 /(rhomax*zeta0*sqrt(piupCdown) ),
						//~ row2 /(rhomax*zeta0*sqrt(piupCdown) ),
						//~ row3 /(rhomax*zeta0*sqrt(piupCdown) ) 
					  //~ };
					  //~ 
		//~ double rho=r[0];
		//~ for(l=1;l<6;l++)
			//~ if(r[l]>rho)
				//~ rho=r[l];
		//~ if(fabs(rho)<1e-5)
			//~ rho=1;
				//~ 
		//~ for(l=0;l<Npi;l++)
			//~ HydroGrid[i][j][k].pi[l] *=  tanh(rho)/rho;
						
			
		HydroGrid[i][j][k].PI += ( HydroGrid[i][j][k].Result[VARN+Npi] + KSPI*HydroGrid[i][j][k].Source[VARN+Npi]*ts);
	}
}


#ifdef KT
void CalcCentreFlux(GRID HydroGrid, double tau)
{
	int i,j,k;
	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		DECLp5u4;
		DECLePPIa;

#if defined CON			
		HydroGrid[i][j][k].Fx[0] =  tau*( HydroGrid[i][j][k].T00*u1/u0    );
		HydroGrid[i][j][k].Fx[1] =  tau*( HydroGrid[i][j][k].T10*u1/u0    );
		HydroGrid[i][j][k].Fx[2] =  tau*( HydroGrid[i][j][k].T20*u1/u0    );
		HydroGrid[i][j][k].Fx[3] =  tau*( HydroGrid[i][j][k].T30*u1/u0    );
#else
		HydroGrid[i][j][k].Fx[0] =  tau*(  A2 + u1*u0*(e+P-PI)            );
		HydroGrid[i][j][k].Fx[1] =  tau*(  p1 + (P-PI) + u1*u1*(e+P-PI)   );
		HydroGrid[i][j][k].Fx[2] =  tau*(  p3 + u1*u2*(e+P-PI)            );
		HydroGrid[i][j][k].Fx[3] =  tau*(  p4 + u1*u3*(e+P-PI)            );
#endif

        HydroGrid[i][j][k].Fx[4]=  (p1*u1)/u0;
        HydroGrid[i][j][k].Fx[5]=  (p2*u1)/u0;
        HydroGrid[i][j][k].Fx[6]=  (p3*u1)/u0;
        HydroGrid[i][j][k].Fx[7]=  (p4*u1)/u0;
        HydroGrid[i][j][k].Fx[8]=  (p5*u1)/u0;
        HydroGrid[i][j][k].Fx[9]=  (PI*u1)/u0;
#if  defined CON
		HydroGrid[i][j][k].Fy[0] =  tau*(  HydroGrid[i][j][k].T00*u2/u0      );
		HydroGrid[i][j][k].Fy[1] =  tau*(  HydroGrid[i][j][k].T10*u2/u0      );
		HydroGrid[i][j][k].Fy[2] =  tau*(  HydroGrid[i][j][k].T20*u2/u0      );
		HydroGrid[i][j][k].Fy[3] =  tau*(  HydroGrid[i][j][k].T30*u2/u0      );
#else
		HydroGrid[i][j][k].Fy[0] =  tau*(  A3 + u2*u0*(e+P-PI)            );
		HydroGrid[i][j][k].Fy[1] =  tau*(  p3 + u2*u1*(e+P-PI)            );
		HydroGrid[i][j][k].Fy[2] =  tau*(  p2 + (P-PI) + u2*u2*(e+P-PI)   );
		HydroGrid[i][j][k].Fy[3] =  tau*(  p5 + u2*u3*(e+P-PI)            );
#endif

		HydroGrid[i][j][k].Fy[4]=  (p1*u2)/u0;
        HydroGrid[i][j][k].Fy[5]=  (p2*u2)/u0;
        HydroGrid[i][j][k].Fy[6]=  (p3*u2)/u0;
        HydroGrid[i][j][k].Fy[7]=  (p4*u2)/u0;
        HydroGrid[i][j][k].Fy[8]=  (p5*u2)/u0;
        HydroGrid[i][j][k].Fy[9]=  (PI*u2)/u0;


#if defined CON			
		HydroGrid[i][j][k].Fz[0] =  tau*(  HydroGrid[i][j][k].T00*u3/u0       );
		HydroGrid[i][j][k].Fz[1] =  tau*(  HydroGrid[i][j][k].T10*u3/u0       );
		HydroGrid[i][j][k].Fz[2] =  tau*(  HydroGrid[i][j][k].T20*u3/u0       );
		HydroGrid[i][j][k].Fz[3] =  tau*(  HydroGrid[i][j][k].T30*u3/u0       );	
#else	
		
  		HydroGrid[i][j][k].Fz[0] =  tau*(  A4 + u3*u0*(e+P-PI)                     );
		HydroGrid[i][j][k].Fz[1] =  tau*(  p4 + u3*u1*(e+P-PI)                     );
		HydroGrid[i][j][k].Fz[2] =  tau*(  p5 + u3*u2*(e+P-PI)                     );
		HydroGrid[i][j][k].Fz[3] =  tau*(  A5 + ((P-PI)/(tau*tau)) + u3*u3*(e+P-PI)     );	
#endif
		HydroGrid[i][j][k].Fz[4]=  (p1*u3)/u0;
        HydroGrid[i][j][k].Fz[5]=  (p2*u3)/u0;
        HydroGrid[i][j][k].Fz[6]=  (p3*u3)/u0;
        HydroGrid[i][j][k].Fz[7]=  (p4*u3)/u0;
        HydroGrid[i][j][k].Fz[8]=  (p5*u3)/u0;
        HydroGrid[i][j][k].Fz[9]=  (PI*u3)/u0;
	}	
}

#endif

void CalcDer4Vel(GRID HydroGrid, double tau, double tstep)
{
	int i,j,k;
	int l;

	for(l=0; l<VARN; l++)
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
		HydroGrid[i][j][k].du[l][0]  =  (HydroGrid[i][j][k].u[l] - HydroGrid[i][j][k].prevu[l])/tstep; 
	
 	
	for(l=0;l<VARN;l++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
	{
		for(i=2;i<XCM-2;i++)
		{
			double  qm2,qm1, qc,qp1,qp2;
			
			
			qc=HydroGrid[i][j][k].u[l];
			qm1=HydroGrid[i-1][j][k].u[l];
			qm2=HydroGrid[i-2][j][k].u[l];
			qp1=HydroGrid[i+1][j][k].u[l];
			qp2=HydroGrid[i+2][j][k].u[l];
			
			//~ HydroGrid[i][j][k].du[l][1] = genminmod(qm1,qc,qp1,XS,GMINV);
			HydroGrid[i][j][k].du[l][1] = genWENOder( qm2, qm1,  qc, qp1, qp2, XS);
		}
	}	 
	
	for(l=0;l<VARN;l++)
	for(i=il;i<ir;i++)
	for(k=kl;k<kr;k++)
	{
		for(j=2;j<YCM-2;j++)
		{
			double  qm2,qm1, qc,qp1,qp2;
			
			qc=HydroGrid[i][j][k].u[l];
			qm1=HydroGrid[i][j-1][k].u[l];
			qm2=HydroGrid[i][j-2][k].u[l];
			qp1=HydroGrid[i][j+1][k].u[l];
			qp2=HydroGrid[i][j+2][k].u[l];
	
			
			//~ HydroGrid[i][j][k].du[l][2] = genminmod(qm1,qc,qp1,YS,GMINV);
			HydroGrid[i][j][k].du[l][2] = genWENOder( qm2, qm1,  qc, qp1, qp2, YS);
		}
	}

			
#if !defined LBI			
	for(l=0;l<VARN;l++)
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	{
		for(k=2;k<ZCM-2;k++)
		{
			double  qm2,qm1, qc,qp1,qp2;
			
			qc=HydroGrid[i][j][k].u[l];
			qm1=HydroGrid[i][j][k-1].u[l];
			qm2=HydroGrid[i][j][k-2].u[l];
			qp1=HydroGrid[i][j][k+1].u[l];
			qp2=HydroGrid[i][j][k+2].u[l];
			

			//~ HydroGrid[i][j][k].du[l][3] = genminmod(qm1,qc,qp1,ZS,GMINV);
			HydroGrid[i][j][k].du[l][3] = genWENOder( qm2, qm1,  qc, qp1, qp2, ZS);
		}
	}
	
#endif
}


void CalcSource(GRID HydroGrid, double tau, double ts)
{
	int i,j,k;
	double max1=0,min1=0;

	double tau2 = tau*tau;	
	double tau3 = tau*tau2;	
	double tau4 = tau2*tau2;	
	double tau5 = tau2*tau3;	
	double tau6 = tau3*tau3;	
		   
		   
//Geometical Source terms
 
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{

		DECLePPIa;
		DECLp5u4;		
		
		HydroGrid[i][j][k].Source[0] =  -(P-PI) - tau2*(A5 + (e + P - PI)*u3*u3);
		HydroGrid[i][j][k].Source[1] = 0;
		HydroGrid[i][j][k].Source[2] = 0;
		HydroGrid[i][j][k].Source[3] = ( -2.0* ( A4 + (e + P - PI)*u0*u3) );
	
	}
	
#ifdef CON
/*Xderivatives*/
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	for(i=il;i<ir;i++)
	{
		double q0[5];
		double q1[5];
		double q2[5];
		double q3[5];
		
		int t=i;
		for(i=t-2;i<=t+2;i++)
		{
			DECLePPIa;
			DECLp10u4;	
				
			q0[i-t+2] = -tau*(       (u1/u0)*(P-PI-p1)  + A2  );
			q1[i-t+2] = -tau*(P-PI - (u1/u0)*A2         + p1  );
			q2[i-t+2] = -tau*(     - (u1/u0)*A3         + p3  );
			q3[i-t+2] = -tau*(     - (u1/u0)*A4         + p4  );
		}
		i=t;
		
		//~ double temp0 = genminmod(q0[1],q0[2],q0[3],XS,GMINV);
		//~ double temp1 = genminmod(q1[1],q1[2],q1[3],XS,GMINV);
		//~ double temp2 = genminmod(q2[1],q2[2],q2[3],XS,GMINV);
		//~ double temp3 = genminmod(q3[1],q3[2],q3[3],XS,GMINV);
		
		double temp0 = genWENOder(q0[0],q0[1],q0[2],q0[3],q0[4],XS);
		double temp1 = genWENOder(q1[0],q1[1],q1[2],q1[3],q1[4],XS);
		double temp2 = genWENOder(q2[0],q2[1],q2[2],q2[3],q2[4],XS);
		double temp3 = genWENOder(q3[0],q3[1],q3[2],q3[3],q3[4],XS);
		
		HydroGrid[i][j][k].Source[0] += temp0;
		HydroGrid[i][j][k].Source[1] += temp1;
		HydroGrid[i][j][k].Source[2] += temp2;
		HydroGrid[i][j][k].Source[3] += temp3;	
	}

/* Y derivatives*/
	for(k=kl;k<kr;k++)
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	{		
		double q0[5];
		double q1[5];
		double q2[5];
		double q3[5];
		
		int t=j;
		for(j=t-2;j<=t+2;j++)
		{
			DECLePPIa;
			DECLp10u4;	
				
			q0[j-t+2] = -tau*(       (u2/u0)*(P-PI-p1)  + A3  );
			q1[j-t+2] = -tau*(	   - (u2/u0)*A2         + p3  );
			q2[j-t+2] = -tau*(P-PI - (u2/u0)*A3         + p2  );
			q3[j-t+2] = -tau*(     - (u2/u0)*A4         + p5  );
		}
		j=t;
		
		//~ double temp0 = genminmod(q0[1],q0[2],q0[3],YS,GMINV);
		//~ double temp1 = genminmod(q1[1],q1[2],q1[3],YS,GMINV);
		//~ double temp2 = genminmod(q2[1],q2[2],q2[3],YS,GMINV);
		//~ double temp3 = genminmod(q3[1],q3[2],q3[3],YS,GMINV);
		
		double temp0 = genWENOder(q0[0],q0[1],q0[2],q0[3],q0[4],YS);
		double temp1 = genWENOder(q1[0],q1[1],q1[2],q1[3],q1[4],YS);
		double temp2 = genWENOder(q2[0],q2[1],q2[2],q2[3],q2[4],YS);
		double temp3 = genWENOder(q3[0],q3[1],q3[2],q3[3],q3[4],YS);
		
		HydroGrid[i][j][k].Source[0] += temp0;
		HydroGrid[i][j][k].Source[1] += temp1;
		HydroGrid[i][j][k].Source[2] += temp2;
		HydroGrid[i][j][k].Source[3] += temp3;	
	}

#endif

	CalcDer4Vel(HydroGrid, tau, ts);   
	
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		
		DECLp5u4;
		DECLcoord;
		DECLePPIa;
		

		double T = FT(e,r);	
		double s = FS(e, P, T);
		double SH =   Feta(  s, e, r);
		double tpi= Ftaupi(  SH, P, e, r);	
		double BU =  FZeta(  s, e, r);
		double tPI= FtauPI(  BU, P, e, r);	
						
		 HydroGrid[i][j][k].Source[4]=  (  ( -(pow(tau,-1)*pow(tpi,-1)*pow(u0,-1)*(p1*(3*tau + 4*tpi*u0) + 6*tau2*tpi*u1*u3*(-2*p4*u0 + A2*u3) - 2*SH*u0*(1 + tau2*pow(u3,2) + pow(u1,2)*(1 - 2*tau2*pow(u3,2)))))/3. )
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-2*SH*u0 + 4*p1*tpi*u0 + 6*A2*tpi*u1*pow(u0,2) + 2*SH*pow(u0,3) - 2*SH*u0*pow(u1,2) - 4*SH*pow(u0,3)*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*p1*tpi*u1 - 4*SH*u1*pow(u0,2) + 6*A2*tpi*u0*pow(u1,2) - 4*SH*pow(u0,2)*pow(u1,3)))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*p1*tpi*u2 + 6*A2*tpi*u0*u1*u2 + 2*SH*u2*pow(u0,2) - 4*SH*u2*pow(u0,2)*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[0][2]) + (-(u3*pow(tpi,-1)*pow(u0,-2)*(3*p1*tpi + 2*u0*(SH*u0 + 3*A2*tpi*u1 - 2*SH*u0*pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u1*pow(u0,2) - 6*p1*tpi*u1*pow(u0,2) + 4*SH*pow(u0,2)*pow(u1,3)))/3.)*(HydroGrid[i][j][k].du[1][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u0 + p1*tpi*u0 + 8*SH*u0*pow(u1,2) - 6*p1*tpi*u0*pow(u1,2) + 4*SH*u0*pow(u1,4)))/3.)*(HydroGrid[i][j][k].du[1][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u0*u1*u2 - 6*p1*tpi*u0*u1*u2 + 4*SH*u0*u2*pow(u1,3)))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((-2*u1*u3*pow(tpi,-1)*pow(u0,-1)*(-3*p1*tpi + 2*SH*(1 + pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-6*p3*tpi*u1*pow(u0,2) - 2*SH*u2*pow(u0,2) + 4*SH*u2*pow(u0,2)*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[2][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u0*u1*u2 - 6*p3*tpi*u0*pow(u1,2) + 4*SH*u0*u2*pow(u1,3)))/3.)*(HydroGrid[i][j][k].du[2][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(-2*SH*u0 + p1*tpi*u0 - 6*p3*tpi*u0*u1*u2 - 2*SH*u0*pow(u1,2) - 2*SH*u0*pow(u2,2) + 4*SH*u0*pow(u1,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(-6*p3*tpi*u0*u1*u3 - 2*SH*u0*u2*u3 + 4*SH*u0*u2*u3*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-6*p4*tau2*tpi*u1*pow(u0,2) - 2*SH*tau2*u3*pow(u0,2) + 4*SH*tau2*u3*pow(u0,2)*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[3][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*tau2*u0*u1*u3 - 6*p4*tau2*tpi*u0*pow(u1,2) + 4*SH*tau2*u0*u3*pow(u1,3)))/3.)*(HydroGrid[i][j][k].du[3][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(-6*p4*tau2*tpi*u0*u1*u2 - 2*SH*tau2*u0*u2*u3 + 4*SH*tau2*u0*u2*u3*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[3][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(-2*SH*u0 + p1*tpi*u0 - 6*p4*tau2*tpi*u0*u1*u3 - 2*SH*u0*pow(u1,2) - 2*SH*tau2*u0*pow(u3,2) + 4*SH*tau2*u0*pow(u1,2)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
		 HydroGrid[i][j][k].Source[5]=  (  ( -(pow(tau,-1)*pow(tpi,-1)*pow(u0,-1)*(p2*(3*tau + 4*tpi*u0) + 6*tau2*tpi*u2*u3*(-2*p5*u0 + A3*u3) - 2*SH*u0*(1 + tau2*pow(u3,2) + pow(u2,2)*(1 - 2*tau2*pow(u3,2)))))/3. )
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-2*SH*u0 + 4*p2*tpi*u0 + 6*A3*tpi*u2*pow(u0,2) + 2*SH*pow(u0,3) - 2*SH*u0*pow(u2,2) - 4*SH*pow(u0,3)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*p2*tpi*u1 + 6*A3*tpi*u0*u1*u2 + 2*SH*u1*pow(u0,2) - 4*SH*u1*pow(u0,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*p2*tpi*u2 - 4*SH*u2*pow(u0,2) + 6*A3*tpi*u0*pow(u2,2) - 4*SH*pow(u0,2)*pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[0][2]) + (-(u3*pow(tpi,-1)*pow(u0,-2)*(3*p2*tpi + 2*u0*(SH*u0 + 3*A3*tpi*u2 - 2*SH*u0*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-2*SH*u1*pow(u0,2) - 6*p3*tpi*u2*pow(u0,2) + 4*SH*u1*pow(u0,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[1][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(-2*SH*u0 + p2*tpi*u0 - 6*p3*tpi*u0*u1*u2 - 2*SH*u0*pow(u1,2) - 2*SH*u0*pow(u2,2) + 4*SH*u0*pow(u1,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[1][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u0*u1*u2 - 6*p3*tpi*u0*pow(u2,2) + 4*SH*u0*u1*pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((2*u3*pow(tpi,-1)*pow(u0,-1)*(3*p3*tpi*u2 + SH*(u1 - 2*u1*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u2*pow(u0,2) - 6*p2*tpi*u2*pow(u0,2) + 4*SH*pow(u0,2)*pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[2][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u0*u1*u2 - 6*p2*tpi*u0*u1*u2 + 4*SH*u0*u1*pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[2][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u0 + p2*tpi*u0 + 8*SH*u0*pow(u2,2) - 6*p2*tpi*u0*pow(u2,2) + 4*SH*u0*pow(u2,4)))/3.)*(HydroGrid[i][j][k].du[2][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u0*u2*u3 - 6*p2*tpi*u0*u2*u3 + 4*SH*u0*u3*pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-6*p5*tau2*tpi*u2*pow(u0,2) - 2*SH*tau2*u3*pow(u0,2) + 4*SH*tau2*u3*pow(u0,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[3][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(-6*p5*tau2*tpi*u0*u1*u2 - 2*SH*tau2*u0*u1*u3 + 4*SH*tau2*u0*u1*u3*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[3][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*tau2*u0*u2*u3 - 6*p5*tau2*tpi*u0*pow(u2,2) + 4*SH*tau2*u0*u3*pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[3][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(-2*SH*u0 + p2*tpi*u0 - 6*p5*tau2*tpi*u0*u2*u3 - 2*SH*u0*pow(u2,2) - 2*SH*tau2*u0*pow(u3,2) + 4*SH*tau2*u0*pow(u2,2)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
		 HydroGrid[i][j][k].Source[6]=  (  ( -(pow(tau,-1)*pow(tpi,-1)*pow(u0,-1)*(p3*(3*tau + 4*tpi*u0) + 3*tau2*tpi*u3*(-2*p5*u0*u1 - 2*p4*u0*u2 + A3*u1*u3 + A2*u2*u3) + 2*SH*u0*u1*u2*(-1 + 2*tau2*pow(u3,2))))/3. )
		+ (-(pow(tpi,-1)*pow(u0,-2)*(4*p3*tpi*u0 - 2*SH*u0*u1*u2 + 3*A3*tpi*u1*pow(u0,2) + 3*A2*tpi*u2*pow(u0,2) - 4*SH*u1*u2*pow(u0,3)))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*p3*tpi*u1 + 3*A2*tpi*u0*u1*u2 - 3*SH*u2*pow(u0,2) + 3*A3*tpi*u0*pow(u1,2) - 4*SH*u2*pow(u0,2)*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*p3*tpi*u2 + 3*A3*tpi*u0*u1*u2 - 3*SH*u1*pow(u0,2) + 3*A2*tpi*u0*pow(u2,2) - 4*SH*u1*pow(u0,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[0][2]) + (-((3*p3*tpi + u0*(3*A3*tpi*u1 + 3*A2*tpi*u2 - 4*SH*u0*u1*u2))*u3*pow(tpi,-1)*pow(u0,-2))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-3*p3*tpi*u1*pow(u0,2) + 3*SH*u2*pow(u0,2) - 3*p1*tpi*u2*pow(u0,2) + 4*SH*u2*pow(u0,2)*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(p3*tpi*u0 + 4*SH*u0*u1*u2 - 3*p1*tpi*u0*u1*u2 - 3*p3*tpi*u0*pow(u1,2) + 4*SH*u0*u2*pow(u1,3)))/3.)*(HydroGrid[i][j][k].du[1][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*SH*u0 - 3*p3*tpi*u0*u1*u2 + 3*SH*u0*pow(u1,2) + 3*SH*u0*pow(u2,2) - 3*p1*tpi*u0*pow(u2,2) + 4*SH*u0*pow(u1,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[1][2]) + (-(u3*pow(tpi,-1)*pow(u0,-1)*(-3*p3*tpi*u1 + 3*SH*u2 - 3*p1*tpi*u2 + 4*SH*u2*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(3*SH*u1*pow(u0,2) - 3*p2*tpi*u1*pow(u0,2) - 3*p3*tpi*u2*pow(u0,2) + 4*SH*u1*pow(u0,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*SH*u0 - 3*p3*tpi*u0*u1*u2 + 3*SH*u0*pow(u1,2) - 3*p2*tpi*u0*pow(u1,2) + 3*SH*u0*pow(u2,2) + 4*SH*u0*pow(u1,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(p3*tpi*u0 + 4*SH*u0*u1*u2 - 3*p2*tpi*u0*u1*u2 - 3*p3*tpi*u0*pow(u2,2) + 4*SH*u0*u1*pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[2][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*SH*u0*u1*u3 - 3*p2*tpi*u0*u1*u3 - 3*p3*tpi*u0*u2*u3 + 4*SH*u0*u1*u3*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-3*p5*tau2*tpi*u1*pow(u0,2) - 3*p4*tau2*tpi*u2*pow(u0,2) + 4*SH*tau2*u1*u2*u3*pow(u0,2)))/3.)*(HydroGrid[i][j][k].du[3][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(-3*p4*tau2*tpi*u0*u1*u2 + 3*SH*tau2*u0*u2*u3 - 3*p5*tau2*tpi*u0*pow(u1,2) + 4*SH*tau2*u0*u2*u3*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[3][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(-3*p5*tau2*tpi*u0*u1*u2 + 3*SH*tau2*u0*u1*u3 - 3*p4*tau2*tpi*u0*pow(u2,2) + 4*SH*tau2*u0*u1*u3*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[3][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(p3*tpi*u0 - 2*SH*u0*u1*u2 - 3*p5*tau2*tpi*u0*u1*u3 - 3*p4*tau2*tpi*u0*u2*u3 + 4*SH*tau2*u0*u1*u2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
		 HydroGrid[i][j][k].Source[7]=  (  ( (pow(tau,-1)*pow(tpi,-1)*pow(u0,-1)*(p4*(-3*tau + tpi*u0*(-7 + 6*tau2*pow(u3,2))) - u3*(3*A2*(tpi + tau2*tpi*pow(u3,2)) + u1*(3*tau2*tpi*(-2*A5*u0 + A4*u3) + 4*SH*(u0 + tau2*u0*pow(u3,2))))))/3. )
		+ (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(4*p4*tau2*tpi*u0 - 2*SH*tau2*u0*u1*u3 + 3*A4*tau2*tpi*u1*pow(u0,2) + 3*A2*tau2*tpi*u3*pow(u0,2) - 4*SH*tau2*u1*u3*pow(u0,3)))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(3*p4*tau2*tpi*u1 + 3*A2*tau2*tpi*u0*u1*u3 - 3*SH*tau2*u3*pow(u0,2) + 3*A4*tau2*tpi*u0*pow(u1,2) - 4*SH*tau2*u3*pow(u0,2)*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(3*p4*tau2*tpi*u2 + 3*A4*tau2*tpi*u0*u1*u2 + 3*A2*tau2*tpi*u0*u2*u3 - 4*SH*tau2*u1*u2*u3*pow(u0,2)))/3.)*(HydroGrid[i][j][k].du[0][2]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(3*tau2*tpi*u3*(p4 + A4*u0*u1 + A2*u0*u3) - SH*u1*pow(u0,2)*(3 + 4*tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(-3*p4*tau2*tpi*u1*pow(u0,2) + 3*SH*tau2*u3*pow(u0,2) - 3*p1*tau2*tpi*u3*pow(u0,2) + 4*SH*tau2*u3*pow(u0,2)*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][0]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(p4*tau2*tpi*u0 + 4*SH*tau2*u0*u1*u3 - 3*p1*tau2*tpi*u0*u1*u3 - 3*p4*tau2*tpi*u0*pow(u1,2) + 4*SH*tau2*u0*u3*pow(u1,3)))/3.)*(HydroGrid[i][j][k].du[1][1]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(-3*p4*tau2*tpi*u0*u1*u2 + 3*SH*tau2*u0*u2*u3 - 3*p1*tau2*tpi*u0*u2*u3 + 4*SH*tau2*u0*u2*u3*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][2]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-1)*(-3*tau2*tpi*u3*(p4*u1 + p1*u3) + SH*(3 + 3*tau2*pow(u3,2) + pow(u1,2)*(3 + 4*tau2*pow(u3,2)))))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(-3*p5*tau2*tpi*u1*pow(u0,2) - 3*p3*tau2*tpi*u3*pow(u0,2) + 4*SH*tau2*u1*u2*u3*pow(u0,2)))/3.)*(HydroGrid[i][j][k].du[2][0]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(-3*p3*tau2*tpi*u0*u1*u3 + 3*SH*tau2*u0*u2*u3 - 3*p5*tau2*tpi*u0*pow(u1,2) + 4*SH*tau2*u0*u2*u3*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[2][1]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(p4*tau2*tpi*u0 - 3*p5*tau2*tpi*u0*u1*u2 - 2*SH*tau2*u0*u1*u3 - 3*p3*tau2*tpi*u0*u2*u3 + 4*SH*tau2*u0*u1*u3*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][2]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(3*SH*u0*u1*u2 - 3*p5*tau2*tpi*u0*u1*u3 - 3*p3*tau2*tpi*u0*pow(u3,2) + 4*SH*tau2*u0*u1*u2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(3*SH*tau2*u1*pow(u0,2) - 3*A5*tau4*tpi*u1*pow(u0,2) - 3*p4*tau4*tpi*u3*pow(u0,2) + 4*SH*tau4*u1*pow(u0,2)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][0]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(3*SH*tau2*u0 - 3*p4*tau4*tpi*u0*u1*u3 + 3*SH*tau2*u0*pow(u1,2) - 3*A5*tau4*tpi*u0*pow(u1,2) + 3*SH*tau4*u0*pow(u3,2) + 4*SH*tau4*u0*pow(u1,2)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][1]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(3*SH*tau2*u0*u1*u2 - 3*A5*tau4*tpi*u0*u1*u2 - 3*p4*tau4*tpi*u0*u2*u3 + 4*SH*tau4*u0*u1*u2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][2]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(p4*tau2*tpi*u0 + 4*SH*tau2*u0*u1*u3 - 3*A5*tau4*tpi*u0*u1*u3 - 3*p4*tau4*tpi*u0*pow(u3,2) + 4*SH*tau4*u0*u1*pow(u3,3)))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
		 HydroGrid[i][j][k].Source[8]=  (  ( (pow(tau,-1)*pow(tpi,-1)*pow(u0,-1)*(p5*(-3*tau + tpi*u0*(-7 + 6*tau2*pow(u3,2))) - u3*(3*A3*(tpi + tau2*tpi*pow(u3,2)) + u2*(3*tau2*tpi*(-2*A5*u0 + A4*u3) + 4*SH*(u0 + tau2*u0*pow(u3,2))))))/3. )
		+ (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(4*p5*tau2*tpi*u0 - 2*SH*tau2*u0*u2*u3 + 3*A4*tau2*tpi*u2*pow(u0,2) + 3*A3*tau2*tpi*u3*pow(u0,2) - 4*SH*tau2*u2*u3*pow(u0,3)))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(3*p5*tau2*tpi*u1 + 3*A4*tau2*tpi*u0*u1*u2 + 3*A3*tau2*tpi*u0*u1*u3 - 4*SH*tau2*u1*u2*u3*pow(u0,2)))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(3*p5*tau2*tpi*u2 + 3*A3*tau2*tpi*u0*u2*u3 - 3*SH*tau2*u3*pow(u0,2) + 3*A4*tau2*tpi*u0*pow(u2,2) - 4*SH*tau2*u3*pow(u0,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[0][2]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(3*tau2*tpi*u3*(p5 + A4*u0*u2 + A3*u0*u3) - SH*u2*pow(u0,2)*(3 + 4*tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(-3*p4*tau2*tpi*u2*pow(u0,2) - 3*p3*tau2*tpi*u3*pow(u0,2) + 4*SH*tau2*u1*u2*u3*pow(u0,2)))/3.)*(HydroGrid[i][j][k].du[1][0]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(p5*tau2*tpi*u0 - 3*p4*tau2*tpi*u0*u1*u2 - 3*p3*tau2*tpi*u0*u1*u3 - 2*SH*tau2*u0*u2*u3 + 4*SH*tau2*u0*u2*u3*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][1]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(3*SH*tau2*u0*u1*u3 - 3*p3*tau2*tpi*u0*u2*u3 - 3*p4*tau2*tpi*u0*pow(u2,2) + 4*SH*tau2*u0*u1*u3*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[1][2]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-1)*(-3*tau2*tpi*u3*(p4*u2 + p3*u3) + SH*u1*u2*(3 + 4*tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(-3*p5*tau2*tpi*u2*pow(u0,2) + 3*SH*tau2*u3*pow(u0,2) - 3*p2*tau2*tpi*u3*pow(u0,2) + 4*SH*tau2*u3*pow(u0,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][0]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(-3*p5*tau2*tpi*u0*u1*u2 + 3*SH*tau2*u0*u1*u3 - 3*p2*tau2*tpi*u0*u1*u3 + 4*SH*tau2*u0*u1*u3*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][1]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(p5*tau2*tpi*u0 + 4*SH*tau2*u0*u2*u3 - 3*p2*tau2*tpi*u0*u2*u3 - 3*p5*tau2*tpi*u0*pow(u2,2) + 4*SH*tau2*u0*u3*pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[2][2]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(3*SH*u0 - 3*p5*tau2*tpi*u0*u2*u3 + 3*SH*u0*pow(u2,2) + 3*SH*tau2*u0*pow(u3,2) - 3*p2*tau2*tpi*u0*pow(u3,2) + 4*SH*tau2*u0*pow(u2,2)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(3*SH*tau2*u2*pow(u0,2) - 3*A5*tau4*tpi*u2*pow(u0,2) - 3*p5*tau4*tpi*u3*pow(u0,2) + 4*SH*tau4*u2*pow(u0,2)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][0]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(3*SH*tau2*u0*u1*u2 - 3*A5*tau4*tpi*u0*u1*u2 - 3*p5*tau4*tpi*u0*u1*u3 + 4*SH*tau4*u0*u1*u2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][1]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(3*SH*tau2*u0 - 3*p5*tau4*tpi*u0*u2*u3 + 3*SH*tau2*u0*pow(u2,2) - 3*A5*tau4*tpi*u0*pow(u2,2) + 3*SH*tau4*u0*pow(u3,2) + 4*SH*tau4*u0*pow(u2,2)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][2]) + (-(pow(tau,-2)*pow(tpi,-1)*pow(u0,-2)*(p5*tau2*tpi*u0 + 4*SH*tau2*u0*u2*u3 - 3*A5*tau4*tpi*u0*u2*u3 - 3*p5*tau4*tpi*u0*pow(u3,2) + 4*SH*tau4*u0*u2*pow(u3,3)))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
		
		
		//~ HydroGrid[i][j][k].Source[4]= 0;
		//~ HydroGrid[i][j][k].Source[5]= 0;
		//~ HydroGrid[i][j][k].Source[6]= 0;
		//~ HydroGrid[i][j][k].Source[7]= 0;
		//~ HydroGrid[i][j][k].Source[8]= 0;
		 
		HydroGrid[i][j][k].Source[9]= 0;
		
	}
}

 
void CalcNS(GRID HydroGrid, double tau, double ts)
{
	int i,j,k;
	
	double tau2 = tau*tau;	
	double tau3 = tau*tau2;	
	double tau4 = tau2*tau2;	
	double tau5 = tau2*tau3;	
		   
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		double r = HydroGrid[i][j][k].r; 
		double PI = HydroGrid[i][j][k].PI; 		
		double e = HydroGrid[i][j][k].En;
		double P =  HydroGrid[i][j][k].P;
				 
		DECLp5u4;
		
		double T =  FT(e,r);	
		double s = FS(e, P, T);
		
		double SH =   Feta(  s, e, r);
		double tpi= Ftaupi(  SH, P, e, r);	
			
		double BU =  FZeta(  s, e, r);
		double tPI= FtauPI(  BU, P, e, r);	
		       
		       
		       
		 HydroGrid[i][j][k].nspi[0]=  (  ( (-2*SH*u0*pow(tau,-1)*(-1 - tau2*pow(u3,2) + pow(u1,2)*(-1 + 2*tau2*pow(u3,2))))/3. )
		+ ((2*SH*(1 + pow(u1,2) + pow(u0,2)*(-1 + 2*pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[0][0]) + ((4*SH*u0*u1*(1 + pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[0][1]) + ((2*SH*u0*u2*(-1 + 2*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[0][2]) + ((2*SH*(u0*u3*(1 + tau*u3)*(-1 + 2*pow(u1,2)) - tau*u0*(-1 + 2*pow(u1,2))*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ ((-4*SH*u0*u1*(1 + pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((-4*SH*pow(1 + pow(u1,2),2))/3.)*(HydroGrid[i][j][k].du[1][1]) + ((-4*SH*u1*u2*(1 + pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((-4*SH*u1*u3*(1 + pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ ((-2*SH*u0*u2*(-1 + 2*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[2][0]) + ((-4*SH*u1*u2*(1 + pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((2*SH*(1 + pow(u2,2) - pow(u1,2)*(-1 + 2*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[2][2]) + ((-2*SH*u2*u3*(-1 + 2*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ ((2*SH*(-(tau*u0*u3*(tau + u3)*(-1 + 2*pow(u1,2))) + tau*u0*(-1 + 2*pow(u1,2))*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][0]) + ((-4*SH*tau2*u1*u3*(1 + pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[3][1]) + ((-2*SH*tau2*u2*u3*(-1 + 2*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((2*SH*(u0*pow(tau,-1)*(-1 - tau2*pow(u3,2) + pow(u1,2)*(-1 + 2*tau2*pow(u3,2))) - (tau + u0)*pow(tau,-1)*(-1 - tau2*pow(u3,2) + pow(u1,2)*(-1 + 2*tau2*pow(u3,2)))))/3.)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].nspi[1]=  (  ( (-2*SH*u0*pow(tau,-1)*(-1 - tau2*pow(u3,2) + pow(u2,2)*(-1 + 2*tau2*pow(u3,2))))/3. )
		+ ((2*SH*(1 + pow(u2,2) + pow(u0,2)*(-1 + 2*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[0][0]) + ((2*SH*u0*u1*(-1 + 2*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[0][1]) + ((4*SH*u0*u2*(1 + pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[0][2]) + ((2*SH*(u0*u3*(1 + tau*u3)*(-1 + 2*pow(u2,2)) - tau*u0*(-1 + 2*pow(u2,2))*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ ((-2*SH*u0*u1*(-1 + 2*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((2*SH*(1 + pow(u2,2) - pow(u1,2)*(-1 + 2*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[1][1]) + ((-4*SH*u1*u2*(1 + pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((-2*SH*u1*u3*(-1 + 2*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ ((-4*SH*u0*u2*(1 + pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][0]) + ((-4*SH*u1*u2*(1 + pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((-4*SH*pow(1 + pow(u2,2),2))/3.)*(HydroGrid[i][j][k].du[2][2]) + ((-4*SH*u2*u3*(1 + pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ ((2*SH*(-(tau*u0*u3*(tau + u3)*(-1 + 2*pow(u2,2))) + tau*u0*(-1 + 2*pow(u2,2))*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][0]) + ((-2*SH*tau2*u1*u3*(-1 + 2*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[3][1]) + ((-4*SH*tau2*u2*u3*(1 + pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((2*SH*(u0*pow(tau,-1)*(-1 - tau2*pow(u3,2) + pow(u2,2)*(-1 + 2*tau2*pow(u3,2))) - (tau + u0)*pow(tau,-1)*(-1 - tau2*pow(u3,2) + pow(u2,2)*(-1 + 2*tau2*pow(u3,2)))))/3.)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].nspi[2]=  (  ( (SH*u0*u1*u2*pow(tau,-1)*(2 - 4*tau2*pow(u3,2)))/3. )
		+ ((SH*pow(tau,-1)*(2*tau*u1*u2 + 4*tau*u1*u2*pow(u0,2) + 2*u0*u1*u2*(-1 + 2*tau2*pow(u3,2)) + u0*u1*(2*u2 - 4*tau2*u2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[0][0]) + ((SH*pow(tau,-1)*(2*u0*u1*u2*(-1 + 2*tau2*pow(u3,2)) + u0*(3*tau*u2 + 4*tau*u2*pow(u1,2) + u1*(2*u2 - 4*tau2*u2*pow(u3,2)))))/3.)*(HydroGrid[i][j][k].du[0][1]) + ((SH*pow(tau,-1)*(2*u0*u1*u2*(-1 + 2*tau2*pow(u3,2)) + u0*u1*(2*u2 + tau*(3 + 4*pow(u2,2)) - 4*tau2*u2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[0][2]) + ((SH*pow(tau,-1)*(2*u0*u1*u2*(-1 + 2*tau2*pow(u3,2)) + u0*u1*(2*u2 + 4*tau*u2*u3 - 4*tau2*u2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ ((SH*pow(tau,-1)*(2*u0*u1*u2*(-1 + 2*tau2*pow(u3,2)) + u0*(-3*tau*u2 - 4*tau*u2*pow(u1,2) + u1*(2*u2 - 4*tau2*u2*pow(u3,2)))))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((SH*pow(tau,-1)*(-2*tau*u1*u2*(2 + 2*pow(u1,2)) + 2*u0*u1*u2*(-1 + 2*tau2*pow(u3,2)) + u0*u1*(2*u2 - 4*tau2*u2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[1][1]) + ((SH*pow(tau,-1)*(-(tau*(3*(1 + pow(u1,2)) + (3 + 4*pow(u1,2))*pow(u2,2))) + 2*u0*u1*u2*(-1 + 2*tau2*pow(u3,2)) + u0*u1*(2*u2 - 4*tau2*u2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((SH*pow(tau,-1)*(-(tau*u2*u3*(3 + 4*pow(u1,2))) + 2*u0*u1*u2*(-1 + 2*tau2*pow(u3,2)) + u0*u1*(2*u2 - 4*tau2*u2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ ((SH*pow(tau,-1)*(2*u0*u1*u2*(-1 + 2*tau2*pow(u3,2)) + u0*u1*(-3*tau + 2*u2 - 4*tau*pow(u2,2) - 4*tau2*u2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[2][0]) + ((SH*pow(tau,-1)*(-(tau*(3*(1 + pow(u1,2)) + (3 + 4*pow(u1,2))*pow(u2,2))) + 2*u0*u1*u2*(-1 + 2*tau2*pow(u3,2)) + u0*u1*(2*u2 - 4*tau2*u2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((SH*pow(tau,-1)*(-(tau*(4*u1*u2 + 4*u1*pow(u2,3))) + 2*u0*u1*u2*(-1 + 2*tau2*pow(u3,2)) + u0*u1*(2*u2 - 4*tau2*u2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[2][2]) + ((SH*pow(tau,-1)*(-(tau*(3*u1*u3 + 4*u1*u3*pow(u2,2))) + 2*u0*u1*u2*(-1 + 2*tau2*pow(u3,2)) + u0*u1*(2*u2 - 4*tau2*u2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ ((SH*pow(tau,-1)*(2*u0*u1*u2*(-1 + 2*tau2*pow(u3,2)) + u0*u1*(2*u2 - 4*tau3*u2*u3 - 4*tau2*u2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[3][0]) + ((SH*pow(tau,-1)*(-(tau3*u2*u3*(3 + 4*pow(u1,2))) + 2*u0*u1*u2*(-1 + 2*tau2*pow(u3,2)) + u0*u1*(2*u2 - 4*tau2*u2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[3][1]) + ((SH*pow(tau,-1)*(-(tau*(3*tau2*u1*u3 + 4*tau2*u1*u3*pow(u2,2))) + 2*u0*u1*u2*(-1 + 2*tau2*pow(u3,2)) + u0*u1*(2*u2 - 4*tau2*u2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((SH*pow(tau,-1)*(2*u0*u1*u2*(-1 + 2*tau2*pow(u3,2)) - tau*u2*(-2*u1 + 4*tau2*u1*pow(u3,2)) + u0*u1*(2*u2 - 4*tau2*u2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].nspi[3]=  (  ( (-4*SH*u0*u1*u3*pow(tau,-1)*(1 + tau2*pow(u3,2)))/3. )
		+ ((2*SH*u1*u3*(1 + 2*pow(u0,2)))/3.)*(HydroGrid[i][j][k].du[0][0]) + ((SH*u0*u3*(3 + 4*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[0][1]) + ((4*SH*u0*u1*u2*u3)/3.)*(HydroGrid[i][j][k].du[0][2]) + ((SH*(u0*u1*(1 + tau*u3)*pow(tau,-2)*(3 + 4*tau2*pow(u3,2)) - u0*u1*u3*pow(tau,-1)*(3 + 4*tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(SH*u0*u3*(3 + 4*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((-4*SH*u1*u3*(1 + pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][1]) + (-(SH*u2*u3*(3 + 4*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][2]) + (-(SH*pow(tau,-2)*(3 + 3*tau2*pow(u3,2) + pow(u1,2)*(3 + 4*tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ ((-4*SH*u0*u1*u2*u3)/3.)*(HydroGrid[i][j][k].du[2][0]) + (-(SH*u2*u3*(3 + 4*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((-2*SH*u1*u3*(-1 + 2*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][2]) + (-(SH*u1*u2*pow(tau,-2)*(3 + 4*tau2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ ((SH*(u0*u1*u3*pow(tau,-1)*(3 + 4*tau2*pow(u3,2)) - u0*u1*(tau + u3)*pow(tau,-1)*(3 + 4*tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[3][0]) + ((SH*(-3 - 3*tau2*pow(u3,2) - pow(u1,2)*(3 + 4*tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[3][1]) + (-(SH*u1*u2*(3 + 4*tau2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((SH*(4*u0*u1*u3*pow(tau,-1)*(1 + tau2*pow(u3,2)) - 4*(tau + u0)*u1*u3*pow(tau,-1)*(1 + tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].nspi[4]=  (  ( (-4*SH*u0*u2*u3*pow(tau,-1)*(1 + tau2*pow(u3,2)))/3. )
		+ ((2*SH*u2*u3*(1 + 2*pow(u0,2)))/3.)*(HydroGrid[i][j][k].du[0][0]) + ((4*SH*u0*u1*u2*u3)/3.)*(HydroGrid[i][j][k].du[0][1]) + ((SH*u0*u3*(3 + 4*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[0][2]) + ((SH*(u0*u2*(1 + tau*u3)*pow(tau,-2)*(3 + 4*tau2*pow(u3,2)) - u0*u2*u3*pow(tau,-1)*(3 + 4*tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ ((-4*SH*u0*u1*u2*u3)/3.)*(HydroGrid[i][j][k].du[1][0]) + ((-2*SH*u2*u3*(-1 + 2*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][1]) + (-(SH*u1*u3*(3 + 4*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[1][2]) + (-(SH*u1*u2*pow(tau,-2)*(3 + 4*tau2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(SH*u0*u3*(3 + 4*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][0]) + (-(SH*u1*u3*(3 + 4*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((-4*SH*u2*u3*(1 + pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][2]) + (-(SH*pow(tau,-2)*(3 + 3*tau2*pow(u3,2) + pow(u2,2)*(3 + 4*tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ ((SH*(u0*u2*u3*pow(tau,-1)*(3 + 4*tau2*pow(u3,2)) - u0*u2*(tau + u3)*pow(tau,-1)*(3 + 4*tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[3][0]) + (-(SH*u1*u2*(3 + 4*tau2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][1]) + ((SH*(-3 - 3*tau2*pow(u3,2) - pow(u2,2)*(3 + 4*tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((SH*(4*u0*u2*u3*pow(tau,-1)*(1 + tau2*pow(u3,2)) - 4*(tau + u0)*u2*u3*pow(tau,-1)*(1 + tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[3][3])  );       
		
		HydroGrid[i][j][k].nsPI=  (  ( -(BU*u0*pow(tau,-1)) )
		+ (-BU)*(HydroGrid[i][j][k].du[0][0]) + (0)*(HydroGrid[i][j][k].du[0][1]) + (0)*(HydroGrid[i][j][k].du[0][2]) + (0)*(HydroGrid[i][j][k].du[0][3])
		+ (0)*(HydroGrid[i][j][k].du[1][0]) + (-BU)*(HydroGrid[i][j][k].du[1][1]) + (0)*(HydroGrid[i][j][k].du[1][2]) + (0)*(HydroGrid[i][j][k].du[1][3])
		+ (0)*(HydroGrid[i][j][k].du[2][0]) + (0)*(HydroGrid[i][j][k].du[2][1]) + (-BU)*(HydroGrid[i][j][k].du[2][2]) + (0)*(HydroGrid[i][j][k].du[2][3])
		+ (0)*(HydroGrid[i][j][k].du[3][0]) + (0)*(HydroGrid[i][j][k].du[3][1]) + (0)*(HydroGrid[i][j][k].du[3][2]) + (-BU)*(HydroGrid[i][j][k].du[3][3]) );
	}
}
/*
void Calcrempi(GRID HydroGrid, int i, int j, int k, double tau)
{	
	double u0 = HydroGrid[i][j][k].u[0];
	double u1 = HydroGrid[i][j][k].u[1];
	double u2 = HydroGrid[i][j][k].u[2];
	double u3 = HydroGrid[i][j][k].u[3];	 
	
	double p1 = HydroGrid[i][j][k].pi[0];
	double p2 = HydroGrid[i][j][k].pi[1];
	double p3 = HydroGrid[i][j][k].pi[2];
	double p4 = HydroGrid[i][j][k].pi[3];
	double p5 = HydroGrid[i][j][k].pi[4];
	
//----------- \[Pi]^\[Mu]\[Nu]- Matrix ---------------------------------------------
	HydroGrid[i][j][k].rempi[0]= -((2*p3*u1*u2 + p1*pow(u1,2) + p2*pow(u2,2) + pow(tau,2)*(2*(p4*u1 + p5*u2)*u3 - (p1 + p2)*pow(u3,2)))*pow(-pow(u0,2) + pow(tau,2)*pow(u3,2),-1));
	HydroGrid[i][j][k].rempi[1]= (p1*u1 + p3*u2 + p4*u3*pow(tau,2))*pow(u0,-1);
	HydroGrid[i][j][k].rempi[2]= (p3*u1 + p2*u2 + p5*u3*pow(tau,2))*pow(u0,-1);
	HydroGrid[i][j][k].rempi[3]= pow(u0,-1)*(2*p3*u1*u2*u3 + (p4*u1 + p5*u2 - (p1 + p2)*u3)*pow(u0,2) + p1*u3*pow(u1,2) + p2*u3*pow(u2,2) + p4*u1*pow(tau,2)*pow(u3,2) + p5*u2*pow(tau,2)*pow(u3,2))*pow(pow(u0,2) - pow(tau,2)*pow(u3,2),-1);
	HydroGrid[i][j][k].rempi[4]= -(pow(tau,-2)*(2*p3*u1*u2 + 2*(p4*u1 + p5*u2)*u3*pow(tau,2) - (p1 + p2)*pow(u0,2) + p1*pow(u1,2) + p2*pow(u2,2))*pow(-pow(u0,2) + pow(tau,2)*pow(u3,2),-1));	
}
*/


