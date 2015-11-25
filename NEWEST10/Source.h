

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

#ifdef VORT

void AddVorticity(GRID HydroGrid, double tau)
{
	int i,j,k,l;
	
	double tau2 = tau*tau;	
	double tau3 = tau*tau2;	
	double tau4 = tau2*tau2;	
	double tau5 = tau2*tau3;	
	double tau6 = tau3*tau3;	
		   
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		
		DECLp10u4; 
		DECLePPIa;
	

		 HydroGrid[i][j][k].Vort[0]=  (  ( -2*p4*tau*u0*u3 + 2*p4*tau*u3*pow(u0,-1) + 2*p2*tau*u1*pow(u0,-1)*pow(u3,2) + 2*p3*tau*u2*pow(u0,-1)*pow(u3,2) + 2*p4*tau3*pow(u0,-1)*pow(u3,3) )
		+ (0)*(HydroGrid[i][j][k].du[0][0]) + (-(p2*u0) + p2*pow(u0,-1) + p3*u1*u2*pow(u0,-1) + p4*tau2*u1*u3*pow(u0,-1) + p2*pow(u0,-1)*pow(u1,2))*(HydroGrid[i][j][k].du[0][1]) + (-(p3*u0) + p3*pow(u0,-1) + p2*u1*u2*pow(u0,-1) + p4*tau2*u2*u3*pow(u0,-1) + p3*pow(u0,-1)*pow(u2,2))*(HydroGrid[i][j][k].du[0][2]) + (-(p4*u0) + p4*pow(u0,-1) + p2*u1*u3*pow(u0,-1) + p3*u2*u3*pow(u0,-1) + p4*tau2*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[0][3])
		+ (-(p2*u0) + p2*pow(u0,-1) + p3*u1*u2*pow(u0,-1) + p4*tau2*u1*u3*pow(u0,-1) + p2*pow(u0,-1)*pow(u1,2))*(HydroGrid[i][j][k].du[1][0]) + (0)*(HydroGrid[i][j][k].du[1][1]) + (p3*u1 - p2*u2)*(HydroGrid[i][j][k].du[1][2]) + (p4*u1 - p2*u3)*(HydroGrid[i][j][k].du[1][3])
		+ (-(p3*u0) + p3*pow(u0,-1) + p2*u1*u2*pow(u0,-1) + p4*tau2*u2*u3*pow(u0,-1) + p3*pow(u0,-1)*pow(u2,2))*(HydroGrid[i][j][k].du[2][0]) + (-(p3*u1) + p2*u2)*(HydroGrid[i][j][k].du[2][1]) + (0)*(HydroGrid[i][j][k].du[2][2]) + (p4*u2 - p3*u3)*(HydroGrid[i][j][k].du[2][3])
		+ (-(p4*tau2*u0) + p4*tau2*pow(u0,-1) + p2*tau2*u1*u3*pow(u0,-1) + p3*tau2*u2*u3*pow(u0,-1) + p4*tau4*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[3][0]) + (-(p4*tau2*u1) + p2*tau2*u3)*(HydroGrid[i][j][k].du[3][1]) + (-(p4*tau2*u2) + p3*tau2*u3)*(HydroGrid[i][j][k].du[3][2]) + (0)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].Vort[1]=  (  ( -(p7*tau*u0*u3) - p4*tau*u1*u3 + p7*tau*u3*pow(u0,-1) + p1*tau*u1*pow(u0,-1)*pow(u3,2) + p5*tau*u1*pow(u0,-1)*pow(u3,2) + p6*tau*u2*pow(u0,-1)*pow(u3,2) + p7*tau3*pow(u0,-1)*pow(u3,3) )
		+ (0)*(HydroGrid[i][j][k].du[0][0]) + (-(p1*u0)/2. - (p5*u0)/2. + (p3*u2)/2. + (p4*tau2*u3)/2. + (p1*pow(u0,-1))/2. + (p5*pow(u0,-1))/2. + (p6*u1*u2*pow(u0,-1))/2. + (p7*tau2*u1*u3*pow(u0,-1))/2. + (p1*pow(u0,-1)*pow(u1,2))/2. + (p5*pow(u0,-1)*pow(u1,2))/2.)*(HydroGrid[i][j][k].du[0][1]) + (-(p6*u0)/2. - (p3*u1)/2. + (p6*pow(u0,-1))/2. + (p1*u1*u2*pow(u0,-1))/2. + (p5*u1*u2*pow(u0,-1))/2. + (p7*tau2*u2*u3*pow(u0,-1))/2. + (p6*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[0][2]) + (-(p7*u0)/2. - (p4*u1)/2. + (p7*pow(u0,-1))/2. + (p1*u1*u3*pow(u0,-1))/2. + (p5*u1*u3*pow(u0,-1))/2. + (p6*u2*u3*pow(u0,-1))/2. + (p7*tau2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(p1*u0)/2. - (p5*u0)/2. + (p3*u2)/2. + (p4*tau2*u3)/2. + (p1*pow(u0,-1))/2. + (p5*pow(u0,-1))/2. + (p6*u1*u2*pow(u0,-1))/2. + (p7*tau2*u1*u3*pow(u0,-1))/2. + (p1*pow(u0,-1)*pow(u1,2))/2. + (p5*pow(u0,-1)*pow(u1,2))/2.)*(HydroGrid[i][j][k].du[1][0]) + (0)*(HydroGrid[i][j][k].du[1][1]) + ((p6*u1)/2. - (p1*u2)/2. - (p5*u2)/2. + (p3*pow(u0,-1))/2. + (p4*tau2*u2*u3*pow(u0,-1))/2. + (p3*pow(u0,-1)*pow(u1,2))/2. + (p3*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[1][2]) + ((p7*u1)/2. - (p1*u3)/2. - (p5*u3)/2. + (p4*pow(u0,-1))/2. + (p3*u2*u3*pow(u0,-1))/2. + (p4*pow(u0,-1)*pow(u1,2))/2. + (p4*tau2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(p6*u0)/2. - (p3*u1)/2. + (p6*pow(u0,-1))/2. + (p1*u1*u2*pow(u0,-1))/2. + (p5*u1*u2*pow(u0,-1))/2. + (p7*tau2*u2*u3*pow(u0,-1))/2. + (p6*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[2][0]) + (-(p6*u1)/2. + (p1*u2)/2. + (p5*u2)/2. - (p3*pow(u0,-1))/2. - (p4*tau2*u2*u3*pow(u0,-1))/2. - (p3*pow(u0,-1)*pow(u1,2))/2. - (p3*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[2][1]) + (0)*(HydroGrid[i][j][k].du[2][2]) + ((p7*u2)/2. - (p6*u3)/2. + (p4*u1*u2*pow(u0,-1))/2. - (p3*u1*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(p7*tau2*u0)/2. - (p4*tau2*u1)/2. + (p7*tau2*pow(u0,-1))/2. + (p1*tau2*u1*u3*pow(u0,-1))/2. + (p5*tau2*u1*u3*pow(u0,-1))/2. + (p6*tau2*u2*u3*pow(u0,-1))/2. + (p7*tau4*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][0]) + (-(p7*tau2*u1)/2. + (p1*tau2*u3)/2. + (p5*tau2*u3)/2. - (p4*tau2*pow(u0,-1))/2. - (p3*tau2*u2*u3*pow(u0,-1))/2. - (p4*tau2*pow(u0,-1)*pow(u1,2))/2. - (p4*tau4*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][1]) + (-(p7*tau2*u2)/2. + (p6*tau2*u3)/2. - (p4*tau2*u1*u2*pow(u0,-1))/2. + (p3*tau2*u1*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[3][2]) + (0)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].Vort[2]=  (  ( -(p9*tau*u0*u3) - p4*tau*u2*u3 + p9*tau*u3*pow(u0,-1) + p6*tau*u1*pow(u0,-1)*pow(u3,2) + p1*tau*u2*pow(u0,-1)*pow(u3,2) + p8*tau*u2*pow(u0,-1)*pow(u3,2) + p9*tau3*pow(u0,-1)*pow(u3,3) )
		+ (0)*(HydroGrid[i][j][k].du[0][0]) + (-(p6*u0)/2. - (p2*u2)/2. + (p6*pow(u0,-1))/2. + (p1*u1*u2*pow(u0,-1))/2. + (p8*u1*u2*pow(u0,-1))/2. + (p9*tau2*u1*u3*pow(u0,-1))/2. + (p6*pow(u0,-1)*pow(u1,2))/2.)*(HydroGrid[i][j][k].du[0][1]) + (-(p1*u0)/2. - (p8*u0)/2. + (p2*u1)/2. + (p4*tau2*u3)/2. + (p1*pow(u0,-1))/2. + (p8*pow(u0,-1))/2. + (p6*u1*u2*pow(u0,-1))/2. + (p9*tau2*u2*u3*pow(u0,-1))/2. + (p1*pow(u0,-1)*pow(u2,2))/2. + (p8*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[0][2]) + (-(p9*u0)/2. - (p4*u2)/2. + (p9*pow(u0,-1))/2. + (p6*u1*u3*pow(u0,-1))/2. + (p1*u2*u3*pow(u0,-1))/2. + (p8*u2*u3*pow(u0,-1))/2. + (p9*tau2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(p6*u0)/2. - (p2*u2)/2. + (p6*pow(u0,-1))/2. + (p1*u1*u2*pow(u0,-1))/2. + (p8*u1*u2*pow(u0,-1))/2. + (p9*tau2*u1*u3*pow(u0,-1))/2. + (p6*pow(u0,-1)*pow(u1,2))/2.)*(HydroGrid[i][j][k].du[1][0]) + (0)*(HydroGrid[i][j][k].du[1][1]) + ((p1*u1)/2. + (p8*u1)/2. - (p6*u2)/2. - (p2*pow(u0,-1))/2. - (p4*tau2*u1*u3*pow(u0,-1))/2. - (p2*pow(u0,-1)*pow(u1,2))/2. - (p2*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[1][2]) + ((p9*u1)/2. - (p6*u3)/2. + (p4*u1*u2*pow(u0,-1))/2. - (p2*u2*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(p1*u0)/2. - (p8*u0)/2. + (p2*u1)/2. + (p4*tau2*u3)/2. + (p1*pow(u0,-1))/2. + (p8*pow(u0,-1))/2. + (p6*u1*u2*pow(u0,-1))/2. + (p9*tau2*u2*u3*pow(u0,-1))/2. + (p1*pow(u0,-1)*pow(u2,2))/2. + (p8*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[2][0]) + (-(p1*u1)/2. - (p8*u1)/2. + (p6*u2)/2. + (p2*pow(u0,-1))/2. + (p4*tau2*u1*u3*pow(u0,-1))/2. + (p2*pow(u0,-1)*pow(u1,2))/2. + (p2*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[2][1]) + (0)*(HydroGrid[i][j][k].du[2][2]) + ((p9*u2)/2. - (p1*u3)/2. - (p8*u3)/2. + (p4*pow(u0,-1))/2. + (p2*u1*u3*pow(u0,-1))/2. + (p4*pow(u0,-1)*pow(u2,2))/2. + (p4*tau2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(p9*tau2*u0)/2. - (p4*tau2*u2)/2. + (p9*tau2*pow(u0,-1))/2. + (p6*tau2*u1*u3*pow(u0,-1))/2. + (p1*tau2*u2*u3*pow(u0,-1))/2. + (p8*tau2*u2*u3*pow(u0,-1))/2. + (p9*tau4*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][0]) + (-(p9*tau2*u1)/2. + (p6*tau2*u3)/2. - (p4*tau2*u1*u2*pow(u0,-1))/2. + (p2*tau2*u2*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[3][1]) + (-(p9*tau2*u2)/2. + (p1*tau2*u3)/2. + (p8*tau2*u3)/2. - (p4*tau2*pow(u0,-1))/2. - (p2*tau2*u1*u3*pow(u0,-1))/2. - (p4*tau2*pow(u0,-1)*pow(u2,2))/2. - (p4*tau4*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][2]) + (0)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].Vort[3]=  (  ( -(p10*tau*u0*u3) - p1*u0*u3*pow(tau,-1) + p2*u1*u3*pow(tau,-1) + p3*u2*u3*pow(tau,-1) + p10*tau*u3*pow(u0,-1) + p1*u3*pow(tau,-1)*pow(u0,-1) + p7*tau*u1*pow(u0,-1)*pow(u3,2) + p9*tau*u2*pow(u0,-1)*pow(u3,2) + p1*tau*pow(u0,-1)*pow(u3,3) + p10*tau3*pow(u0,-1)*pow(u3,3) )
		+ (0)*(HydroGrid[i][j][k].du[0][0]) + (-(p7*u0)/2. - (p2*u3)/2. + (p7*pow(u0,-1))/2. + (p9*u1*u2*pow(u0,-1))/2. + (p1*u1*u3*pow(u0,-1))/2. + (p10*tau2*u1*u3*pow(u0,-1))/2. + (p7*pow(u0,-1)*pow(u1,2))/2.)*(HydroGrid[i][j][k].du[0][1]) + (-(p9*u0)/2. - (p3*u3)/2. + (p9*pow(u0,-1))/2. + (p7*u1*u2*pow(u0,-1))/2. + (p1*u2*u3*pow(u0,-1))/2. + (p10*tau2*u2*u3*pow(u0,-1))/2. + (p9*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[0][2]) + (-(p10*u0)/2. - (p1*u0*pow(tau,-2))/2. + (p2*u1*pow(tau,-2))/2. + (p3*u2*pow(tau,-2))/2. + (p10*pow(u0,-1))/2. + (p7*u1*u3*pow(u0,-1))/2. + (p9*u2*u3*pow(u0,-1))/2. + (p1*pow(tau,-2)*pow(u0,-1))/2. + (p1*pow(u0,-1)*pow(u3,2))/2. + (p10*tau2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(p7*u0)/2. - (p2*u3)/2. + (p7*pow(u0,-1))/2. + (p9*u1*u2*pow(u0,-1))/2. + (p1*u1*u3*pow(u0,-1))/2. + (p10*tau2*u1*u3*pow(u0,-1))/2. + (p7*pow(u0,-1)*pow(u1,2))/2.)*(HydroGrid[i][j][k].du[1][0]) + (0)*(HydroGrid[i][j][k].du[1][1]) + ((p9*u1)/2. - (p7*u2)/2. + (p3*u1*u3*pow(u0,-1))/2. - (p2*u2*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[1][2]) + ((p10*u1)/2. - (p7*u3)/2. + (p1*u1*pow(tau,-2))/2. - (p2*pow(tau,-2)*pow(u0,-1))/2. - (p3*u1*u2*pow(tau,-2)*pow(u0,-1))/2. - (p2*pow(tau,-2)*pow(u0,-1)*pow(u1,2))/2. - (p2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(p9*u0)/2. - (p3*u3)/2. + (p9*pow(u0,-1))/2. + (p7*u1*u2*pow(u0,-1))/2. + (p1*u2*u3*pow(u0,-1))/2. + (p10*tau2*u2*u3*pow(u0,-1))/2. + (p9*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[2][0]) + (-(p9*u1)/2. + (p7*u2)/2. - (p3*u1*u3*pow(u0,-1))/2. + (p2*u2*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[2][1]) + (0)*(HydroGrid[i][j][k].du[2][2]) + ((p10*u2)/2. - (p9*u3)/2. + (p1*u2*pow(tau,-2))/2. - (p3*pow(tau,-2)*pow(u0,-1))/2. - (p2*u1*u2*pow(tau,-2)*pow(u0,-1))/2. - (p3*pow(tau,-2)*pow(u0,-1)*pow(u2,2))/2. - (p3*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(p1*u0)/2. - (p10*tau2*u0)/2. + (p2*u1)/2. + (p3*u2)/2. + (p1*pow(u0,-1))/2. + (p10*tau2*pow(u0,-1))/2. + (p7*tau2*u1*u3*pow(u0,-1))/2. + (p9*tau2*u2*u3*pow(u0,-1))/2. + (p1*tau2*pow(u0,-1)*pow(u3,2))/2. + (p10*tau4*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][0]) + (-(p1*u1)/2. - (p10*tau2*u1)/2. + (p7*tau2*u3)/2. + (p2*pow(u0,-1))/2. + (p3*u1*u2*pow(u0,-1))/2. + (p2*pow(u0,-1)*pow(u1,2))/2. + (p2*tau2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][1]) + (-(p1*u2)/2. - (p10*tau2*u2)/2. + (p9*tau2*u3)/2. + (p3*pow(u0,-1))/2. + (p2*u1*u2*pow(u0,-1))/2. + (p3*pow(u0,-1)*pow(u2,2))/2. + (p3*tau2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][2]) + (0)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].Vort[4]=  (  ( -2*p7*tau*u1*u3 + 2*p2*tau*u1*pow(u0,-1)*pow(u3,2) )
		+ (0)*(HydroGrid[i][j][k].du[0][0]) + (-(p2*u0) + p6*u2 + p7*tau2*u3 + p2*pow(u0,-1) + p2*pow(u0,-1)*pow(u1,2))*(HydroGrid[i][j][k].du[0][1]) + (-(p6*u1) + p2*u1*u2*pow(u0,-1))*(HydroGrid[i][j][k].du[0][2]) + (-(p7*u1) + p2*u1*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[0][3])
		+ (-(p2*u0) + p6*u2 + p7*tau2*u3 + p2*pow(u0,-1) + p2*pow(u0,-1)*pow(u1,2))*(HydroGrid[i][j][k].du[1][0]) + (0)*(HydroGrid[i][j][k].du[1][1]) + (-(p2*u2) + p6*pow(u0,-1) + p7*tau2*u2*u3*pow(u0,-1) + p6*pow(u0,-1)*pow(u1,2) + p6*pow(u0,-1)*pow(u2,2))*(HydroGrid[i][j][k].du[1][2]) + (-(p2*u3) + p7*pow(u0,-1) + p6*u2*u3*pow(u0,-1) + p7*pow(u0,-1)*pow(u1,2) + p7*tau2*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[1][3])
		+ (-(p6*u1) + p2*u1*u2*pow(u0,-1))*(HydroGrid[i][j][k].du[2][0]) + (p2*u2 - p6*pow(u0,-1) - p7*tau2*u2*u3*pow(u0,-1) - p6*pow(u0,-1)*pow(u1,2) - p6*pow(u0,-1)*pow(u2,2))*(HydroGrid[i][j][k].du[2][1]) + (0)*(HydroGrid[i][j][k].du[2][2]) + (p7*u1*u2*pow(u0,-1) - p6*u1*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[2][3])
		+ (-(p7*tau2*u1) + p2*tau2*u1*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[3][0]) + (p2*tau2*u3 - p7*tau2*pow(u0,-1) - p6*tau2*u2*u3*pow(u0,-1) - p7*tau2*pow(u0,-1)*pow(u1,2) - p7*tau4*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[3][1]) + (-(p7*tau2*u1*u2*pow(u0,-1)) + p6*tau2*u1*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[3][2]) + (0)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].Vort[5]=  (  ( -(p9*tau*u1*u3) - p7*tau*u2*u3 + p3*tau*u1*pow(u0,-1)*pow(u3,2) + p2*tau*u2*pow(u0,-1)*pow(u3,2) )
		+ (0)*(HydroGrid[i][j][k].du[0][0]) + (-(p3*u0)/2. - (p5*u2)/2. + (p8*u2)/2. + (p9*tau2*u3)/2. + (p3*pow(u0,-1))/2. + (p2*u1*u2*pow(u0,-1))/2. + (p3*pow(u0,-1)*pow(u1,2))/2.)*(HydroGrid[i][j][k].du[0][1]) + (-(p2*u0)/2. + (p5*u1)/2. - (p8*u1)/2. + (p7*tau2*u3)/2. + (p2*pow(u0,-1))/2. + (p3*u1*u2*pow(u0,-1))/2. + (p2*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[0][2]) + (-(p9*u1)/2. - (p7*u2)/2. + (p3*u1*u3*pow(u0,-1))/2. + (p2*u2*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(p3*u0)/2. - (p5*u2)/2. + (p8*u2)/2. + (p9*tau2*u3)/2. + (p3*pow(u0,-1))/2. + (p2*u1*u2*pow(u0,-1))/2. + (p3*pow(u0,-1)*pow(u1,2))/2.)*(HydroGrid[i][j][k].du[1][0]) + (0)*(HydroGrid[i][j][k].du[1][1]) + ((p2*u1)/2. - (p3*u2)/2. - (p5*pow(u0,-1))/2. + (p8*pow(u0,-1))/2. - (p7*tau2*u1*u3*pow(u0,-1))/2. + (p9*tau2*u2*u3*pow(u0,-1))/2. - (p5*pow(u0,-1)*pow(u1,2))/2. + (p8*pow(u0,-1)*pow(u1,2))/2. - (p5*pow(u0,-1)*pow(u2,2))/2. + (p8*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[1][2]) + (-(p3*u3)/2. + (p9*pow(u0,-1))/2. + (p7*u1*u2*pow(u0,-1))/2. - (p5*u2*u3*pow(u0,-1))/2. + (p8*u2*u3*pow(u0,-1))/2. + (p9*pow(u0,-1)*pow(u1,2))/2. + (p9*tau2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(p2*u0)/2. + (p5*u1)/2. - (p8*u1)/2. + (p7*tau2*u3)/2. + (p2*pow(u0,-1))/2. + (p3*u1*u2*pow(u0,-1))/2. + (p2*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[2][0]) + (-(p2*u1)/2. + (p3*u2)/2. + (p5*pow(u0,-1))/2. - (p8*pow(u0,-1))/2. + (p7*tau2*u1*u3*pow(u0,-1))/2. - (p9*tau2*u2*u3*pow(u0,-1))/2. + (p5*pow(u0,-1)*pow(u1,2))/2. - (p8*pow(u0,-1)*pow(u1,2))/2. + (p5*pow(u0,-1)*pow(u2,2))/2. - (p8*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[2][1]) + (0)*(HydroGrid[i][j][k].du[2][2]) + (-(p2*u3)/2. + (p7*pow(u0,-1))/2. + (p9*u1*u2*pow(u0,-1))/2. + (p5*u1*u3*pow(u0,-1))/2. - (p8*u1*u3*pow(u0,-1))/2. + (p7*pow(u0,-1)*pow(u2,2))/2. + (p7*tau2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(p9*tau2*u1)/2. - (p7*tau2*u2)/2. + (p3*tau2*u1*u3*pow(u0,-1))/2. + (p2*tau2*u2*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[3][0]) + ((p3*tau2*u3)/2. - (p9*tau2*pow(u0,-1))/2. - (p7*tau2*u1*u2*pow(u0,-1))/2. + (p5*tau2*u2*u3*pow(u0,-1))/2. - (p8*tau2*u2*u3*pow(u0,-1))/2. - (p9*tau2*pow(u0,-1)*pow(u1,2))/2. - (p9*tau4*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][1]) + ((p2*tau2*u3)/2. - (p7*tau2*pow(u0,-1))/2. - (p9*tau2*u1*u2*pow(u0,-1))/2. - (p5*tau2*u1*u3*pow(u0,-1))/2. + (p8*tau2*u1*u3*pow(u0,-1))/2. - (p7*tau2*pow(u0,-1)*pow(u2,2))/2. - (p7*tau4*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][2]) + (0)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].Vort[6]=  (  ( -(p10*tau*u1*u3) - p2*u0*u3*pow(tau,-1) + p5*u1*u3*pow(tau,-1) + p6*u2*u3*pow(tau,-1) + p2*u3*pow(tau,-1)*pow(u0,-1) + p4*tau*u1*pow(u0,-1)*pow(u3,2) + p2*tau*pow(u0,-1)*pow(u3,3) )
		+ (0)*(HydroGrid[i][j][k].du[0][0]) + (-(p4*u0)/2. + (p9*u2)/2. - (p5*u3)/2. + (p10*tau2*u3)/2. + (p4*pow(u0,-1))/2. + (p2*u1*u3*pow(u0,-1))/2. + (p4*pow(u0,-1)*pow(u1,2))/2.)*(HydroGrid[i][j][k].du[0][1]) + (-(p9*u1)/2. - (p6*u3)/2. + (p4*u1*u2*pow(u0,-1))/2. + (p2*u2*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[0][2]) + (-(p10*u1)/2. - (p2*u0*pow(tau,-2))/2. + (p5*u1*pow(tau,-2))/2. + (p6*u2*pow(tau,-2))/2. + (p4*u1*u3*pow(u0,-1))/2. + (p2*pow(tau,-2)*pow(u0,-1))/2. + (p2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(p4*u0)/2. + (p9*u2)/2. - (p5*u3)/2. + (p10*tau2*u3)/2. + (p4*pow(u0,-1))/2. + (p2*u1*u3*pow(u0,-1))/2. + (p4*pow(u0,-1)*pow(u1,2))/2.)*(HydroGrid[i][j][k].du[1][0]) + (0)*(HydroGrid[i][j][k].du[1][1]) + (-(p4*u2)/2. + (p9*pow(u0,-1))/2. + (p6*u1*u3*pow(u0,-1))/2. - (p5*u2*u3*pow(u0,-1))/2. + (p10*tau2*u2*u3*pow(u0,-1))/2. + (p9*pow(u0,-1)*pow(u1,2))/2. + (p9*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[1][2]) + (-(p4*u3)/2. + (p2*u1*pow(tau,-2))/2. + (p10*pow(u0,-1))/2. + (p9*u2*u3*pow(u0,-1))/2. - (p5*pow(tau,-2)*pow(u0,-1))/2. - (p6*u1*u2*pow(tau,-2)*pow(u0,-1))/2. + (p10*pow(u0,-1)*pow(u1,2))/2. - (p5*pow(tau,-2)*pow(u0,-1)*pow(u1,2))/2. - (p5*pow(u0,-1)*pow(u3,2))/2. + (p10*tau2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(p9*u1)/2. - (p6*u3)/2. + (p4*u1*u2*pow(u0,-1))/2. + (p2*u2*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[2][0]) + ((p4*u2)/2. - (p9*pow(u0,-1))/2. - (p6*u1*u3*pow(u0,-1))/2. + (p5*u2*u3*pow(u0,-1))/2. - (p10*tau2*u2*u3*pow(u0,-1))/2. - (p9*pow(u0,-1)*pow(u1,2))/2. - (p9*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[2][1]) + (0)*(HydroGrid[i][j][k].du[2][2]) + ((p2*u2*pow(tau,-2))/2. + (p10*u1*u2*pow(u0,-1))/2. - (p9*u1*u3*pow(u0,-1))/2. - (p6*pow(tau,-2)*pow(u0,-1))/2. - (p5*u1*u2*pow(tau,-2)*pow(u0,-1))/2. - (p6*pow(tau,-2)*pow(u0,-1)*pow(u2,2))/2. - (p6*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(p2*u0)/2. + (p5*u1)/2. - (p10*tau2*u1)/2. + (p6*u2)/2. + (p2*pow(u0,-1))/2. + (p4*tau2*u1*u3*pow(u0,-1))/2. + (p2*tau2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][0]) + (-(p2*u1)/2. + (p4*tau2*u3)/2. + (p5*pow(u0,-1))/2. - (p10*tau2*pow(u0,-1))/2. + (p6*u1*u2*pow(u0,-1))/2. - (p9*tau2*u2*u3*pow(u0,-1))/2. + (p5*pow(u0,-1)*pow(u1,2))/2. - (p10*tau2*pow(u0,-1)*pow(u1,2))/2. + (p5*tau2*pow(u0,-1)*pow(u3,2))/2. - (p10*tau4*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][1]) + (-(p2*u2)/2. + (p6*pow(u0,-1))/2. + (p5*u1*u2*pow(u0,-1))/2. - (p10*tau2*u1*u2*pow(u0,-1))/2. + (p9*tau2*u1*u3*pow(u0,-1))/2. + (p6*pow(u0,-1)*pow(u2,2))/2. + (p6*tau2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][2]) + (0)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].Vort[7]=  (  ( -2*p9*tau*u2*u3 + 2*p3*tau*u2*pow(u0,-1)*pow(u3,2) )
		+ (0)*(HydroGrid[i][j][k].du[0][0]) + (-(p6*u2) + p3*u1*u2*pow(u0,-1))*(HydroGrid[i][j][k].du[0][1]) + (-(p3*u0) + p6*u1 + p9*tau2*u3 + p3*pow(u0,-1) + p3*pow(u0,-1)*pow(u2,2))*(HydroGrid[i][j][k].du[0][2]) + (-(p9*u2) + p3*u2*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[0][3])
		+ (-(p6*u2) + p3*u1*u2*pow(u0,-1))*(HydroGrid[i][j][k].du[1][0]) + (0)*(HydroGrid[i][j][k].du[1][1]) + (p3*u1 - p6*pow(u0,-1) - p9*tau2*u1*u3*pow(u0,-1) - p6*pow(u0,-1)*pow(u1,2) - p6*pow(u0,-1)*pow(u2,2))*(HydroGrid[i][j][k].du[1][2]) + (p9*u1*u2*pow(u0,-1) - p6*u2*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[1][3])
		+ (-(p3*u0) + p6*u1 + p9*tau2*u3 + p3*pow(u0,-1) + p3*pow(u0,-1)*pow(u2,2))*(HydroGrid[i][j][k].du[2][0]) + (-(p3*u1) + p6*pow(u0,-1) + p9*tau2*u1*u3*pow(u0,-1) + p6*pow(u0,-1)*pow(u1,2) + p6*pow(u0,-1)*pow(u2,2))*(HydroGrid[i][j][k].du[2][1]) + (0)*(HydroGrid[i][j][k].du[2][2]) + (-(p3*u3) + p9*pow(u0,-1) + p6*u1*u3*pow(u0,-1) + p9*pow(u0,-1)*pow(u2,2) + p9*tau2*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[2][3])
		+ (-(p9*tau2*u2) + p3*tau2*u2*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[3][0]) + (-(p9*tau2*u1*u2*pow(u0,-1)) + p6*tau2*u2*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[3][1]) + (p3*tau2*u3 - p9*tau2*pow(u0,-1) - p6*tau2*u1*u3*pow(u0,-1) - p9*tau2*pow(u0,-1)*pow(u2,2) - p9*tau4*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[3][2]) + (0)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].Vort[8]=  (  ( -(p10*tau*u2*u3) - p3*u0*u3*pow(tau,-1) + p6*u1*u3*pow(tau,-1) + p8*u2*u3*pow(tau,-1) + p3*u3*pow(tau,-1)*pow(u0,-1) + p4*tau*u2*pow(u0,-1)*pow(u3,2) + p3*tau*pow(u0,-1)*pow(u3,3) )
		+ (0)*(HydroGrid[i][j][k].du[0][0]) + (-(p7*u2)/2. - (p6*u3)/2. + (p4*u1*u2*pow(u0,-1))/2. + (p3*u1*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[0][1]) + (-(p4*u0)/2. + (p7*u1)/2. - (p8*u3)/2. + (p10*tau2*u3)/2. + (p4*pow(u0,-1))/2. + (p3*u2*u3*pow(u0,-1))/2. + (p4*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[0][2]) + (-(p10*u2)/2. - (p3*u0*pow(tau,-2))/2. + (p6*u1*pow(tau,-2))/2. + (p8*u2*pow(tau,-2))/2. + (p4*u2*u3*pow(u0,-1))/2. + (p3*pow(tau,-2)*pow(u0,-1))/2. + (p3*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(p7*u2)/2. - (p6*u3)/2. + (p4*u1*u2*pow(u0,-1))/2. + (p3*u1*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[1][0]) + (0)*(HydroGrid[i][j][k].du[1][1]) + ((p4*u1)/2. - (p7*pow(u0,-1))/2. + (p8*u1*u3*pow(u0,-1))/2. - (p10*tau2*u1*u3*pow(u0,-1))/2. - (p6*u2*u3*pow(u0,-1))/2. - (p7*pow(u0,-1)*pow(u1,2))/2. - (p7*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[1][2]) + ((p3*u1*pow(tau,-2))/2. + (p10*u1*u2*pow(u0,-1))/2. - (p7*u2*u3*pow(u0,-1))/2. - (p6*pow(tau,-2)*pow(u0,-1))/2. - (p8*u1*u2*pow(tau,-2)*pow(u0,-1))/2. - (p6*pow(tau,-2)*pow(u0,-1)*pow(u1,2))/2. - (p6*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(p4*u0)/2. + (p7*u1)/2. - (p8*u3)/2. + (p10*tau2*u3)/2. + (p4*pow(u0,-1))/2. + (p3*u2*u3*pow(u0,-1))/2. + (p4*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[2][0]) + (-(p4*u1)/2. + (p7*pow(u0,-1))/2. - (p8*u1*u3*pow(u0,-1))/2. + (p10*tau2*u1*u3*pow(u0,-1))/2. + (p6*u2*u3*pow(u0,-1))/2. + (p7*pow(u0,-1)*pow(u1,2))/2. + (p7*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[2][1]) + (0)*(HydroGrid[i][j][k].du[2][2]) + (-(p4*u3)/2. + (p3*u2*pow(tau,-2))/2. + (p10*pow(u0,-1))/2. + (p7*u1*u3*pow(u0,-1))/2. - (p8*pow(tau,-2)*pow(u0,-1))/2. - (p6*u1*u2*pow(tau,-2)*pow(u0,-1))/2. + (p10*pow(u0,-1)*pow(u2,2))/2. - (p8*pow(tau,-2)*pow(u0,-1)*pow(u2,2))/2. - (p8*pow(u0,-1)*pow(u3,2))/2. + (p10*tau2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(p3*u0)/2. + (p6*u1)/2. + (p8*u2)/2. - (p10*tau2*u2)/2. + (p3*pow(u0,-1))/2. + (p4*tau2*u2*u3*pow(u0,-1))/2. + (p3*tau2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][0]) + (-(p3*u1)/2. + (p6*pow(u0,-1))/2. + (p8*u1*u2*pow(u0,-1))/2. - (p10*tau2*u1*u2*pow(u0,-1))/2. + (p7*tau2*u2*u3*pow(u0,-1))/2. + (p6*pow(u0,-1)*pow(u1,2))/2. + (p6*tau2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][1]) + (-(p3*u2)/2. + (p4*tau2*u3)/2. + (p8*pow(u0,-1))/2. - (p10*tau2*pow(u0,-1))/2. + (p6*u1*u2*pow(u0,-1))/2. - (p7*tau2*u1*u3*pow(u0,-1))/2. + (p8*pow(u0,-1)*pow(u2,2))/2. - (p10*tau2*pow(u0,-1)*pow(u2,2))/2. + (p8*tau2*pow(u0,-1)*pow(u3,2))/2. - (p10*tau4*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][2]) + (0)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].Vort[9]=  (  ( -2*p4*u0*u3*pow(tau,-1) + 2*p7*u1*u3*pow(tau,-1) + 2*p9*u2*u3*pow(tau,-1) + 2*p4*u3*pow(tau,-1)*pow(u0,-1) + 2*p4*tau*pow(u0,-1)*pow(u3,3) )
		+ (0)*(HydroGrid[i][j][k].du[0][0]) + (-(p7*u3) + p4*u1*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[0][1]) + (-(p9*u3) + p4*u2*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[0][2]) + (-(p4*u0*pow(tau,-2)) + p7*u1*pow(tau,-2) + p9*u2*pow(tau,-2) + p4*pow(tau,-2)*pow(u0,-1) + p4*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[0][3])
		+ (-(p7*u3) + p4*u1*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[1][0]) + (0)*(HydroGrid[i][j][k].du[1][1]) + (p9*u1*u3*pow(u0,-1) - p7*u2*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[1][2]) + (p4*u1*pow(tau,-2) - p7*pow(tau,-2)*pow(u0,-1) - p9*u1*u2*pow(tau,-2)*pow(u0,-1) - p7*pow(tau,-2)*pow(u0,-1)*pow(u1,2) - p7*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[1][3])
		+ (-(p9*u3) + p4*u2*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[2][0]) + (-(p9*u1*u3*pow(u0,-1)) + p7*u2*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[2][1]) + (0)*(HydroGrid[i][j][k].du[2][2]) + (p4*u2*pow(tau,-2) - p9*pow(tau,-2)*pow(u0,-1) - p7*u1*u2*pow(tau,-2)*pow(u0,-1) - p9*pow(tau,-2)*pow(u0,-1)*pow(u2,2) - p9*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[2][3])
		+ (-(p4*u0) + p7*u1 + p9*u2 + p4*pow(u0,-1) + p4*tau2*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[3][0]) + (-(p4*u1) + p7*pow(u0,-1) + p9*u1*u2*pow(u0,-1) + p7*pow(u0,-1)*pow(u1,2) + p7*tau2*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[3][1]) + (-(p4*u2) + p9*pow(u0,-1) + p7*u1*u2*pow(u0,-1) + p9*pow(u0,-1)*pow(u2,2) + p9*tau2*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[3][2]) + (0)*(HydroGrid[i][j][k].du[3][3])  );
	}
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	for(l=0;l<Npi;l++)
		  HydroGrid[i][j][k].Source[VARN+l] += HydroGrid[i][j][k].Vort[l];
		  
}

#endif

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
		DECLp10u4;		
		
		HydroGrid[i][j][k].Source[0] =  -(P-PI) - tau2*(p10 + (e + P - PI)*u3*u3);
		HydroGrid[i][j][k].Source[1] = 0;
		HydroGrid[i][j][k].Source[2] = 0;
		HydroGrid[i][j][k].Source[3] = ( -2.0*( p4 + (e + P - PI)*u0*u3) );
	
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
				
			q0[i-t+2] = -tau*(    		  0	    			  );
			q1[i-t+2] = -tau*(P-PI - (u1/u0)*p2         + p5  );
			q2[i-t+2] = -tau*(     - (u1/u0)*p3         + p6  );
			q3[i-t+2] = -tau*(     - (u1/u0)*p4         + p7  );
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
				
			q0[j-t+2] = -tau*(       0                        );
			q1[j-t+2] = -tau*(	   - (u2/u0)*p2         + p6  );
			q2[j-t+2] = -tau*(P-PI - (u2/u0)*p3         + p8  );
			q3[j-t+2] = -tau*(     - (u2/u0)*p4         + p9  );
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
		
		DECLp10u4;
		DECLcoord;
		DECLePPIa;
		

		double T = FT(e,r);	
		double s = FS(e, P, T);
		double SH =   Feta(  s, e, r);
		double tpi= Ftaupi(  SH, P, e, r);	
		double BU =  FZeta(  s, e, r);
		double tPI= FtauPI(  BU, P, e, r);	
		
		 HydroGrid[i][j][k].Source[4]=  (  ( (-3*(2*p4*tau*u3 + p1*pow(tpi,-1))*pow(u0,-1) - 2*pow(tau,-1)*(2*p1 + pow(tpi,-1)*(SH - SH*pow(u0,2))) - 2*tau*(-6*p4*u0*u3 + (3*p1 + pow(tpi,-1)*(SH + 2*SH*pow(u0,2)))*pow(u3,2)))/3. )
		+ ((-2*p1*(3*u0 + 2*pow(u0,-1)) + 4*SH*pow(tpi,-1)*(-2*u0 + pow(u0,-1) + pow(u0,3)))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(u1*(3*p1*(2 + pow(u0,-2)) - 4*SH*pow(tpi,-1)*(-1 + pow(u0,2))))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(u2*(3*p1*(2 + pow(u0,-2)) - 4*SH*pow(tpi,-1)*(-1 + pow(u0,2))))/3.)*(HydroGrid[i][j][k].du[0][2]) + (-(u3*(3*p1*(2 + pow(u0,-2)) - 4*SH*pow(tpi,-1)*(-1 + pow(u0,2))))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ ((2*(3*p2*u0 - 2*SH*u1*pow(tpi,-1)*(-1 + pow(u0,2))))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((6*p2*u1 - p1*pow(u0,-1) - 2*SH*pow(tpi,-1)*(pow(u0,-1)*(1 + pow(u1,2)) + u0*(-1 + 2*pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[1][1]) + ((2*u2*(3*p2 - SH*u1*pow(tpi,-1)*(2*u0 + pow(u0,-1))))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((2*u3*(3*p2 - SH*u1*pow(tpi,-1)*(2*u0 + pow(u0,-1))))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ ((2*(3*p3*u0 - 2*SH*u2*pow(tpi,-1)*(-1 + pow(u0,2))))/3.)*(HydroGrid[i][j][k].du[2][0]) + ((2*u1*(3*p3 - SH*u2*pow(tpi,-1)*(2*u0 + pow(u0,-1))))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((6*p3*u2 - p1*pow(u0,-1) - 2*SH*pow(tpi,-1)*(pow(u0,-1)*(1 + pow(u2,2)) + u0*(-1 + 2*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[2][2]) + ((2*u3*(3*p3 - SH*u2*pow(tpi,-1)*(2*u0 + pow(u0,-1))))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ ((2*tau2*(3*p4*u0 - 2*SH*u3*pow(tpi,-1)*(-1 + pow(u0,2))))/3.)*(HydroGrid[i][j][k].du[3][0]) + ((2*tau2*u1*(3*p4 - SH*u3*pow(tpi,-1)*(2*u0 + pow(u0,-1))))/3.)*(HydroGrid[i][j][k].du[3][1]) + ((2*tau2*u2*(3*p4 - SH*u3*pow(tpi,-1)*(2*u0 + pow(u0,-1))))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((6*p4*tau2*u3 - p1*pow(u0,-1) - 2*SH*pow(tpi,-1)*(pow(u0,-1)*(1 + tau2*pow(u3,2)) + u0*(-1 + 2*tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
		 HydroGrid[i][j][k].Source[5]=  (  ( (pow(tau,-1)*(-4*p2 + 2*SH*u0*u1*pow(tpi,-1)) - 3*pow(u0,-1)*(p7*tau*u3 + p2*pow(tpi,-1) + p1*tau*u1*pow(u3,2)) + tau*(6*(p7*u0 + p4*u1)*u3 - (3*p2 + 4*SH*u0*u1*pow(tpi,-1))*pow(u3,2)))/3. )
		+ ((-3*(p2*u0 + p1*u1) - 4*p2*pow(u0,-1) + 4*SH*u1*pow(tpi,-1)*(-1 + pow(u0,2)))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(p2*u1) - p2*u1*pow(u0,-2) - p1*pow(u0,-1)*pow(u1,2) + (SH*pow(tpi,-1)*(-3*pow(u0,-1)*(1 + pow(u1,2)) + u0*(3 + 4*pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(u2*(SH*u1*pow(tpi,-1)*(-4*u0 + 3*pow(u0,-1)) + 3*(p2 + p2*pow(u0,-2) + p1*u1*pow(u0,-1))))/3.)*(HydroGrid[i][j][k].du[0][2]) + (-(u3*(SH*u1*pow(tpi,-1)*(-4*u0 + 3*pow(u0,-1)) + 3*(p2 + p2*pow(u0,-2) + p1*u1*pow(u0,-1))))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (p5*u0 + p2*u1 - (SH*pow(tpi,-1)*(-3*pow(u0,-1)*(1 + pow(u1,2)) + u0*(3 + 4*pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((3*p5*u1 + p2*pow(u0,-1)*(-1 + 3*pow(u1,2)) - 4*SH*pow(tpi,-1)*(u1 + pow(u1,3)))/3.)*(HydroGrid[i][j][k].du[1][1]) + (u2*(p5 + p2*u1*pow(u0,-1)) - (SH*u2*pow(tpi,-1)*(3 + 4*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][2]) + (u3*(p5 + p2*u1*pow(u0,-1)) - (SH*u3*pow(tpi,-1)*(3 + 4*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (p6*u0 + p3*u1 + pow(tpi,-1)*((-4*SH*u0*u1*u2)/3. + SH*u1*u2*pow(u0,-1)))*(HydroGrid[i][j][k].du[2][0]) + (p6*u1 + p3*pow(u0,-1)*pow(u1,2) - (SH*u2*pow(tpi,-1)*(3 + 4*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((3*p6*u2 - (p2 - 3*p3*u1*u2)*pow(u0,-1) + 2*SH*u1*pow(tpi,-1)*(1 - 2*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][2]) + (p6*u3 - (4*SH*u1*u2*u3*pow(tpi,-1))/3. + p3*u1*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[2][3])
		+ (tau2*(p7*u0 + p4*u1) - (SH*tau2*u1*u3*pow(tpi,-1)*(4*u0 - 3*pow(u0,-1)))/3.)*(HydroGrid[i][j][k].du[3][0]) + (-(SH*tau2*u3*pow(tpi,-1)*(3 + 4*pow(u1,2)))/3. + tau2*(p7*u1 + p4*pow(u0,-1)*pow(u1,2)))*(HydroGrid[i][j][k].du[3][1]) + (p7*tau2*u2 - (4*SH*tau2*u1*u2*u3*pow(tpi,-1))/3. + p4*tau2*u1*u2*pow(u0,-1))*(HydroGrid[i][j][k].du[3][2]) + ((3*p7*tau2*u3 - (p2 - 3*p4*tau2*u1*u3)*pow(u0,-1) + 2*SH*u1*pow(tpi,-1)*(1 - 2*tau2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
		 HydroGrid[i][j][k].Source[6]=  (  ( (pow(tau,-1)*(-4*p3 + 2*SH*u0*u2*pow(tpi,-1)) - 3*pow(u0,-1)*(p9*tau*u3 + p3*pow(tpi,-1) + p1*tau*u2*pow(u3,2)) + tau*(6*(p9*u0 + p4*u2)*u3 - (3*p3 + 4*SH*u0*u2*pow(tpi,-1))*pow(u3,2)))/3. )
		+ ((-3*(p3*u0 + p1*u2) - 4*p3*pow(u0,-1) + 4*SH*u2*pow(tpi,-1)*(-1 + pow(u0,2)))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(u1*(SH*u2*pow(tpi,-1)*(-4*u0 + 3*pow(u0,-1)) + 3*(p3 + p3*pow(u0,-2) + p1*u2*pow(u0,-1))))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(p3*u2) - p3*u2*pow(u0,-2) - p1*pow(u0,-1)*pow(u2,2) + (SH*pow(tpi,-1)*(-3*pow(u0,-1)*(1 + pow(u2,2)) + u0*(3 + 4*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[0][2]) + (-(u3*(SH*u2*pow(tpi,-1)*(-4*u0 + 3*pow(u0,-1)) + 3*(p3 + p3*pow(u0,-2) + p1*u2*pow(u0,-1))))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (p6*u0 + p2*u2 + pow(tpi,-1)*((-4*SH*u0*u1*u2)/3. + SH*u1*u2*pow(u0,-1)))*(HydroGrid[i][j][k].du[1][0]) + ((3*p6*u1 - (p3 - 3*p2*u1*u2)*pow(u0,-1) + 2*SH*u2*pow(tpi,-1)*(1 - 2*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][1]) + (p6*u2 + p2*pow(u0,-1)*pow(u2,2) - (SH*u1*pow(tpi,-1)*(3 + 4*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[1][2]) + (p6*u3 - (4*SH*u1*u2*u3*pow(tpi,-1))/3. + p2*u2*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[1][3])
		+ (p8*u0 + p3*u2 - (SH*pow(tpi,-1)*(-3*pow(u0,-1)*(1 + pow(u2,2)) + u0*(3 + 4*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[2][0]) + (u1*(p8 + p3*u2*pow(u0,-1)) - (SH*u1*pow(tpi,-1)*(3 + 4*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((3*p8*u2 + p3*pow(u0,-1)*(-1 + 3*pow(u2,2)) - 4*SH*pow(tpi,-1)*(u2 + pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[2][2]) + (u3*(p8 + p3*u2*pow(u0,-1)) - (SH*u3*pow(tpi,-1)*(3 + 4*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (tau2*(p9*u0 + p4*u2) - (SH*tau2*u2*u3*pow(tpi,-1)*(4*u0 - 3*pow(u0,-1)))/3.)*(HydroGrid[i][j][k].du[3][0]) + (p9*tau2*u1 - (4*SH*tau2*u1*u2*u3*pow(tpi,-1))/3. + p4*tau2*u1*u2*pow(u0,-1))*(HydroGrid[i][j][k].du[3][1]) + (-(SH*tau2*u3*pow(tpi,-1)*(3 + 4*pow(u2,2)))/3. + tau2*(p9*u2 + p4*pow(u0,-1)*pow(u2,2)))*(HydroGrid[i][j][k].du[3][2]) + ((3*p9*tau2*u3 - (p3 - 3*p4*tau2*u2*u3)*pow(u0,-1) + 2*SH*u2*pow(tpi,-1)*(1 - 2*tau2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
		 HydroGrid[i][j][k].Source[7]=  (  ( (-(pow(tau,-1)*(7*p4 + 4*SH*u0*u3*pow(tpi,-1) + 3*p1*u3*pow(u0,-1))) - 3*pow(u0,-1)*(p10*tau*u3 + p4*pow(tpi,-1) + p1*tau*pow(u3,3)) + tau*(6*p10*u0*u3 + 3*p4*pow(u3,2) - 4*SH*u0*pow(tpi,-1)*pow(u3,3)))/3. )
		+ ((-3*(p4*u0 + p1*u3) - 4*p4*pow(u0,-1) + 4*SH*u3*pow(tpi,-1)*(-1 + pow(u0,2)))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(u1*(SH*u3*pow(tpi,-1)*(-4*u0 + 3*pow(u0,-1)) + 3*(p4 + p4*pow(u0,-2) + p1*u3*pow(u0,-1))))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(u2*(SH*u3*pow(tpi,-1)*(-4*u0 + 3*pow(u0,-1)) + 3*(p4 + p4*pow(u0,-2) + p1*u3*pow(u0,-1))))/3.)*(HydroGrid[i][j][k].du[0][2]) + (-(p4*u3) - p4*u3*pow(u0,-2) + SH*pow(tau,-2)*pow(tpi,-1)*(u0 - pow(u0,-1)) + (4*SH*u0*pow(tpi,-1)*pow(u3,2))/3. - p1*pow(u0,-1)*pow(u3,2) - SH*pow(tpi,-1)*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[0][3])
		+ (p7*u0 + p2*u3 + pow(tpi,-1)*((-4*SH*u0*u1*u3)/3. + SH*u1*u3*pow(u0,-1)))*(HydroGrid[i][j][k].du[1][0]) + ((3*p7*u1 - (p4 - 3*p2*u1*u3)*pow(u0,-1) + 2*SH*u3*pow(tpi,-1)*(1 - 2*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][1]) + (p7*u2 - (4*SH*u1*u2*u3*pow(tpi,-1))/3. + p2*u2*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[1][2]) + (p7*u3 - SH*u1*pow(tau,-2)*pow(tpi,-1) - (4*SH*u1*pow(tpi,-1)*pow(u3,2))/3. + p2*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[1][3])
		+ (p9*u0 + p3*u3 + pow(tpi,-1)*((-4*SH*u0*u2*u3)/3. + SH*u2*u3*pow(u0,-1)))*(HydroGrid[i][j][k].du[2][0]) + (p9*u1 - (4*SH*u1*u2*u3*pow(tpi,-1))/3. + p3*u1*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[2][1]) + ((3*p9*u2 - (p4 - 3*p3*u2*u3)*pow(u0,-1) + 2*SH*u3*pow(tpi,-1)*(1 - 2*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][2]) + (p9*u3 - SH*u2*pow(tau,-2)*pow(tpi,-1) - (4*SH*u2*pow(tpi,-1)*pow(u3,2))/3. + p3*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[2][3])
		+ (tau2*(p10*u0 + p4*u3) - (SH*pow(tpi,-1)*(-3*pow(u0,-1)*(1 + tau2*pow(u3,2)) + u0*(3 + 4*tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[3][0]) + (tau2*u1*(p10 + p4*u3*pow(u0,-1)) - (SH*u1*pow(tpi,-1)*(3 + 4*tau2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][1]) + (tau2*u2*(p10 + p4*u3*pow(u0,-1)) - (SH*u2*pow(tpi,-1)*(3 + 4*tau2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((3*p10*tau2*u3 + p4*pow(u0,-1)*(-1 + 3*tau2*pow(u3,2)) - 4*SH*pow(tpi,-1)*(u3 + tau2*pow(u3,3)))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
		 HydroGrid[i][j][k].Source[8]=  (  ( (pow(tau,-1)*(-4*p5 + 2*SH*pow(tpi,-1)*(1 + pow(u1,2))) + 6*tau*u1*(2*p7*u3 - p2*pow(u0,-1)*pow(u3,2)) + pow(tpi,-1)*(-3*p5*pow(u0,-1) + 2*SH*tau*(1 - 2*pow(u1,2))*pow(u3,2)))/3. )
		+ ((-2*(3*p2*u1 + 2*p5*pow(u0,-1) - SH*pow(tpi,-1)*(pow(u0,-1)*(1 + pow(u1,2)) + u0*(-1 + 2*pow(u1,2)))))/3.)*(HydroGrid[i][j][k].du[0][0]) + ((-3*(p5*u1*pow(u0,-2) + 2*p2*pow(u0,-1)*pow(u1,2)) + 4*SH*pow(tpi,-1)*(u1 + pow(u1,3)))/3.)*(HydroGrid[i][j][k].du[0][1]) + ((u2*(-3*(p5*pow(u0,-2) + 2*p2*u1*pow(u0,-1)) + 2*SH*pow(tpi,-1)*(-1 + 2*pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[0][2]) + ((u3*(-3*(p5*pow(u0,-2) + 2*p2*u1*pow(u0,-1)) + 2*SH*pow(tpi,-1)*(-1 + 2*pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ ((2*(3*p5*u1 - 2*SH*pow(tpi,-1)*(u1 + pow(u1,3))))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((pow(u0,-1)*(p5*(-1 + 6*pow(u1,2)) - 4*SH*pow(tpi,-1)*(1 + 2*pow(u1,2) + pow(u1,4))))/3.)*(HydroGrid[i][j][k].du[1][1]) + ((2*u2*pow(u0,-1)*(3*p5*u1 - 2*SH*pow(tpi,-1)*(u1 + pow(u1,3))))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((2*u3*pow(u0,-1)*(3*p5*u1 - 2*SH*pow(tpi,-1)*(u1 + pow(u1,3))))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ ((2*(3*p6*u1 + SH*u2*pow(tpi,-1)*(1 - 2*pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[2][0]) + ((-2*pow(u0,-1)*(-3*p6*pow(u1,2) + 2*SH*u2*pow(tpi,-1)*(u1 + pow(u1,3))))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((pow(u0,-1)*(-p5 + 6*p6*u1*u2 + 2*SH*pow(tpi,-1)*(1 + pow(u1,2)*(1 - 2*pow(u2,2)) + pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[2][2]) + ((2*u3*pow(u0,-1)*(3*p6*u1 + SH*u2*pow(tpi,-1)*(1 - 2*pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ ((2*tau2*(3*p7*u1 + SH*u3*pow(tpi,-1)*(1 - 2*pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[3][0]) + ((-2*tau2*pow(u0,-1)*(-3*p7*pow(u1,2) + 2*SH*u3*pow(tpi,-1)*(u1 + pow(u1,3))))/3.)*(HydroGrid[i][j][k].du[3][1]) + ((2*tau2*u2*pow(u0,-1)*(3*p7*u1 + SH*u3*pow(tpi,-1)*(1 - 2*pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((pow(u0,-1)*(-p5 + 6*p7*tau2*u1*u3 + 2*SH*pow(tpi,-1)*(1 + tau2*pow(u3,2) + pow(u1,2)*(1 - 2*tau2*pow(u3,2)))))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
		 HydroGrid[i][j][k].Source[9]=  (  ( (pow(tau,-1)*(-4*p6 + 2*SH*u1*u2*pow(tpi,-1)) - pow(tpi,-1)*(3*p6*pow(u0,-1) + 4*SH*tau*u1*u2*pow(u3,2)) + 3*tau*(2*(p9*u1 + p7*u2)*u3 - (p3*u1 + p2*u2)*pow(u0,-1)*pow(u3,2)))/3. )
		+ ((-3*(p3*u1 + p2*u2) - 4*p6*pow(u0,-1) + 2*SH*u1*u2*pow(tpi,-1)*(2*u0 + pow(u0,-1)))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(p6*u1*pow(u0,-2)) + (SH*u2*pow(tpi,-1)*(3 + 4*pow(u1,2)))/3. - pow(u0,-1)*(p2*u1*u2 + p3*pow(u1,2)))*(HydroGrid[i][j][k].du[0][1]) + (-(p6*u2*pow(u0,-2)) + (SH*u1*pow(tpi,-1)*(3 + 4*pow(u2,2)))/3. - pow(u0,-1)*(p3*u1*u2 + p2*pow(u2,2)))*(HydroGrid[i][j][k].du[0][2]) + ((u3*(4*SH*u1*u2*pow(tpi,-1) - 3*(p6*pow(u0,-2) + (p3*u1 + p2*u2)*pow(u0,-1))))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (p6*u1 + p5*u2 - (SH*u2*pow(tpi,-1)*(3 + 4*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((pow(u0,-1)*(-p6 + 3*p5*u1*u2 + 3*p6*pow(u1,2) - 4*SH*u2*pow(tpi,-1)*(u1 + pow(u1,3))))/3.)*(HydroGrid[i][j][k].du[1][1]) + ((pow(u0,-1)*(3*(p6*u1*u2 + p5*pow(u2,2)) - SH*pow(tpi,-1)*(3*(1 + pow(u2,2)) + pow(u1,2)*(3 + 4*pow(u2,2)))))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((u3*pow(u0,-1)*(3*(p6*u1 + p5*u2) - SH*u2*pow(tpi,-1)*(3 + 4*pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (p8*u1 + p6*u2 - (SH*u1*pow(tpi,-1)*(3 + 4*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][0]) + ((pow(u0,-1)*(3*(p6*u1*u2 + p8*pow(u1,2)) - SH*pow(tpi,-1)*(3*(1 + pow(u2,2)) + pow(u1,2)*(3 + 4*pow(u2,2)))))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((pow(u0,-1)*(-p6 + 3*p8*u1*u2 + 3*p6*pow(u2,2) - 4*SH*u1*pow(tpi,-1)*(u2 + pow(u2,3))))/3.)*(HydroGrid[i][j][k].du[2][2]) + ((u3*pow(u0,-1)*(3*(p8*u1 + p6*u2) - SH*u1*pow(tpi,-1)*(3 + 4*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (p9*tau2*u1 + p7*tau2*u2 - (4*SH*tau2*u1*u2*u3*pow(tpi,-1))/3.)*(HydroGrid[i][j][k].du[3][0]) + ((tau2*pow(u0,-1)*(-(SH*u2*u3*pow(tpi,-1)*(3 + 4*pow(u1,2))) + 3*(p7*u1*u2 + p9*pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[3][1]) + ((tau2*pow(u0,-1)*(-(SH*u1*u3*pow(tpi,-1)*(3 + 4*pow(u2,2))) + 3*(p9*u1*u2 + p7*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[3][2]) + (-(pow(u0,-1)*(p6 - 3*tau2*(p9*u1 + p7*u2)*u3 + 2*SH*u1*u2*pow(tpi,-1)*(-1 + 2*tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
		 HydroGrid[i][j][k].Source[10]=  (  ( (-(pow(tau,-1)*(7*p7 + 4*SH*u1*u3*pow(tpi,-1) + 3*p2*u3*pow(u0,-1))) - pow(tpi,-1)*(3*p7*pow(u0,-1) + 4*SH*tau*u1*pow(u3,3)) + 3*tau*(2*p10*u1*u3 + (2*p7 - p4*u1*pow(u0,-1))*pow(u3,2) - p2*pow(u0,-1)*pow(u3,3)))/3. )
		+ ((-3*(p4*u1 + p2*u3) - 4*p7*pow(u0,-1) + 2*SH*u1*u3*pow(tpi,-1)*(2*u0 + pow(u0,-1)))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(p7*u1*pow(u0,-2)) + (SH*u3*pow(tpi,-1)*(3 + 4*pow(u1,2)))/3. - pow(u0,-1)*(p2*u1*u3 + p4*pow(u1,2)))*(HydroGrid[i][j][k].du[0][1]) + ((u2*(4*SH*u1*u3*pow(tpi,-1) - 3*(p7*pow(u0,-2) + (p4*u1 + p2*u3)*pow(u0,-1))))/3.)*(HydroGrid[i][j][k].du[0][2]) + (SH*u1*pow(tau,-2)*pow(tpi,-1) - p7*u3*pow(u0,-2) - p4*u1*u3*pow(u0,-1) + (4*SH*u1*pow(tpi,-1)*pow(u3,2))/3. - p2*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[0][3])
		+ (p7*u1 + p5*u3 - (SH*u3*pow(tpi,-1)*(3 + 4*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((pow(u0,-1)*(-p7 + 3*p5*u1*u3 + 3*p7*pow(u1,2) - 4*SH*u3*pow(tpi,-1)*(u1 + pow(u1,3))))/3.)*(HydroGrid[i][j][k].du[1][1]) + ((u2*pow(u0,-1)*(3*(p7*u1 + p5*u3) - SH*u3*pow(tpi,-1)*(3 + 4*pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((pow(u0,-1)*(3*p7*u1*u3 - 3*SH*pow(tau,-2)*pow(tpi,-1)*(1 + pow(u1,2)) + (3*p5 - SH*pow(tpi,-1)*(3 + 4*pow(u1,2)))*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (p9*u1 + p6*u3 - (4*SH*u1*u2*u3*pow(tpi,-1))/3.)*(HydroGrid[i][j][k].du[2][0]) + ((pow(u0,-1)*(-(SH*u2*u3*pow(tpi,-1)*(3 + 4*pow(u1,2))) + 3*(p6*u1*u3 + p9*pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[2][1]) + (-(pow(u0,-1)*(p7 - 3*u2*(p9*u1 + p6*u3) + 2*SH*u1*u3*pow(tpi,-1)*(-1 + 2*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[2][2]) + ((pow(u0,-1)*(3*p9*u1*u3 - 3*SH*u1*u2*pow(tau,-2)*pow(tpi,-1) + (3*p6 - 4*SH*u1*u2*pow(tpi,-1))*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (tau2*(p10*u1 + p7*u3) - (SH*u1*pow(tpi,-1)*(3 + 4*tau2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][0]) + ((pow(u0,-1)*(3*tau2*(p7*u1*u3 + p10*pow(u1,2)) - SH*pow(tpi,-1)*(3 + 3*tau2*pow(u3,2) + pow(u1,2)*(3 + 4*tau2*pow(u3,2)))))/3.)*(HydroGrid[i][j][k].du[3][1]) + ((u2*pow(u0,-1)*(3*tau2*(p10*u1 + p7*u3) - SH*u1*pow(tpi,-1)*(3 + 4*tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((pow(u0,-1)*(-p7 + 3*p10*tau2*u1*u3 + 3*p7*tau2*pow(u3,2) - 4*SH*u1*pow(tpi,-1)*(u3 + tau2*pow(u3,3))))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
		 HydroGrid[i][j][k].Source[11]=  (  ( (pow(tau,-1)*(-4*p8 + 2*SH*pow(tpi,-1)*(1 + pow(u2,2))) + 6*tau*u2*(2*p9*u3 - p3*pow(u0,-1)*pow(u3,2)) + pow(tpi,-1)*(-3*p8*pow(u0,-1) + 2*SH*tau*(1 - 2*pow(u2,2))*pow(u3,2)))/3. )
		+ ((-2*(3*p3*u2 + 2*p8*pow(u0,-1) - SH*pow(tpi,-1)*(pow(u0,-1)*(1 + pow(u2,2)) + u0*(-1 + 2*pow(u2,2)))))/3.)*(HydroGrid[i][j][k].du[0][0]) + ((u1*(-3*(p8*pow(u0,-2) + 2*p3*u2*pow(u0,-1)) + 2*SH*pow(tpi,-1)*(-1 + 2*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[0][1]) + ((-3*(p8*u2*pow(u0,-2) + 2*p3*pow(u0,-1)*pow(u2,2)) + 4*SH*pow(tpi,-1)*(u2 + pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[0][2]) + ((u3*(-3*(p8*pow(u0,-2) + 2*p3*u2*pow(u0,-1)) + 2*SH*pow(tpi,-1)*(-1 + 2*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ ((2*(3*p6*u2 + SH*u1*pow(tpi,-1)*(1 - 2*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((pow(u0,-1)*(-p8 + 6*p6*u1*u2 + 2*SH*pow(tpi,-1)*(1 + pow(u1,2)*(1 - 2*pow(u2,2)) + pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[1][1]) + ((-2*pow(u0,-1)*(-3*p6*pow(u2,2) + 2*SH*u1*pow(tpi,-1)*(u2 + pow(u2,3))))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((2*u3*pow(u0,-1)*(3*p6*u2 + SH*u1*pow(tpi,-1)*(1 - 2*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ ((2*(3*p8*u2 - 2*SH*pow(tpi,-1)*(u2 + pow(u2,3))))/3.)*(HydroGrid[i][j][k].du[2][0]) + ((2*u1*pow(u0,-1)*(3*p8*u2 - 2*SH*pow(tpi,-1)*(u2 + pow(u2,3))))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((pow(u0,-1)*(p8*(-1 + 6*pow(u2,2)) - 4*SH*pow(tpi,-1)*(1 + 2*pow(u2,2) + pow(u2,4))))/3.)*(HydroGrid[i][j][k].du[2][2]) + ((2*u3*pow(u0,-1)*(3*p8*u2 - 2*SH*pow(tpi,-1)*(u2 + pow(u2,3))))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ ((2*tau2*(3*p9*u2 + SH*u3*pow(tpi,-1)*(1 - 2*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[3][0]) + ((2*tau2*u1*pow(u0,-1)*(3*p9*u2 + SH*u3*pow(tpi,-1)*(1 - 2*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[3][1]) + ((-2*tau2*pow(u0,-1)*(-3*p9*pow(u2,2) + 2*SH*u3*pow(tpi,-1)*(u2 + pow(u2,3))))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((pow(u0,-1)*(-p8 + 6*p9*tau2*u2*u3 + 2*SH*pow(tpi,-1)*(1 + tau2*pow(u3,2) + pow(u2,2)*(1 - 2*tau2*pow(u3,2)))))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
		 HydroGrid[i][j][k].Source[12]=  (  ( (-(pow(tau,-1)*(7*p9 + 4*SH*u2*u3*pow(tpi,-1) + 3*p3*u3*pow(u0,-1))) - pow(tpi,-1)*(3*p9*pow(u0,-1) + 4*SH*tau*u2*pow(u3,3)) + 3*tau*(2*p10*u2*u3 + (2*p9 - p4*u2*pow(u0,-1))*pow(u3,2) - p3*pow(u0,-1)*pow(u3,3)))/3. )
		+ ((-3*(p4*u2 + p3*u3) - 4*p9*pow(u0,-1) + 2*SH*u2*u3*pow(tpi,-1)*(2*u0 + pow(u0,-1)))/3.)*(HydroGrid[i][j][k].du[0][0]) + ((u1*(4*SH*u2*u3*pow(tpi,-1) - 3*(p9*pow(u0,-2) + (p4*u2 + p3*u3)*pow(u0,-1))))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(p9*u2*pow(u0,-2)) + (SH*u3*pow(tpi,-1)*(3 + 4*pow(u2,2)))/3. - pow(u0,-1)*(p3*u2*u3 + p4*pow(u2,2)))*(HydroGrid[i][j][k].du[0][2]) + (SH*u2*pow(tau,-2)*pow(tpi,-1) - p9*u3*pow(u0,-2) - p4*u2*u3*pow(u0,-1) + (4*SH*u2*pow(tpi,-1)*pow(u3,2))/3. - p3*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[0][3])
		+ (p7*u2 + p6*u3 - (4*SH*u1*u2*u3*pow(tpi,-1))/3.)*(HydroGrid[i][j][k].du[1][0]) + (-(pow(u0,-1)*(p9 - 3*u1*(p7*u2 + p6*u3) + 2*SH*u2*u3*pow(tpi,-1)*(-1 + 2*pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[1][1]) + ((pow(u0,-1)*(-(SH*u1*u3*pow(tpi,-1)*(3 + 4*pow(u2,2))) + 3*(p6*u2*u3 + p7*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((pow(u0,-1)*(3*p7*u2*u3 - 3*SH*u1*u2*pow(tau,-2)*pow(tpi,-1) + (3*p6 - 4*SH*u1*u2*pow(tpi,-1))*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (p9*u2 + p8*u3 - (SH*u3*pow(tpi,-1)*(3 + 4*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][0]) + ((u1*pow(u0,-1)*(3*(p9*u2 + p8*u3) - SH*u3*pow(tpi,-1)*(3 + 4*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((pow(u0,-1)*(-p9 + 3*p8*u2*u3 + 3*p9*pow(u2,2) - 4*SH*u3*pow(tpi,-1)*(u2 + pow(u2,3))))/3.)*(HydroGrid[i][j][k].du[2][2]) + ((pow(u0,-1)*(3*p9*u2*u3 - 3*SH*pow(tau,-2)*pow(tpi,-1)*(1 + pow(u2,2)) + (3*p8 - SH*pow(tpi,-1)*(3 + 4*pow(u2,2)))*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (tau2*(p10*u2 + p9*u3) - (SH*u2*pow(tpi,-1)*(3 + 4*tau2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][0]) + ((u1*pow(u0,-1)*(3*tau2*(p10*u2 + p9*u3) - SH*u2*pow(tpi,-1)*(3 + 4*tau2*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[3][1]) + ((pow(u0,-1)*(3*tau2*(p9*u2*u3 + p10*pow(u2,2)) - SH*pow(tpi,-1)*(3 + 3*tau2*pow(u3,2) + pow(u2,2)*(3 + 4*tau2*pow(u3,2)))))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((pow(u0,-1)*(-p9 + 3*p10*tau2*u2*u3 + 3*p9*tau2*pow(u3,2) - 4*SH*u2*pow(tpi,-1)*(u3 + tau2*pow(u3,3))))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
		 HydroGrid[i][j][k].Source[13]=  (  ( (-4*SH*pow(tau,-3)*pow(tpi,-1) - 3*p10*pow(tpi,-1)*pow(u0,-1) + 12*p10*tau*pow(u3,2) - 2*pow(tau,-1)*(5*p10 + 3*p4*u3*pow(u0,-1) + 4*SH*pow(tpi,-1)*pow(u3,2)) - 6*p4*tau*pow(u0,-1)*pow(u3,3) - 4*SH*tau*pow(tpi,-1)*pow(u3,4))/3. )
		+ ((-2*(3*p4*u3 + SH*pow(tau,-2)*pow(tpi,-1)*(u0 - pow(u0,-1)) - 2*SH*u0*pow(tpi,-1)*pow(u3,2) + pow(u0,-1)*(2*p10 - SH*pow(tpi,-1)*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(u1*(2*SH*pow(tau,-2)*pow(tpi,-1) + 3*p10*pow(u0,-2) + 6*p4*u3*pow(u0,-1) - 4*SH*pow(tpi,-1)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(u2*(2*SH*pow(tau,-2)*pow(tpi,-1) + 3*p10*pow(u0,-2) + 6*p4*u3*pow(u0,-1) - 4*SH*pow(tpi,-1)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[0][2]) + ((4*SH*u3*pow(tau,-2)*pow(tpi,-1))/3. - p10*u3*pow(u0,-2) - 2*p4*pow(u0,-1)*pow(u3,2) + (4*SH*pow(tpi,-1)*pow(u3,3))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ ((2*(3*p7*u3 + SH*u1*pow(tau,-2)*pow(tpi,-1) - 2*SH*u1*pow(tpi,-1)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((pow(u0,-1)*(-p10 + 6*p7*u1*u3 + 2*SH*pow(tau,-2)*pow(tpi,-1)*(1 + pow(u1,2)) + 2*SH*pow(tpi,-1)*(1 - 2*pow(u1,2))*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[1][1]) + ((2*u2*pow(u0,-1)*(3*p7*u3 + SH*u1*pow(tau,-2)*pow(tpi,-1) - 2*SH*u1*pow(tpi,-1)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((-2*pow(u0,-1)*(2*SH*u1*u3*pow(tau,-2)*pow(tpi,-1) - 3*p7*pow(u3,2) + 2*SH*u1*pow(tpi,-1)*pow(u3,3)))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ ((2*(3*p9*u3 + SH*u2*pow(tau,-2)*pow(tpi,-1) - 2*SH*u2*pow(tpi,-1)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[2][0]) + ((2*u1*pow(u0,-1)*(3*p9*u3 + SH*u2*pow(tau,-2)*pow(tpi,-1) - 2*SH*u2*pow(tpi,-1)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((pow(u0,-1)*(-p10 + 6*p9*u2*u3 + 2*SH*pow(tau,-2)*pow(tpi,-1)*(1 + pow(u2,2)) + 2*SH*pow(tpi,-1)*(1 - 2*pow(u2,2))*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[2][2]) + ((-2*pow(u0,-1)*(2*SH*u2*u3*pow(tau,-2)*pow(tpi,-1) - 3*p9*pow(u3,2) + 2*SH*u2*pow(tpi,-1)*pow(u3,3)))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ ((2*(3*p10*tau2*u3 - 2*SH*pow(tpi,-1)*(u3 + tau2*pow(u3,3))))/3.)*(HydroGrid[i][j][k].du[3][0]) + ((2*u1*pow(u0,-1)*(3*p10*tau2*u3 - 2*SH*pow(tpi,-1)*(u3 + tau2*pow(u3,3))))/3.)*(HydroGrid[i][j][k].du[3][1]) + ((2*u2*pow(u0,-1)*(3*p10*tau2*u3 - 2*SH*pow(tpi,-1)*(u3 + tau2*pow(u3,3))))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((pow(u0,-1)*(-p10 - 4*SH*pow(tau,-2)*pow(tpi,-1) + (6*p10*tau2 - 8*SH*pow(tpi,-1))*pow(u3,2) - 4*SH*tau2*pow(tpi,-1)*pow(u3,4)))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
#ifdef BULK
		 HydroGrid[i][j][k].Source[14]=  ( (-((3*PI*tau + 3*BU*u0 + 4*PI*tPI*u0)*pow(tau,-1)*pow(tPI,-1)*pow(u0,-1))/3. )
		+ (-((3*BU + 4*PI*tPI)*pow(tPI,-1)*pow(u0,-1))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(PI*u1*pow(u0,-2)))*(HydroGrid[i][j][k].du[0][1]) + (-(PI*u2*pow(u0,-2)))*(HydroGrid[i][j][k].du[0][2]) + (-(PI*u3*pow(u0,-2)))*(HydroGrid[i][j][k].du[0][3])
		+ (0)*(HydroGrid[i][j][k].du[1][0]) + (-((3*BU + PI*tPI)*pow(tPI,-1)*pow(u0,-1))/3.)*(HydroGrid[i][j][k].du[1][1]) + (0)*(HydroGrid[i][j][k].du[1][2]) + (0)*(HydroGrid[i][j][k].du[1][3])
		+ (0)*(HydroGrid[i][j][k].du[2][0]) + (0)*(HydroGrid[i][j][k].du[2][1]) + (-((3*BU + PI*tPI)*pow(tPI,-1)*pow(u0,-1))/3.)*(HydroGrid[i][j][k].du[2][2]) + (0)*(HydroGrid[i][j][k].du[2][3])
		+ (0)*(HydroGrid[i][j][k].du[3][0]) + (0)*(HydroGrid[i][j][k].du[3][1]) + (0)*(HydroGrid[i][j][k].du[3][2]) + (-((3*BU + PI*tPI)*pow(tPI,-1)*pow(u0,-1))/3.)*(HydroGrid[i][j][k].du[3][3]));
#else	        
		HydroGrid[i][j][k].Source[14]= 0;
#endif


	}
		
#ifdef VORT
	AddVorticity(HydroGrid,tau);
#endif 
}

 
void CalcNS(GRID HydroGrid, double tau, double ts)
{
	int i,j,k;
	double max1=0,min1=0;
	double tau2 = tau*tau;	

	CalcDer4Vel(HydroGrid, tau, ts); //du[4][4] is computed
		

	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		double r = HydroGrid[i][j][k].r;
		double u0 = HydroGrid[i][j][k].u[0];
		double u1 = HydroGrid[i][j][k].u[1];
		double u2 = HydroGrid[i][j][k].u[2];
		double u3 = HydroGrid[i][j][k].u[3];	
		double e = HydroGrid[i][j][k].En;	
		double P =  HydroGrid[i][j][k].P;
		double T =  FT(e,r);	
		double s = FS(e, P, T);
		
		double SH = Feta(s, e, r);
		double tpi = Ftaupi( SH, P, e, r);	
		
		double BU = FZeta(s, e, r);
		double tPI = FtauPI(BU,P, e, r);		
						   
       
		            
		 HydroGrid[i][j][k].nspi[0]=  (  ( (-2*SH*u0*pow(tau,-1))/3. + (2*SH*pow(tau,-1)*pow(u0,3))/3. - (2*SH*tau*u0*pow(u3,2))/3. - (4*SH*tau*pow(u0,3)*pow(u3,2))/3. )
		+ ((4*SH)/3. - (8*SH*pow(u0,2))/3. + (4*SH*pow(u0,4))/3.)*(HydroGrid[i][j][k].du[0][0]) + ((-4*SH*u0*u1)/3. + (4*SH*u1*pow(u0,3))/3.)*(HydroGrid[i][j][k].du[0][1]) + ((-4*SH*u0*u2)/3. + (4*SH*u2*pow(u0,3))/3.)*(HydroGrid[i][j][k].du[0][2]) + ((-4*SH*u0*u3)/3. + (4*SH*u3*pow(u0,3))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ ((4*SH*u0*u1)/3. - (4*SH*u1*pow(u0,3))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((-2*SH)/3. + (2*SH*pow(u0,2))/3. - (2*SH*pow(u1,2))/3. - (4*SH*pow(u0,2)*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[1][1]) + ((-2*SH*u1*u2)/3. - (4*SH*u1*u2*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((-2*SH*u1*u3)/3. - (4*SH*u1*u3*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ ((4*SH*u0*u2)/3. - (4*SH*u2*pow(u0,3))/3.)*(HydroGrid[i][j][k].du[2][0]) + ((-2*SH*u1*u2)/3. - (4*SH*u1*u2*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((-2*SH)/3. + (2*SH*pow(u0,2))/3. - (2*SH*pow(u2,2))/3. - (4*SH*pow(u0,2)*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[2][2]) + ((-2*SH*u2*u3)/3. - (4*SH*u2*u3*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ ((4*SH*tau2*u0*u3)/3. - (4*SH*tau2*u3*pow(u0,3))/3.)*(HydroGrid[i][j][k].du[3][0]) + ((-2*SH*tau2*u1*u3)/3. - (4*SH*tau2*u1*u3*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[3][1]) + ((-2*SH*tau2*u2*u3)/3. - (4*SH*tau2*u2*u3*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((-2*SH)/3. + (2*SH*pow(u0,2))/3. - (2*SH*tau2*pow(u3,2))/3. - (4*SH*tau2*pow(u0,2)*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].nspi[1]=  (  ( (2*SH*u1*pow(tau,-1)*pow(u0,2))/3. - (4*SH*tau*u1*pow(u0,2)*pow(u3,2))/3. )
		+ ((-4*SH*u0*u1)/3. + (4*SH*u1*pow(u0,3))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-SH + SH*pow(u0,2) - SH*pow(u1,2) + (4*SH*pow(u0,2)*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(SH*u1*u2) + (4*SH*u1*u2*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[0][2]) + (-(SH*u1*u3) + (4*SH*u1*u3*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (SH - SH*pow(u0,2) + SH*pow(u1,2) - (4*SH*pow(u0,2)*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((-4*SH*u0*u1)/3. - (4*SH*u0*pow(u1,3))/3.)*(HydroGrid[i][j][k].du[1][1]) + (-(SH*u0*u2) - (4*SH*u0*u2*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[1][2]) + (-(SH*u0*u3) - (4*SH*u0*u3*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (SH*u1*u2 - (4*SH*u1*u2*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[2][0]) + (-(SH*u0*u2) - (4*SH*u0*u2*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((2*SH*u0*u1)/3. - (4*SH*u0*u1*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[2][2]) + ((-4*SH*u0*u1*u2*u3)/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (SH*tau2*u1*u3 - (4*SH*tau2*u1*u3*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[3][0]) + (-(SH*tau2*u0*u3) - (4*SH*tau2*u0*u3*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[3][1]) + ((-4*SH*tau2*u0*u1*u2*u3)/3.)*(HydroGrid[i][j][k].du[3][2]) + ((2*SH*u0*u1)/3. - (4*SH*tau2*u0*u1*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].nspi[2]=  (  ( (2*SH*u2*pow(tau,-1)*pow(u0,2))/3. - (4*SH*tau*u2*pow(u0,2)*pow(u3,2))/3. )
		+ ((-4*SH*u0*u2)/3. + (4*SH*u2*pow(u0,3))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(SH*u1*u2) + (4*SH*u1*u2*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-SH + SH*pow(u0,2) - SH*pow(u2,2) + (4*SH*pow(u0,2)*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[0][2]) + (-(SH*u2*u3) + (4*SH*u2*u3*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (SH*u1*u2 - (4*SH*u1*u2*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((2*SH*u0*u2)/3. - (4*SH*u0*u2*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[1][1]) + (-(SH*u0*u1) - (4*SH*u0*u1*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((-4*SH*u0*u1*u2*u3)/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (SH - SH*pow(u0,2) + SH*pow(u2,2) - (4*SH*pow(u0,2)*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[2][0]) + (-(SH*u0*u1) - (4*SH*u0*u1*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((-4*SH*u0*u2)/3. - (4*SH*u0*pow(u2,3))/3.)*(HydroGrid[i][j][k].du[2][2]) + (-(SH*u0*u3) - (4*SH*u0*u3*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (SH*tau2*u2*u3 - (4*SH*tau2*u2*u3*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[3][0]) + ((-4*SH*tau2*u0*u1*u2*u3)/3.)*(HydroGrid[i][j][k].du[3][1]) + (-(SH*tau2*u0*u3) - (4*SH*tau2*u0*u3*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((2*SH*u0*u2)/3. - (4*SH*tau2*u0*u2*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].nspi[3]=  (  ( (-4*SH*u3*pow(tau,-1)*pow(u0,2))/3. - (4*SH*tau*pow(u0,2)*pow(u3,3))/3. )
		+ ((-4*SH*u0*u3)/3. + (4*SH*u3*pow(u0,3))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(SH*u1*u3) + (4*SH*u1*u3*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(SH*u2*u3) + (4*SH*u2*u3*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[0][2]) + (-(SH*pow(tau,-2)) + SH*pow(tau,-2)*pow(u0,2) - SH*pow(u3,2) + (4*SH*pow(u0,2)*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (SH*u1*u3 - (4*SH*u1*u3*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((2*SH*u0*u3)/3. - (4*SH*u0*u3*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[1][1]) + ((-4*SH*u0*u1*u2*u3)/3.)*(HydroGrid[i][j][k].du[1][2]) + (-(SH*u0*u1*pow(tau,-2)) - (4*SH*u0*u1*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (SH*u2*u3 - (4*SH*u2*u3*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[2][0]) + ((-4*SH*u0*u1*u2*u3)/3.)*(HydroGrid[i][j][k].du[2][1]) + ((2*SH*u0*u3)/3. - (4*SH*u0*u3*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[2][2]) + (-(SH*u0*u2*pow(tau,-2)) - (4*SH*u0*u2*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (SH - SH*pow(u0,2) + SH*tau2*pow(u3,2) - (4*SH*tau2*pow(u0,2)*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[3][0]) + (-(SH*u0*u1) - (4*SH*tau2*u0*u1*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[3][1]) + (-(SH*u0*u2) - (4*SH*tau2*u0*u2*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((-4*SH*u0*u3)/3. - (4*SH*tau2*u0*pow(u3,3))/3.)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].nspi[4]=  (  ( (2*SH*u0*pow(tau,-1))/3. + (2*SH*u0*pow(tau,-1)*pow(u1,2))/3. + (2*SH*tau*u0*pow(u3,2))/3. - (4*SH*tau*u0*pow(u1,2)*pow(u3,2))/3. )
		+ ((2*SH)/3. - (2*SH*pow(u0,2))/3. + (2*SH*pow(u1,2))/3. + (4*SH*pow(u0,2)*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[0][0]) + ((4*SH*u0*u1)/3. + (4*SH*u0*pow(u1,3))/3.)*(HydroGrid[i][j][k].du[0][1]) + ((-2*SH*u0*u2)/3. + (4*SH*u0*u2*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[0][2]) + ((-2*SH*u0*u3)/3. + (4*SH*u0*u3*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ ((-4*SH*u0*u1)/3. - (4*SH*u0*pow(u1,3))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((-4*SH)/3. - (8*SH*pow(u1,2))/3. - (4*SH*pow(u1,4))/3.)*(HydroGrid[i][j][k].du[1][1]) + ((-4*SH*u1*u2)/3. - (4*SH*u2*pow(u1,3))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((-4*SH*u1*u3)/3. - (4*SH*u3*pow(u1,3))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ ((2*SH*u0*u2)/3. - (4*SH*u0*u2*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[2][0]) + ((-4*SH*u1*u2)/3. - (4*SH*u2*pow(u1,3))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((2*SH)/3. + (2*SH*pow(u1,2))/3. + (2*SH*pow(u2,2))/3. - (4*SH*pow(u1,2)*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[2][2]) + ((2*SH*u2*u3)/3. - (4*SH*u2*u3*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ ((2*SH*tau2*u0*u3)/3. - (4*SH*tau2*u0*u3*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[3][0]) + ((-4*SH*tau2*u1*u3)/3. - (4*SH*tau2*u3*pow(u1,3))/3.)*(HydroGrid[i][j][k].du[3][1]) + ((2*SH*tau2*u2*u3)/3. - (4*SH*tau2*u2*u3*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((2*SH)/3. + (2*SH*pow(u1,2))/3. + (2*SH*tau2*pow(u3,2))/3. - (4*SH*tau2*pow(u1,2)*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].nspi[5]=  (  ( (2*SH*u0*u1*u2*pow(tau,-1))/3. - (4*SH*tau*u0*u1*u2*pow(u3,2))/3. )
		+ ((2*SH*u1*u2)/3. + (4*SH*u1*u2*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[0][0]) + (SH*u0*u2 + (4*SH*u0*u2*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[0][1]) + (SH*u0*u1 + (4*SH*u0*u1*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[0][2]) + ((4*SH*u0*u1*u2*u3)/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(SH*u0*u2) - (4*SH*u0*u2*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((-4*SH*u1*u2)/3. - (4*SH*u2*pow(u1,3))/3.)*(HydroGrid[i][j][k].du[1][1]) + (-SH - SH*pow(u1,2) - SH*pow(u2,2) - (4*SH*pow(u1,2)*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[1][2]) + (-(SH*u2*u3) - (4*SH*u2*u3*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(SH*u0*u1) - (4*SH*u0*u1*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[2][0]) + (-SH - SH*pow(u1,2) - SH*pow(u2,2) - (4*SH*pow(u1,2)*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((-4*SH*u1*u2)/3. - (4*SH*u1*pow(u2,3))/3.)*(HydroGrid[i][j][k].du[2][2]) + (-(SH*u1*u3) - (4*SH*u1*u3*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ ((-4*SH*tau2*u0*u1*u2*u3)/3.)*(HydroGrid[i][j][k].du[3][0]) + (-(SH*tau2*u2*u3) - (4*SH*tau2*u2*u3*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[3][1]) + (-(SH*tau2*u1*u3) - (4*SH*tau2*u1*u3*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((2*SH*u1*u2)/3. - (4*SH*tau2*u1*u2*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].nspi[6]=  (  ( (-4*SH*u0*u1*u3*pow(tau,-1))/3. - (4*SH*tau*u0*u1*pow(u3,3))/3. )
		+ ((2*SH*u1*u3)/3. + (4*SH*u1*u3*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[0][0]) + (SH*u0*u3 + (4*SH*u0*u3*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[0][1]) + ((4*SH*u0*u1*u2*u3)/3.)*(HydroGrid[i][j][k].du[0][2]) + (SH*u0*u1*pow(tau,-2) + (4*SH*u0*u1*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(SH*u0*u3) - (4*SH*u0*u3*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((-4*SH*u1*u3)/3. - (4*SH*u3*pow(u1,3))/3.)*(HydroGrid[i][j][k].du[1][1]) + (-(SH*u2*u3) - (4*SH*u2*u3*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[1][2]) + (-(SH*pow(tau,-2)) - SH*pow(tau,-2)*pow(u1,2) - SH*pow(u3,2) - (4*SH*pow(u1,2)*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ ((-4*SH*u0*u1*u2*u3)/3.)*(HydroGrid[i][j][k].du[2][0]) + (-(SH*u2*u3) - (4*SH*u2*u3*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((2*SH*u1*u3)/3. - (4*SH*u1*u3*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[2][2]) + (-(SH*u1*u2*pow(tau,-2)) - (4*SH*u1*u2*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(SH*u0*u1) - (4*SH*tau2*u0*u1*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[3][0]) + (-SH - SH*pow(u1,2) - SH*tau2*pow(u3,2) - (4*SH*tau2*pow(u1,2)*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[3][1]) + (-(SH*u1*u2) - (4*SH*tau2*u1*u2*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((-4*SH*u1*u3)/3. - (4*SH*tau2*u1*pow(u3,3))/3.)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].nspi[7]=  (  ( (2*SH*u0*pow(tau,-1))/3. + (2*SH*u0*pow(tau,-1)*pow(u2,2))/3. + (2*SH*tau*u0*pow(u3,2))/3. - (4*SH*tau*u0*pow(u2,2)*pow(u3,2))/3. )
		+ ((2*SH)/3. - (2*SH*pow(u0,2))/3. + (2*SH*pow(u2,2))/3. + (4*SH*pow(u0,2)*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[0][0]) + ((-2*SH*u0*u1)/3. + (4*SH*u0*u1*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[0][1]) + ((4*SH*u0*u2)/3. + (4*SH*u0*pow(u2,3))/3.)*(HydroGrid[i][j][k].du[0][2]) + ((-2*SH*u0*u3)/3. + (4*SH*u0*u3*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ ((2*SH*u0*u1)/3. - (4*SH*u0*u1*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((2*SH)/3. + (2*SH*pow(u1,2))/3. + (2*SH*pow(u2,2))/3. - (4*SH*pow(u1,2)*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[1][1]) + ((-4*SH*u1*u2)/3. - (4*SH*u1*pow(u2,3))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((2*SH*u1*u3)/3. - (4*SH*u1*u3*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ ((-4*SH*u0*u2)/3. - (4*SH*u0*pow(u2,3))/3.)*(HydroGrid[i][j][k].du[2][0]) + ((-4*SH*u1*u2)/3. - (4*SH*u1*pow(u2,3))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((-4*SH)/3. - (8*SH*pow(u2,2))/3. - (4*SH*pow(u2,4))/3.)*(HydroGrid[i][j][k].du[2][2]) + ((-4*SH*u2*u3)/3. - (4*SH*u3*pow(u2,3))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ ((2*SH*tau2*u0*u3)/3. - (4*SH*tau2*u0*u3*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[3][0]) + ((2*SH*tau2*u1*u3)/3. - (4*SH*tau2*u1*u3*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[3][1]) + ((-4*SH*tau2*u2*u3)/3. - (4*SH*tau2*u3*pow(u2,3))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((2*SH)/3. + (2*SH*pow(u2,2))/3. + (2*SH*tau2*pow(u3,2))/3. - (4*SH*tau2*pow(u2,2)*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].nspi[8]=  (  ( (-4*SH*u0*u2*u3*pow(tau,-1))/3. - (4*SH*tau*u0*u2*pow(u3,3))/3. )
		+ ((2*SH*u2*u3)/3. + (4*SH*u2*u3*pow(u0,2))/3.)*(HydroGrid[i][j][k].du[0][0]) + ((4*SH*u0*u1*u2*u3)/3.)*(HydroGrid[i][j][k].du[0][1]) + (SH*u0*u3 + (4*SH*u0*u3*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[0][2]) + (SH*u0*u2*pow(tau,-2) + (4*SH*u0*u2*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ ((-4*SH*u0*u1*u2*u3)/3.)*(HydroGrid[i][j][k].du[1][0]) + ((2*SH*u2*u3)/3. - (4*SH*u2*u3*pow(u1,2))/3.)*(HydroGrid[i][j][k].du[1][1]) + (-(SH*u1*u3) - (4*SH*u1*u3*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[1][2]) + (-(SH*u1*u2*pow(tau,-2)) - (4*SH*u1*u2*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(SH*u0*u3) - (4*SH*u0*u3*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[2][0]) + (-(SH*u1*u3) - (4*SH*u1*u3*pow(u2,2))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((-4*SH*u2*u3)/3. - (4*SH*u3*pow(u2,3))/3.)*(HydroGrid[i][j][k].du[2][2]) + (-(SH*pow(tau,-2)) - SH*pow(tau,-2)*pow(u2,2) - SH*pow(u3,2) - (4*SH*pow(u2,2)*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(SH*u0*u2) - (4*SH*tau2*u0*u2*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[3][0]) + (-(SH*u1*u2) - (4*SH*tau2*u1*u2*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[3][1]) + (-SH - SH*pow(u2,2) - SH*tau2*pow(u3,2) - (4*SH*tau2*pow(u2,2)*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((-4*SH*u2*u3)/3. - (4*SH*tau2*u2*pow(u3,3))/3.)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].nspi[9]=  (  ( (-4*SH*u0*pow(tau,-3))/3. - (8*SH*u0*pow(tau,-1)*pow(u3,2))/3. - (4*SH*tau*u0*pow(u3,4))/3. )
		+ ((2*SH*pow(tau,-2))/3. - (2*SH*pow(tau,-2)*pow(u0,2))/3. + (2*SH*pow(u3,2))/3. + (4*SH*pow(u0,2)*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[0][0]) + ((-2*SH*u0*u1*pow(tau,-2))/3. + (4*SH*u0*u1*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[0][1]) + ((-2*SH*u0*u2*pow(tau,-2))/3. + (4*SH*u0*u2*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[0][2]) + ((4*SH*u0*u3*pow(tau,-2))/3. + (4*SH*u0*pow(u3,3))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ ((2*SH*u0*u1*pow(tau,-2))/3. - (4*SH*u0*u1*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[1][0]) + ((2*SH*pow(tau,-2))/3. + (2*SH*pow(tau,-2)*pow(u1,2))/3. + (2*SH*pow(u3,2))/3. - (4*SH*pow(u1,2)*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[1][1]) + ((2*SH*u1*u2*pow(tau,-2))/3. - (4*SH*u1*u2*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((-4*SH*u1*u3*pow(tau,-2))/3. - (4*SH*u1*pow(u3,3))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ ((2*SH*u0*u2*pow(tau,-2))/3. - (4*SH*u0*u2*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[2][0]) + ((2*SH*u1*u2*pow(tau,-2))/3. - (4*SH*u1*u2*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[2][1]) + ((2*SH*pow(tau,-2))/3. + (2*SH*pow(tau,-2)*pow(u2,2))/3. + (2*SH*pow(u3,2))/3. - (4*SH*pow(u2,2)*pow(u3,2))/3.)*(HydroGrid[i][j][k].du[2][2]) + ((-4*SH*u2*u3*pow(tau,-2))/3. - (4*SH*u2*pow(u3,3))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ ((-4*SH*u0*u3)/3. - (4*SH*tau2*u0*pow(u3,3))/3.)*(HydroGrid[i][j][k].du[3][0]) + ((-4*SH*u1*u3)/3. - (4*SH*tau2*u1*pow(u3,3))/3.)*(HydroGrid[i][j][k].du[3][1]) + ((-4*SH*u2*u3)/3. - (4*SH*tau2*u2*pow(u3,3))/3.)*(HydroGrid[i][j][k].du[3][2]) + ((-4*SH*pow(tau,-2))/3. - (8*SH*pow(u3,2))/3. - (4*SH*tau2*pow(u3,4))/3.)*(HydroGrid[i][j][k].du[3][3])  );	
	}
}
