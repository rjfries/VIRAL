

void CalcDer4Vel(GRID HydroGrid, double tau, double tstep)
{
	int i,j,k;
	int l;

	for(l=0;l<VARN;l++)
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
		 
			
			HydroGrid[i][j][k].du[l][3] = genWENOder( qm2, qm1,  qc, qp1, qp2, ZS);
		}
	}	
#endif


}

#ifdef VORT

void AddVorticity(GRID HydroGrid, double tau)
{
	int i,j,k,l;
	 	
		   
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		
		DECLp5u4; 
		DECLePPIa;
	
       
		 HydroGrid[i][j][k].Vort[0]=  (  ( 0 )
		+ (0)*(HydroGrid[i][j][k].du[0][0]) + (-(A2*u0) + p3*u2 + p4*u3 + A2*pow(u0,-1) + A2*pow(u0,-1)*pow(u1,2))*(HydroGrid[i][j][k].du[0][1]) + (-(p3*u1) + A2*u1*u2*pow(u0,-1))*(HydroGrid[i][j][k].du[0][2]) + (-(p4*u1) + A2*u1*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[0][3])
		+ (-(A2*u0) + p3*u2 + p4*u3 + A2*pow(u0,-1) + A2*pow(u0,-1)*pow(u1,2))*(HydroGrid[i][j][k].du[1][0]) + (0)*(HydroGrid[i][j][k].du[1][1]) + (-(A2*u2) + p3*pow(u0,-1) + p4*u2*u3*pow(u0,-1) + p3*pow(u0,-1)*pow(u1,2) + p3*pow(u0,-1)*pow(u2,2))*(HydroGrid[i][j][k].du[1][2]) + (-(A2*u3) + p4*pow(u0,-1) + p3*u2*u3*pow(u0,-1) + p4*pow(u0,-1)*pow(u1,2) + p4*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[1][3])
		+ (-(p3*u1) + A2*u1*u2*pow(u0,-1))*(HydroGrid[i][j][k].du[2][0]) + (A2*u2 - p3*pow(u0,-1) - p4*u2*u3*pow(u0,-1) - p3*pow(u0,-1)*pow(u1,2) - p3*pow(u0,-1)*pow(u2,2))*(HydroGrid[i][j][k].du[2][1]) + (0)*(HydroGrid[i][j][k].du[2][2]) + (p4*u1*u2*pow(u0,-1) - p3*u1*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[2][3])
		+ (-(p4*u1) + A2*u1*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[3][0]) + (A2*u3 - p4*pow(u0,-1) - p3*u2*u3*pow(u0,-1) - p4*pow(u0,-1)*pow(u1,2) - p4*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[3][1]) + (-(p4*u1*u2*pow(u0,-1)) + p3*u1*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[3][2]) + (0)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].Vort[1]=  (  ( 0 )
		+ (0)*(HydroGrid[i][j][k].du[0][0]) + (-(p3*u2) + A3*u1*u2*pow(u0,-1))*(HydroGrid[i][j][k].du[0][1]) + (-(A3*u0) + p3*u1 + p5*u3 + A3*pow(u0,-1) + A3*pow(u0,-1)*pow(u2,2))*(HydroGrid[i][j][k].du[0][2]) + (-(p5*u2) + A3*u2*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[0][3])
		+ (-(p3*u2) + A3*u1*u2*pow(u0,-1))*(HydroGrid[i][j][k].du[1][0]) + (0)*(HydroGrid[i][j][k].du[1][1]) + (A3*u1 - p3*pow(u0,-1) - p5*u1*u3*pow(u0,-1) - p3*pow(u0,-1)*pow(u1,2) - p3*pow(u0,-1)*pow(u2,2))*(HydroGrid[i][j][k].du[1][2]) + (p5*u1*u2*pow(u0,-1) - p3*u2*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[1][3])
		+ (-(A3*u0) + p3*u1 + p5*u3 + A3*pow(u0,-1) + A3*pow(u0,-1)*pow(u2,2))*(HydroGrid[i][j][k].du[2][0]) + (-(A3*u1) + p3*pow(u0,-1) + p5*u1*u3*pow(u0,-1) + p3*pow(u0,-1)*pow(u1,2) + p3*pow(u0,-1)*pow(u2,2))*(HydroGrid[i][j][k].du[2][1]) + (0)*(HydroGrid[i][j][k].du[2][2]) + (-(A3*u3) + p5*pow(u0,-1) + p3*u1*u3*pow(u0,-1) + p5*pow(u0,-1)*pow(u2,2) + p5*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[2][3])
		+ (-(p5*u2) + A3*u2*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[3][0]) + (-(p5*u1*u2*pow(u0,-1)) + p3*u2*u3*pow(u0,-1))*(HydroGrid[i][j][k].du[3][1]) + (A3*u3 - p5*pow(u0,-1) - p3*u1*u3*pow(u0,-1) - p5*pow(u0,-1)*pow(u2,2) - p5*pow(u0,-1)*pow(u3,2))*(HydroGrid[i][j][k].du[3][2]) + (0)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].Vort[2]=  (  ( 0 )
		+ (0)*(HydroGrid[i][j][k].du[0][0]) + (-(A3*u0)/2. - (p1*u2)/2. + (p2*u2)/2. + (p5*u3)/2. + (A3*pow(u0,-1))/2. + (A2*u1*u2*pow(u0,-1))/2. + (A3*pow(u0,-1)*pow(u1,2))/2.)*(HydroGrid[i][j][k].du[0][1]) + (-(A2*u0)/2. + (p1*u1)/2. - (p2*u1)/2. + (p4*u3)/2. + (A2*pow(u0,-1))/2. + (A3*u1*u2*pow(u0,-1))/2. + (A2*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[0][2]) + (-(p5*u1)/2. - (p4*u2)/2. + (A3*u1*u3*pow(u0,-1))/2. + (A2*u2*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(A3*u0)/2. - (p1*u2)/2. + (p2*u2)/2. + (p5*u3)/2. + (A3*pow(u0,-1))/2. + (A2*u1*u2*pow(u0,-1))/2. + (A3*pow(u0,-1)*pow(u1,2))/2.)*(HydroGrid[i][j][k].du[1][0]) + (0)*(HydroGrid[i][j][k].du[1][1]) + ((A2*u1)/2. - (A3*u2)/2. - (p1*pow(u0,-1))/2. + (p2*pow(u0,-1))/2. - (p4*u1*u3*pow(u0,-1))/2. + (p5*u2*u3*pow(u0,-1))/2. - (p1*pow(u0,-1)*pow(u1,2))/2. + (p2*pow(u0,-1)*pow(u1,2))/2. - (p1*pow(u0,-1)*pow(u2,2))/2. + (p2*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[1][2]) + (-(A3*u3)/2. + (p5*pow(u0,-1))/2. + (p4*u1*u2*pow(u0,-1))/2. - (p1*u2*u3*pow(u0,-1))/2. + (p2*u2*u3*pow(u0,-1))/2. + (p5*pow(u0,-1)*pow(u1,2))/2. + (p5*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(A2*u0)/2. + (p1*u1)/2. - (p2*u1)/2. + (p4*u3)/2. + (A2*pow(u0,-1))/2. + (A3*u1*u2*pow(u0,-1))/2. + (A2*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[2][0]) + (-(A2*u1)/2. + (A3*u2)/2. + (p1*pow(u0,-1))/2. - (p2*pow(u0,-1))/2. + (p4*u1*u3*pow(u0,-1))/2. - (p5*u2*u3*pow(u0,-1))/2. + (p1*pow(u0,-1)*pow(u1,2))/2. - (p2*pow(u0,-1)*pow(u1,2))/2. + (p1*pow(u0,-1)*pow(u2,2))/2. - (p2*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[2][1]) + (0)*(HydroGrid[i][j][k].du[2][2]) + (-(A2*u3)/2. + (p4*pow(u0,-1))/2. + (p5*u1*u2*pow(u0,-1))/2. + (p1*u1*u3*pow(u0,-1))/2. - (p2*u1*u3*pow(u0,-1))/2. + (p4*pow(u0,-1)*pow(u2,2))/2. + (p4*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(p5*u1)/2. - (p4*u2)/2. + (A3*u1*u3*pow(u0,-1))/2. + (A2*u2*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[3][0]) + ((A3*u3)/2. - (p5*pow(u0,-1))/2. - (p4*u1*u2*pow(u0,-1))/2. + (p1*u2*u3*pow(u0,-1))/2. - (p2*u2*u3*pow(u0,-1))/2. - (p5*pow(u0,-1)*pow(u1,2))/2. - (p5*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][1]) + ((A2*u3)/2. - (p4*pow(u0,-1))/2. - (p5*u1*u2*pow(u0,-1))/2. - (p1*u1*u3*pow(u0,-1))/2. + (p2*u1*u3*pow(u0,-1))/2. - (p4*pow(u0,-1)*pow(u2,2))/2. - (p4*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][2]) + (0)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].Vort[3]=  (  ( 0 )
		+ (0)*(HydroGrid[i][j][k].du[0][0]) + (-(A4*u0)/2. + (p5*u2)/2. + (A5*u3)/2. - (p1*u3)/2. + (A4*pow(u0,-1))/2. + (A2*u1*u3*pow(u0,-1))/2. + (A4*pow(u0,-1)*pow(u1,2))/2.)*(HydroGrid[i][j][k].du[0][1]) + (-(p5*u1)/2. - (p3*u3)/2. + (A4*u1*u2*pow(u0,-1))/2. + (A2*u2*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[0][2]) + (-(A2*u0)/2. - (A5*u1)/2. + (p1*u1)/2. + (p3*u2)/2. + (A2*pow(u0,-1))/2. + (A4*u1*u3*pow(u0,-1))/2. + (A2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(A4*u0)/2. + (p5*u2)/2. + (A5*u3)/2. - (p1*u3)/2. + (A4*pow(u0,-1))/2. + (A2*u1*u3*pow(u0,-1))/2. + (A4*pow(u0,-1)*pow(u1,2))/2.)*(HydroGrid[i][j][k].du[1][0]) + (0)*(HydroGrid[i][j][k].du[1][1]) + (-(A4*u2)/2. + (p5*pow(u0,-1))/2. + (p3*u1*u3*pow(u0,-1))/2. + (A5*u2*u3*pow(u0,-1))/2. - (p1*u2*u3*pow(u0,-1))/2. + (p5*pow(u0,-1)*pow(u1,2))/2. + (p5*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[1][2]) + ((A2*u1)/2. - (A4*u3)/2. + (A5*pow(u0,-1))/2. - (p1*pow(u0,-1))/2. - (p3*u1*u2*pow(u0,-1))/2. + (p5*u2*u3*pow(u0,-1))/2. + (A5*pow(u0,-1)*pow(u1,2))/2. - (p1*pow(u0,-1)*pow(u1,2))/2. + (A5*pow(u0,-1)*pow(u3,2))/2. - (p1*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(p5*u1)/2. - (p3*u3)/2. + (A4*u1*u2*pow(u0,-1))/2. + (A2*u2*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[2][0]) + ((A4*u2)/2. - (p5*pow(u0,-1))/2. - (p3*u1*u3*pow(u0,-1))/2. - (A5*u2*u3*pow(u0,-1))/2. + (p1*u2*u3*pow(u0,-1))/2. - (p5*pow(u0,-1)*pow(u1,2))/2. - (p5*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[2][1]) + (0)*(HydroGrid[i][j][k].du[2][2]) + ((A2*u2)/2. - (p3*pow(u0,-1))/2. + (A5*u1*u2*pow(u0,-1))/2. - (p1*u1*u2*pow(u0,-1))/2. - (p5*u1*u3*pow(u0,-1))/2. - (p3*pow(u0,-1)*pow(u2,2))/2. - (p3*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(A2*u0)/2. - (A5*u1)/2. + (p1*u1)/2. + (p3*u2)/2. + (A2*pow(u0,-1))/2. + (A4*u1*u3*pow(u0,-1))/2. + (A2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][0]) + (-(A2*u1)/2. + (A4*u3)/2. - (A5*pow(u0,-1))/2. + (p1*pow(u0,-1))/2. + (p3*u1*u2*pow(u0,-1))/2. - (p5*u2*u3*pow(u0,-1))/2. - (A5*pow(u0,-1)*pow(u1,2))/2. + (p1*pow(u0,-1)*pow(u1,2))/2. - (A5*pow(u0,-1)*pow(u3,2))/2. + (p1*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][1]) + (-(A2*u2)/2. + (p3*pow(u0,-1))/2. - (A5*u1*u2*pow(u0,-1))/2. + (p1*u1*u2*pow(u0,-1))/2. + (p5*u1*u3*pow(u0,-1))/2. + (p3*pow(u0,-1)*pow(u2,2))/2. + (p3*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][2]) + (0)*(HydroGrid[i][j][k].du[3][3])  );
			   
		 HydroGrid[i][j][k].Vort[4]=  (  ( 0 )
		+ (0)*(HydroGrid[i][j][k].du[0][0]) + (-(p4*u2)/2. - (p3*u3)/2. + (A4*u1*u2*pow(u0,-1))/2. + (A3*u1*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[0][1]) + (-(A4*u0)/2. + (p4*u1)/2. + (A5*u3)/2. - (p2*u3)/2. + (A4*pow(u0,-1))/2. + (A3*u2*u3*pow(u0,-1))/2. + (A4*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[0][2]) + (-(A3*u0)/2. + (p3*u1)/2. - (A5*u2)/2. + (p2*u2)/2. + (A3*pow(u0,-1))/2. + (A4*u2*u3*pow(u0,-1))/2. + (A3*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(p4*u2)/2. - (p3*u3)/2. + (A4*u1*u2*pow(u0,-1))/2. + (A3*u1*u3*pow(u0,-1))/2.)*(HydroGrid[i][j][k].du[1][0]) + (0)*(HydroGrid[i][j][k].du[1][1]) + ((A4*u1)/2. - (p4*pow(u0,-1))/2. - (A5*u1*u3*pow(u0,-1))/2. + (p2*u1*u3*pow(u0,-1))/2. - (p3*u2*u3*pow(u0,-1))/2. - (p4*pow(u0,-1)*pow(u1,2))/2. - (p4*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[1][2]) + ((A3*u1)/2. - (p3*pow(u0,-1))/2. + (A5*u1*u2*pow(u0,-1))/2. - (p2*u1*u2*pow(u0,-1))/2. - (p4*u2*u3*pow(u0,-1))/2. - (p3*pow(u0,-1)*pow(u1,2))/2. - (p3*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(A4*u0)/2. + (p4*u1)/2. + (A5*u3)/2. - (p2*u3)/2. + (A4*pow(u0,-1))/2. + (A3*u2*u3*pow(u0,-1))/2. + (A4*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[2][0]) + (-(A4*u1)/2. + (p4*pow(u0,-1))/2. + (A5*u1*u3*pow(u0,-1))/2. - (p2*u1*u3*pow(u0,-1))/2. + (p3*u2*u3*pow(u0,-1))/2. + (p4*pow(u0,-1)*pow(u1,2))/2. + (p4*pow(u0,-1)*pow(u2,2))/2.)*(HydroGrid[i][j][k].du[2][1]) + (0)*(HydroGrid[i][j][k].du[2][2]) + ((A3*u2)/2. - (A4*u3)/2. + (A5*pow(u0,-1))/2. - (p2*pow(u0,-1))/2. - (p3*u1*u2*pow(u0,-1))/2. + (p4*u1*u3*pow(u0,-1))/2. + (A5*pow(u0,-1)*pow(u2,2))/2. - (p2*pow(u0,-1)*pow(u2,2))/2. + (A5*pow(u0,-1)*pow(u3,2))/2. - (p2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(A3*u0)/2. + (p3*u1)/2. - (A5*u2)/2. + (p2*u2)/2. + (A3*pow(u0,-1))/2. + (A4*u2*u3*pow(u0,-1))/2. + (A3*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][0]) + (-(A3*u1)/2. + (p3*pow(u0,-1))/2. - (A5*u1*u2*pow(u0,-1))/2. + (p2*u1*u2*pow(u0,-1))/2. + (p4*u2*u3*pow(u0,-1))/2. + (p3*pow(u0,-1)*pow(u1,2))/2. + (p3*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][1]) + (-(A3*u2)/2. + (A4*u3)/2. - (A5*pow(u0,-1))/2. + (p2*pow(u0,-1))/2. + (p3*u1*u2*pow(u0,-1))/2. - (p4*u1*u3*pow(u0,-1))/2. - (A5*pow(u0,-1)*pow(u2,2))/2. + (p2*pow(u0,-1)*pow(u2,2))/2. - (A5*pow(u0,-1)*pow(u3,2))/2. + (p2*pow(u0,-1)*pow(u3,2))/2.)*(HydroGrid[i][j][k].du[3][2]) + (0)*(HydroGrid[i][j][k].du[3][3])  );

	
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
		   
		   
//Geometical Source terms
 
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{

		DECLePPIa;
		DECLp5u4;		
		
		HydroGrid[i][j][k].Source[0] = 0;
		HydroGrid[i][j][k].Source[1] = 0;
		HydroGrid[i][j][k].Source[2] = 0;
		HydroGrid[i][j][k].Source[3] = 0;
	
	}
	
#ifdef SHAS
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
				
			q0[i-t+2] = -(           0					  );
			q1[i-t+2] = -(P-PI - (u1/u0)*A2         + p1  );
			q2[i-t+2] = -(     - (u1/u0)*A3         + p3  );
			q3[i-t+2] = -(     - (u1/u0)*A4         + p4  );
		}
		i=t;
		 
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
				
			q0[j-t+2] = -(       0					 	  );
			q1[j-t+2] = -(	   - (u2/u0)*A2         + p3  );
			q2[j-t+2] = -(P-PI - (u2/u0)*A3         + p2  );
			q3[j-t+2] = -(     - (u2/u0)*A4         + p5  );
		}
		j=t;
		 
		
		double temp0 = genWENOder(q0[0],q0[1],q0[2],q0[3],q0[4],YS);
		double temp1 = genWENOder(q1[0],q1[1],q1[2],q1[3],q1[4],YS);
		double temp2 = genWENOder(q2[0],q2[1],q2[2],q2[3],q2[4],YS);
		double temp3 = genWENOder(q3[0],q3[1],q3[2],q3[3],q3[4],YS);
		
		HydroGrid[i][j][k].Source[0] += temp0;
		HydroGrid[i][j][k].Source[1] += temp1;
		HydroGrid[i][j][k].Source[2] += temp2;
		HydroGrid[i][j][k].Source[3] += temp3;	
	}
	
#if !defined LBI	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{		
		double q0[5];
		double q1[5];
		double q2[5]; 
		double q3[5];
		
		int t=k;
		for(k=t-2;k<=t+2;k++)
		{ 
			DECLePPIa;
			DECLp10u4;	
				
			q0[k-t+2] = -(       0					 	  );
			q1[k-t+2] = -(	   - (u3/u0)*A2         + p4  );
			q2[k-t+2] = -(	   - (u3/u0)*A3         + p5  );
			q3[k-t+2] = -(P-PI - (u3/u0)*A4         + A5 );
		}
		k=t;
		double temp0 = genWENOder(q0[0],q0[1],q0[2],q0[3],q0[4],ZS);
		double temp1 = genWENOder(q1[0],q1[1],q1[2],q1[3],q1[4],ZS);
		double temp2 = genWENOder(q2[0],q2[1],q2[2],q2[3],q2[4],ZS);
		double temp3 = genWENOder(q3[0],q3[1],q3[2],q3[3],q3[4],ZS);
		
		HydroGrid[i][j][k].Source[0] += temp0;
		HydroGrid[i][j][k].Source[1] += temp1;
		HydroGrid[i][j][k].Source[2] += temp2;
		HydroGrid[i][j][k].Source[3] += temp3;	
	}
#endif


#endif


	CalcDer4Vel(HydroGrid, tau, ts);   
	
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=0;k<ZCMA;k++)
	{
		
		DECLp5u4;
		DECLcoord;
		DECLePPIa; 
			
		

		double T = FT(e );	
		double s = FS(e, P, T);		
		double SH =   Feta(  s, e );			
		double tpi= Ftaupi(  SH, P, e );	
   
#ifdef BULK
		double BU =  FZeta(  s, e );
		double tPI= FtauPI(  BU, P, e );	
#endif
		       
		 HydroGrid[i][j][k].Source[4]=  (  ( -(p1*pow(tpi,-1)*pow(u0,-1)) )
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-2*SH*u0 + 4*p1*tpi*u0 + 6*A2*tpi*u1*pow(u0,2) + 2*SH*pow(u0,3) - 2*SH*u0*pow(u1,2) - 4*SH*pow(u0,3)*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*p1*tpi*u1 - 4*SH*u1*pow(u0,2) + 6*A2*tpi*u0*pow(u1,2) - 4*SH*pow(u0,2)*pow(u1,3)))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*p1*tpi*u2 + 6*A2*tpi*u0*u1*u2 + 2*SH*u2*pow(u0,2) - 4*SH*u2*pow(u0,2)*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[0][2]) + (-(u3*pow(tpi,-1)*pow(u0,-2)*(3*p1*tpi + 2*u0*(SH*u0 + 3*A2*tpi*u1 - 2*SH*u0*pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u1*pow(u0,2) - 6*p1*tpi*u1*pow(u0,2) + 4*SH*pow(u0,2)*pow(u1,3)))/3.)*(HydroGrid[i][j][k].du[1][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u0 + p1*tpi*u0 + 8*SH*u0*pow(u1,2) - 6*p1*tpi*u0*pow(u1,2) + 4*SH*u0*pow(u1,4)))/3.)*(HydroGrid[i][j][k].du[1][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u0*u1*u2 - 6*p1*tpi*u0*u1*u2 + 4*SH*u0*u2*pow(u1,3)))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((-2*u1*u3*pow(tpi,-1)*pow(u0,-1)*(-3*p1*tpi + 2*SH*(1 + pow(u1,2))))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-6*p3*tpi*u1*pow(u0,2) - 2*SH*u2*pow(u0,2) + 4*SH*u2*pow(u0,2)*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[2][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u0*u1*u2 - 6*p3*tpi*u0*pow(u1,2) + 4*SH*u0*u2*pow(u1,3)))/3.)*(HydroGrid[i][j][k].du[2][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(-2*SH*u0 + p1*tpi*u0 - 6*p3*tpi*u0*u1*u2 - 2*SH*u0*pow(u1,2) - 2*SH*u0*pow(u2,2) + 4*SH*u0*pow(u1,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(-6*p3*tpi*u0*u1*u3 - 2*SH*u0*u2*u3 + 4*SH*u0*u2*u3*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-6*p4*tpi*u1*pow(u0,2) - 2*SH*u3*pow(u0,2) + 4*SH*u3*pow(u0,2)*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[3][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u0*u1*u3 - 6*p4*tpi*u0*pow(u1,2) + 4*SH*u0*u3*pow(u1,3)))/3.)*(HydroGrid[i][j][k].du[3][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(-6*p4*tpi*u0*u1*u2 - 2*SH*u0*u2*u3 + 4*SH*u0*u2*u3*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[3][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(-2*SH*u0 + p1*tpi*u0 - 6*p4*tpi*u0*u1*u3 - 2*SH*u0*pow(u1,2) - 2*SH*u0*pow(u3,2) + 4*SH*u0*pow(u1,2)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
		 HydroGrid[i][j][k].Source[5]=  (  ( -(p2*pow(tpi,-1)*pow(u0,-1)) )
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-2*SH*u0 + 4*p2*tpi*u0 + 6*A3*tpi*u2*pow(u0,2) + 2*SH*pow(u0,3) - 2*SH*u0*pow(u2,2) - 4*SH*pow(u0,3)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*p2*tpi*u1 + 6*A3*tpi*u0*u1*u2 + 2*SH*u1*pow(u0,2) - 4*SH*u1*pow(u0,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*p2*tpi*u2 - 4*SH*u2*pow(u0,2) + 6*A3*tpi*u0*pow(u2,2) - 4*SH*pow(u0,2)*pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[0][2]) + (-(u3*pow(tpi,-1)*pow(u0,-2)*(3*p2*tpi + 2*u0*(SH*u0 + 3*A3*tpi*u2 - 2*SH*u0*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-2*SH*u1*pow(u0,2) - 6*p3*tpi*u2*pow(u0,2) + 4*SH*u1*pow(u0,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[1][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(-2*SH*u0 + p2*tpi*u0 - 6*p3*tpi*u0*u1*u2 - 2*SH*u0*pow(u1,2) - 2*SH*u0*pow(u2,2) + 4*SH*u0*pow(u1,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[1][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u0*u1*u2 - 6*p3*tpi*u0*pow(u2,2) + 4*SH*u0*u1*pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[1][2]) + ((2*u3*pow(tpi,-1)*pow(u0,-1)*(3*p3*tpi*u2 + SH*(u1 - 2*u1*pow(u2,2))))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u2*pow(u0,2) - 6*p2*tpi*u2*pow(u0,2) + 4*SH*pow(u0,2)*pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[2][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u0*u1*u2 - 6*p2*tpi*u0*u1*u2 + 4*SH*u0*u1*pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[2][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u0 + p2*tpi*u0 + 8*SH*u0*pow(u2,2) - 6*p2*tpi*u0*pow(u2,2) + 4*SH*u0*pow(u2,4)))/3.)*(HydroGrid[i][j][k].du[2][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u0*u2*u3 - 6*p2*tpi*u0*u2*u3 + 4*SH*u0*u3*pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-6*p5*tpi*u2*pow(u0,2) - 2*SH*u3*pow(u0,2) + 4*SH*u3*pow(u0,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[3][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(-6*p5*tpi*u0*u1*u2 - 2*SH*u0*u1*u3 + 4*SH*u0*u1*u3*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[3][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(4*SH*u0*u2*u3 - 6*p5*tpi*u0*pow(u2,2) + 4*SH*u0*u3*pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[3][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(-2*SH*u0 + p2*tpi*u0 - 6*p5*tpi*u0*u2*u3 - 2*SH*u0*pow(u2,2) - 2*SH*u0*pow(u3,2) + 4*SH*u0*pow(u2,2)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
		 HydroGrid[i][j][k].Source[6]=  (  ( -(p3*pow(tpi,-1)*pow(u0,-1)) )
		+ (-(pow(tpi,-1)*pow(u0,-2)*(4*p3*tpi*u0 - 2*SH*u0*u1*u2 + 3*A3*tpi*u1*pow(u0,2) + 3*A2*tpi*u2*pow(u0,2) - 4*SH*u1*u2*pow(u0,3)))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*p3*tpi*u1 + 3*A2*tpi*u0*u1*u2 - 3*SH*u2*pow(u0,2) + 3*A3*tpi*u0*pow(u1,2) - 4*SH*u2*pow(u0,2)*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*p3*tpi*u2 + 3*A3*tpi*u0*u1*u2 - 3*SH*u1*pow(u0,2) + 3*A2*tpi*u0*pow(u2,2) - 4*SH*u1*pow(u0,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[0][2]) + (-((3*p3*tpi + u0*(3*A3*tpi*u1 + 3*A2*tpi*u2 - 4*SH*u0*u1*u2))*u3*pow(tpi,-1)*pow(u0,-2))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-3*p3*tpi*u1*pow(u0,2) + 3*SH*u2*pow(u0,2) - 3*p1*tpi*u2*pow(u0,2) + 4*SH*u2*pow(u0,2)*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(p3*tpi*u0 + 4*SH*u0*u1*u2 - 3*p1*tpi*u0*u1*u2 - 3*p3*tpi*u0*pow(u1,2) + 4*SH*u0*u2*pow(u1,3)))/3.)*(HydroGrid[i][j][k].du[1][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*SH*u0 - 3*p3*tpi*u0*u1*u2 + 3*SH*u0*pow(u1,2) + 3*SH*u0*pow(u2,2) - 3*p1*tpi*u0*pow(u2,2) + 4*SH*u0*pow(u1,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[1][2]) + (-(u3*pow(tpi,-1)*pow(u0,-1)*(-3*p3*tpi*u1 + 3*SH*u2 - 3*p1*tpi*u2 + 4*SH*u2*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(3*SH*u1*pow(u0,2) - 3*p2*tpi*u1*pow(u0,2) - 3*p3*tpi*u2*pow(u0,2) + 4*SH*u1*pow(u0,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*SH*u0 - 3*p3*tpi*u0*u1*u2 + 3*SH*u0*pow(u1,2) - 3*p2*tpi*u0*pow(u1,2) + 3*SH*u0*pow(u2,2) + 4*SH*u0*pow(u1,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(p3*tpi*u0 + 4*SH*u0*u1*u2 - 3*p2*tpi*u0*u1*u2 - 3*p3*tpi*u0*pow(u2,2) + 4*SH*u0*u1*pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[2][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*SH*u0*u1*u3 - 3*p2*tpi*u0*u1*u3 - 3*p3*tpi*u0*u2*u3 + 4*SH*u0*u1*u3*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-3*p5*tpi*u1*pow(u0,2) - 3*p4*tpi*u2*pow(u0,2) + 4*SH*u1*u2*u3*pow(u0,2)))/3.)*(HydroGrid[i][j][k].du[3][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(-3*p4*tpi*u0*u1*u2 + 3*SH*u0*u2*u3 - 3*p5*tpi*u0*pow(u1,2) + 4*SH*u0*u2*u3*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[3][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(-3*p5*tpi*u0*u1*u2 + 3*SH*u0*u1*u3 - 3*p4*tpi*u0*pow(u2,2) + 4*SH*u0*u1*u3*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[3][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(p3*tpi*u0 - 2*SH*u0*u1*u2 - 3*p5*tpi*u0*u1*u3 - 3*p4*tpi*u0*u2*u3 + 4*SH*u0*u1*u2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
		 HydroGrid[i][j][k].Source[7]=  (  ( -(p4*pow(tpi,-1)*pow(u0,-1)) )
		+ (-(pow(tpi,-1)*pow(u0,-2)*(4*p4*tpi*u0 - 2*SH*u0*u1*u3 + 3*A4*tpi*u1*pow(u0,2) + 3*A2*tpi*u3*pow(u0,2) - 4*SH*u1*u3*pow(u0,3)))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*p4*tpi*u1 + 3*A2*tpi*u0*u1*u3 - 3*SH*u3*pow(u0,2) + 3*A4*tpi*u0*pow(u1,2) - 4*SH*u3*pow(u0,2)*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*p4*tpi*u2 + 3*A4*tpi*u0*u1*u2 + 3*A2*tpi*u0*u2*u3 - 4*SH*u1*u2*u3*pow(u0,2)))/3.)*(HydroGrid[i][j][k].du[0][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*tpi*u3*(p4 + A4*u0*u1 + A2*u0*u3) - SH*u1*pow(u0,2)*(3 + 4*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-3*p4*tpi*u1*pow(u0,2) + 3*SH*u3*pow(u0,2) - 3*p1*tpi*u3*pow(u0,2) + 4*SH*u3*pow(u0,2)*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(p4*tpi*u0 + 4*SH*u0*u1*u3 - 3*p1*tpi*u0*u1*u3 - 3*p4*tpi*u0*pow(u1,2) + 4*SH*u0*u3*pow(u1,3)))/3.)*(HydroGrid[i][j][k].du[1][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(-3*p4*tpi*u0*u1*u2 + 3*SH*u0*u2*u3 - 3*p1*tpi*u0*u2*u3 + 4*SH*u0*u2*u3*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][2]) + (-(pow(tpi,-1)*pow(u0,-1)*(-3*tpi*u3*(p4*u1 + p1*u3) + SH*(3*(1 + pow(u3,2)) + pow(u1,2)*(3 + 4*pow(u3,2)))))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-3*p5*tpi*u1*pow(u0,2) - 3*p3*tpi*u3*pow(u0,2) + 4*SH*u1*u2*u3*pow(u0,2)))/3.)*(HydroGrid[i][j][k].du[2][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(-3*p3*tpi*u0*u1*u3 + 3*SH*u0*u2*u3 - 3*p5*tpi*u0*pow(u1,2) + 4*SH*u0*u2*u3*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[2][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(p4*tpi*u0 - 3*p5*tpi*u0*u1*u2 - 2*SH*u0*u1*u3 - 3*p3*tpi*u0*u2*u3 + 4*SH*u0*u1*u3*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*SH*u0*u1*u2 - 3*p5*tpi*u0*u1*u3 - 3*p3*tpi*u0*pow(u3,2) + 4*SH*u0*u1*u2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(3*SH*u1*pow(u0,2) - 3*A5*tpi*u1*pow(u0,2) - 3*p4*tpi*u3*pow(u0,2) + 4*SH*u1*pow(u0,2)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*SH*u0 - 3*p4*tpi*u0*u1*u3 + 3*SH*u0*pow(u1,2) - 3*A5*tpi*u0*pow(u1,2) + 3*SH*u0*pow(u3,2) + 4*SH*u0*pow(u1,2)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*SH*u0*u1*u2 - 3*A5*tpi*u0*u1*u2 - 3*p4*tpi*u0*u2*u3 + 4*SH*u0*u1*u2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(p4*tpi*u0 + 4*SH*u0*u1*u3 - 3*A5*tpi*u0*u1*u3 - 3*p4*tpi*u0*pow(u3,2) + 4*SH*u0*u1*pow(u3,3)))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);
			   
		 HydroGrid[i][j][k].Source[8]=  (  ( -(p5*pow(tpi,-1)*pow(u0,-1)) )
		+ (-(pow(tpi,-1)*pow(u0,-2)*(4*p5*tpi*u0 - 2*SH*u0*u2*u3 + 3*A4*tpi*u2*pow(u0,2) + 3*A3*tpi*u3*pow(u0,2) - 4*SH*u2*u3*pow(u0,3)))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*p5*tpi*u1 + 3*A4*tpi*u0*u1*u2 + 3*A3*tpi*u0*u1*u3 - 4*SH*u1*u2*u3*pow(u0,2)))/3.)*(HydroGrid[i][j][k].du[0][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*p5*tpi*u2 + 3*A3*tpi*u0*u2*u3 - 3*SH*u3*pow(u0,2) + 3*A4*tpi*u0*pow(u2,2) - 4*SH*u3*pow(u0,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[0][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*tpi*u3*(p5 + A4*u0*u2 + A3*u0*u3) - SH*u2*pow(u0,2)*(3 + 4*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[0][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-3*p4*tpi*u2*pow(u0,2) - 3*p3*tpi*u3*pow(u0,2) + 4*SH*u1*u2*u3*pow(u0,2)))/3.)*(HydroGrid[i][j][k].du[1][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(p5*tpi*u0 - 3*p4*tpi*u0*u1*u2 - 3*p3*tpi*u0*u1*u3 - 2*SH*u0*u2*u3 + 4*SH*u0*u2*u3*pow(u1,2)))/3.)*(HydroGrid[i][j][k].du[1][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*SH*u0*u1*u3 - 3*p3*tpi*u0*u2*u3 - 3*p4*tpi*u0*pow(u2,2) + 4*SH*u0*u1*u3*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[1][2]) + (-(pow(tpi,-1)*pow(u0,-1)*(-3*tpi*u3*(p4*u2 + p3*u3) + SH*u1*u2*(3 + 4*pow(u3,2))))/3.)*(HydroGrid[i][j][k].du[1][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(-3*p5*tpi*u2*pow(u0,2) + 3*SH*u3*pow(u0,2) - 3*p2*tpi*u3*pow(u0,2) + 4*SH*u3*pow(u0,2)*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(-3*p5*tpi*u0*u1*u2 + 3*SH*u0*u1*u3 - 3*p2*tpi*u0*u1*u3 + 4*SH*u0*u1*u3*pow(u2,2)))/3.)*(HydroGrid[i][j][k].du[2][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(p5*tpi*u0 + 4*SH*u0*u2*u3 - 3*p2*tpi*u0*u2*u3 - 3*p5*tpi*u0*pow(u2,2) + 4*SH*u0*u3*pow(u2,3)))/3.)*(HydroGrid[i][j][k].du[2][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*SH*u0 - 3*p5*tpi*u0*u2*u3 + 3*SH*u0*pow(u2,2) + 3*SH*u0*pow(u3,2) - 3*p2*tpi*u0*pow(u3,2) + 4*SH*u0*pow(u2,2)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[2][3])
		+ (-(pow(tpi,-1)*pow(u0,-2)*(3*SH*u2*pow(u0,2) - 3*A5*tpi*u2*pow(u0,2) - 3*p5*tpi*u3*pow(u0,2) + 4*SH*u2*pow(u0,2)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][0]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*SH*u0*u1*u2 - 3*A5*tpi*u0*u1*u2 - 3*p5*tpi*u0*u1*u3 + 4*SH*u0*u1*u2*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][1]) + (-(pow(tpi,-1)*pow(u0,-2)*(3*SH*u0 - 3*p5*tpi*u0*u2*u3 + 3*SH*u0*pow(u2,2) - 3*A5*tpi*u0*pow(u2,2) + 3*SH*u0*pow(u3,2) + 4*SH*u0*pow(u2,2)*pow(u3,2)))/3.)*(HydroGrid[i][j][k].du[3][2]) + (-(pow(tpi,-1)*pow(u0,-2)*(p5*tpi*u0 + 4*SH*u0*u2*u3 - 3*A5*tpi*u0*u2*u3 - 3*p5*tpi*u0*pow(u3,2) + 4*SH*u0*u2*pow(u3,3)))/3.)*(HydroGrid[i][j][k].du[3][3])  
		);


#ifdef BULK		        
	 HydroGrid[i][j][k].Source[9]=  (  ( -(PI*pow(tPI,-1)*pow(u0,-1)) )
	+ (-((3*BU + 4*PI*tPI)*pow(tPI,-1)*pow(u0,-1))/3.)*(HydroGrid[i][j][k].du[0][0]) + (-(PI*u1*pow(u0,-2)))*(HydroGrid[i][j][k].du[0][1]) + (-(PI*u2*pow(u0,-2)))*(HydroGrid[i][j][k].du[0][2]) + (-(PI*u3*pow(u0,-2)))*(HydroGrid[i][j][k].du[0][3])
	+ (0)*(HydroGrid[i][j][k].du[1][0]) + (-((3*BU + PI*tPI)*pow(tPI,-1)*pow(u0,-1))/3.)*(HydroGrid[i][j][k].du[1][1]) + (0)*(HydroGrid[i][j][k].du[1][2]) + (0)*(HydroGrid[i][j][k].du[1][3])
	+ (0)*(HydroGrid[i][j][k].du[2][0]) + (0)*(HydroGrid[i][j][k].du[2][1]) + (-((3*BU + PI*tPI)*pow(tPI,-1)*pow(u0,-1))/3.)*(HydroGrid[i][j][k].du[2][2]) + (0)*(HydroGrid[i][j][k].du[2][3])
	+ (0)*(HydroGrid[i][j][k].du[3][0]) + (0)*(HydroGrid[i][j][k].du[3][1]) + (0)*(HydroGrid[i][j][k].du[3][2]) + (-((3*BU + PI*tPI)*pow(tPI,-1)*pow(u0,-1))/3.)*(HydroGrid[i][j][k].du[3][3]));
#else	        
		HydroGrid[i][j][k].Source[9]= 0;
#endif


	}
	
		
#ifdef VORT
	AddVorticity(HydroGrid,tau);
#endif 
}

 
  
