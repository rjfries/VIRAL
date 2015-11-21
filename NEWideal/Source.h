void UpdatePrimaryVariablesAtEndOfTimeStep(GRID HydroGrid, double ts)
{
	int i,j,k,l;

	double tweaksource=1;
	double tweakresult=1;
		
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		double X = HydroGrid[i][j][k].X;
		double Y = HydroGrid[i][j][k].Y;
		double r = HydroGrid[i][j][k].r;
		
		HydroGrid[i][j][k].T00 =  (tau*HydroGrid[i][j][k].T00 
								+ HydroGrid[i][j][k].Result[0]
								+ HydroGrid[i][j][k].Source[0]*ts)/(tau+ts);
								
		HydroGrid[i][j][k].T10 =  (tau*HydroGrid[i][j][k].T10 
								+ HydroGrid[i][j][k].Result[1])/(tau+ts);
								
		HydroGrid[i][j][k].T20 =  (tau*HydroGrid[i][j][k].T20 
								+ HydroGrid[i][j][k].Result[2])/(tau+ts);
								
		HydroGrid[i][j][k].T30 =  (tau*HydroGrid[i][j][k].T30 
								+ HydroGrid[i][j][k].Result[3]
								+ HydroGrid[i][j][k].Source[3]*ts)/(tau+ts);	

	}
	
		//~ int temp=0;
		//~ 
		//~ for(i=il;i<ir;i++)
		//~ for(j=jl;j<jr;j++)
		//~ for(k=kl;k<kr;k++)
		//~ {
			//~ if(fabs(HydroGrid[i][j][k].T10) > 1e-17)
			//~ { 
	//~ 
				//~ double e = HydroGrid[i][j][k].En;
				//~ double P = HydroGrid[i][j][k].P;
				//~ double Vx = HydroGrid[i][j][k].Vx;
				//~ double Vy = HydroGrid[i][j][k].Vy;
				//~ double Ve = HydroGrid[i][j][k].Ve;
//~ 
				//~ double u0 = pow(1-Vx*Vx-Vy*Vy-Ve*Ve, -0.5);
				//~ double u1 = u0*Vx;
				//~ double u2 = u0*Vy;
				//~ double u3 = u0*Ve; 
				//~ 
			//~ //	cout<<std::scientific<<"i j k  "<<i<< " "<<j<<" "<<k<<"  "<<HydroGrid[i][j][k].T10<<endl;
				//~ temp=1;
			//~ }
		//~ }
		//~ 
		//~ if (temp)
			//~ exit(1);
}



void CalcCentreFlux(GRID HydroGrid, double tau, double ts)
{
	int i,j,k;
	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		double e = HydroGrid[i][j][k].En;
		double P = HydroGrid[i][j][k].P;
		double Vx = HydroGrid[i][j][k].Vx;
		double Vy = HydroGrid[i][j][k].Vy;
		double Ve = HydroGrid[i][j][k].Ve;
		
		double u0 = pow(1-Vx*Vx-Vy*Vy-tau*tau*Ve*Ve, -0.5);
		double u1 = u0*Vx;
		double u2 = u0*Vy;
		double u3 = u0*Ve;
		
		HydroGrid[i][j][k].Fx[0] =  tau*(   u1*u0*(e+P)       );
		HydroGrid[i][j][k].Fx[1] =  tau*(   P + u1*u1*(e+P)   );
		HydroGrid[i][j][k].Fx[2] =  tau*(   u1*u2*(e+P)       );
		HydroGrid[i][j][k].Fx[3] =  tau*(   u1*u3*(e+P)       );
		
		HydroGrid[i][j][k].Fy[0] =  tau*(   u2*u0*(e+P)       );
		HydroGrid[i][j][k].Fy[1] =  tau*(   u2*u1*(e+P)       );
		HydroGrid[i][j][k].Fy[2] =  tau*(   P + u2*u2*(e+P)   );
		HydroGrid[i][j][k].Fy[3] =  tau*(   u2*u3*(e+P)       );
		
		HydroGrid[i][j][k].Fe[0] =  tau*(   u3*u0*(e+P)                     );
		HydroGrid[i][j][k].Fe[1] =  tau*(   u3*u1*(e+P)                     );
		HydroGrid[i][j][k].Fe[2] =  tau*(   u3*u2*(e+P)                     );
		HydroGrid[i][j][k].Fe[3] =  tau*(   (P/(tau*tau)) + u3*u3*(e+P)     );		
					
	}		
}



void CalcSource(GRID HydroGrid, double tau, double ts)
{
	int i,j,k;
	double max1=0,min1=0;
	double tau2 = tau*tau;
		

//Source Term for T00	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		double En = HydroGrid[i][j][k].En;
		double P = HydroGrid[i][j][k].P;
		double Vx = HydroGrid[i][j][k].Vx;
		double Vy = HydroGrid[i][j][k].Vy;
		double Ve = HydroGrid[i][j][k].Ve;
		
		double u0 = pow(1 - pow(Vx,2) - pow(Vy,2) - pow(tau*Ve,2), -0.5);
		double u3 = u0*Ve;
		
		HydroGrid[i][j][k].Source[0] =  -(P) - tau2*( (En + P )*u3*u3);			
	}
	
//Source Term for T30	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{

		double En = HydroGrid[i][j][k].En;			
		double P = HydroGrid[i][j][k].P;
		double Vx = HydroGrid[i][j][k].Vx;
		double Vy = HydroGrid[i][j][k].Vy;
		double Ve = HydroGrid[i][j][k].Ve;
			
		double u0 = pow(1 - pow(Vx,2) - pow(Vy,2) - pow(Ve,2), -0.5);
		double u3 = u0*Ve;
		
		HydroGrid[i][j][k].Source[3] = ( -2.0*(  (En + P )*u3*u3) );
	}
}
