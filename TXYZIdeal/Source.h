
void CalcSource(GRID HydroGrid, double tau, double ts)
{
	int i,j,k;
	double max1=0,min1=0; 	
		   
		   
//Geometical Source terms
 
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{

		DECLePa;
		DECLu4;		
		
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
			DECLePa;
			DECLu4;	
				
			q0[i-t+2] = -(           0					  );
			q1[i-t+2] = -(           P 					  );
			q2[i-t+2] = -(           0					  );
			q3[i-t+2] = -(           0					  );
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
			DECLePa;
			DECLu4;	
				
			q0[j-t+2] = -(           0		   		 	  );
			q1[j-t+2] = -(           0					  );
			q2[j-t+2] = -(         	 P					  );
			q3[j-t+2] = -(           0					  );
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
			DECLePa;
			DECLu4;	
				
			q0[k-t+2] = -(           0   			 	  );
			q1[k-t+2] = -(           0					  );
			q2[k-t+2] = -(           0					  );
			q3[k-t+2] = -(           P					  );
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


}
  
