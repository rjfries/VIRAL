void XYBoundaryCopy(GRID HydroGrid)
{

	
	int i,j,k,l;
	double x[2],y[2];
	
	if(MYXLEFT == BOUNDARY)
	{
		for(j=jl;j<jr;j++)
		for(k=kl;k<kr;k++)
		for(i=0;i<BORDER;i++)
		{

			HydroGrid[i][j][k].T00 = HydroGrid[il][j][k].T00;
			HydroGrid[i][j][k].T10 = HydroGrid[il][j][k].T10;
			HydroGrid[i][j][k].T20 = HydroGrid[il][j][k].T20;
			HydroGrid[i][j][k].T30 = HydroGrid[il][j][k].T30;

			
			HydroGrid[i][j][k].En = HydroGrid[il][j][k].En;
			HydroGrid[i][j][k].Vx = HydroGrid[il][j][k].Vx;
			HydroGrid[i][j][k].Vy = HydroGrid[il][j][k].Vy;
			HydroGrid[i][j][k].Ve = HydroGrid[il][j][k].Ve;
			HydroGrid[i][j][k].P = HydroGrid[il][j][k].P;
			
			HydroGrid[i][j][k].u[0] = HydroGrid[il][j][k].u[0];
			HydroGrid[i][j][k].u[1] = HydroGrid[il][j][k].u[1];
			HydroGrid[i][j][k].u[2] = HydroGrid[il][j][k].u[2];
			HydroGrid[i][j][k].u[3] = HydroGrid[il][j][k].u[3];
			
#ifdef GUBSER

			HydroGrid[i][j][k].u[1] = HydroGrid[il+BORDER-i-1][j][k].u[1];
			HydroGrid[i][j][k].u[2] = HydroGrid[il+BORDER-i-1][j][k].u[2];

#endif	 
			HydroGrid[i][j][k].prevu[0] = HydroGrid[il][j][k].prevu[0];
			HydroGrid[i][j][k].prevu[1] = HydroGrid[il][j][k].prevu[1];
			HydroGrid[i][j][k].prevu[2] = HydroGrid[il][j][k].prevu[2];
			HydroGrid[i][j][k].prevu[3] = HydroGrid[il][j][k].prevu[3];
			
			for(l=0;l<Npi;l++)
				HydroGrid[i][j][k].pi[l] = HydroGrid[il][j][k].pi[l];
			
			HydroGrid[i][j][k].PI = HydroGrid[il][j][k].PI;
			

		}
	}
	
	if(MYXRIGHT == BOUNDARY)
	{	
		for(j=jl;j<jr;j++)
		for(k=kl;k<kr;k++)
		for(i=0;i<BORDER;i++)
		{
			HydroGrid[ir+i][j][k].T00 = HydroGrid[ir-1][j][k].T00;
			HydroGrid[ir+i][j][k].T10 = HydroGrid[ir-1][j][k].T10;
			HydroGrid[ir+i][j][k].T20 = HydroGrid[ir-1][j][k].T20;
			HydroGrid[ir+i][j][k].T30 = HydroGrid[ir-1][j][k].T30;



			
			HydroGrid[ir+i][j][k].En = HydroGrid[ir-1][j][k].En;
			HydroGrid[ir+i][j][k].Vx = HydroGrid[ir-1][j][k].Vx;
			HydroGrid[ir+i][j][k].Vy = HydroGrid[ir-1][j][k].Vy;
			HydroGrid[ir+i][j][k].Ve = HydroGrid[ir-1][j][k].Ve;
			HydroGrid[ir+i][j][k].P = HydroGrid[ir-1][j][k].P;
			
			HydroGrid[ir+i][j][k].u[0] = HydroGrid[ir-1][j][k].u[0];
			HydroGrid[ir+i][j][k].u[1] = HydroGrid[ir-1][j][k].u[1];
			HydroGrid[ir+i][j][k].u[2] = HydroGrid[ir-1][j][k].u[2];
			HydroGrid[ir+i][j][k].u[3] = HydroGrid[ir-1][j][k].u[3];	
			
#ifdef GUBSER
			HydroGrid[ir+i][j][k].u[1] = HydroGrid[ir-i-1][j][k].u[1];
			HydroGrid[ir+i][j][k].u[2] = HydroGrid[ir-i-1][j][k].u[2];

#endif
			HydroGrid[ir+i][j][k].prevu[0] = HydroGrid[ir-1][j][k].prevu[0];
			HydroGrid[ir+i][j][k].prevu[1] = HydroGrid[ir-1][j][k].prevu[1];
			HydroGrid[ir+i][j][k].prevu[2] = HydroGrid[ir-1][j][k].prevu[2];
			HydroGrid[ir+i][j][k].prevu[3] = HydroGrid[ir-1][j][k].prevu[3];	

			for(l=0;l<Npi;l++)
				HydroGrid[ir+i][j][k].pi[l] = HydroGrid[ir-1][j][k].pi[l];
			
			HydroGrid[ir+i][j][k].PI = HydroGrid[ir-1][j][k].PI;
		}
	}

	if(MYYLEFT == BOUNDARY)
	{	
		for(i=il;i<ir;i++)
		for(k=kl;k<kr;k++)
		for(j=0;j<BORDER;j++)
		{
			HydroGrid[i][j][k].T00 = HydroGrid[i][jl][k].T00;
			HydroGrid[i][j][k].T10 = HydroGrid[i][jl][k].T10;
			HydroGrid[i][j][k].T20 = HydroGrid[i][jl][k].T20;
			HydroGrid[i][j][k].T30 = HydroGrid[i][jl][k].T30;
			

			
			HydroGrid[i][j][k].En = HydroGrid[i][jl][k].En;
			HydroGrid[i][j][k].Vx = HydroGrid[i][jl][k].Vx;
			HydroGrid[i][j][k].Vy = HydroGrid[i][jl][k].Vy;
			HydroGrid[i][j][k].Ve = HydroGrid[i][jl][k].Ve;
			HydroGrid[i][j][k].P  = HydroGrid[i][jl][k].P;
			
			HydroGrid[i][j][k].u[0] = HydroGrid[i][jl][k].u[0];
			HydroGrid[i][j][k].u[1] = HydroGrid[i][jl][k].u[1];
			HydroGrid[i][j][k].u[2] = HydroGrid[i][jl][k].u[2];
			HydroGrid[i][j][k].u[3] = HydroGrid[i][jl][k].u[3];	
						
#ifdef GUBSER
			HydroGrid[i][j][k].u[1] = HydroGrid[i][jl+BORDER-j-1][k].u[1];
			HydroGrid[i][j][k].u[2] = HydroGrid[i][jl+BORDER-j-1][k].u[2];
			
#endif	 
		
			HydroGrid[i][j][k].prevu[0] = HydroGrid[i][jl][k].prevu[0];
			HydroGrid[i][j][k].prevu[1] = HydroGrid[i][jl][k].prevu[1];
			HydroGrid[i][j][k].prevu[2] = HydroGrid[i][jl][k].prevu[2];
			HydroGrid[i][j][k].prevu[3] = HydroGrid[i][jl][k].prevu[3];	

			for(l=0;l<Npi;l++)
				HydroGrid[i][j][k].pi[l] = HydroGrid[i][jl][k].pi[l];
			
			HydroGrid[i][j][k].PI = HydroGrid[i][jl][k].PI;
		}
	}
	
	if(MYYRIGHT == BOUNDARY)
	{	
		for(i=il;i<ir;i++)
		for(k=kl;k<kr;k++)
		for(j=0;j<BORDER;j++)
		{
			HydroGrid[i][jr+j][k].T00 = HydroGrid[i][jr-1][k].T00;
			HydroGrid[i][jr+j][k].T10 = HydroGrid[i][jr-1][k].T10;
			HydroGrid[i][jr+j][k].T20 = HydroGrid[i][jr-1][k].T20;
			HydroGrid[i][jr+j][k].T30 = HydroGrid[i][jr-1][k].T30;


			
			HydroGrid[i][jr+j][k].En = HydroGrid[i][jr-1][k].En;
			HydroGrid[i][jr+j][k].Vx = HydroGrid[i][jr-1][k].Vx;
			HydroGrid[i][jr+j][k].Vy = HydroGrid[i][jr-1][k].Vy;
			HydroGrid[i][jr+j][k].Ve = HydroGrid[i][jr-1][k].Ve;
			HydroGrid[i][jr+j][k].P  = HydroGrid[i][jr-1][k].P;
			
			HydroGrid[i][jr+j][k].u[0] = HydroGrid[i][jr-1][k].u[0];
			HydroGrid[i][jr+j][k].u[1] = HydroGrid[i][jr-1][k].u[1];
			HydroGrid[i][jr+j][k].u[2] = HydroGrid[i][jr-1][k].u[2];
			HydroGrid[i][jr+j][k].u[3] = HydroGrid[i][jr-1][k].u[3];
		
#ifdef GUBSER

			HydroGrid[i][jr+j][k].u[1] = HydroGrid[i][jr-j-1][k].u[1];
			HydroGrid[i][jr+j][k].u[2] = HydroGrid[i][jr-j-1][k].u[2];
			
#endif			
						
			HydroGrid[i][jr+j][k].prevu[0] = HydroGrid[i][jr-1][k].prevu[0];
			HydroGrid[i][jr+j][k].prevu[1] = HydroGrid[i][jr-1][k].prevu[1];
			HydroGrid[i][jr+j][k].prevu[2] = HydroGrid[i][jr-1][k].prevu[2];
			HydroGrid[i][jr+j][k].prevu[3] = HydroGrid[i][jr-1][k].prevu[3];
			
			for(l=0;l<Npi;l++)
				HydroGrid[i][jr+j][k].pi[l] = HydroGrid[i][jr-1][k].pi[l];
			
			HydroGrid[i][jr+j][k].PI = HydroGrid[i][jr-1][k].PI;
		}
	}
	
	
}




#if !defined LBI

void ZBoundaryCopy(GRID HydroGrid)
{

	int i,j,k,l;	
	
	if(MYZLEFT == BOUNDARY)
	{	
		for(i=il;i<ir;i++)
		for(j=jl;j<jr;j++)
		for(k=0;k<BORDER;k++)
		{
			HydroGrid[i][j][k].T00 = HydroGrid[i][j][kl].T00;
			HydroGrid[i][j][k].T10 = HydroGrid[i][j][kl].T10;
			HydroGrid[i][j][k].T20 = HydroGrid[i][j][kl].T20;
			HydroGrid[i][j][k].T30 = HydroGrid[i][j][kl].T30;
			
			for(l=0;l<Npi;l++)
				HydroGrid[i][j][k].pi[l] = HydroGrid[i][j][kl].pi[l];
			
			HydroGrid[i][j][k].PI = HydroGrid[i][j][kl].PI;
			
			HydroGrid[i][j][k].En = HydroGrid[i][j][kl].En;
			HydroGrid[i][j][k].Vx = HydroGrid[i][j][kl].Vx;
			HydroGrid[i][j][k].Vy = HydroGrid[i][j][kl].Vy;
			HydroGrid[i][j][k].Ve = HydroGrid[i][j][kl].Ve;
			HydroGrid[i][j][k].P  = HydroGrid[i][j][kl].P;
			
			HydroGrid[i][j][k].u[0] = HydroGrid[i][j][kl].u[0];
			HydroGrid[i][j][k].u[1] = HydroGrid[i][j][kl].u[1];
			HydroGrid[i][j][k].u[2] = HydroGrid[i][j][kl].u[2];
			HydroGrid[i][j][k].u[3] = HydroGrid[i][j][kl].u[3];	
			
			HydroGrid[i][j][k].prevu[0] = HydroGrid[i][j][kl].prevu[0];
			HydroGrid[i][j][k].prevu[1] = HydroGrid[i][j][kl].prevu[1];
			HydroGrid[i][j][k].prevu[2] = HydroGrid[i][j][kl].prevu[2];
			HydroGrid[i][j][k].prevu[3] = HydroGrid[i][j][kl].prevu[3];	
		}
	}
	
	if(MYZRIGHT == BOUNDARY)
	{	
		for(i=il;i<ir;i++)
		for(j=jl;j<jr;j++)
		for(k=0;k<BORDER;k++)
		{
			HydroGrid[i][j][kr+k].T00 = HydroGrid[i][j][kr-1].T00;
			HydroGrid[i][j][kr+k].T10 = HydroGrid[i][j][kr-1].T10;
			HydroGrid[i][j][kr+k].T20 = HydroGrid[i][j][kr-1].T20;
			HydroGrid[i][j][kr+k].T30 = HydroGrid[i][j][kr-1].T30;

			for(l=0;l<Npi;l++)
				HydroGrid[i][j][kr+k].pi[l] = HydroGrid[i][j][kr-1].pi[l];
			
			HydroGrid[i][j][kr+k].PI = HydroGrid[i][j][kr-1].PI;
			
			HydroGrid[i][j][kr+k].En = HydroGrid[i][j][kr-1].En;
			HydroGrid[i][j][kr+k].Vx = HydroGrid[i][j][kr-1].Vx;
			HydroGrid[i][j][kr+k].Vy = HydroGrid[i][j][kr-1].Vy;
			HydroGrid[i][j][kr+k].Ve = HydroGrid[i][j][kr-1].Ve;
			HydroGrid[i][j][kr+k].P  = HydroGrid[i][j][kr-1].P;
			
			HydroGrid[i][j][kr+k].u[0] = HydroGrid[i][j][kr-1].u[0];
			HydroGrid[i][j][kr+k].u[1] = HydroGrid[i][j][kr-1].u[1];
			HydroGrid[i][j][kr+k].u[2] = HydroGrid[i][j][kr-1].u[2];
			HydroGrid[i][j][kr+k].u[3] = HydroGrid[i][j][kr-1].u[3];	
			
			HydroGrid[i][j][kr+k].prevu[0] = HydroGrid[i][j][kr-1].prevu[0];
			HydroGrid[i][j][kr+k].prevu[1] = HydroGrid[i][j][kr-1].prevu[1];
			HydroGrid[i][j][kr+k].prevu[2] = HydroGrid[i][j][kr-1].prevu[2];
			HydroGrid[i][j][kr+k].prevu[3] = HydroGrid[i][j][kr-1].prevu[3];	
			 
		}
	}
}

#endif



void boundary(GRID HydroGrid)
{
	XYBoundaryCopy(HydroGrid);
	
	
#if !defined LBI
	ZBoundaryCopy(HydroGrid);	
#endif

}

