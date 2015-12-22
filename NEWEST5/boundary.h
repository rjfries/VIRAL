
inline double line(double y2, double y1, double x2, double x1, double x3)
{
	double m,c;

	m=(y2-y1)/(x2-x1);
	c=y1-m*x1;

	return(m*x3+c);
}

inline double gaussian(double y[3], double x[3], double X)
{
	double a,b,c;
	double y1=log(y[0]);
	double y2=log(y[1]);
	double y3=log(y[2]);
	double x1=x[0];
	double x2=x[1];
	double x3=x[2];
	double D= 1.0/((x1-x2)*(x1-x3)*(x2-x3));
	
	a=(x3*(y2-y1) + x2*(y1-y3) + x1*(y3-y2))*D;
	b=(x3*x3*(y1-y2) + x1*x1*(y2-y3) + x2*x2*(y3-y1))*D;
	c=(x1*x3*(x3-x1)*y2 + x2*x2*(x3*y1-x1*y3) + x2*(x1*x1*y3-x3*x3*y1))*D;

	return(exp(a*X*X + b*X+c));

}


void XYBoundaryGaussian(GRID HydroGrid)
{	
	int i,j,k;
	double x[3],y[3];
	
#define var HydroGrid
	if(MYXLEFT == BOUNDARY)
	{
		for(j=jl;j<jr;j++)
		for(k=0;k<ZCMA;k++)
		for(i=0;i<BORDER;i++)
		{
			y[2]=(var[il+2][j][k].En);
			y[1]=(var[il+1][j][k].En);
			y[0]=(var[il][j][k].En);
			x[2]=HydroGrid[il+2][j][k].X;
			x[1]=HydroGrid[il+1][j][k].X;
			x[0]=HydroGrid[il][j][k].X;
				var[i][j][k].En  =  gaussian(y,x,HydroGrid[i][j][k].X);
		}
			
	}
	
	if(MYXRIGHT == BOUNDARY)
	{	
		for(j=jl;j<jr;j++)
		for(k=0;k<ZCMA;k++)
		for(i=0;i<BORDER;i++)	
		{
			y[2]=(var[ir-3][j][k].En);
			y[1]=(var[ir-2][j][k].En);
			y[0]=(var[ir-1][j][k].En);
			x[2]=HydroGrid[ir-3][j][k].X;
			x[1]=HydroGrid[ir-2][j][k].X;
			x[0]=HydroGrid[ir-1][j][k].X;
				var[ir+i][j][k].En  =   gaussian(y,x,HydroGrid[ir+i][j][k].X);
		}
	}

	if(MYYLEFT == BOUNDARY)
	{	
		for(i=il;i<ir;i++)
		for(k=0;k<ZCMA;k++)
		for(j=0;j<BORDER;j++)		
		{
			y[2]=var[i][jl+2][k].En;
			y[1]=var[i][jl+1][k].En;
			y[0]=var[i][jl][k].En;
			x[2]=HydroGrid[i][jl+2][k].Y;
			x[1]=HydroGrid[i][jl+1][k].Y;
			x[0]=HydroGrid[i][jl][k].Y;
				var[i][j][k].En  =  gaussian(y,x,HydroGrid[i][j][k].Y);
		}
	}
	
	if(MYYRIGHT == BOUNDARY)
	{	
		for(i=il;i<ir;i++)
		for(k=0;k<ZCMA;k++)
		for(j=0;j<BORDER;j++)
		{
			y[2]=var[i][jr-3][k].En;
			y[1]=var[i][jr-2][k].En;
			y[0]=var[i][jr-1][k].En;
			x[2]=HydroGrid[i][jr-3][k].Y;
			x[1]=HydroGrid[i][jr-2][k].Y;
			x[0]=HydroGrid[i][jr-1][k].Y;
				var[i][jr+j][k].En  =  gaussian(y,x,HydroGrid[i][jr+j][k].Y);
		}
	}
	
#undef var
}



void XYBoundaryCopy(GRID HydroGrid)
{

	
	int i,j,k,l;
	double x[2],y[2];
	
	if(MYXLEFT == BOUNDARY)
	{
		for(j=jl;j<jr;j++)
		for(k=0;k<ZCMA;k++)
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
		for(k=0;k<ZCMA;k++)
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
		for(k=0;k<ZCMA;k++)
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
		for(k=0;k<ZCMA;k++)
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



 

void boundary(GRID HydroGrid)
{
	XYBoundaryCopy(HydroGrid);

}

