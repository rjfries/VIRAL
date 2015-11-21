
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
		for(k=kl;k<kr;k++)
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
		for(k=kl;k<kr;k++)
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
		for(k=kl;k<kr;k++)
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
		for(k=kl;k<kr;k++)
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
			HydroGrid[i][j][k].P = HydroGrid[i][jl][k].P;
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
			HydroGrid[i][jr+j][k].P = HydroGrid[i][jr-1][k].P;
		}
	}
}

void ZBoundaryCopy(GRID HydroGrid)
{

	int i,j,k;	
	
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
			
			HydroGrid[i][j][k].En = HydroGrid[i][j][kl].En;
			HydroGrid[i][j][k].Vx = HydroGrid[i][j][kl].Vx;
			HydroGrid[i][j][k].Vy = HydroGrid[i][j][kl].Vy;
			HydroGrid[i][j][k].Ve = HydroGrid[i][j][kl].Ve;
			HydroGrid[i][j][k].P = HydroGrid[i][j][kl].P;
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

			HydroGrid[i][j][kr+k].En = HydroGrid[i][j][kr-1].En;
			HydroGrid[i][j][kr+k].Vx = HydroGrid[i][j][kr-1].Vx;
			HydroGrid[i][j][kr+k].Vy = HydroGrid[i][j][kr-1].Vy;
			HydroGrid[i][j][kr+k].Ve = HydroGrid[i][j][kr-1].Ve;
			HydroGrid[i][j][kr+k].P = HydroGrid[i][j][kr-1].P;
		}
	}
}



void boundary(GRID HydroGrid)
{

	XYBoundaryCopy(HydroGrid);
	ZBoundaryCopy(HydroGrid);	
}





/*
 * 
 * 
 */
 
 
void XYBoundaryCopyF(GRID HydroGrid)
{

	
	int i,j,k;
	
	if(MYXLEFT == BOUNDARY)
	{
		for(j=jl;j<jr;j++)
		for(k=kl;k<kr;k++)
		for(i=0;i<BORDER;i++)
		{

			HydroGrid[i][j][k].fluxL = HydroGrid[il][j][k].fluxL;
			HydroGrid[i][j][k].fluxR = HydroGrid[il][j][k].fluxR;
			HydroGrid[i][j][k].velL= HydroGrid[il][j][k].velL;
			HydroGrid[i][j][k].velR = HydroGrid[il][j][k].velR;


			HydroGrid[i][j][k].varL = HydroGrid[il][j][k].varL;
			HydroGrid[i][j][k].varR = HydroGrid[il][j][k].varR;
			HydroGrid[i][j][k].u0L = HydroGrid[il][j][k].u0L;
			HydroGrid[i][j][k].u0R = HydroGrid[il][j][k].u0R;


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
			HydroGrid[i][j][k].P = HydroGrid[i][jl][k].P;
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
			HydroGrid[i][jr+j][k].P = HydroGrid[i][jr-1][k].P;
		}
	}
}

void ZBoundaryCopyF(GRID HydroGrid)
{

	int i,j,k;	
	
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
			
			HydroGrid[i][j][k].En = HydroGrid[i][j][kl].En;
			HydroGrid[i][j][k].Vx = HydroGrid[i][j][kl].Vx;
			HydroGrid[i][j][k].Vy = HydroGrid[i][j][kl].Vy;
			HydroGrid[i][j][k].Ve = HydroGrid[i][j][kl].Ve;
			HydroGrid[i][j][k].P = HydroGrid[i][j][kl].P;
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

			HydroGrid[i][j][kr+k].En = HydroGrid[i][j][kr-1].En;
			HydroGrid[i][j][kr+k].Vx = HydroGrid[i][j][kr-1].Vx;
			HydroGrid[i][j][kr+k].Vy = HydroGrid[i][j][kr-1].Vy;
			HydroGrid[i][j][kr+k].Ve = HydroGrid[i][j][kr-1].Ve;
			HydroGrid[i][j][kr+k].P = HydroGrid[i][j][kr-1].P;
		}
	}
}



void boundaryF(GRID HydroGrid)
{

	XYBoundaryCopyF(HydroGrid);
	ZBoundaryCopyF(HydroGrid);	
}
