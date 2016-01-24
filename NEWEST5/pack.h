void pack(GRID HydroGrid);

#define CHUNKX (PACKVAR*BORDER*YCMA*ZCMA) //in double's
#define CHUNKY (PACKVAR*BORDER*XCMA*ZCMA) //in double's
#define CHUNKZ (PACKVAR*BORDER*XCMA*YCMA) //in double's

#define CHUNKXY (PACKVAR*BORDER*BORDER*ZCMA) //in double's



void PackToRight(GRID HydroGrid);
void PackToLeft(GRID HydroGrid);
void UnPackFromRight(GRID HydroGrid);
void UnPackFromLeft(GRID HydroGrid);



void PackToRight(GRID HydroGrid)
{
	int i,j,k,l,m;
	int ii,jj,kk;

	
	if(MYXRIGHT != BOUNDARY || MYYRIGHT != BOUNDARY || MYZRIGHT != BOUNDARY)
	{
		if(MYXRIGHT != BOUNDARY)
		{
		
			for(m=0;m<BORDER;m++)
			for(j=jl;j<jr;j++)
			for(k=kl;k<kr;k++)
			{
				jj=j-jl;
				kk=k-kl;

				rightX[0][m][jj][kk] =	HydroGrid[ir-m-1][j][k].T00;
				rightX[1][m][jj][kk] =	HydroGrid[ir-m-1][j][k].T10;
				rightX[2][m][jj][kk] =	HydroGrid[ir-m-1][j][k].T20;
				rightX[3][m][jj][kk] =	HydroGrid[ir-m-1][j][k].T30;
				
				rightX[4][m][jj][kk] =	HydroGrid[ir-m-1][j][k].En;
				rightX[5][m][jj][kk] =	HydroGrid[ir-m-1][j][k].Vx;
				rightX[6][m][jj][kk] =	HydroGrid[ir-m-1][j][k].Vy;
				rightX[7][m][jj][kk] =	HydroGrid[ir-m-1][j][k].Ve;
				
				for(l=0;l<VARN;l++)
				{
					rightX[2*VARN+l][m][jj][kk] =	HydroGrid[ir-m-1][j][k].u[l];
					rightX[3*VARN+l][m][jj][kk] =	HydroGrid[ir-m-1][j][k].prevu[l];
				}
				
				for(l=0;l<Npi;l++)
					rightX[4*VARN+l][m][jj][kk] =	HydroGrid[ir-m-1][j][k].pi[l];
					
				rightX[4*VARN+Npi][m][jj][kk] =	HydroGrid[ir-m-1][j][k].PI;
				rightX[4*VARN+Npi+NPI][m][jj][kk] =	HydroGrid[ir-m-1][j][k].P;
			}
				
		}
	

		if(MYYRIGHT != BOUNDARY)
		{
		
			for(m=0;m<BORDER;m++)
			for(i=il;i<ir;i++)
			for(k=kl;k<kr;k++)
			{
				ii=i-il;
				kk=k-kl;

				rightY[0][m][ii][kk] =	HydroGrid[i][jr-m-1][k].T00;
				rightY[1][m][ii][kk] =	HydroGrid[i][jr-m-1][k].T10;
				rightY[2][m][ii][kk] =	HydroGrid[i][jr-m-1][k].T20;
				rightY[3][m][ii][kk] =	HydroGrid[i][jr-m-1][k].T30;
				
				rightY[4][m][ii][kk] =	HydroGrid[i][jr-m-1][k].En;
				rightY[5][m][ii][kk] =	HydroGrid[i][jr-m-1][k].Vx;
				rightY[6][m][ii][kk] =	HydroGrid[i][jr-m-1][k].Vy;
				rightY[7][m][ii][kk] =	HydroGrid[i][jr-m-1][k].Ve;
		
				for(l=0;l<VARN;l++)
				{
					rightY[2*VARN+l][m][ii][kk] =	HydroGrid[i][jr-m-1][k].u[l];
					rightY[3*VARN+l][m][ii][kk] =	HydroGrid[i][jr-m-1][k].prevu[l];
				}
				
				for(l=0;l<Npi;l++)
					rightY[4*VARN+l][m][ii][kk] =  HydroGrid[i][jr-m-1][k].pi[l];
					
				rightY[4*VARN+Npi][m][ii][kk] =	HydroGrid[i][jr-m-1][k].PI;
				rightY[4*VARN+Npi+NPI][m][ii][kk] =	HydroGrid[i][jr-m-1][k].P;
			}
		}
 
	}
}



void UnPackFromLeft(GRID HydroGrid)
{


	int i,j,k,l,m;
	int ii,jj,kk;
	if(MYXLEFT!= BOUNDARY || MYYLEFT!= BOUNDARY || MYZLEFT!= BOUNDARY)
	{
		if(MYXLEFT!= BOUNDARY)
		{
			for(m=0;m<BORDER;m++)
			for(j=jl;j<jr;j++)
			for(k=kl;k<kr;k++)
			{
				jj=j-jl;
				kk=k-kl;
				i=il-1-m;
				
				HydroGrid[i][j][k].T00 = leftX[0][m][jj][kk];
				HydroGrid[i][j][k].T10 = leftX[1][m][jj][kk];
				HydroGrid[i][j][k].T20 = leftX[2][m][jj][kk];
				HydroGrid[i][j][k].T30 = leftX[3][m][jj][kk];
				
				
				HydroGrid[i][j][k].En = leftX[4][m][jj][kk];
				HydroGrid[i][j][k].Vx = leftX[5][m][jj][kk];
				HydroGrid[i][j][k].Vy = leftX[6][m][jj][kk];
				HydroGrid[i][j][k].Ve = leftX[7][m][jj][kk];
				
				for(l=0;l<VARN;l++)
				{
					HydroGrid[i][j][k].u[l] = leftX[2*VARN+l][m][jj][kk];
					HydroGrid[i][j][k].prevu[l] = leftX[3*VARN+l][m][jj][kk];
				}
				
				for(l=0;l<Npi;l++)
					HydroGrid[i][j][k].pi[l]= leftX[4*VARN+l][m][jj][kk];
				
				HydroGrid[i][j][k].PI= leftX[4*VARN+Npi][m][jj][kk];	
				HydroGrid[i][j][k].P= leftX[4*VARN+Npi+NPI][m][jj][kk];	
			}
		}
		 
		if(MYYLEFT!= BOUNDARY)
		{
			for(m=0;m<BORDER;m++)
			for(i=il;i<ir;i++)
			for(k=kl;k<kr;k++)
			{
				
				ii=i-il;
				kk=k-kl;
				j=jl-1-m;
				
				HydroGrid[i][j][k].T00 = leftY[0][m][ii][kk];  
				HydroGrid[i][j][k].T10 = leftY[1][m][ii][kk];
				HydroGrid[i][j][k].T20 = leftY[2][m][ii][kk];
				HydroGrid[i][j][k].T30 = leftY[3][m][ii][kk];
			
				HydroGrid[i][j][k].En = leftY[4][m][ii][kk]; 
				HydroGrid[i][j][k].Vx = leftY[5][m][ii][kk];
				HydroGrid[i][j][k].Vy = leftY[6][m][ii][kk];
				HydroGrid[i][j][k].Ve = leftY[7][m][ii][kk];
				
				for(l=0;l<VARN;l++)
				{
					HydroGrid[i][j][k].u[l] = leftY[2*VARN+l][m][ii][kk];
					HydroGrid[i][j][k].prevu[l] = leftY[3*VARN+l][m][ii][kk];
				}
				
				for(l=0;l<Npi;l++)
					HydroGrid[i][j][k].pi[l] = leftY[4*VARN+l][m][ii][kk];
				
				HydroGrid[i][j][k].PI = leftY[4*VARN+Npi][m][ii][kk];	
				HydroGrid[i][j][k].P = leftY[4*VARN+Npi+NPI][m][ii][kk];	
				
				
								
			}
		}
 
	}
}

void PackToLeft(GRID HydroGrid)
{
	int i,j,k,l,m;
	int ii,jj,kk;

	if(MYXLEFT!= BOUNDARY || MYYLEFT!= BOUNDARY || MYZLEFT!= BOUNDARY)
	{

		if(MYXLEFT != BOUNDARY)
		{
		
			for(m=0;m<BORDER;m++)
			for(j=jl;j<jr;j++)
			for(k=kl;k<kr;k++)
			{
				jj=j-jl;
				kk=k-kl;
				
				
				leftX[0][m][jj][kk] =	HydroGrid[il+m][j][k].T00;
				leftX[1][m][jj][kk] =	HydroGrid[il+m][j][k].T10;
				leftX[2][m][jj][kk] =	HydroGrid[il+m][j][k].T20;
				leftX[3][m][jj][kk] =	HydroGrid[il+m][j][k].T30;
				
				leftX[4][m][jj][kk] =	HydroGrid[il+m][j][k].En;
				leftX[5][m][jj][kk] =	HydroGrid[il+m][j][k].Vx;
				leftX[6][m][jj][kk] =	HydroGrid[il+m][j][k].Vy;
				leftX[7][m][jj][kk] =	HydroGrid[il+m][j][k].Ve;
				
				for(l=0;l<VARN;l++)
				{
					leftX[2*VARN+l][m][jj][kk] =	HydroGrid[il+m][j][k].u[l];
					leftX[3*VARN+l][m][jj][kk] =	HydroGrid[il+m][j][k].prevu[l];
				}
				
				for(l=0;l<Npi;l++)
					leftX[4*VARN+l][m][jj][kk] =	HydroGrid[il+m][j][k].pi[l];
					
				leftX[4*VARN+Npi][m][jj][kk] =	HydroGrid[il+m][j][k].PI;
				leftX[4*VARN+Npi+NPI][m][jj][kk] =	HydroGrid[il+m][j][k].P;
			}
		}
	

		if(MYYLEFT != BOUNDARY)
		{
			for(m=0;m<BORDER;m++)
			for(i=il;i<ir;i++)
			for(k=kl;k<kr;k++)
			{
				ii=i-il;
				kk=k-kl;
				
				leftY[0][m][ii][kk] =	HydroGrid[i][jl+m][k].T00;
				leftY[1][m][ii][kk] =	HydroGrid[i][jl+m][k].T10;
				leftY[2][m][ii][kk] =	HydroGrid[i][jl+m][k].T20;
				leftY[3][m][ii][kk] =	HydroGrid[i][jl+m][k].T30;
				
				leftY[4][m][ii][kk] =	HydroGrid[i][jl+m][k].En;
				leftY[5][m][ii][kk] =	HydroGrid[i][jl+m][k].Vx;
				leftY[6][m][ii][kk] =	HydroGrid[i][jl+m][k].Vy;
				leftY[7][m][ii][kk] =	HydroGrid[i][jl+m][k].Ve;
		
				for(l=0;l<VARN;l++)
				{
					leftY[2*VARN+l][m][ii][kk] =	HydroGrid[i][jl+m][k].u[l];
					leftY[3*VARN+l][m][ii][kk] =	HydroGrid[i][jl+m][k].prevu[l];
				}
				
				for(l=0;l<Npi;l++)
					leftY[4*VARN+l][m][ii][kk] =  HydroGrid[i][jl+m][k].pi[l];
					
				leftY[4*VARN+Npi][m][ii][kk] =	HydroGrid[i][jl+m][k].PI;
				leftY[4*VARN+Npi+NPI][m][ii][kk] =	HydroGrid[i][jl+m][k].P;
			}
		}
 
	}
}

void UnPackFromRight(GRID HydroGrid)
{


	int i,j,k,l,m;
	int ii,jj,kk;

	
	if(MYXRIGHT != BOUNDARY || MYYRIGHT != BOUNDARY || MYZRIGHT != BOUNDARY)
	{
		if(MYXRIGHT != BOUNDARY)
		{
			for(m=0;m<BORDER;m++)
			for(j=jl;j<jr;j++)
			for(k=kl;k<kr;k++)
			{
				jj=j-jl;
				kk=k-kl;
				i = ir+m;
					
				HydroGrid[i][j][k].T00 = rightX[0][m][jj][kk];
				HydroGrid[i][j][k].T10 = rightX[1][m][jj][kk];
				HydroGrid[i][j][k].T20 = rightX[2][m][jj][kk];
				HydroGrid[i][j][k].T30 = rightX[3][m][jj][kk];
				
				
				HydroGrid[i][j][k].En = rightX[4][m][jj][kk]; 
				HydroGrid[i][j][k].Vx = rightX[5][m][jj][kk];
				HydroGrid[i][j][k].Vy = rightX[6][m][jj][kk];
				HydroGrid[i][j][k].Ve = rightX[7][m][jj][kk];
				
				for(l=0;l<VARN;l++)
				{
					HydroGrid[i][j][k].u[l] = rightX[2*VARN+l][m][jj][kk];
					HydroGrid[i][j][k].prevu[l] = rightX[3*VARN+l][m][jj][kk];
				}
				
				for(l=0;l<Npi;l++)
					HydroGrid[i][j][k].pi[l]= rightX[4*VARN+l][m][jj][kk];
				
				HydroGrid[i][j][k].PI= rightX[4*VARN+Npi][m][jj][kk];	
				HydroGrid[i][j][k].P = rightX[4*VARN+Npi+NPI][m][jj][kk];	
				
						
			}
		}
		
		if(MYYRIGHT!= BOUNDARY)
		{
			for(m=0;m<BORDER;m++)
			for(i=il;i<ir;i++)
			for(k=kl;k<kr;k++)
			{
				ii=i-il;
				kk=k-kl;
				j=jr+m;
				
				HydroGrid[i][j][k].T00 = rightY[0][m][ii][kk];  
				HydroGrid[i][j][k].T10 = rightY[1][m][ii][kk];
				HydroGrid[i][j][k].T20 = rightY[2][m][ii][kk];
				HydroGrid[i][j][k].T30 = rightY[3][m][ii][kk];
			
				HydroGrid[i][j][k].En = rightY[4][m][ii][kk];  
				HydroGrid[i][j][k].Vx = rightY[5][m][ii][kk];
				HydroGrid[i][j][k].Vy = rightY[6][m][ii][kk];
				HydroGrid[i][j][k].Ve = rightY[7][m][ii][kk];
				
				for(l=0;l<VARN;l++)
				{
					HydroGrid[i][j][k].u[l] = rightY[2*VARN+l][m][ii][kk];
					HydroGrid[i][j][k].prevu[l] = rightY[3*VARN+l][m][ii][kk];
				}
				
				for(l=0;l<Npi;l++)
					HydroGrid[i][j][k].pi[l] = rightY[4*VARN+l][m][ii][kk];
				
				HydroGrid[i][j][k].PI = rightY[4*VARN+Npi][m][ii][kk];				
				HydroGrid[i][j][k].P  = rightY[4*VARN+Npi+NPI][m][ii][kk];				
			}
		}
 
	}
}



/*
 * 
 * 
 * 
 * Diagonal Communication
 * 
 *
*/

void PackXL( GRID HydroGrid )
{
	int k,l,m,n;

	
	if(MYXLYL != BOUNDARY)
	{
		for(m=0;m<BORDER;m++)
		for(n=0;n<BORDER;n++)
		for(k=kl;k<kr;k++)
		{  
			int ioff = il;
			int joff = jl;
			
			xlyl[0][m][n][k] = HydroGrid[ioff + m][joff + n][k].T00;
			xlyl[1][m][n][k] = HydroGrid[ioff + m][joff + n][k].T10;
			xlyl[2][m][n][k] = HydroGrid[ioff + m][joff + n][k].T20;
			xlyl[3][m][n][k] = HydroGrid[ioff + m][joff + n][k].T30;
			
			xlyl[4][m][n][k] = HydroGrid[ioff + m][joff + n][k].En;
			xlyl[5][m][n][k] = HydroGrid[ioff + m][joff + n][k].Vx;
			xlyl[6][m][n][k] = HydroGrid[ioff + m][joff + n][k].Vy;
			xlyl[7][m][n][k] = HydroGrid[ioff + m][joff + n][k].Ve;
			
			for(l=0;l<VARN;l++)
			{
				xlyl[2*VARN+l][m][n][k] = HydroGrid[ioff + m][joff + n][k].u[l];
				xlyl[3*VARN+l][m][n][k] = HydroGrid[ioff + m][joff + n][k].prevu[l];
			}
			
			for(l=0;l<Npi;l++)
				xlyl[4*VARN+l][m][n][k] = HydroGrid[ioff + m][joff + n][k].pi[l];
			
			xlyl[4*VARN+Npi][m][n][k] = HydroGrid[ioff + m][joff + n][k].PI;
			xlyl[4*VARN+Npi+NPI][m][n][k] = HydroGrid[ioff + m][joff + n][k].P;
		}
	}		
	 
	
	
	
	if(  MYXLYR != BOUNDARY)
	{
		if(MYXLYR != BOUNDARY)
		{
			for(m=0;m<BORDER;m++)
			for(n=0;n<BORDER;n++)
			for(k=kl;k<kr;k++)
			{
			
				int ioff = il;
				int joff = jr-BORDER;
					
				xlyr[0][m][n][k] = HydroGrid[ioff + m][joff + n][k].T00;
				xlyr[1][m][n][k] = HydroGrid[ioff + m][joff + n][k].T10;
				xlyr[2][m][n][k] = HydroGrid[ioff + m][joff + n][k].T20;
				xlyr[3][m][n][k] = HydroGrid[ioff + m][joff + n][k].T30;
				
				
				xlyr[4][m][n][k] = HydroGrid[ioff + m][joff + n][k].En;
				xlyr[5][m][n][k] = HydroGrid[ioff + m][joff + n][k].Vx;
				xlyr[6][m][n][k] = HydroGrid[ioff + m][joff + n][k].Vy;
				xlyr[7][m][n][k] = HydroGrid[ioff + m][joff + n][k].Ve;
				
				for(l=0;l<VARN;l++)
				{
					xlyr[2*VARN+l][m][n][k] = HydroGrid[ioff + m][joff + n][k].u[l];
					xlyr[3*VARN+l][m][n][k] = HydroGrid[ioff + m][joff + n][k].prevu[l];
				}
				
				for(l=0;l<Npi;l++)
					xlyr[4*VARN+l][m][n][k] = HydroGrid[ioff + m][joff + n][k].pi[l];
				
				xlyr[4*VARN+Npi][m][n][k] = HydroGrid[ioff + m][joff + n][k].PI;
				xlyr[4*VARN+Npi+NPI][m][n][k] = HydroGrid[ioff + m][joff + n][k].P;
			}
		}		
	}
}






void PackXR( GRID HydroGrid )
{
	int  k,l,m,n; 

	if(MYXRYL != BOUNDARY || MYXRYR != BOUNDARY)
	{
		if(MYXRYL != BOUNDARY)
		{
			for(m=0;m<BORDER;m++)
			for(n=0;n<BORDER;n++)
			for(k=kl;k<kr;k++)
			{ 
				
				int ioff = ir - BORDER;
				int joff = jl;
					
				xryl[0][m][n][k] = HydroGrid[ioff + m][joff + n][k].T00;
				xryl[1][m][n][k] = HydroGrid[ioff + m][joff + n][k].T10;
				xryl[2][m][n][k] = HydroGrid[ioff + m][joff + n][k].T20;
				xryl[3][m][n][k] = HydroGrid[ioff + m][joff + n][k].T30;
				
				
				xryl[4][m][n][k] = HydroGrid[ioff + m][joff + n][k].En;
				xryl[5][m][n][k] = HydroGrid[ioff + m][joff + n][k].Vx;
				xryl[6][m][n][k] = HydroGrid[ioff + m][joff + n][k].Vy;
				xryl[7][m][n][k] = HydroGrid[ioff + m][joff + n][k].Ve;
				
				for(l=0;l<VARN;l++)
				{
					xryl[2*VARN+l][m][n][k] = HydroGrid[ioff + m][joff + n][k].u[l];
					xryl[3*VARN+l][m][n][k] = HydroGrid[ioff + m][joff + n][k].prevu[l];
				}
				
				for(l=0;l<Npi;l++)
					xryl[4*VARN+l][m][n][k] = HydroGrid[ioff + m][joff + n][k].pi[l];
				
				xryl[4*VARN+Npi][m][n][k] = HydroGrid[ioff + m][joff + n][k].PI;
				xryl[4*VARN+Npi+NPI][m][n][k] = HydroGrid[ioff + m][joff + n][k].P;
			}
		}	
	
	
		if(MYXRYR != BOUNDARY)
		{
			for(m=0;m<BORDER;m++)
			for(n=0;n<BORDER;n++)
			for(k=kl;k<kr;k++)
			{
				
				int ioff = ir-BORDER;					
				int joff = jr-BORDER;
					
				xryr[0][m][n][k] = HydroGrid[ioff + m][joff + n][k].T00;
				xryr[1][m][n][k] = HydroGrid[ioff + m][joff + n][k].T10;
				xryr[2][m][n][k] = HydroGrid[ioff + m][joff + n][k].T20;
				xryr[3][m][n][k] = HydroGrid[ioff + m][joff + n][k].T30;
				
				
				xryr[4][m][n][k] = HydroGrid[ioff + m][joff + n][k].En;
				xryr[5][m][n][k] = HydroGrid[ioff + m][joff + n][k].Vx;
				xryr[6][m][n][k] = HydroGrid[ioff + m][joff + n][k].Vy;
				xryr[7][m][n][k] = HydroGrid[ioff + m][joff + n][k].Ve;
				
				for(l=0;l<VARN;l++)
				{
					xryr[2*VARN+l][m][n][k] = HydroGrid[ioff + m][joff + n][k].u[l];
					xryr[3*VARN+l][m][n][k] = HydroGrid[ioff + m][joff + n][k].prevu[l];
				}
				
				for(l=0;l<Npi;l++)
					xryr[4*VARN+l][m][n][k] = HydroGrid[ioff + m][joff + n][k].pi[l];
				
				xryr[4*VARN+Npi][m][n][k] = HydroGrid[ioff + m][joff + n][k].PI;
				xryr[4*VARN+Npi+NPI][m][n][k] = HydroGrid[ioff + m][joff + n][k].P;
			}
		}	
	}
}





void UnpackXR(GRID  HydroGrid )
{
	int  k,l,m,n; 

	if(MYXRYL != BOUNDARY || MYXRYR != BOUNDARY)
	{
		if(MYXRYL != BOUNDARY)
		{
			for(m=0;m<BORDER;m++)
			for(n=0;n<BORDER;n++)
			for(k=kl;k<kr;k++)
			{ 
				
				int ioff = ir;
				int joff = 0;
					
				HydroGrid[ioff + m][joff + n][k].T00 = xryl[0][m][n][k] ;
				HydroGrid[ioff + m][joff + n][k].T10 = xryl[1][m][n][k] ;
				HydroGrid[ioff + m][joff + n][k].T20 = xryl[2][m][n][k] ;
				HydroGrid[ioff + m][joff + n][k].T30 = xryl[3][m][n][k] ;
				
				
				HydroGrid[ioff + m][joff + n][k].En = xryl[4][m][n][k] ;
				HydroGrid[ioff + m][joff + n][k].Vx = xryl[5][m][n][k] ;
				HydroGrid[ioff + m][joff + n][k].Vy = xryl[6][m][n][k] ;
				HydroGrid[ioff + m][joff + n][k].Ve = xryl[7][m][n][k] ;
				
				for(l=0;l<VARN;l++)
				{
					HydroGrid[ioff + m][joff + n][k].u[l] = xryl[2*VARN+l][m][n][k];
					HydroGrid[ioff + m][joff + n][k].prevu[l] = xryl[3*VARN+l][m][n][k];
				}
				
				for(l=0;l<Npi;l++)
					HydroGrid[ioff + m][joff + n][k].pi[l] = xryl[4*VARN+l][m][n][k];
				
				HydroGrid[ioff + m][joff + n][k].PI = xryl[4*VARN+Npi][m][n][k];
				HydroGrid[ioff + m][joff + n][k].P = xryl[4*VARN+Npi+NPI][m][n][k];
			}
		}	
	
	
		if(MYXRYR != BOUNDARY)
		{
			for(m=0;m<BORDER;m++)
			for(n=0;n<BORDER;n++)
			for(k=kl;k<kr;k++)
			{
				
				int ioff = ir;					
				int joff = jr;
					
				HydroGrid[ioff + m][joff + n][k].T00 = xryr[0][m][n][k];
				HydroGrid[ioff + m][joff + n][k].T10 = xryr[1][m][n][k];
				HydroGrid[ioff + m][joff + n][k].T20 = xryr[2][m][n][k];
				HydroGrid[ioff + m][joff + n][k].T30 = xryr[3][m][n][k];
				
				
				HydroGrid[ioff + m][joff + n][k].En = xryr[4][m][n][k];
				HydroGrid[ioff + m][joff + n][k].Vx = xryr[5][m][n][k];
				HydroGrid[ioff + m][joff + n][k].Vy = xryr[6][m][n][k];
				HydroGrid[ioff + m][joff + n][k].Ve = xryr[7][m][n][k];
				
				for(l=0;l<VARN;l++)
				{
					HydroGrid[ioff + m][joff + n][k].u[l] = xryr[2*VARN+l][m][n][k];
					HydroGrid[ioff + m][joff + n][k].prevu[l] = xryr[3*VARN+l][m][n][k];
				}
				
				for(l=0;l<Npi;l++)
					HydroGrid[ioff + m][joff + n][k].pi[l] = xryr[4*VARN+l][m][n][k];
				
				HydroGrid[ioff + m][joff + n][k].PI = xryr[4*VARN+Npi][m][n][k];
				HydroGrid[ioff + m][joff + n][k].P = xryr[4*VARN+Npi+NPI][m][n][k];
			}
		}	
	}	
}


void UnpackXL( GRID HydroGrid )
{
	int  k,l,m,n; 

	if(MYXLYL != BOUNDARY || MYXLYR != BOUNDARY)
	{
		if(MYXLYL != BOUNDARY)
		{
			for(m=0;m<BORDER;m++)
			for(n=0;n<BORDER;n++)
			for(k=kl;k<kr;k++)
			{  
					
				int ioff = 0;					
				int joff = 0;
					
				HydroGrid[ioff + m][joff + n][k].T00 = xlyl[0][m][n][k] ;
				HydroGrid[ioff + m][joff + n][k].T10 = xlyl[1][m][n][k] ;
				HydroGrid[ioff + m][joff + n][k].T20 = xlyl[2][m][n][k] ;
				HydroGrid[ioff + m][joff + n][k].T30 = xlyl[3][m][n][k] ;
				
				
				HydroGrid[ioff + m][joff + n][k].En = xlyl[4][m][n][k] ;
				HydroGrid[ioff + m][joff + n][k].Vx = xlyl[5][m][n][k] ;
				HydroGrid[ioff + m][joff + n][k].Vy = xlyl[6][m][n][k] ;
				HydroGrid[ioff + m][joff + n][k].Ve = xlyl[7][m][n][k] ;
				
				for(l=0;l<VARN;l++)
				{
					HydroGrid[ioff + m][joff + n][k].u[l] = xlyl[2*VARN+l][m][n][k];
					HydroGrid[ioff + m][joff + n][k].prevu[l] = xlyl[3*VARN+l][m][n][k];
				}
				
				for(l=0;l<Npi;l++)
					HydroGrid[ioff + m][joff + n][k].pi[l] = xlyl[4*VARN+l][m][n][k];
				
				HydroGrid[ioff + m][joff + n][k].PI = xlyl[4*VARN+Npi][m][n][k];
				HydroGrid[ioff + m][joff + n][k].P = xlyl[4*VARN+Npi+NPI][m][n][k];
			}
		}	
	
	
		if(MYXLYR != BOUNDARY)
		{
			for(m=0;m<BORDER;m++)
			for(n=0;n<BORDER;n++)
			for(k=kl;k<kr;k++)
			{
				
				int ioff = 0;
				int joff = jr;
					
				HydroGrid[ioff + m][joff + n][k].T00 = xlyr[0][m][n][k];
				HydroGrid[ioff + m][joff + n][k].T10 = xlyr[1][m][n][k];
				HydroGrid[ioff + m][joff + n][k].T20 = xlyr[2][m][n][k];
				HydroGrid[ioff + m][joff + n][k].T30 = xlyr[3][m][n][k];
				
				
				HydroGrid[ioff + m][joff + n][k].En = xlyr[4][m][n][k];
				HydroGrid[ioff + m][joff + n][k].Vx = xlyr[5][m][n][k];
				HydroGrid[ioff + m][joff + n][k].Vy = xlyr[6][m][n][k];
				HydroGrid[ioff + m][joff + n][k].Ve = xlyr[7][m][n][k];
				
				for(l=0;l<VARN;l++)
				{
					HydroGrid[ioff + m][joff + n][k].u[l] = xlyr[2*VARN+l][m][n][k];
					HydroGrid[ioff + m][joff + n][k].prevu[l] = xlyr[3*VARN+l][m][n][k];
				}
				
				for(l=0;l<Npi;l++)
					HydroGrid[ioff + m][joff + n][k].pi[l] = xlyr[4*VARN+l][m][n][k];
				
				HydroGrid[ioff + m][joff + n][k].PI = xlyr[4*VARN+Npi][m][n][k];
				HydroGrid[ioff + m][joff + n][k].P = xlyr[4*VARN+Npi+NPI][m][n][k];
			}
		}	
	}	
}


int tag=1;
void pack(GRID HydroGrid )  
{

	MPI_Request req;

    MPI_Request reqLtoR[3],reqRtoL[3];
    MPI_Status status[3];

	for(int i=0;i<3;i++)
	{
		reqLtoR[i]=MPI_REQUEST_NULL ;
		reqRtoL[i]=MPI_REQUEST_NULL ;
	}

	
//
//If a process to my right , then pack and send data
//
	PackToRight(HydroGrid);
	
	if(MYXRIGHT != BOUNDARY)
	{
		MPI_Isend(rightX,CHUNKX,MPI_DOUBLE,MYXRIGHT,tag,mpi_grid,&req );
	}

	if(MYYRIGHT != BOUNDARY)
	{
		MPI_Isend(rightY,CHUNKY,MPI_DOUBLE,MYYRIGHT,tag,mpi_grid,&req );		
	}
	
	 

//If a process to my left, then it must have sent so
// receive from left and unpack

	
	if(MYXLEFT!= BOUNDARY)
	{
		MPI_Irecv(leftX,CHUNKX,MPI_DOUBLE,MYXLEFT,tag,mpi_grid,&reqLtoR[0]);
	}
	if(MYYLEFT != BOUNDARY)
	{
		MPI_Irecv(leftY,CHUNKY,MPI_DOUBLE,MYYLEFT,tag,mpi_grid,&reqLtoR[1]);		
	}
 


	 


	MPI_Waitall(3,reqLtoR,status);
	UnPackFromLeft( HydroGrid);


	//~ MPI_Barrier(mpi_grid);
	
	
//
//If a process to my left i pack and send data
//
	PackToLeft( HydroGrid);
	if(MYXLEFT != BOUNDARY)
	{
		MPI_Isend(leftX,CHUNKX,MPI_DOUBLE,MYXLEFT,tag,mpi_grid,&req );
	}
	if(MYYLEFT != BOUNDARY)
	{
		MPI_Isend(leftY,CHUNKY,MPI_DOUBLE,MYYLEFT,tag,mpi_grid,&req );		
	}
 
 
// If a process to my right, then it must have sent so
// receive from right and unpack
 
	if(MYXRIGHT!= BOUNDARY)
	{
		MPI_Irecv(rightX,CHUNKX,MPI_DOUBLE,MYXRIGHT,tag,mpi_grid,&reqRtoL[0]);
	}
	if(MYYRIGHT!= BOUNDARY)
	{
		MPI_Irecv(rightY,CHUNKY,MPI_DOUBLE,MYYRIGHT,tag,mpi_grid,&reqRtoL[1]);		
	}
 
	

	MPI_Waitall(3,reqRtoL,status);
	UnPackFromRight(HydroGrid);	
	


	//~ MPI_Barrier(mpi_grid);




	
	
	
	/*
	 * Diagonal Communication
	 */
	
	//~ MPI_Request reqXL2XR[2],reqXR2XL[2];
  
	//~ for(int i=0;i<2;i++)
	//~ {
		//~ reqXL2XR[i]=MPI_REQUEST_NULL ;
		//~ reqXR2XL[i]=MPI_REQUEST_NULL ;
	//~ }



	
//~ // Send from XL towards XR

	//~ PackXL( HydroGrid );

	//~ if(MYXLYL != BOUNDARY)
	//~ {
		//~ MPI_Isend(xlyl,CHUNKXY,MPI_DOUBLE,MYXLYL,tag,mpi_grid,&req );
	//~ }

	//~ if(MYXLYR != BOUNDARY)
	//~ {
		//~ MPI_Isend(xlyr,CHUNKXY,MPI_DOUBLE,MYXLYR,tag,mpi_grid,&req );		
	//~ }

//~ //Receive data into XR and unpack
	//~ if(MYXRYL!= BOUNDARY)
	//~ {
		//~ MPI_Irecv(xryl,CHUNKXY,MPI_DOUBLE,MYXRYL,tag,mpi_grid,&reqXL2XR[0]);
	//~ }
	//~ if(MYXRYR != BOUNDARY)
	//~ {
		
		//~ MPI_Irecv(xryr,CHUNKXY,MPI_DOUBLE,MYXRYR,tag,mpi_grid,&reqXL2XR[1]);		
	//~ }
 
	//~ MPI_Waitall(2,reqXL2XR,status);
	//~ UnpackXR(HydroGrid);
	//~ MPI_Barrier(mpi_grid);
//~ //Second half of diagonal communication
	
	
//~ // Send from XR towards XL

	//~ PackXR( HydroGrid );

	//~ if(MYXRYL != BOUNDARY)
	//~ {
		//~ MPI_Isend(xryl,CHUNKXY,MPI_DOUBLE,MYXRYL,tag,mpi_grid,&req );
	//~ }

	//~ if(MYXRYR != BOUNDARY)
	//~ {
		
		//~ MPI_Isend(xryr,CHUNKXY,MPI_DOUBLE,MYXRYR,tag,mpi_grid,&req );		
	//~ }

//~ //Receive data into XL and unpack
	//~ MPI_Barrier(mpi_grid);
	//~ if(MYXLYL!= BOUNDARY)
	//~ {
		//~ MPI_Irecv(xlyl,CHUNKXY,MPI_DOUBLE,MYXLYL,tag,mpi_grid,&reqXR2XL[0]);
	//~ }
	
	//~ if(MYXLYR != BOUNDARY)
	//~ {
		//~ MPI_Irecv(xlyr,CHUNKXY,MPI_DOUBLE,MYXLYR,tag,mpi_grid,&reqXR2XL[1]);		
	//~ }
 
	//~ MPI_Waitall(2,reqXR2XL,status);
	//~ UnpackXL(HydroGrid); 
	
	
}

