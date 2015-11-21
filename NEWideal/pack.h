void pack(GRID HydroGrid);

#define CHUNKX (PACKVAR*BORDER*YCMA*ZCMA) //in double's
#define CHUNKY (PACKVAR*BORDER*XCMA*ZCMA) //in double's
#define CHUNKZ (PACKVAR*BORDER*XCMA*YCMA) //in double's



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

				rightX[8][m][jj][kk] =	HydroGrid[ir-m-1][j][k].P;
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

				rightY[8][m][ii][kk] =	HydroGrid[i][jr-m-1][k].P;
		
			}
		}

		if(MYZRIGHT != BOUNDARY)
		{
			for(m=0;m<BORDER;m++)
			for(i=il;i<ir;i++)
			for(j=jl;j<jr;j++)
			{
				ii=i-il;
				jj=j-jl;


				rightZ[0][m][ii][jj] =	HydroGrid[i][j][kr-m-1].T00;
				rightZ[1][m][ii][jj] =	HydroGrid[i][j][kr-m-1].T10;
				rightZ[2][m][ii][jj] =	HydroGrid[i][j][kr-m-1].T20;
				rightZ[3][m][ii][jj] =	HydroGrid[i][j][kr-m-1].T30;
				
				rightZ[4][m][ii][jj] =	HydroGrid[i][j][kr-m-1].En;
				rightZ[5][m][ii][jj] =	HydroGrid[i][j][kr-m-1].Vx;
				rightZ[6][m][ii][jj] =	HydroGrid[i][j][kr-m-1].Vy;
				rightZ[7][m][ii][jj] =	HydroGrid[i][j][kr-m-1].Ve;

				rightZ[8][m][ii][jj] =	HydroGrid[i][j][kr-m-1].P;
				
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

				HydroGrid[i][j][k].P = leftX[8][m][jj][kk];


				
		
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


				HydroGrid[i][j][k].P = leftY[8][m][ii][kk];
								
			}
		}

		if(MYZLEFT!= BOUNDARY)
		{

			for(m=0;m<BORDER;m++)
			for(i=il;i<ir;i++)
			for(j=jl;j<jr;j++)
			{
				ii = i-il;
				jj = j-jl;
				k = kl-1-m;
				
				HydroGrid[i][j][k].T00 = leftZ[0][m][ii][jj];
				HydroGrid[i][j][k].T10 = leftZ[1][m][ii][jj];
				HydroGrid[i][j][k].T20 = leftZ[2][m][ii][jj];
				HydroGrid[i][j][k].T30 = leftZ[3][m][ii][jj];
				
				HydroGrid[i][j][k].En = leftZ[4][m][ii][jj];
				HydroGrid[i][j][k].Vx = leftZ[5][m][ii][jj];
				HydroGrid[i][j][k].Vy = leftZ[6][m][ii][jj];
				HydroGrid[i][j][k].Ve = leftZ[7][m][ii][jj];

				HydroGrid[i][j][k].P =   leftZ[8][m][ii][jj];
				
				
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

				leftX[8][m][jj][kk] =	HydroGrid[il+m][j][k].P;
				
			
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

				leftY[8][m][ii][kk] =	HydroGrid[i][jl+m][k].P;
		
			}
		}

		if(MYZLEFT != BOUNDARY)
		{
			for(m=0;m<BORDER;m++)
			for(i=il;i<ir;i++)
			for(j=jl;j<jr;j++)
			{
				ii=i-il;
				jj=j-jl;

				leftZ[0][m][ii][jj] =	HydroGrid[i][j][kr-m-1].T00;
				leftZ[1][m][ii][jj] =	HydroGrid[i][j][kr-m-1].T10;
				leftZ[2][m][ii][jj] =	HydroGrid[i][j][kr-m-1].T20;
				leftZ[3][m][ii][jj] =	HydroGrid[i][j][kr-m-1].T30;
				
				leftZ[4][m][ii][jj] =	HydroGrid[i][j][kr-m-1].En;
				leftZ[5][m][ii][jj] =	HydroGrid[i][j][kr-m-1].Vx;
				leftZ[6][m][ii][jj] =	HydroGrid[i][j][kr-m-1].Vy;
				leftZ[7][m][ii][jj] =	HydroGrid[i][j][kr-m-1].Ve;


				leftZ[8][m][ii][jj] =	HydroGrid[i][j][kr-m-1].P;
				
			
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
				
				HydroGrid[i][j][k].P = rightX[8][m][jj][kk];
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

				HydroGrid[i][j][k].P = rightY[8][m][ii][kk];
				
		
			}
		}

		if(MYZRIGHT!= BOUNDARY)
		{
			for(m=0;m<BORDER;m++)
			for(i=il;i<ir;i++)
			for(j=jl;j<jr;j++)
			{
				ii=i-il;
				jj=j-jl;
				k = kr+m;

				HydroGrid[i][j][k].T00 = rightZ[0][m][ii][jj];
				HydroGrid[i][j][k].T10 = rightZ[1][m][ii][jj];
				HydroGrid[i][j][k].T20 = rightZ[2][m][ii][jj];
				HydroGrid[i][j][k].T30 = rightZ[3][m][ii][jj];
				
				HydroGrid[i][j][k].En = rightZ[4][m][ii][jj];  HydroGrid[i][j][k].P =   EOS(HydroGrid[i][j][k].En , HydroGrid[i][j][k].r);
				HydroGrid[i][j][k].Vx = rightZ[5][m][ii][jj];
				HydroGrid[i][j][k].Vy = rightZ[6][m][ii][jj];
				HydroGrid[i][j][k].Ve = rightZ[7][m][ii][jj];
				HydroGrid[i][j][k].P = rightZ[8][m][ii][jj];
				
			
			}
		}
	}
}


int tag=1;
void pack(GRID HydroGrid) //does not exchange Tnu0 .. last line added for that functionality
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
	
	if(MYZRIGHT != BOUNDARY)
	{
		MPI_Isend(rightZ,CHUNKZ,MPI_DOUBLE,MYZRIGHT,tag,mpi_grid,&req );		
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

	if(MYZLEFT != BOUNDARY)
	{
		MPI_Irecv(leftZ,CHUNKZ,MPI_DOUBLE,MYZLEFT,tag,mpi_grid,&reqLtoR[2]);		
	}


	 


	MPI_Waitall(3,reqLtoR,status);
	UnPackFromLeft( HydroGrid);


	MPI_Barrier(mpi_grid);
/*************************Second Half Of communication ***************************/
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
	if(MYZLEFT != BOUNDARY)
	{
		MPI_Isend(leftZ,CHUNKZ,MPI_DOUBLE,MYZLEFT,tag,mpi_grid,&req );		
	}
//
//If a process to my right, then it must have sent so
// receive from right and unpack
//
	if(MYXRIGHT!= BOUNDARY)
	{
		MPI_Irecv(rightX,CHUNKX,MPI_DOUBLE,MYXRIGHT,tag,mpi_grid,&reqRtoL[0]);
	}
	if(MYYRIGHT!= BOUNDARY)
	{
		MPI_Irecv(rightY,CHUNKY,MPI_DOUBLE,MYYRIGHT,tag,mpi_grid,&reqRtoL[1]);		
	}
	if(MYZRIGHT!= BOUNDARY)
	{
		MPI_Irecv(rightZ,CHUNKZ,MPI_DOUBLE,MYZRIGHT,tag,mpi_grid,&reqRtoL[2]);		
	}
	

	MPI_Waitall(3,reqRtoL,status);
	UnPackFromRight(HydroGrid);	
}

