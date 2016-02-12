#define PRINTVAR (6) 

/*********Write out only XY plane*********/
inline int OFFSETXYALLVARS(int gi)  
{
	double sizeline=sizeof(double)*PRINTVAR;
	 
	
	if( (gi%FREQ)==0   )
	{	
		int ypoints = GRIDYPOINTS/FREQ;
		gi /= FREQ;
		return((gi*ypoints + (jSTART/FREQ)) * sizeline);
	}	
	else
	{
		cout<<"problemo2"<<endl;
		return -1;
	}
} 


inline int OFFSETXYALLVARS(int gi , int size)  
{
	double sizeline=sizeof(double)*size;
	
	if(gi<0)
		cout<<"negative index"<<endl;
		
	if( (gi%FREQ)==0   )
	{	
		int ypoints = (NPY*YCMA+2*BORDER)/FREQ;
		gi /= FREQ;
		
		if( MYYLEFT==BOUNDARY)
			return( (gi*ypoints + (jSTART/FREQ)) * sizeline);
		else
			return( (gi*ypoints + BORDER+ (jSTART/FREQ)) * sizeline);
			
	}	
	else
	{
		cout<<"problemo3"<<endl;
		return -1;
	}
}

 


void WriteResultsXY(double tau, GRID HydroGrid)
{

	int len = snprintf(0,0,"res/tau%2.3ffm.bin",tau);
	char *str = new char[len+3];
	sprintf(str,"res/tau%2.3ffm.bin",tau);

	MPI_Offset offset;
	int i,j;
	int gi;
	MPI_Status ierr;
	MPI_File fh;	
	
	MPI_File_open(mpi_grid, str,  MPI_MODE_WRONLY |MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	double *buf;

	int f = FREQ;
	
	int chunk = PRINTVAR*YCMA/f;
	
	buf = new double [chunk];
	int off=0;
	
	for(i=il;i<ir;i=i+f)
	{
		gi = MYGLOBALiWB(i);		
		offset = OFFSETXYALLVARS(gi); //(in bytes)
	
		MPI_File_seek(fh,offset,MPI_SEEK_SET);	//takes offset in bytes

		
		for(j=jl;j<jr;j=j+f)
		{
			int jj=j-jl;
			jj=jj/f;
		
			buf[jj*PRINTVAR+0]=HydroGrid[i][j][k0+off].En;
			buf[jj*PRINTVAR+1]=HydroGrid[i][j][k0+off].Temp;
			buf[jj*PRINTVAR+2]=HydroGrid[i][j][k0+off].P;
			buf[jj*PRINTVAR+3]=HydroGrid[i][j][k0+off].Vx;
			buf[jj*PRINTVAR+4]=HydroGrid[i][j][k0+off].Vy;
			buf[jj*PRINTVAR+5]=HydroGrid[i][j][k0+off].Ve;	
		}


		int error = MPI_File_write(fh, (void*)buf, chunk, MPI_DOUBLE, &ierr);
		
		if(error!=MPI_SUCCESS)
		{
			cout<<"hell0"<<endl;
			MPI_Finalize();
			exit(0);
		}
	}
	
	MPI_File_close(&fh);
	delete buf;

}
 
/*************Write out entire grid***********/

inline int OFFSET(int gi, int gj)  
{
	
	double sizeline=sizeof(double)*PRINTVAR;
	
	if( (gi%FREQ)==0  && (gj%FREQ)==0 )
	{	
		int ypoints = GRIDYPOINTS/FREQ;
		int zpoints = GRIDZPOINTS/FREQZ;
		gi /= FREQ;
		gj /= FREQ;
		return((gi*ypoints*zpoints + gj*zpoints) * sizeline);
	}	
	else
	{
		cout<<"problemo5"<<endl;
		return -1;
	}
}

 





void WriteResults(double tau, GRID HydroGrid)
{


	MPI_Offset offset;

	int i,j,k;
	int gi,gj;
	MPI_Status ierr;
	
	MPI_File fh;	

	int len = snprintf(0,0,"res/tau%2.3ffm.bin",tau);
	char *str = new char[len+3];
	sprintf(str,"res/tau%2.3ffm.bin",tau);
	
	MPI_File_open(mpi_grid, str,  MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	delete str;
	
  

	int  chunk = ZCMA*PRINTVAR/FREQZ; 
	double *buf = new double [chunk];
	
	for(i=il;i<ir;i=i+FREQ)
	for(j=jl;j<jr;j=j+FREQ)
	{	
		gi = MYGLOBALiWB(i);
		gj = MYGLOBALjWB(j);
		
		offset=OFFSET(gi,gj);

		for(k=kl;k<kr;k=k+FREQZ)
		{ 
			int kk = (k-kl)/FREQZ;
			buf[kk*PRINTVAR+0]=HydroGrid[i][j][k].En;	
			buf[kk*PRINTVAR+1]=HydroGrid[i][j][k].Temp;
			buf[kk*PRINTVAR+2]=HydroGrid[i][j][k].P;						
			buf[kk*PRINTVAR+3]=HydroGrid[i][j][k].Vx;
			buf[kk*PRINTVAR+4]=HydroGrid[i][j][k].Vy;
			buf[kk*PRINTVAR+5]=HydroGrid[i][j][k].Ve;
						
		}

		MPI_File_seek(fh,offset,MPI_SEEK_SET);

		int error=MPI_File_write(fh, (void*)buf, chunk, MPI_DOUBLE, &ierr);
		
		if(error!=MPI_SUCCESS)
		{
			cout<<"hell2"<<endl;
			MPI_Finalize();
			exit(0);
		}
	}

	MPI_File_close(&fh);
	delete buf;
}


 
