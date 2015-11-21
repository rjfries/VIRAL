

#ifdef GUBSER
	#define PRINTVAR (4+10+1)
#else
	#define PRINTVAR (4+10+1)
#endif

#define SOURCEVAR 4



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



inline int OFFSETXYALLVARSSOURCE(int gi)  
{
	double sizeline=sizeof(double)*4;
	
	if( (gi%FREQ)==0   )
	{	
		int ypoints = GRIDYPOINTS/FREQ;
		gi /= FREQ;
		return((gi*ypoints + (jSTART/FREQ)) * sizeline);
	}	
	else
	{
		cout<<"problemo4"<<endl;
		return -1;
	}
}



void WriteResultsXY(double tau, GRID HydroGrid)
{

	int len = snprintf(0,0,"res/tau%2.2ffm.bin",tau);
	char *str = new char[len+3];
	sprintf(str,"res/tau%2.2ffm.bin",tau);

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
			buf[jj*PRINTVAR+1]=HydroGrid[i][j][k0+off].Vx;
			buf[jj*PRINTVAR+2]=HydroGrid[i][j][k0+off].Vy;
			buf[jj*PRINTVAR+3]=HydroGrid[i][j][k0+off].Ve;
			
			buf[jj*PRINTVAR+4]=HydroGrid[i][j][k0+off].pi[0];
			buf[jj*PRINTVAR+5]=HydroGrid[i][j][k0+off].pi[1];
			buf[jj*PRINTVAR+6]=HydroGrid[i][j][k0+off].pi[2];
			buf[jj*PRINTVAR+7]=HydroGrid[i][j][k0+off].pi[3];
			buf[jj*PRINTVAR+8]=HydroGrid[i][j][k0+off].pi[4];
			buf[jj*PRINTVAR+9]=HydroGrid[i][j][k0+off].pi[5];
			buf[jj*PRINTVAR+10]=HydroGrid[i][j][k0+off].pi[6];
			buf[jj*PRINTVAR+11]=HydroGrid[i][j][k0+off].pi[7];
			buf[jj*PRINTVAR+12]=HydroGrid[i][j][k0+off].pi[8];
			buf[jj*PRINTVAR+13]=HydroGrid[i][j][k0+off].pi[9];
			buf[jj*PRINTVAR+14]=HydroGrid[i][j][k0+off].PI;
			
			
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


void WriteResultsXYCom(double tau, GRID HydroGrid)
{

	int len = snprintf(0,0,"res/tau%2.2ffm.bin",tau);
	char *str = new char[len+3];
	sprintf(str,"res/tau%2.2ffm.bin",tau);

	MPI_Offset offset;
	int i,j;
	int gi;
	MPI_Status ierr;
	MPI_File fh;	
	
	MPI_File_open(mpi_grid, str,  MPI_MODE_WRONLY |MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	double *buf;

	int IL,IR,JL,JR;
	
	IL=il;IR=ir;JL=jl;JR=jr;
	
	if(MYXLEFT == BOUNDARY)
		IL=0;
	if(MYXRIGHT == BOUNDARY)
		IR=XCM;
	if(MYYLEFT == BOUNDARY)
		JL=0;
	if(MYYRIGHT == BOUNDARY)
		JR=YCM;
	
		
	int f = FREQ;
	
	int chunk = PRINTVAR*(JR-JL)/f;
	
	buf = new double [chunk];
	
	int off=0;
	
	for(i=IL;i<IR;i=i+f)
	{
		gi = MYGLOBALi(i);		
		offset = OFFSETXYALLVARS(gi, PRINTVAR); //(in bytes)
	
		MPI_File_seek(fh,offset,MPI_SEEK_SET);	//takes offset in bytes

		
		for(j=JL;j<JR;j=j+f)
		{
			int jj=j-JL;
			jj=jj/f;
		
 
			buf[jj*PRINTVAR+0]=HydroGrid[i][j][k0+off].En;	
			buf[jj*PRINTVAR+1]=HydroGrid[i][j][k0+off].Vx;
			buf[jj*PRINTVAR+2]=HydroGrid[i][j][k0+off].Vy;
			buf[jj*PRINTVAR+3]=HydroGrid[i][j][k0+off].Ve;

					
			buf[jj*PRINTVAR+4]=HydroGrid[i][j][k0+off].pi[0];
			buf[jj*PRINTVAR+5]=HydroGrid[i][j][k0+off].pi[1];
			buf[jj*PRINTVAR+6]=HydroGrid[i][j][k0+off].pi[2];
			buf[jj*PRINTVAR+7]=HydroGrid[i][j][k0+off].pi[3];
			buf[jj*PRINTVAR+8]=HydroGrid[i][j][k0+off].pi[4];
			buf[jj*PRINTVAR+9]=HydroGrid[i][j][k0+off].pi[5];
			buf[jj*PRINTVAR+10]=HydroGrid[i][j][k0+off].pi[6];
			buf[jj*PRINTVAR+11]=HydroGrid[i][j][k0+off].pi[7];
			buf[jj*PRINTVAR+12]=HydroGrid[i][j][k0+off].pi[8];
			buf[jj*PRINTVAR+13]=HydroGrid[i][j][k0+off].pi[9];
			buf[jj*PRINTVAR+14]=HydroGrid[i][j][k0+off].PI; 
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




void WriteTempXYCom(double tau, GRID HydroGrid)
{

	int len = snprintf(0,0,"res/tau%2.2ffm.bin",tau);
	char *str = new char[len+3];
	sprintf(str,"res/tau%2.2ffm.bin",tau);

	MPI_Offset offset;
	int i,j;
	int gi;
	MPI_Status ierr;
	MPI_File fh;	
	
	MPI_File_open(mpi_grid, str,  MPI_MODE_WRONLY |MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	double *buf;

	int IL,IR,JL,JR;
	
	IL=il;IR=ir;JL=jl;JR=jr;
	
	if(MYXLEFT == BOUNDARY)
		IL=0;
	if(MYXRIGHT == BOUNDARY)
		IR=XCM;
	if(MYYLEFT == BOUNDARY)
		JL=0;
	if(MYYRIGHT == BOUNDARY)
		JR=YCM;
	
		
	int f = FREQ;
	int numofvar = TEMPVAR;
	int chunk = numofvar*(JR-JL)/f;
	
	buf = new double [chunk];
	int off=0;
	
	for(i=IL;i<IR;i=i+f)
	{
		gi = MYGLOBALi(i);		
		offset = OFFSETXYALLVARS(gi, numofvar); //(in bytes)
	
		MPI_File_seek(fh,offset,MPI_SEEK_SET);	//takes offset in bytes

		
		for(j=JL;j<JR;j=j+f)
		{
			int jj=j-JL;
			jj=jj/f;
		
			for(int z=0;z<numofvar;z++)
				buf[jj*numofvar+z]=HydroGrid[i][j][k0+off].temp[z];
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





void WriteSourceXY(double tau, GRID HydroGrid)
{

	int len = snprintf(0,0,"res/source%2.2ffm.bin",tau);
	char *str = new char[len+3];
	sprintf(str,"res/source%2.2ffm.bin",tau);

	MPI_Offset offset;
	int i,j;
	int gi;
	MPI_Status ierr;
	MPI_File fh;	
	
	MPI_File_open(mpi_grid, str,  MPI_MODE_WRONLY |MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	double *buf;

	int f = FREQ;
	
	int chunk = SOURCEVAR*YCMA/f;
	
	buf = new double [chunk];
	
	for(i=il;i<ir;i=i+f)
	{
		gi = MYGLOBALiWB(i);		
		offset = OFFSETXYALLVARSSOURCE(gi); //(in bytes)
	
		MPI_File_seek(fh,offset,MPI_SEEK_SET);	//takes offset in bytes

		
		for(j=jl;j<jr;j=j+f)
		{
			int jj=j-jl;
			jj=jj/f;


			for(int l=0;l<SOURCEVAR;l++)
				buf[jj*SOURCEVAR+l] = HydroGrid[i][j][k0].Source[l];
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
		return((gi*ypoints*zpoints + gj*zpoints + (kSTART/FREQZ)) * sizeline);
	}	
	else
	{
		cout<<"problemo5"<<endl;
		return -1;
	}
}




inline int OFFSETSOURCE(int gi, int gj)  
{
	
	double sizeline=sizeof(double)*SOURCEVAR;
	
	if( (gi%FREQ)==0  && (gj%FREQ)==0 )
	{	
		int ypoints = GRIDYPOINTS/FREQ;
		int zpoints = GRIDZPOINTS/FREQZ;
		gi /= FREQ;
		gj /= FREQ;
		return((gi*ypoints*zpoints + gj*zpoints + (kSTART/FREQZ)) * sizeline);
	}	
	else
	{
		cout<<"problemo6"<<endl;
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

	int len = snprintf(0,0,"res/tau%2.2ffm.bin",tau);
	char *str = new char[len+3];
	sprintf(str,"res/tau%2.2ffm.bin",tau);
	
	MPI_File_open(mpi_grid, str,  MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	double *buf;
  

	int  chunk = ZCMA*PRINTVAR/FREQZ; 
	buf = new double [chunk];
	
	for(i=il;i<ir;i=i+FREQ)
	for(j=jl;j<jr;j=j+FREQ)
	{	
		gi = MYGLOBALiWB(i);
		gj = MYGLOBALjWB(j);
		
		offset=OFFSET(gi,gj);

		for(k=kl;k<kr;k=k+FREQZ)
		{
			int kk = k-kl;
			kk=kk/FREQZ;
			
			buf[kk*PRINTVAR+0]=HydroGrid[i][j][k].En;
			buf[kk*PRINTVAR+1]=HydroGrid[i][j][k].Vx;
			buf[kk*PRINTVAR+2]=HydroGrid[i][j][k].Vy;
			buf[kk*PRINTVAR+3]=HydroGrid[i][j][k].Ve;
			
			buf[kk*PRINTVAR+4]=HydroGrid[i][j][k].T00;
			buf[kk*PRINTVAR+5]=HydroGrid[i][j][k].T10;
			buf[kk*PRINTVAR+6]=HydroGrid[i][j][k].T20;
			buf[kk*PRINTVAR+7]=HydroGrid[i][j][k].T30;			
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




void WriteSource(double tau, GRID HydroGrid)
{


	MPI_Offset offset;

	int i,j,k;
	int gi,gj;
	MPI_Status ierr;
	
	MPI_File fh;	

	int len = snprintf(0,0,"res/source%2.2ffm.bin",tau);
	char *str = new char[len+3];
	sprintf(str,"res/source%2.2ffm.bin",tau);


	
	MPI_File_open(mpi_grid, str,  MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	double *buf;
  
  	int chunk = ZCMA*SOURCEVAR/FREQZ; 
	buf = new double [chunk];


	for(i=il;i<ir;i=i+FREQ)
	for(j=jl;j<jr;j=j+FREQ)
	{	
		gi = MYGLOBALiWB(i);
		gj = MYGLOBALjWB(j);
		
	
		offset=OFFSETSOURCE(gi,gj);
				

		for(k=kl;k<kr;k=k+FREQZ)
		{
			int kk = k-kl;
			kk=kk/FREQ;
			
			for(int l=0;l<SOURCEVAR;l++)
				buf[kk*SOURCEVAR+l] = HydroGrid[i][j][k].Source[l];
		}
		

		MPI_File_seek(fh,offset,MPI_SEEK_SET);

		int error=MPI_File_write(fh, (void*)buf, chunk, MPI_DOUBLE, &ierr);
		
		if(error!=MPI_SUCCESS)
		{
			cout<<"hell3"<<endl;
			MPI_Finalize();
			exit(0);
		}
	

	}

	MPI_File_close(&fh);
	delete buf;
}


void WriteResultsXYSteven(double tau, GRID HydroGrid)
{

	int len = snprintf(0,0,"res/tau%2.2ffm.bin",tau);
	char *str = new char[len+3];
	sprintf(str,"res/tau%2.2ffm.bin",tau);

	MPI_Offset offset;
	int i,j;
	int gi;
	MPI_Status ierr;
	MPI_File fh;	
	
	MPI_File_open(mpi_grid, str,  MPI_MODE_WRONLY |MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	double *buf;

	int f = FREQ;
	
	int chunk = (SOURCEVAR)*YCMA/f;
	
	buf = new double [chunk];
	int off=0;
	
	for(i=il;i<ir;i=i+f)
	{
		gi = MYGLOBALiWB(i);		
		offset = OFFSETXYALLVARSSOURCE(gi); //(in bytes)
	
		MPI_File_seek(fh,offset,MPI_SEEK_SET);	//takes offset in bytes

		
		for(j=jl;j<jr;j=j+f)
		{
			int jj=j-jl;
			jj=jj/f;
		
			buf[jj*SOURCEVAR+0]=HydroGrid[i][j][k0+off].En;
			buf[jj*SOURCEVAR+1]=HydroGrid[i][j][k0+off].Vx;
			buf[jj*SOURCEVAR+2]=HydroGrid[i][j][k0+off].Vy;
			buf[jj*SOURCEVAR+3]=HydroGrid[i][j][k0+off].Ve;
			
			
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
