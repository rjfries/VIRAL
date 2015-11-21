

#define PRINTVAR  8
#define DEBUG  9
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
		cout<<"problemo"<<endl;
		return -1;
	}
}


inline int OFFSETXYALLVARS(int gi , int js)  
{
	double sizeline=sizeof(double)*PRINTVAR;
	
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
		cout<<"problemo"<<endl;
		return -1;
	}
}



inline int OFFSETXYALLVARSDebug(int gi , int js)  
{
	double sizeline=sizeof(double)*DEBUG;
	
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
		cout<<"problemo"<<endl;
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
		cout<<"problemo"<<endl;
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
			
			buf[jj*PRINTVAR+4]=HydroGrid[i][j][k0+off].T00;
			buf[jj*PRINTVAR+5]=HydroGrid[i][j][k0+off].T10;
			buf[jj*PRINTVAR+6]=HydroGrid[i][j][k0+off].T20;
			buf[jj*PRINTVAR+7]=HydroGrid[i][j][k0+off].T30;
			
			
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
		offset = OFFSETXYALLVARS(gi, JR-JL); //(in bytes)
	
		MPI_File_seek(fh,offset,MPI_SEEK_SET);	//takes offset in bytes

		
		for(j=JL;j<JR;j=j+f)
		{
			int jj=j-JL;
			jj=jj/f;
		
			buf[jj*PRINTVAR+0]=HydroGrid[i][j][k0+off].En;
			buf[jj*PRINTVAR+1]=HydroGrid[i][j][k0+off].Vx;
			buf[jj*PRINTVAR+2]=HydroGrid[i][j][k0+off].Vy;
			buf[jj*PRINTVAR+3]=HydroGrid[i][j][k0+off].Ve;
			
			buf[jj*PRINTVAR+4]=HydroGrid[i][j][k0+off].T00;
			buf[jj*PRINTVAR+5]=HydroGrid[i][j][k0+off].T10;
			buf[jj*PRINTVAR+6]=HydroGrid[i][j][k0+off].T20;
			buf[jj*PRINTVAR+7]=HydroGrid[i][j][k0+off].T30;
			
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




void WriteResultsXYComDebug(double tau, GRID HydroGrid)
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
	
	int chunk = DEBUG*(JR-JL)/f;
	
	buf = new double [chunk];
	int off=0;
	
	for(i=IL;i<IR;i=i+f)
	{
		gi = MYGLOBALi(i);		
		offset = OFFSETXYALLVARSDebug(gi, JR-JL); //(in bytes)
	
		MPI_File_seek(fh,offset,MPI_SEEK_SET);	//takes offset in bytes

		
		for(j=JL;j<JR;j=j+f)
		{
			int jj=j-JL;
			jj=jj/f;
		
			buf[jj*DEBUG+0]=HydroGrid[i][j][k0+off].fluxT;
			buf[jj*DEBUG+1]=HydroGrid[i][j][k0+off].PartialResult[0];
			if(i!=0)
				buf[jj*DEBUG+2]=(HydroGrid[i][j][k0+off].T00 - HydroGrid[i-1][j][k0+off].T00);
			else
				buf[jj*DEBUG+2]=0;
			
			buf[jj*DEBUG+3]=HydroGrid[i][j][k0+off].Vx;			
			buf[jj*DEBUG+4]=HydroGrid[i][j][k0+off].velL;
			buf[jj*DEBUG+5]=HydroGrid[i][j][k0+off].velR;
			
			buf[jj*DEBUG+6]=HydroGrid[i][j][k0+off].BufA[1];
			buf[jj*DEBUG+7]=HydroGrid[i][j][k0+off].varL;
			buf[jj*DEBUG+8]=HydroGrid[i][j][k0+off].varR;			
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
		cout<<"problemo"<<endl;
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
		cout<<"problemo"<<endl;
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


void WriteResults1(double tau, GRID HydroGrid)
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
			
			buf[kk*PRINTVAR+0]=HydroGrid[i][j][k].fluxL;
			buf[kk*PRINTVAR+1]=HydroGrid[i][j][k].fluxR;
			buf[kk*PRINTVAR+2]=HydroGrid[i][j][k].velL;
			buf[kk*PRINTVAR+3]=HydroGrid[i][j][k].velR;
			
			buf[kk*PRINTVAR+4]=HydroGrid[i][j][k].varL;
			buf[kk*PRINTVAR+5]=HydroGrid[i][j][k].varR;
			buf[kk*PRINTVAR+6]=HydroGrid[i][j][k].u0L;
			buf[kk*PRINTVAR+7]=HydroGrid[i][j][k].u0R;	
									
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


