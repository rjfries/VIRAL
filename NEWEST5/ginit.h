
#define XPOSNUCLEI 3
#define IMPACTPARAMETER (2*XPOSNUCLEI)
 inline double genWENOder(double qm2, double qm1, double qc, double qp1, double qp2, double dx  );

#define NM 4  


typedef double (*ARRXY) [YCM + 2*NOS*BORDER];
ARRXY mu1,mu2,mu1mu2;
ARRXY Dxmu1mu2,Dymu1mu2;
ARRXY Dxmu1,Dxmu2;
ARRXY Dymu1,Dymu2;
ARRXY eps,ax,ay,bx,by; 
ARRXY Dxax,Dyay,Dxbx,Dyby;
ARRXY DivA,DivB;
ARRXY enxy,vxxy,vyxy,vzxy;
ARRXY p1xy,p2xy,p3xy,p4xy,p5xy;
ARRXY TField[4][4];


struct PAR{
double x; 
double y; 
double x_off;
};

int xt = XCM + 2*NOS*BORDER;
int yt = YCM + 2*NOS*BORDER;


//c == complete
int xcs = NOS*BORDER;
int xce = NOS*BORDER+XCM;
int ycs = NOS*BORDER;
int yce = NOS*BORDER+YCM;

//i == inner box
int xis = (NOS+1)*BORDER;
int xie = (NOS+1)*BORDER+XCMA;
int yis = (NOS+1)*BORDER;
int yie = (NOS+1)*BORDER+YCMA;

inline double  WoodSaxon(const double Z, void *param)
{
	PAR* P;

	P =  (PAR *)param;
	double X, Y, Xoff;

	
	X = P->x;
	Y = P->y;
	Xoff = P->x_off;

	double r = sqrt( (X+Xoff)*(X+Xoff) + Y*Y + Z*Z );
	double a = 0.546;

	//	GOLD PARAMETERS
	double A = 197;
	static double R = 1.12*pow(A,1.0/3.0) - 0.86*pow(A,-1.0/3.0);
//	static double R = 1.2*pow(A,1.0/3.0);
	static double rho0 = A/((4.0/3.0)*PIE*R*R*R);	
	double ret =(rho0)/( 1.0 + exp(  (r - R)/a ) );
	
	return( ret ); 
}


void RotateTMuNu(double Tmunu[4][4],double eta, double tau)
{
	int i,j;  

	
	double TRANSFORM[16]={cosh(eta),0,0,-sinh(eta),0,1,0,0,0,0,1,0, -sinh(eta)/tau,0,0,cosh(eta)/tau};	
	double INVERSETRANSFORM[16]={cosh(eta),0,0,-sinh(eta)/tau,0,1,0,0,0,0,1,0,-sinh(eta),0,0,cosh(eta)/tau};	
  	gsl_matrix_view ROT = gsl_matrix_view_array (TRANSFORM, 4, 4); 
  	gsl_matrix_view INVROT = gsl_matrix_view_array (INVERSETRANSFORM, 4, 4); 
	
	double data16[16];
	double temp16[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double result16[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	for(i=0;i<4;i++)
	for(j=0;j<4;j++)
		data16[4*i+j] = Tmunu[i][j] ; 

		

  	gsl_matrix_view TMN = gsl_matrix_view_array (data16, 4, 4); 
  	gsl_matrix_view TMNT = gsl_matrix_view_array (temp16, 4, 4);
  	gsl_matrix_view TMNR = gsl_matrix_view_array (result16, 4, 4);

	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &ROT.matrix, &TMN.matrix,
                  0.0, &TMNT.matrix);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &TMNT.matrix, &INVROT.matrix,
                  0.0, &TMNR.matrix);
	
	for(i=0;i<4;i++)
	for(j=0;j<4;j++)
		Tmunu[i][j] = result16[4*i+j]; 
	
}





void gsl_eigen_complex_print(  gsl_vector_complex*  v ,gsl_matrix_complex*  m)
{
	int i,j;
	double Abs,R,I;
	gsl_complex vec[NM]; 

	cout<<"Eigen Values"<<endl;

	for(i=0;i<4;i++)
	{
		gsl_complex eig_val_i = gsl_vector_complex_get (v, i);
		Abs= gsl_complex_abs(eig_val_i);
		R= GSL_REAL(eig_val_i);
		I= GSL_IMAG(eig_val_i);

		cout<<std::fixed<<std::setprecision(7)<<std::scientific<<R<<" + i "<<I<<" & MAG-> "<<Abs<<endl;
	}

	cout<<endl<<endl;
	cout<<"& Correponding Eigen Vectors"<<endl;

	for(i=0;i<4;i++)
	{
		cout<<"Eigven Vector Number "<<i<<endl<<endl;

		if(GSL_IMAG(gsl_vector_complex_get (v, i)) == 0)
		{
	 		for(j=0;j<4;j++)
			{
				vec[j] = gsl_matrix_complex_get(m, j,i);					
	 			cout<<GSL_REAL(vec[j])/GSL_REAL(vec[0])<<" + i "<< GSL_IMAG(vec[j])<<endl;
	 			if(j==3)
	 				cout<<sqrt(pow(GSL_REAL(vec[1])/GSL_REAL(vec[0]),2) + pow(GSL_REAL(vec[2])/GSL_REAL(vec[0]),2) + pow(GSL_REAL(vec[3])/GSL_REAL(vec[0]),2) );
	 		}			
				
	 	}
		else
		{
	 		for(j=0;j<4;j++)
			{
				vec[j] = gsl_matrix_complex_get(m, j,i);			
	 			cout<<GSL_REAL(vec[j])<<" + i "<< GSL_IMAG(vec[j])<<endl;
	 		}
	 	}
 		cout<<endl<<endl;
 	} 
}



void weno_xy_XDER(ARRXY array,ARRXY result)
{
	double w[3], q[3], d[3], alpha[3];
	double wt[3], qt[3], dt[3], alphat[3];

	double   beta[4];

	double p,eps=WENOEPS,sum;
	int i,j,k,c,l;

	d[0]= dt[2]=0.3;
	d[1]= dt[1]=0.6 ;
	d[2]= dt[0]=0.1 ;
	
	for(j=0;j<yt;j++)
	{
		for(i=2;i<xt-2;i++)
		{
			double qm2,qm1,qc,qp1,qp2;
				
			qm2 = array[i-2][j];
			qm1 = array[i-1][j];
			qc  = array  [i][j];
			qp1 = array[i+1][j];
			qp2 = array[i+2][j];
				
			result[i][j] = genWENOder(qm2,qm1,qc,qp1,qp2,XS);		
		}
	}
}





void weno_xy_YDER(ARRXY array,ARRXY result)
{
	double w[3], q[3], d[3], alpha[3];
	double wt[3], qt[3], dt[3], alphat[3];

	double   beta[4];

	double p,eps=WENOEPS,sum;
	int i,j,k,c,l;

	d[0]= dt[2]=0.3;
	d[1]= dt[1]=0.6 ;
	d[2]= dt[0]=0.1 ;
		
	for(i=0; i<xt; i++)
	for(j=2; j<yt-2; j++)
	{
		
		double qm2,qm1,qc,qp1,qp2;
		
		qm2 = array[i][j-2];
		qm1 = array[i][j-1];
		qc  = array[i][j];
		qp1 = array[i][j+1];
		qp2 = array[i][j+2]; 
		 
		result[i][j] = genWENOder(qm2,qm1,qc,qp1,qp2,YS);
	}
}




void weno_xylog_XDER(ARRXY array,ARRXY result)
{
	double w[3], q[3], d[3], alpha[3];
	double wt[3], qt[3], dt[3], alphat[3];

	double   beta[4];

	double p,eps=WENOEPS,sum;
	int i,j,k,c,l;

	d[0]= dt[2]=0.3;
	d[1]= dt[1]=0.6 ;
	d[2]= dt[0]=0.1 ;
	
	for(j=0;j<yt;j++)
	{
		for(i=2;i<xt-2;i++)
		{
			double qm2,qm1,qc,qp1,qp2; 
			
				
			qm2 = log( array[i-2][j] );			
			qm1 = log( array[i-1][j] );
			qc  = log( array  [i][j] );
			qp1 = log( array[i+1][j] );
			qp2 = log( array[i+2][j] );
			
			double var = array[i][j];
			result[i][j] = var*genWENOder(qm2,qm1,qc,qp1,qp2,XS);				
		}
	}
}


void weno_xylog_YDER(ARRXY array,ARRXY result)
{
	double w[3], q[3], d[3], alpha[3];
	double wt[3], qt[3], dt[3], alphat[3];

	double   beta[4];

	double p,eps=WENOEPS,sum;
	int i,j,k,c,l;

	d[0]= dt[2]=0.3;
	d[1]= dt[1]=0.6 ;
	d[2]= dt[0]=0.1 ;
		
	for(i=0; i<xt; i++)
	for(j=2; j<yt-2; j++)
	{
		
		double qm2,qm1,qc,qp1,qp2; 
		
		qm2 = log( array[i][j-2]);
		qm1 = log( array[i][j-1]);
		qc  = log( array[i][j]  );
		qp1 = log( array[i][j+1]);
		qp2 = log( array[i][j+2]);
 
		double var = array[i][j];
		result[i][j] = var*genWENOder(qm2,qm1,qc,qp1,qp2,YS);
	}
}





inline void FillPiMuNuUP(GRID HydroGrid,double pimn[4][4],int i,int j,int k)
{
	double u0 = HydroGrid[i][j][k].u[0];
	double u1 = HydroGrid[i][j][k].u[1];
	double u2 = HydroGrid[i][j][k].u[2];
	double u3 = HydroGrid[i][j][k].u[3];	 
	
	double p1 = HydroGrid[i][j][k].pi[0];
	double p2 = HydroGrid[i][j][k].pi[1];
	double p3 = HydroGrid[i][j][k].pi[2];
	double p4 = HydroGrid[i][j][k].pi[3];
	double p5 = HydroGrid[i][j][k].pi[4];
	
//----------- \[Pi]^\[Mu]\[Nu]- Matrix ---------------------------------------------
	double pi00  = -((2*p3*u1*u2 + p1*pow(u1,2) + p2*pow(u2,2) + pow(tau,2)*(2*(p4*u1 + p5*u2)*u3 - (p1 + p2)*pow(u3,2)))*pow(-pow(u0,2) + pow(tau,2)*pow(u3,2),-1));
	double pi10  = (p1*u1 + p3*u2 + p4*u3*pow(tau,2))*pow(u0,-1);
	double pi20  = (p3*u1 + p2*u2 + p5*u3*pow(tau,2))*pow(u0,-1);
	double pi30  = pow(u0,-1)*(2*p3*u1*u2*u3 + (p4*u1 + p5*u2 - (p1 + p2)*u3)*pow(u0,2) + p1*u3*pow(u1,2) + p2*u3*pow(u2,2) + p4*u1*pow(tau,2)*pow(u3,2) + p5*u2*pow(tau,2)*pow(u3,2))*pow(pow(u0,2) - pow(tau,2)*pow(u3,2),-1);
	double pi33  = -(pow(tau,-2)*(2*p3*u1*u2 + 2*(p4*u1 + p5*u2)*u3*pow(tau,2) - (p1 + p2)*pow(u0,2) + p1*pow(u1,2) + p2*pow(u2,2))*pow(-pow(u0,2) + pow(tau,2)*pow(u3,2),-1));	



//----------- \[Pi]^\[Mu]\[Nu]- Matrix    
    pimn[0][0]=  pi00;
    pimn[0][1]=  pi10;
    pimn[0][2]=  pi20;
    pimn[0][3]=  pi30;
	pimn[1][0]=  pi10;
    pimn[1][1]=  p1;
    pimn[1][2]=  p3;
    pimn[1][3]=  p4;
    pimn[2][0]=  pi20;
    pimn[2][1]=  p3;
    pimn[2][2]=  p2;
    pimn[2][3]=  p5;
    pimn[3][0]=  pi30;
	pimn[3][1]=  p4;
    pimn[3][2]=  p5;
    pimn[3][3]=  pi33;


}


void AllocateXYArray(ARRXY* buf)
{
	*buf = new double[ XCM + 2*NOS*BORDER][YCM+2*NOS*BORDER];

	for(int i=0;i<xt;i++)
	for(int j=0;j<yt;j++)
			(*buf)[i][j]= 0; 
}


void ReleaseXYArray(ARRXY* buf)
{
	delete [] (*buf);
}


void StoreVars()
{

	for(int i=0;i<4;i++)
	for(int j=0;j<4;j++)
		AllocateXYArray(&TField[i][j]);
		
	
	AllocateXYArray(&mu1);
	AllocateXYArray(&mu2);
	AllocateXYArray(&mu1mu2);
	AllocateXYArray(&Dxmu1mu2);
	AllocateXYArray(&Dymu1mu2);
	AllocateXYArray(&Dxmu1);
	AllocateXYArray(&Dxmu2);
	AllocateXYArray(&Dymu1);
	AllocateXYArray(&Dymu2);
	AllocateXYArray(&eps);
	AllocateXYArray(&ax);
	AllocateXYArray(&ay);
	AllocateXYArray(&bx);
	AllocateXYArray(&by);
	AllocateXYArray(&Dxax);
	AllocateXYArray(&Dxbx);
	AllocateXYArray(&Dyay);
	AllocateXYArray(&Dyby);
	AllocateXYArray(&DivA);
	AllocateXYArray(&DivB);

	AllocateXYArray(&enxy);
	AllocateXYArray(&vxxy);
	AllocateXYArray(&vyxy);
	AllocateXYArray(&vzxy);
	AllocateXYArray(&p1xy);
	AllocateXYArray(&p2xy);
	AllocateXYArray(&p3xy);
	AllocateXYArray(&p4xy);
	AllocateXYArray(&p5xy);
}


void RemoveVars()
{
	
	for(int i=0;i<4;i++)
	for(int j=0;j<4;j++)
		ReleaseXYArray(&TField[i][j]);
		
	ReleaseXYArray(&mu1);
	ReleaseXYArray(&mu2);
	ReleaseXYArray(&mu1mu2);
	ReleaseXYArray(&Dxmu1mu2);
	ReleaseXYArray(&Dymu1mu2);
	ReleaseXYArray(&Dxmu1);
	ReleaseXYArray(&Dxmu2);
	ReleaseXYArray(&Dymu1);
	ReleaseXYArray(&Dymu2);
	ReleaseXYArray(&eps);
	ReleaseXYArray(&ax);
	ReleaseXYArray(&ay);
	ReleaseXYArray(&bx);
	ReleaseXYArray(&by);
	ReleaseXYArray(&Dxax);
	ReleaseXYArray(&Dxbx);
	ReleaseXYArray(&Dyay);
	ReleaseXYArray(&Dyby);
	ReleaseXYArray(&DivA);
	ReleaseXYArray(&DivB);

	
	ReleaseXYArray(&enxy);
	ReleaseXYArray(&vxxy);
	ReleaseXYArray(&vyxy);
	ReleaseXYArray(&vzxy);
	ReleaseXYArray(&p1xy);
	ReleaseXYArray(&p2xy);
	ReleaseXYArray(&p3xy);
	ReleaseXYArray(&p4xy);
	ReleaseXYArray(&p5xy);

}



/* XY PLANE FILE WRITING*/
inline int OFFSETXY(int gi)  
{

	double sizeline=sizeof(double);
	
	if( (gi%FREQ)==0   )
	{	
		int ypoints = GRIDYPOINTS/FREQ;
		gi /= FREQ;
		return((gi*ypoints + (jSTART/FREQ)) * sizeline);
	}	
	else
	{
		cout<<"problemo11"<<endl;
		return -1;
	}
}

void WriteXY(ARRXY VAR, char name[30])
{

	MPI_Offset offset;
	int i,j;
	int gi;
	MPI_Status ierr;
	char str[100];

	int IL = xis;
	int IR = xie;
	int JL = yis;
	int JR = yie;
	
	strcpy(str,"init/");
	strcat(str,name);
		
	MPI_File fh;	
	MPI_File_open(mpi_grid, str,  MPI_MODE_WRONLY |MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	double *buf;

	int f = FREQ;
	//~ int f = 1;
	
	int chunk = YCMA/f;
	
	buf = new double [chunk];
	
	for(i=IL;i<IR;i=i+f)
	{
 		gi = MYGLOBALiWB(i-NOS*BORDER);
		offset = OFFSETXY(gi); //(in bytes)
	
		MPI_File_seek(fh,offset,MPI_SEEK_SET);	//takes offset in bytes

		
		for(j=JL;j<JR;j=j+f)
		{
			int jj=j-(NOS+1)*BORDER;
			jj=jj/f;
			buf[jj]=VAR[i][j];
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




void FindVar(const double tmn[NM][NM], double EV[NM], double tau)
{

	int i,j;
	double VX,VY,VE; 
	double E;
	double X, Y, ETA;

	X=EV[0];
	Y=EV[1];
	ETA=EV[2];

	double inp[NM*NM];
	
	for(i=0;i<NM;i++)
	for(j=0;j<NM;j++)
	{
	
		if( j==1 || j==2 )
			inp[NM*i+j] = -tmn[i][j];
		else if(j==3)
			inp[NM*i+j] = -tau*tau*tmn[i][j];
		else
			inp[NM*i+j] = tmn[i][j];
	}
		
  	gsl_matrix_view TMN = gsl_matrix_view_array (inp, 4, 4);
	gsl_vector_complex *eig_val = gsl_vector_complex_alloc (4);
	gsl_matrix_complex *eig_vecs = gsl_matrix_complex_alloc (4, 4);
	gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc (4);

  	gsl_eigen_nonsymmv (&TMN.matrix, eig_val, eig_vecs, w);
	gsl_eigen_nonsymmv_sort(eig_val,eig_vecs,GSL_EIGEN_SORT_ABS_ASC );
	gsl_complex vec[4]; 

	int ikey,index;
	ikey=0;
	index=-1;

	for(i=0;i<NM;i++)
	{
		gsl_complex eig_val_i = gsl_vector_complex_get (eig_val, i);
		E= gsl_complex_abs(eig_val_i);
		E= fabs(GSL_REAL(eig_val_i));

		if( GSL_IMAG(eig_val_i)==0 )
		{
			for(j=0;j<4;j++)
				vec[j] = gsl_matrix_complex_get(eig_vecs,j,i);
			
			if( GSL_REAL(vec[0]) != 0
				&&  GSL_IMAG (vec[0]) == 0 
				&&  GSL_IMAG (vec[1]) == 0 
				&&  GSL_IMAG (vec[2]) == 0 
				&&  GSL_IMAG (vec[3]) == 0 
			  )
			{
				VX =  GSL_REAL(vec[1])/GSL_REAL(vec[0]);
				VY =  GSL_REAL(vec[2])/GSL_REAL(vec[0]);
				VE =  GSL_REAL(vec[3])/GSL_REAL(vec[0]);

				if( fabs(VX)<1 
					&& fabs(VY)<1 
					&& fabs(VE)<1 
					&& sqrt(VX*VX+VY*VY+VE*VE)<1	
					)
				{
					ikey++;
					index = i;
				}
		
			}
		}
		
	}

 	if(ikey==1)
 	{
		gsl_complex eig_val_i = gsl_vector_complex_get (eig_val, index);
		E= fabs(GSL_REAL(eig_val_i));

		for(j=0;j<4;j++)
				vec[j] = gsl_matrix_complex_get(eig_vecs, j,index);		


		VX =  GSL_REAL(vec[1])/GSL_REAL(vec[0]);
		VY =  GSL_REAL(vec[2])/GSL_REAL(vec[0]);
		VE =  GSL_REAL(vec[3])/GSL_REAL(vec[0]);
 	}
 	else
 	{
		cout<<"EIGEN VALUE SEARCH FAILED"<<endl;
 		cout<<"ikey --> "<< ikey<<endl;
 		cout<<endl<<"@ (X, Y, ETA) --> ( "<<X<<", "<<Y<<","<<ETA<<" )"<<endl;
 		cout<<"EPS0 at this point -->  "<< EV[3]<<endl<<endl; 
 		gsl_eigen_complex_print(eig_val,eig_vecs);
		exit(1);
 	}

	EV[0] = E;
	EV[1] = VX;
	EV[2] = VY;
	EV[3] = VE;
	 
	gsl_eigen_nonsymmv_free (w);		
	gsl_vector_complex_free(eig_val);
	gsl_matrix_complex_free(eig_vecs);
	
} 




























void ginit(GRID HydroGrid, double tau)
{
	int i,j,k;
	
	double xpos = XPOSNUCLEI;
	double X,Y;
	bool AtStart;
	
	if(fabs(tau-TAUSTART)<1e-8)
		AtStart=true;
	else
		AtStart=false;
	
	StoreVars();
		
	
	for(i=0;i<xt;i++)
	for(j=0;j<yt;j++)
	{
		X=XCORDWB(i-NOS*BORDER);
		Y=YCORDWB(j-NOS*BORDER);
		PAR params;

		double fac = 12.2546;
			 
		gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
		double	 abserr;
		gsl_function f;
		
		params.x=X;
		params.y=Y;
		params.x_off= -XPOSNUCLEI;		
		f.params=&params;
		f.function=&WoodSaxon;	
		gsl_integration_qag(&f, 0, 250, 0, 1e-13, 1000, GSL_INTEG_GAUSS61, w, &mu1[i][j], &abserr);
		mu1[i][j] *= (2*fac);
		//~ mu1[i][j] += (PEDESTAL);
		
		params.x=X;
		params.y=Y;
		params.x_off= XPOSNUCLEI; 	
		f.params=&params;
		f.function=&WoodSaxon;	
		gsl_integration_qag(&f, 0, 250, 0, 1e-13, 1000, GSL_INTEG_GAUSS61, w, &mu2[i][j], &abserr);
		mu2[i][j] *= (2*fac);
		//~ mu2[i][j] += (PEDESTAL);


		mu1mu2[i][j] = mu1[i][j]*mu2[i][j];
		gsl_integration_workspace_free (w);
	}
	
	

	
	double mu1mu2max=0;
	for(i=0;i<xt;i++)
	for(j=0;j<yt;j++)
	{
		if(mu1mu2[i][j]>mu1mu2max)
			mu1mu2max=mu1mu2[i][j];
	}
	double globmu1mu2max;
	MPI_Allreduce( &mu1mu2max, &globmu1mu2max, 1,MPI_DOUBLE, MPI_MAX,mpi_grid); 
	
	
	
	
	double Tmax = 0.4/GEVFM;  //0.4GeV @ 0.6 fm/c converted to fm inverse
	double Emaxguess  =   FEnFromTemp(Tmax); 	
	double e0;	
	while( fabs(  s95p_T(Emaxguess) -  Tmax)>1e-10)
	{
		e0 = Emaxguess;
		Emaxguess = e0 - (s95p_T(e0) - Tmax) / (6.608681233 * pow(e0, -0.25));	//6.608681233 
		//~ if(!rank)
			//~ cout<<std::scientific<<"ITERATIVE:: "<< s95p_T(Emaxguess)*GEVFM<<endl;
	}
	double Emax = Emaxguess*pow( tau/0.6,-4./3.);
	
	//~ if(!rank)
	//~ {
		//~ cout<<Emax<<endl;
		//~ cout<<Emax/GEVFM<<endl;
	//~ }
	//~ exit(1);
	
	double Norm = 0.58*Emax/globmu1mu2max;
	
	if(!rank && AtStart)
	{
		cout<<"@0.6:: emax-> "<<Emaxguess<<" & temp in Gev --> "<<s95p_TGev(Emaxguess)<<endl;
		cout<<"@"<<tau<<":: emax-> "<<Emax<<" & temp in Gev --> "<<s95p_TGev(Emax)<<endl;
		cout<<"Normalisation "<<Norm <<endl <<endl <<endl;
	} 
		
	for(i=0;i<xt;i++)
	for(j=0;j<yt;j++)
	{
		eps[i][j] =  Norm*mu1mu2[i][j];//+ (Emax/1e8);
		//~ if(eps[i][j] < 1e-8)
			//~ eps[i][j]=1e-8;
	} 
	
	
	//~ if(AtStart)
	//~ {		
		//~ WriteXY(mu1,"mu1.bin");	
		//~ WriteXY(mu2,"mu2.bin");
		//~ WriteXY(eps,"eps.bin");
	//~ }

	
	
	weno_xylog_XDER(mu1mu2,Dxmu1mu2);
	weno_xylog_YDER(mu1mu2,Dymu1mu2);
	
	weno_xylog_XDER(mu1,Dxmu1);
	weno_xylog_YDER(mu1,Dymu1);
	weno_xylog_XDER(mu2,Dxmu2);
	weno_xylog_YDER(mu2,Dymu2);

 
	



	//~ if(AtStart)
	//~ {
		//~ WriteXY(Dxmu1mu2,"Dxmu1mu2.bin");	
		//~ WriteXY(Dymu1mu2,"Dymu1mu2.bin");
		//~ WriteXY(Dxmu1,"Dxmu1.bin");	
		//~ WriteXY(Dxmu2,"Dxmu2.bin");
		//~ WriteXY(Dymu1,"Dymu1.bin");	
		//~ WriteXY(Dymu2,"Dymu2.bin");
	//~ }




	for(i=0;i<xt;i++)
	for(j=0;j<yt;j++)
	{
		ax[i][j] = (-Norm*Dxmu1mu2[i][j]);
		ay[i][j] = (-Norm*Dymu1mu2[i][j]);
		bx[i][j] = Norm*(-Dxmu1[i][j]*mu2[i][j] + mu1[i][j]*Dxmu2[i][j]);
		by[i][j] = Norm*(-Dymu1[i][j]*mu2[i][j] + mu1[i][j]*Dymu2[i][j]);
	}

	//~ if(AtStart)
	//~ {
		//~ WriteXY(ax,"ax.bin");
		//~ WriteXY(ay,"ay.bin");
		//~ WriteXY(bx,"bx.bin");
		//~ WriteXY(by,"by.bin");
	//~ }
	
	
	weno_xy_XDER(ax,Dxax);
	weno_xy_XDER(bx,Dxbx);
	weno_xy_YDER(ay,Dyay);
	weno_xy_YDER(by,Dyby);

	//~ if(AtStart)
	//~ {	
		//~ WriteXY(Dxax,"Dxax.bin");
		//~ WriteXY(Dxbx,"Dxbx.bin");
		//~ WriteXY(Dyay,"Dyay.bin");
		//~ WriteXY(Dyby,"Dyby.bin");
	//~ }

	for(i=0;i<xt;i++)
	for(j=0;j<yt;j++)
	{
		DivA[i][j] = Dxax[i][j]+Dyay[i][j];
		DivB[i][j] = Dxbx[i][j]+Dyby[i][j];
	}

	//~ if(AtStart)
	//~ {
		//~ WriteXY(DivA,"DivA.bin");
		//~ WriteXY(DivB,"DivB.bin");
	//~ }
	
	
	


	double tau0 = tau;
	double tau02 = tau0*tau0;
	double Tmunu[NM][NM];
	
	
	
	double maxTrField=0;
	double maxTrIFluid=0;
	double maxTrPIFluid=0;
	double maxTrpiFluid=0;

	 
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCMA;k++)
	{
		double EV[NM];
		X =	HydroGrid[i][j][k].X;
		Y = HydroGrid[i][j][k].Y;
		
		double eta = HydroGrid[i][j][k].eta;
		
		int ii = i+NOS*BORDER;
		int jj = j+NOS*BORDER;
			
		EV[0]=X;
		EV[1]=Y;
		EV[2]=eta;
		EV[3]=eps[ii][jj];
		
		double AbsE = fabs(eta);
		double EtaF = 4; //flat part of eta
		double sig = 1;		
		double cutoff =  exp(-  pow(  ( AbsE - (EtaF/2) ) / ( sqrt(2)*sig ),  2)*HeaviSideTheta(  ( AbsE - (EtaF/2) ) )     ) ;

		double ee = eps[ii][jj];
		double dt = 100*ee;
		double AX = ax[ii][jj];
		double AY = ay[ii][jj];
		double BX = bx[ii][jj];
		double BY = by[ii][jj];
		double Le = -DivA[ii][jj];
		double DB = DivB[ii][jj]; 
		
		
		Tmunu[0][0] =  ee - 0.125*tau02*(-2*Le+dt );
		Tmunu[0][1] =  0.5*tau0*(AX );
		Tmunu[0][2] =  0.5*tau0*(AY);
		Tmunu[0][3] =  0.125*tau0*DB;
		Tmunu[1][0] =  Tmunu[0][1];
		Tmunu[1][1] =  ee - 0.25*tau02*(-Le+dt);
		Tmunu[1][2] =  0;
		Tmunu[1][3] =  0.5*(BX);
		Tmunu[2][0] =  Tmunu[0][2];
		Tmunu[2][1] =  Tmunu[1][2];
		Tmunu[2][2] =  Tmunu[1][1];
		Tmunu[2][3] =  0.5*(BY);
		Tmunu[3][0] =  Tmunu[0][3];
		Tmunu[3][1] =  Tmunu[1][3];
		Tmunu[3][2] =  Tmunu[2][3];
		Tmunu[3][3] =  -ee/(tau02)+ 0.125*(-2*Le+3*dt);
		
		for(int  a=0;a<4;a++)
		for(int  b=0;b<4;b++)	
			Tmunu[a][b] *= cutoff;
			
			
		for(int  a=0;a<4;a++)
		for(int  b=0;b<4;b++)
			TField[a][b][ii][jj] = Tmunu[a][b]; 

		double temp1 = Tmunu[0][0] - Tmunu[1][1] - Tmunu[2][2] - tau02*Tmunu[3][3];		
		if(fabs(temp1)>maxTrField)
			maxTrField= fabs(temp1);
			
			
		FindVar(Tmunu,EV,tau0);

		double eps,VX,VY,VE; 
		
		eps = EV[0]+PEDESTAL;
		VX = EV[1];
		VY = EV[2];
		VE = EV[3];
		
		
		HydroGrid[i][j][k].En = eps;
		HydroGrid[i][j][k].Temp = FT(eps );
		HydroGrid[i][j][k].P = EOS(eps );
		HydroGrid[i][j][k].Vx=VX;
		HydroGrid[i][j][k].Vy=VY;
		HydroGrid[i][j][k].Ve=VE;
		HydroGrid[i][j][k].u[0] = 1.0/sqrt(1.0 - VX*VX - VY*VY - tau02*VE*VE);
		HydroGrid[i][j][k].u[1] = HydroGrid[i][j][k].u[0]*VX;
		HydroGrid[i][j][k].u[2] = HydroGrid[i][j][k].u[0]*VY;
		HydroGrid[i][j][k].u[3] = HydroGrid[i][j][k].u[0]*VE;
		
		
		
		DECLePPIa;
		DECLu4;
		if(!AtStart)
		{
			HydroGrid[i][j][k].prevu[0] = u0;
			HydroGrid[i][j][k].prevu[1] = u1;
			HydroGrid[i][j][k].prevu[2] = u2;
			HydroGrid[i][j][k].prevu[3] = u3;
			continue;
		}
		
		
		
		double tid[4][4];
		double PImat[4][4];
		double guu[4][4]={{1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1.0/(tau02)}};
		double uup[4]={u0,u1,u2,u3};
		
		for(int p=0;p<4;p++)
		for(int q=0;q<4;q++)
			tid[p][q]= (e+P)*uup[p]*uup[q] - P*guu[p][q];
	
		double temp2= tid[0][0] - tid[1][1] - tid[2][2] - tau0*tau0*tid[3][3];	
			
		if(fabs(temp2) > maxTrIFluid)
			maxTrIFluid= fabs(temp2);
			
		HydroGrid[i][j][k].PI = -temp2/3.0;
		PI = HydroGrid[i][j][k].PI;
			
		for(int p=0;p<4;p++)
		for(int q=0;q<4;q++)
			PImat[p][q]= PI*( guu[p][q] - uup[p]*uup[q]);	
			
		double temp3= PImat[0][0] - PImat[1][1] - PImat[2][2] - tau0*tau0*PImat[3][3];	
		if(fabs(temp3) > maxTrPIFluid)
			maxTrPIFluid= fabs(temp3);
			
		
		HydroGrid[i][j][k].pi[0] = Tmunu[1][1] - tid[1][1] - PImat[1][1];
		HydroGrid[i][j][k].pi[1] = Tmunu[2][2] - tid[2][2] - PImat[2][2];
		HydroGrid[i][j][k].pi[2] = Tmunu[1][2] - tid[1][2] - PImat[1][2];
		HydroGrid[i][j][k].pi[3] = Tmunu[1][3] - tid[1][3] - PImat[1][3];
		HydroGrid[i][j][k].pi[4] = Tmunu[2][3] - tid[2][3] - PImat[2][3]; 
		
#ifndef BULK		
		HydroGrid[i][j][k].PI = 0;
		PI=0;
#endif
		
		DECLp5;	
		double temp4=  A1 - p1 - p2 - tau0*tau0*A5;	
			
		if(fabs(temp4)>maxTrpiFluid)
			maxTrpiFluid= fabs(temp4);

		HydroGrid[i][j][k].T00 = -P + PI + (e + P - PI)*pow(u0,2) + A1;
		HydroGrid[i][j][k].T10 = (e + P - PI)*u0*u1 + A2;
		HydroGrid[i][j][k].T20 = (e + P - PI)*u0*u2 + A3;
		HydroGrid[i][j][k].T30 = (e + P - PI)*u0*u3 + A4;
	}
	
	
	
	if(AtStart)
	{
			WriteXY(TField[0][0],"t00.bin");
			WriteXY(TField[0][1],"t01.bin");
			WriteXY(TField[0][2],"t02.bin");
			WriteXY(TField[0][3],"t03.bin");
			WriteXY(TField[1][1],"t11.bin");
			WriteXY(TField[1][2],"t12.bin");
			WriteXY(TField[1][3],"t13.bin");
			WriteXY(TField[2][2],"t22.bin");
			WriteXY(TField[2][3],"t23.bin");
			WriteXY(TField[3][3],"t33.bin");
	}
	
	RemoveVars();
	
	
	if(!AtStart)
	{
		if(!rank)
			cout<<"Eigen Value search done at TAUSTART - TS and we filled prevu[4]"<<endl;
	
		return;
	}
	 
	 
	
	double gmaxTrField,gmaxTrIFluid,gmaxTrPIFluid,gmaxTrpiFluid; 
	MPI_Allreduce( &maxTrField, &gmaxTrField, 1,MPI_DOUBLE, MPI_MAX,mpi_grid); 
	MPI_Allreduce( &maxTrIFluid, &gmaxTrIFluid, 1,MPI_DOUBLE, MPI_MAX,mpi_grid); 
	MPI_Allreduce( &maxTrPIFluid, &gmaxTrPIFluid, 1,MPI_DOUBLE, MPI_MAX,mpi_grid); 
	MPI_Allreduce( &maxTrpiFluid, &gmaxTrpiFluid, 1,MPI_DOUBLE, MPI_MAX,mpi_grid); 
	
	if(!rank)
	{
		cout<<"Max TRACE"<<endl;
		cout<<std::scientific<<"Field --> "<<gmaxTrField <<"  IdealFluid --> "<< gmaxTrIFluid <<"  PI Matrix --> "<< gmaxTrPIFluid<<"  pi Matrix --> "<< gmaxTrpiFluid<<endl;
	}

	if(!rank)
		cout<<endl<<"done with eigen value problem at TAUSTART"<<endl;
	
}
