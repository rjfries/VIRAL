
#define XPOSNUCLEI 2.5
#define IMPACTPARAMETER (2*XPOSNUCLEI)


#define NM 4


typedef double (*ARRXY) [YCM + 2*NOS*BORDER];
ARRXY mu1,mu2,mmu1mu2;
ARRXY Dxmu1mu2,Dymu1mu2;
ARRXY Dxmu1,Dxmu2;
ARRXY Dymu1,Dymu2;
ARRXY eps,ax,ay,bx,by; 
ARRXY Dxax,Dyay,Dxbx,Dyby;
ARRXY DivA,DivB;
ARRXY enxy,vxxy,vyxy,vzxy;
ARRXY p1xy,p2xy,p3xy,p4xy,p5xy;



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
			double left,right;
			
				
			qm2 = array[i-2][j];
			qm1 = array[i-1][j];
			qc  = array  [i][j];
			qp1 = array[i+1][j];
			qp2 = array[i+2][j];
			
			
			beta[0]=(13./12.)*pow(qc - 2*qp1 + qp2 , 2) + 0.25*pow(3*qc - 4*qp1 + qp2 , 2);
			beta[1]=(13./12.)*pow(qm1 - 2*qc + qp1 , 2) + 0.25*pow(qm1 - qp1 , 2);
			beta[2]=(13./12.)*pow(qm2 - 2*qm1 + qc , 2) + 0.25*pow(qm2 - 4*qm1 + 3*qc , 2);


			q[0]= (2*qc  + 5*qp1 - 1*qp2)/6.0;
			q[1]= (-1*qm1 + 5*qc + 2*qp1)/6.0;
			q[2]= (2*qm2 - 7*qm1 + 11*qc)/6.0;

			for(c=0;c<3;c++)
				alpha[c]= d[c]/pow(eps + beta[c],WENOP);

			sum=0;
			for(c=0;c<3;c++)
				sum += alpha[c];

			for(c=0;c<3;c++)
				w[c]= alpha[c]/sum;
			
			left=0;
			for(c=0;c<3;c++)
				left += w[c]*q[c];
			
			qt[0]= (11*qc - 7*qp1 + 2*qp2)/6.0;
			qt[1]= (2*qm1 + 5*qc - 1*qp1)/6.0;
			qt[2]= (-1*qm2 + 5*qm1 + 2*qc)/6.0;

			for(c=0;c<3;c++)
				alphat[c]= dt[c]/pow(eps + beta[c],WENOP);
			
			sum=0;
			for(c=0;c<3;c++)
				sum += alphat[c];
			
			for(c=0;c<3;c++)
				wt[c]= alphat[c]/sum;
			
			right=0;
			for(c=0;c<3;c++)
				right += wt[c]*qt[c];
				
			result[i][j] = (right-left)/XS;			
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
		double left,right;
		
		
		qm2 = array[i][j-2];
		qm1 = array[i][j-1];
		qc  = array[i][j];
		qp1 = array[i][j+1];
		qp2 = array[i][j+2];

		
		beta[0]=(13./12.)*pow(qc - 2*qp1 + qp2 , 2) + 0.25*pow(3*qc - 4*qp1 + qp2 , 2);
		beta[1]=(13./12.)*pow(qm1 - 2*qc + qp1 , 2) + 0.25*pow(qm1 - qp1 , 2);
		beta[2]=(13./12.)*pow(qm2 - 2*qm1 + qc , 2) + 0.25*pow(qm2 - 4*qm1 + 3*qc , 2);
		
		q[0]= (2*qc + 5*qp1 - 1*qp2)/6.0;
		q[1]= (-1*qm1 + 5*qc + 2*qp1)/6.0;
		q[2]= (2*qm2 - 7*qm1 + 11*qc)/6.0;

		for(c=0;c<3;c++)
			alpha[c]= d[c]/pow(eps + beta[c],WENOP);

		sum=0;
		for(c=0;c<3;c++)
			sum += alpha[c];

		for(c=0;c<3;c++)
			w[c]= alpha[c]/sum;
		
		left=0;
		for(c=0;c<3;c++)
			left+= w[c]*q[c];


		
		qt[0]= (11*qc   - 7 *qp1 + 2*qp2)/6.0;
		qt[1]= (2*qm1 + 5*qc - 1*qp1)/6.0;
		qt[2]= (-1*qm2 + 5*qm1 + 2*qc)/6.0;

		for(c=0;c<3;c++)
			alphat[c]= dt[c]/pow(eps + beta[c],WENOP);
		
		sum=0;
		for(c=0;c<3;c++)
			sum += alphat[c];
		
		for(c=0;c<3;c++)
			wt[c]= alphat[c]/sum;
		
		right=0;
		for(c=0;c<3;c++)
			right += wt[c]*qt[c];
			
		result[i][j] = (right-left)/YS;
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
			double left,right;
			
				
			qm2 = log( array[i-2][j] );			
			qm1 = log( array[i-1][j] );
			qc  = log( array  [i][j] );
			qp1 = log( array[i+1][j] );
			qp2 = log( array[i+2][j] );
			
			double val = array[i][j] ;
			
			beta[0]=(13./12.)*pow(qc - 2*qp1 + qp2 , 2) + 0.25*pow(3*qc - 4*qp1 + qp2 , 2);
			beta[1]=(13./12.)*pow(qm1 - 2*qc + qp1 , 2) + 0.25*pow(qm1 - qp1 , 2);
			beta[2]=(13./12.)*pow(qm2 - 2*qm1 + qc , 2) + 0.25*pow(qm2 - 4*qm1 + 3*qc , 2);


			q[0]= (2*qc  + 5*qp1 - 1*qp2)/6.0;
			q[1]= (-1*qm1 + 5*qc + 2*qp1)/6.0;
			q[2]= (2*qm2 - 7*qm1 + 11*qc)/6.0;

			for(c=0;c<3;c++)
				alpha[c]= d[c]/pow(eps + beta[c],WENOP);

			sum=0;
			for(c=0;c<3;c++)
				sum += alpha[c];

			for(c=0;c<3;c++)
				w[c]= alpha[c]/sum;
			
			left=0;
			for(c=0;c<3;c++)
				left += w[c]*q[c];
			
			qt[0]= (11*qc - 7*qp1 + 2*qp2)/6.0;
			qt[1]= (2*qm1 + 5*qc - 1*qp1)/6.0;
			qt[2]= (-1*qm2 + 5*qm1 + 2*qc)/6.0;

			for(c=0;c<3;c++)
				alphat[c]= dt[c]/pow(eps + beta[c],WENOP);
			
			sum=0;
			for(c=0;c<3;c++)
				sum += alphat[c];
			
			for(c=0;c<3;c++)
				wt[c]= alphat[c]/sum;
			
			right=0;
			for(c=0;c<3;c++)
				right += wt[c]*qt[c];
				
			result[i][j] = val*(right-left)/XS;			
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
		double left,right;
		
		
		qm2 = log( array[i][j-2]);
		qm1 = log( array[i][j-1]);
		qc  = log( array[i][j]  );
		qp1 = log( array[i][j+1]);
		qp2 = log( array[i][j+2]);

		double val = array[i][j] ;
		beta[0]=(13./12.)*pow(qc - 2*qp1 + qp2 , 2) + 0.25*pow(3*qc - 4*qp1 + qp2 , 2);
		beta[1]=(13./12.)*pow(qm1 - 2*qc + qp1 , 2) + 0.25*pow(qm1 - qp1 , 2);
		beta[2]=(13./12.)*pow(qm2 - 2*qm1 + qc , 2) + 0.25*pow(qm2 - 4*qm1 + 3*qc , 2);
		
		q[0]= (2*qc + 5*qp1 - 1*qp2)/6.0;
		q[1]= (-1*qm1 + 5*qc + 2*qp1)/6.0;
		q[2]= (2*qm2 - 7*qm1 + 11*qc)/6.0;

		for(c=0;c<3;c++)
			alpha[c]= d[c]/pow(eps + beta[c],WENOP);

		sum=0;
		for(c=0;c<3;c++)
			sum += alpha[c];

		for(c=0;c<3;c++)
			w[c]= alpha[c]/sum;
		
		left=0;
		for(c=0;c<3;c++)
			left+= w[c]*q[c];


		
		qt[0]= (11*qc   - 7 *qp1 + 2*qp2)/6.0;
		qt[1]= (2*qm1 + 5*qc - 1*qp1)/6.0;
		qt[2]= (-1*qm2 + 5*qm1 + 2*qc)/6.0;

		for(c=0;c<3;c++)
			alphat[c]= dt[c]/pow(eps + beta[c],WENOP);
		
		sum=0;
		for(c=0;c<3;c++)
			sum += alphat[c];
		
		for(c=0;c<3;c++)
			wt[c]= alphat[c]/sum;
		
		right=0;
		for(c=0;c<3;c++)
			right += wt[c]*qt[c];
			
		result[i][j] = val*(right-left)/YS;
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

	AllocateXYArray(&mu1);
	AllocateXYArray(&mu2);
	AllocateXYArray(&mmu1mu2);
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
	ReleaseXYArray(&mu1);
	ReleaseXYArray(&mu2);
	ReleaseXYArray(&mmu1mu2);
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
		
		params.x=X;
		params.y=Y;
		params.x_off= XPOSNUCLEI; 	
		f.params=&params;
		f.function=&WoodSaxon;	
		gsl_integration_qag(&f, 0, 250, 0, 1e-13, 1000, GSL_INTEG_GAUSS61, w, &mu2[i][j], &abserr);
		mu2[i][j] *= (2*fac);

		eps[i][j] = mu1[i][j]*mu2[i][j];
		gsl_integration_workspace_free (w);
	}
	
	

	
	if(AtStart)
	{		
		WriteXY(mu1,"mu1.bin");	
		WriteXY(mu2,"mu2.bin");
		WriteXY(eps,"eps.bin");
	}

	
	
	weno_xylog_XDER(eps,Dxmu1mu2);
	weno_xylog_YDER(eps,Dymu1mu2);
	
	weno_xylog_XDER(mu1,Dxmu1);
	weno_xylog_YDER(mu1,Dymu1);
	weno_xylog_XDER(mu2,Dxmu2);
	weno_xylog_YDER(mu2,Dymu2);

 
	



	if(AtStart)
	{
		WriteXY(Dxmu1mu2,"Dxmu1mu2.bin");	
		WriteXY(Dymu1mu2,"Dymu1mu2.bin");
		WriteXY(Dxmu1,"Dxmu1.bin");	
		WriteXY(Dxmu2,"Dxmu2.bin");
		WriteXY(Dymu1,"Dymu1.bin");	
		WriteXY(Dymu2,"Dymu2.bin");
	}




	for(i=0;i<xt;i++)
	for(j=0;j<yt;j++)
	{
		ax[i][j] = (-Dxmu1mu2[i][j]);
		ay[i][j] = (-Dymu1mu2[i][j]);
		bx[i][j] = (-Dxmu1[i][j]*mu2[i][j] + mu1[i][j]*Dxmu2[i][j]);
		by[i][j] = (-Dymu1[i][j]*mu2[i][j] + mu1[i][j]*Dymu2[i][j]);
	}

	if(AtStart)
	{
		WriteXY(ax,"ax.bin");
		WriteXY(ay,"ay.bin");
		WriteXY(bx,"bx.bin");
		WriteXY(by,"by.bin");
	}
	
	
	weno_xy_XDER(ax,Dxax);
	weno_xy_XDER(bx,Dxbx);
	weno_xy_YDER(ay,Dyay);
	weno_xy_YDER(by,Dyby);

	if(AtStart)
	{	
		WriteXY(Dxax,"Dxax.bin");
		WriteXY(Dxbx,"Dxbx.bin");
		WriteXY(Dyay,"Dyay.bin");
		WriteXY(Dyby,"Dyby.bin");
	}

	for(i=0;i<xt;i++)
	for(j=0;j<yt;j++)
	{
		DivA[i][j] = Dxax[i][j]+Dyay[i][j];
		DivB[i][j] = Dxbx[i][j]+Dyby[i][j];
	}

	if(AtStart)
	{
		WriteXY(DivA,"DivA.bin");
		WriteXY(DivB,"DivB.bin");
	}
	


	double tau0 = tau;
	double tau02 = tau0*tau0;
	double Tmunu[NM][NM];
	
	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
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

		double ee = eps[ii][jj];
		double dt = 80*ee;
		double AX = ax[ii][jj];
		double AY = ay[ii][jj];
		double BX = bx[ii][jj];
		double BY = by[ii][jj];
		double DA = DivA[ii][jj];
		double DB = DivB[ii][jj];

		double c1= cosh(eta);
		double s1= sinh(eta);
		double c2= cosh(2*eta);
		double s2= sinh(2*eta);
		
		Tmunu[0][0] =  ee + tau02*(-0.25*(dt+DA) + 0.125*dt*(c2) - 0.125*DB*(s2));
		Tmunu[0][1] =  tau0*(AX*(c1) + BX*(s1));
		Tmunu[0][2] =  tau0*(AY*(c1) + BY*(s1));
		Tmunu[0][3] =  tau02*(0.125*dt*(s2) - 0.125*DB*(c2));
		Tmunu[1][0] =  Tmunu[0][1];
		Tmunu[1][1] =  ee - tau02*0.25*(dt+DA);
		Tmunu[1][2] =  0;
		Tmunu[1][3] =  tau0*(AX*(s1) + BX*(c1));
		Tmunu[2][0] =  Tmunu[0][2];
		Tmunu[2][1] =  0;
		Tmunu[2][2] =  Tmunu[1][1];
		Tmunu[2][3] =  tau0*(AY*(s1) + BY*(c1));
		Tmunu[3][0] =  Tmunu[0][3];
		Tmunu[3][1] =  Tmunu[1][3];
		Tmunu[3][2] =  Tmunu[2][3];
		Tmunu[3][3] =  -ee + tau02*( 0.25*(DA+dt) - 0.125*DB*(s2) + 0.125*dt*(c2) );


		RotateTMuNu( Tmunu, eta, tau0);

		FindVar(Tmunu,EV,tau0);

		double eps,VX,VY,VE; 
		
		eps = EV[0];
		VX = EV[1];
		VY = EV[2];
		VE = EV[3];
		

		HydroGrid[i][j][k].En = eps;
		HydroGrid[i][j][k].P = EOS(eps, HydroGrid[i][j][k].r);
		HydroGrid[i][j][k].Vx=VX;
		HydroGrid[i][j][k].Vy=VY;
		HydroGrid[i][j][k].Ve=VE;
		HydroGrid[i][j][k].u[0] = 1.0/sqrt(1.0 - VX*VX - VY*VY - tau0*tau0*VE*VE);
		HydroGrid[i][j][k].u[1] = u0*VX;
		HydroGrid[i][j][k].u[2] = u0*VY;
		HydroGrid[i][j][k].u[3] = u0*VE;
		
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
		
		
		double tid00,tid01,tid02,tid03;
		double tid11,tid12,tid13;
		double tid22,tid23;
		double tid33;
 
		tid00 = -P + ((e + P)*u0*u0);
		tid01 = ((e + P)*u0*u1);
		tid02 = ((e + P)*u0*u2);
		tid03 = ((e + P)*u0*u3);
		
		tid11 = P + ((e + P)*u1*u1);
		tid22 = P + ((e + P)*u2*u2);
		tid12 = ((e + P)*u1*u2);
		tid13 = ((e + P)*u1*u3);
		
		tid23 = ((e + P)*u2*u3); 
		
		tid33 = ((e + P)*u3*u3 + P/(tau0*tau0) ); 

		HydroGrid[i][j][k].pi[0] = Tmunu[1][1] - tid11;
		HydroGrid[i][j][k].pi[1] = Tmunu[2][2] - tid22;
		HydroGrid[i][j][k].pi[2] = Tmunu[1][2] - tid12;
		HydroGrid[i][j][k].pi[3] = Tmunu[1][3] - tid13;
		HydroGrid[i][j][k].pi[4] = Tmunu[2][3] - tid23;
		DECLp5;	

		HydroGrid[i][j][k].T00 = -P + PI + (e + P - PI)*pow(u0,2) + A1;
		HydroGrid[i][j][k].T10 = (e + P - PI)*u0*u1 + A2;
		HydroGrid[i][j][k].T20 = (e + P - PI)*u0*u2 + A3;
		HydroGrid[i][j][k].T30 = (e + P - PI)*u0*u3 + A4;
	}
	
	
	RemoveVars();
	
	
	if(!AtStart)
	{
		if(!rank)
			cout<<"Eigen Value search done at TAUSTART - TS and we filled prevu[4]"<<endl;
	
		return;
	}
	
	double maxT=0;
	
 	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		DECLp5u4;
		double trace =  (A1 - p1 - p2 - A5*(tau0*tau0));
		
		if(maxT>trace)
			maxT=trace;		
	}
	
	if(!rank)
		cout<<"Max TRACE"<<endl;

	MPI_Barrier(mpi_grid);
	
	cout<<maxT<<" --  rank ---"<<rank<<endl;

	MPI_Barrier(mpi_grid);
	
	if(!rank)
		cout<<endl<<"done with eigen value problem at TAUSTART"<<endl;
		
}
