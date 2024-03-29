inline double  MIN(double a,double b){return(fmin(a,b));}
inline double  MIN(double a,double b,double c){return(fmin(fmin(a,b),c) );}
inline double  MAX(double a,double b) {return(fmax(a,b));}
inline double  MAX(double a,double b,double c){return(fmax(fmax(a,b),c) );}
inline double  SGN(double x)  {return((x<0)?-1:1);}


inline double minmod(double r)
{
	double phi=0;

	if(r>0)
		phi = (MAX(0 ,MIN(1,r)));
	
	return (phi);
}


inline double minmod(double a , double b)
{
	double ret;

	if(a*b>0)
		ret = SGN(a)*MIN(fabs(a),fabs(b));
	else 
		ret =0;

	return ret;

}


inline double minmod(double a , double b, double c)
{
	double ret;

	if(a*b>0 && b*c>0)
		ret = SGN(a) * MIN(fabs(a),fabs(b),fabs(c));
	else 
		ret=0;

	return ret;
}

inline double genminmod(double um, double u, double up, double dx, double gminv )
{

	if(gminv<0)
		return ((up-um)/(2*dx));
		
	double diff = gminv;	
	double ret = minmod(diff*((u-um)/dx) , ((up-um)/(2*dx)) ,diff*((up-u)/dx));

	return ret;
}

 inline double genWENOder(double qm2, double qm1, double qc, double qp1, double qp2, double dx  )
{ 	
	double w[3], q[3], d[3], alpha[3];
	double wt[3], qt[3], dt[3], alphat[3];
	double beta[4];
	double p,eps=WENOEPS,sum;
	double left=0,right=0;
	
	d[0]= dt[2]=0.3;
	d[1]= dt[1]=0.6 ;
	d[2]= dt[0]=0.1 ;
	beta[0]=(13./12.)*pow(qc - 2*qp1 + qp2 , 2) + 0.25*pow(3*qc - 4*qp1 + qp2 , 2);
	beta[1]=(13./12.)*pow(qm1 - 2*qc + qp1 , 2) + 0.25*pow(qm1 - qp1 , 2);
	beta[2]=(13./12.)*pow(qm2 - 2*qm1 + qc , 2) + 0.25*pow(qm2 - 4*qm1 + 3*qc , 2);
	q[0]= (2*qc + 5*qp1 - 1*qp2)/6.0;
	q[1]= (-1*qm1 + 5*qc + 2*qp1)/6.0;
	q[2]= (2*qm2 - 7*qm1 + 11*qc)/6.0;
	for(int c=0;c<3;c++)
		alpha[c]= d[c]/pow(eps + beta[c],WENOP);
	sum=0;
	for(int c=0;c<3;c++)
		sum += alpha[c];
	for(int c=0;c<3;c++)
		w[c]= alpha[c]/sum;
	for(int c=0;c<3;c++)
		left+= w[c]*q[c];
		
	qt[0]= (11*qc   - 7 *qp1 + 2*qp2)/6.0;
	qt[1]= (2*qm1 + 5*qc - 1*qp1)/6.0;
	qt[2]= (-1*qm2 + 5*qm1 + 2*qc)/6.0;
	for(int c=0;c<3;c++)
		alphat[c]= dt[c]/pow(eps + beta[c],WENOP);
	sum=0;
	for(int c=0;c<3;c++)
		sum += alphat[c];
	for(int c=0;c<3;c++)
		wt[c]= alphat[c]/sum;
	for(int c=0;c<3;c++)
		right += wt[c]*qt[c];
		
	return((left-right)/dx);
}

  

inline double genWENOR(double qm2, double qm1, double qc, double qp1, double qp2  )
{ 	 
	double wt[3], qt[3], dt[3], alphat[3];
	double beta[4];
	double eps=WENOEPS,sum;
	double right=0;
	
	dt[2]=0.3;
	dt[1]=0.6 ;
	dt[0]=0.1 ;
	
	beta[0]=(13./12.)*pow(qc - 2*qp1 + qp2 , 2) + 0.25*pow(3*qc - 4*qp1 + qp2 , 2);
	beta[1]=(13./12.)*pow(qm1 - 2*qc + qp1 , 2) + 0.25*pow(qm1 - qp1 , 2);
	beta[2]=(13./12.)*pow(qm2 - 2*qm1 + qc , 2) + 0.25*pow(qm2 - 4*qm1 + 3*qc , 2);
 		
	qt[0]= (11*qc   - 7 *qp1 + 2*qp2)/6.0;
	qt[1]= (2*qm1 + 5*qc - 1*qp1)/6.0;
	qt[2]= (-1*qm2 + 5*qm1 + 2*qc)/6.0;
	for(int c=0;c<3;c++)
		alphat[c]= dt[c]/pow(eps + beta[c],WENOP);
	sum=0;
	for(int c=0;c<3;c++)
		sum += alphat[c];
	for(int c=0;c<3;c++)
		wt[c]= alphat[c]/sum;
	for(int c=0;c<3;c++)
		right += wt[c]*qt[c];
		
	return(right);
}
inline double genWENOL(double qm2, double qm1, double qc, double qp1, double qp2 )
{ 	
	double w[3], q[3], d[3], alpha[3]; 
	double beta[4];
	double p,eps=WENOEPS,sum;
	double left=0;
	
	d[0]= 0.3;
	d[1]= 0.6 ;
	d[2]= 0.1 ;
	beta[0]=(13./12.)*pow(qc - 2*qp1 + qp2 , 2) + 0.25*pow(3*qc - 4*qp1 + qp2 , 2);
	beta[1]=(13./12.)*pow(qm1 - 2*qc + qp1 , 2) + 0.25*pow(qm1 - qp1 , 2);
	beta[2]=(13./12.)*pow(qm2 - 2*qm1 + qc , 2) + 0.25*pow(qm2 - 4*qm1 + 3*qc , 2);
	q[0]= (2*qc + 5*qp1 - 1*qp2)/6.0;
	q[1]= (-1*qm1 + 5*qc + 2*qp1)/6.0;
	q[2]= (2*qm2 - 7*qm1 + 11*qc)/6.0;
	for(int c=0;c<3;c++)
		alpha[c]= d[c]/pow(eps + beta[c],WENOP);
	sum=0;
	for(int c=0;c<3;c++)
		sum += alpha[c];
	for(int c=0;c<3;c++)
		w[c]= alpha[c]/sum;
	for(int c=0;c<3;c++)
		left+= w[c]*q[c];
	
		
	return(left);
}

 


inline double SBee(double r)
{
	double phi=0;
	if(r>0)
		phi =( MAX(0 ,MIN(2*r,1), MIN(r,2)) );
	
	return (phi);
}

inline double VLeer(double r)
{
	double phi=0;
	if(r>0)
		phi = ( ( r + fabs(r))/(1+fabs(r)));

	return (phi);
}
//ospre flux limiter
inline double OFLimiter(double r)
{
	double phi=0;
	if(r>0)
		phi = 1.5*(( r*r + r)/( r*r + r + 1));

	return (phi);
}


inline double phi(double r, int method)
{
	double ret;

		
	switch(method)
	{	
		case 0:
			ret = OFLimiter(r);
			break;

		case 1:
			ret = SBee(r);
			break;

		case 2:
			ret = minmod(r);		
			break;

		case 3:
			ret = VLeer(r);		
			break;
			
		case -1:
			ret =1;
			break;
	}
	return(ret);
}









void CheckPhysics(GRID HydroGrid, double tau, int stage )
{
	int i,j,k;
	bool error = false;

	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{

		double vx = HydroGrid[i][j][k].Vx;
		double vy = HydroGrid[i][j][k].Vy;
		double ve = HydroGrid[i][j][k].Ve;
		double b = sqrt(vx*vx+vy*vy+ tau*tau*ve*ve);
		double en = HydroGrid[i][j][k].En;
		if(b >=1)
		{
		 	cout<<"stage "<< stage<< " Causality violated ---> ( "<< vx<< ", "<<vy<< ","<<tau*ve<< ") @ "<<  HydroGrid[i][j][k].X<< ", "<< HydroGrid[i][j][k].Y<< ","<< HydroGrid[i][j][k].eta<<endl;
			 error=true;
		}		
		if(en<=0)
		{
			 cout<<"stage "<< stage<< " Energy Negative "<<en<< " -> "<< HydroGrid[i][j][k].Y<< ", "<<HydroGrid[i][j][k].Y<< ", "<<HydroGrid[i][j][k].eta<<endl;
			 error=true;
		}		
	}
	//~ MPI_Barrier(mpi_grid);
}

 


inline double MaxTempGev(GRID HydroGrid)
{
	int i,j,k;
	
	double tmax=0;
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		DECLePPIa;		
		double temp =   GEVFM* FT(HydroGrid[i][j][k].En );// FT returns in 1/fm 
		
		if(temp>tmax)
			tmax=temp;		
	}	
		
	double globaltmax;
	MPI_Allreduce( &tmax, &globaltmax, 1,MPI_DOUBLE, MPI_MAX,mpi_grid); 	
	
	return globaltmax;
}

inline double FluctuatingEn(double x, double y)
{
	double ret;	
	double width = 1;
	double sharp=1.0/pow(width,2);
	double E =((double)2.7182818284590452353602874713526624977572470937000);

	ret = 7.44478*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(-6.41375 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(-6.41375 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(1.28275 + x,2) - pow(-5.68075 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.64925 + x,2) - pow(-5.68075 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(1.64925 + x,2) - pow(-5.31425 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(2.01575 + x,2) - pow(-5.31425 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(2.38225 + x,2) - pow(-5.31425 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(-4.94775 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(-4.94775 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(-4.94775 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(-4.94775 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.91625 + x,2) - pow(-4.94775 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.28275 + x,2) - pow(-4.94775 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.64925 + x,2) - pow(-4.94775 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(2.38225 + x,2) - pow(-4.94775 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(3.11525 + x,2) - pow(-4.94775 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(-4.58125 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(-4.58125 + y,2))) + 
   44.6687*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(-4.58125 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(-4.58125 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(-4.58125 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(-4.58125 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.18325 + x,2) - pow(-4.58125 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.54975 + x,2) - pow(-4.58125 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(0.91625 + x,2) - pow(-4.58125 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(-4.21475 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(-4.21475 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(-4.21475 + y,2))) + 
   37.2239*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(-4.21475 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(-4.21475 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(-4.21475 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.18325 + x,2) - pow(-4.21475 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.54975 + x,2) - pow(-4.21475 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(0.91625 + x,2) - pow(-4.21475 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(1.28275 + x,2) - pow(-4.21475 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(-3.84825 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(-3.84825 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(-3.84825 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(-3.84825 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(-3.84825 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(-3.84825 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(0.54975 + x,2) - pow(-3.84825 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(-3.48175 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(-3.48175 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(-3.48175 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(-3.48175 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(-3.48175 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.18325 + x,2) - pow(-3.48175 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(0.54975 + x,2) - pow(-3.48175 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.91625 + x,2) - pow(-3.48175 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.64925 + x,2) - pow(-3.48175 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(-3.11525 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(-3.11525 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(-3.11525 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(-3.11525 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(-3.11525 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(0.91625 + x,2) - pow(-3.11525 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(1.28275 + x,2) - pow(-3.11525 + y,2))) + 
   44.6687*pow(E,sharp*(-pow(1.64925 + x,2) - pow(-3.11525 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(3.84825 + x,2) - pow(-3.11525 + y,2))) + 
   52.1135*pow(E,sharp*(-pow(4.21475 + x,2) - pow(-3.11525 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(-2.74875 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(-2.74875 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(-2.74875 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(-2.74875 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(0.18325 + x,2) - pow(-2.74875 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.54975 + x,2) - pow(-2.74875 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.91625 + x,2) - pow(-2.74875 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(4.21475 + x,2) - pow(-2.74875 + y,2))) + 
   52.1135*pow(E,sharp*(-pow(4.58125 + x,2) - pow(-2.74875 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(-2.38225 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(-2.38225 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(-2.38225 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(-2.38225 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(0.54975 + x,2) - pow(-2.38225 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.91625 + x,2) - pow(-2.38225 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(3.48175 + x,2) - pow(-2.38225 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(3.84825 + x,2) - pow(-2.38225 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(4.21475 + x,2) - pow(-2.38225 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(4.58125 + x,2) - pow(-2.38225 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-5.31425 + x,2) - pow(-2.01575 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-4.94775 + x,2) - pow(-2.01575 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-3.11525 + x,2) - pow(-2.01575 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(-2.01575 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(-2.01575 + y,2))) + 
   37.2239*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(-2.01575 + y,2))) + 
   81.8926*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(-2.01575 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(-2.01575 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(-2.01575 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.18325 + x,2) - pow(-2.01575 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.54975 + x,2) - pow(-2.01575 + y,2))) + 
   37.2239*pow(E,sharp*(-pow(0.91625 + x,2) - pow(-2.01575 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.28275 + x,2) - pow(-2.01575 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(3.84825 + x,2) - pow(-2.01575 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(4.21475 + x,2) - pow(-2.01575 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(4.58125 + x,2) - pow(-2.01575 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(4.94775 + x,2) - pow(-2.01575 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-4.94775 + x,2) - pow(-1.64925 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-3.84825 + x,2) - pow(-1.64925 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-3.48175 + x,2) - pow(-1.64925 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-3.11525 + x,2) - pow(-1.64925 + y,2))) + 
   81.8926*pow(E,sharp*(-pow(-2.74875 + x,2) - pow(-1.64925 + y,2))) + 
   81.8926*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(-1.64925 + y,2))) + 
   89.3374*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(-1.64925 + y,2))) + 
   37.2239*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(-1.64925 + y,2))) + 
   67.0031*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(-1.64925 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(-1.64925 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(-1.64925 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(-1.64925 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.54975 + x,2) - pow(-1.64925 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.91625 + x,2) - pow(-1.64925 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(1.28275 + x,2) - pow(-1.64925 + y,2))) + 
   37.2239*pow(E,sharp*(-pow(1.64925 + x,2) - pow(-1.64925 + y,2))) + 
   44.6687*pow(E,sharp*(-pow(4.21475 + x,2) - pow(-1.64925 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-3.84825 + x,2) - pow(-1.28275 + y,2))) + 
   37.2239*pow(E,sharp*(-pow(-3.48175 + x,2) - pow(-1.28275 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-3.11525 + x,2) - pow(-1.28275 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-2.74875 + x,2) - pow(-1.28275 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(-1.28275 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(-1.28275 + y,2))) + 
   44.6687*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(-1.28275 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(-1.28275 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(-1.28275 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(-1.28275 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(-1.28275 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(1.28275 + x,2) - pow(-1.28275 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.64925 + x,2) - pow(-1.28275 + y,2))) + 
   74.4478*pow(E,sharp*(-pow(2.01575 + x,2) - pow(-1.28275 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(3.84825 + x,2) - pow(-1.28275 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(4.21475 + x,2) - pow(-1.28275 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-3.48175 + x,2) - pow(-0.91625 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-3.11525 + x,2) - pow(-0.91625 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(-0.91625 + y,2))) + 
   44.6687*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(-0.91625 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(-0.91625 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(-0.91625 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(-0.91625 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(-0.91625 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.18325 + x,2) - pow(-0.91625 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.54975 + x,2) - pow(-0.91625 + y,2))) + 
   59.5583*pow(E,sharp*(-pow(2.01575 + x,2) - pow(-0.91625 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(2.38225 + x,2) - pow(-0.91625 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(4.21475 + x,2) - pow(-0.91625 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(4.94775 + x,2) - pow(-0.91625 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(-0.54975 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(-0.54975 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(-0.54975 + y,2))) + 
   52.1135*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(-0.54975 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.18325 + x,2) - pow(-0.54975 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.54975 + x,2) - pow(-0.54975 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.64925 + x,2) - pow(-0.54975 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(2.01575 + x,2) - pow(-0.54975 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(4.58125 + x,2) - pow(-0.54975 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(4.94775 + x,2) - pow(-0.54975 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(5.31425 + x,2) - pow(-0.54975 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-3.11525 + x,2) - pow(-0.18325 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-2.74875 + x,2) - pow(-0.18325 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(-0.18325 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(-0.18325 + y,2))) + 
   44.6687*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(-0.18325 + y,2))) + 
   59.5583*pow(E,sharp*(-pow(0.18325 + x,2) - pow(-0.18325 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(0.54975 + x,2) - pow(-0.18325 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.91625 + x,2) - pow(-0.18325 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.64925 + x,2) - pow(-0.18325 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(2.01575 + x,2) - pow(-0.18325 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(4.58125 + x,2) - pow(-0.18325 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(4.94775 + x,2) - pow(-0.18325 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(5.31425 + x,2) - pow(-0.18325 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-3.11525 + x,2) - pow(0.18325 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-2.74875 + x,2) - pow(0.18325 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(0.18325 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(0.18325 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(0.18325 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(0.18325 + y,2))) + 
   52.1135*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(0.18325 + y,2))) + 
   44.6687*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(0.18325 + y,2))) + 
   44.6687*pow(E,sharp*(-pow(0.18325 + x,2) - pow(0.18325 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(0.54975 + x,2) - pow(0.18325 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.91625 + x,2) - pow(0.18325 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(1.28275 + x,2) - pow(0.18325 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(1.64925 + x,2) - pow(0.18325 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(2.01575 + x,2) - pow(0.18325 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(3.48175 + x,2) - pow(0.18325 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(4.21475 + x,2) - pow(0.18325 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(4.58125 + x,2) - pow(0.18325 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(4.94775 + x,2) - pow(0.18325 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(6.41375 + x,2) - pow(0.18325 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-3.84825 + x,2) - pow(0.54975 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-3.11525 + x,2) - pow(0.54975 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.74875 + x,2) - pow(0.54975 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(0.54975 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(0.54975 + y,2))) + 
   52.1135*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(0.54975 + y,2))) + 
   52.1135*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(0.54975 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(0.54975 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(0.18325 + x,2) - pow(0.54975 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(0.54975 + x,2) - pow(0.54975 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.91625 + x,2) - pow(0.54975 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(1.28275 + x,2) - pow(0.54975 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(1.64925 + x,2) - pow(0.54975 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(2.01575 + x,2) - pow(0.54975 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(3.11525 + x,2) - pow(0.54975 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(3.48175 + x,2) - pow(0.54975 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(4.21475 + x,2) - pow(0.54975 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(4.94775 + x,2) - pow(0.54975 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(5.31425 + x,2) - pow(0.54975 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(5.68075 + x,2) - pow(0.54975 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(6.04725 + x,2) - pow(0.54975 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(6.41375 + x,2) - pow(0.54975 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-3.11525 + x,2) - pow(0.91625 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-2.74875 + x,2) - pow(0.91625 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(0.91625 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(0.91625 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(0.91625 + y,2))) + 
   37.2239*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(0.91625 + y,2))) + 
   59.5583*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(0.91625 + y,2))) + 
   37.2239*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(0.91625 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(0.91625 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.18325 + x,2) - pow(0.91625 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.54975 + x,2) - pow(0.91625 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.91625 + x,2) - pow(0.91625 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.28275 + x,2) - pow(0.91625 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.64925 + x,2) - pow(0.91625 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(2.01575 + x,2) - pow(0.91625 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(2.38225 + x,2) - pow(0.91625 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(2.74875 + x,2) - pow(0.91625 + y,2))) + 
   52.1135*pow(E,sharp*(-pow(3.11525 + x,2) - pow(0.91625 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(5.68075 + x,2) - pow(0.91625 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(6.04725 + x,2) - pow(0.91625 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(6.41375 + x,2) - pow(0.91625 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-3.84825 + x,2) - pow(1.28275 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-3.48175 + x,2) - pow(1.28275 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-3.11525 + x,2) - pow(1.28275 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-2.74875 + x,2) - pow(1.28275 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(1.28275 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(1.28275 + y,2))) + 
   44.6687*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(1.28275 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(1.28275 + y,2))) + 
   37.2239*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(1.28275 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(1.28275 + y,2))) + 
   37.2239*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(1.28275 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.18325 + x,2) - pow(1.28275 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.54975 + x,2) - pow(1.28275 + y,2))) + 
   37.2239*pow(E,sharp*(-pow(1.28275 + x,2) - pow(1.28275 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.64925 + x,2) - pow(1.28275 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(2.01575 + x,2) - pow(1.28275 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(2.38225 + x,2) - pow(1.28275 + y,2))) + 
   52.1135*pow(E,sharp*(-pow(2.74875 + x,2) - pow(1.28275 + y,2))) + 
   52.1135*pow(E,sharp*(-pow(3.11525 + x,2) - pow(1.28275 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(3.48175 + x,2) - pow(1.28275 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(3.84825 + x,2) - pow(1.28275 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(5.68075 + x,2) - pow(1.28275 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-3.48175 + x,2) - pow(1.64925 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-3.11525 + x,2) - pow(1.64925 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.74875 + x,2) - pow(1.64925 + y,2))) + 
   52.1135*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(1.64925 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(1.64925 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(1.64925 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(1.64925 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(1.64925 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(1.64925 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.28275 + x,2) - pow(1.64925 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(1.64925 + x,2) - pow(1.64925 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(2.38225 + x,2) - pow(1.64925 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(3.11525 + x,2) - pow(1.64925 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(3.84825 + x,2) - pow(1.64925 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(4.21475 + x,2) - pow(1.64925 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-3.48175 + x,2) - pow(2.01575 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-3.11525 + x,2) - pow(2.01575 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-2.74875 + x,2) - pow(2.01575 + y,2))) + 
   59.5583*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(2.01575 + y,2))) + 
   44.6687*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(2.01575 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(2.01575 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(2.01575 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(2.01575 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.28275 + x,2) - pow(2.01575 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(1.64925 + x,2) - pow(2.01575 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(4.58125 + x,2) - pow(2.01575 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(4.94775 + x,2) - pow(2.01575 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-3.84825 + x,2) - pow(2.38225 + y,2))) + 
   37.2239*pow(E,sharp*(-pow(-3.11525 + x,2) - pow(2.38225 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-2.74875 + x,2) - pow(2.38225 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(2.38225 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(2.38225 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(2.38225 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(2.38225 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(2.38225 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(2.38225 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.28275 + x,2) - pow(2.38225 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.64925 + x,2) - pow(2.38225 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(3.84825 + x,2) - pow(2.38225 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(4.58125 + x,2) - pow(2.38225 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-3.48175 + x,2) - pow(2.74875 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(2.74875 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(2.74875 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(2.74875 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(2.74875 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(2.74875 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(2.74875 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(0.91625 + x,2) - pow(2.74875 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(1.28275 + x,2) - pow(2.74875 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(1.64925 + x,2) - pow(2.74875 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(4.21475 + x,2) - pow(2.74875 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(4.58125 + x,2) - pow(2.74875 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(4.94775 + x,2) - pow(2.74875 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-3.48175 + x,2) - pow(3.11525 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(3.11525 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(3.11525 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(3.11525 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(3.11525 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(0.18325 + x,2) - pow(3.11525 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(0.54975 + x,2) - pow(3.11525 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(0.91625 + x,2) - pow(3.11525 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(1.28275 + x,2) - pow(3.11525 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(1.64925 + x,2) - pow(3.11525 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(2.01575 + x,2) - pow(3.11525 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(4.21475 + x,2) - pow(3.11525 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(4.58125 + x,2) - pow(3.11525 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.74875 + x,2) - pow(3.48175 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(3.48175 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-1.64925 + x,2) - pow(3.48175 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(3.48175 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(3.48175 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.54975 + x,2) - pow(3.48175 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(3.48175 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.18325 + x,2) - pow(3.48175 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.54975 + x,2) - pow(3.48175 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.91625 + x,2) - pow(3.48175 + y,2))) + 
   37.2239*pow(E,sharp*(-pow(1.28275 + x,2) - pow(3.48175 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.64925 + x,2) - pow(3.48175 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(2.01575 + x,2) - pow(3.48175 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(2.38225 + x,2) - pow(3.48175 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.74875 + x,2) - pow(3.84825 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(3.84825 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(3.84825 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-1.28275 + x,2) - pow(3.84825 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-0.91625 + x,2) - pow(3.84825 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-0.18325 + x,2) - pow(3.84825 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.18325 + x,2) - pow(3.84825 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.54975 + x,2) - pow(3.84825 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.91625 + x,2) - pow(3.84825 + y,2))) + 
   29.7791*pow(E,sharp*(-pow(1.28275 + x,2) - pow(3.84825 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.64925 + x,2) - pow(3.84825 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(2.38225 + x,2) - pow(3.84825 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.74875 + x,2) - pow(4.21475 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(4.21475 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.18325 + x,2) - pow(4.21475 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.54975 + x,2) - pow(4.21475 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.91625 + x,2) - pow(4.21475 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(1.28275 + x,2) - pow(4.21475 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(1.64925 + x,2) - pow(4.21475 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(2.01575 + x,2) - pow(4.21475 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(3.11525 + x,2) - pow(4.21475 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-3.11525 + x,2) - pow(4.58125 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(4.58125 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.18325 + x,2) - pow(4.58125 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(0.91625 + x,2) - pow(4.58125 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(1.28275 + x,2) - pow(4.58125 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.64925 + x,2) - pow(4.58125 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(2.38225 + x,2) - pow(4.58125 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(3.11525 + x,2) - pow(4.58125 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(3.48175 + x,2) - pow(4.58125 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.74875 + x,2) - pow(4.94775 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(4.94775 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.28275 + x,2) - pow(4.94775 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(1.64925 + x,2) - pow(4.94775 + y,2))) + 
   22.3344*pow(E,sharp*(-pow(2.01575 + x,2) - pow(4.94775 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(2.38225 + x,2) - pow(4.94775 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(2.74875 + x,2) - pow(4.94775 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(3.48175 + x,2) - pow(4.94775 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.74875 + x,2) - pow(5.31425 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(5.31425 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(5.31425 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(3.11525 + x,2) - pow(5.31425 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(4.21475 + x,2) - pow(5.31425 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.74875 + x,2) - pow(5.68075 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(5.68075 + y,2))) + 
   14.8896*pow(E,sharp*(-pow(-2.01575 + x,2) - pow(5.68075 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(2.01575 + x,2) - pow(5.68075 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(3.11525 + x,2) - pow(5.68075 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(-2.38225 + x,2) - pow(6.04725 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(2.38225 + x,2) - pow(6.04725 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(3.48175 + x,2) - pow(6.04725 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(1.64925 + x,2) - pow(6.41375 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.54975 + x,2) - pow(6.78025 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.91625 + x,2) - pow(6.78025 + y,2))) + 
   7.44478*pow(E,sharp*(-pow(0.54975 + x,2) - pow(7.14675 + y,2)));


	
	return( ret);

}
