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



void CheckPhysics(GRID HydroGrid, int stage )
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
		double b = sqrt(vx*vx+vy*vy+ve*ve);
		double en = HydroGrid[i][j][k].En;
		if(b >=1)
		{
		 	cout<<"stage "<< stage<< " Causality violated ---> ( "<< vx<< ", "<<vy<< ","<<ve<< ") @ "<<  HydroGrid[i][j][k].X<< ", "<< HydroGrid[i][j][k].Y<< ","<< HydroGrid[i][j][k].eta<<endl;
			 error=true;
		}
		
		if(en<=0)
		{
			 cout<<"stage "<< stage<< " Energy Negative "<<en<< " -> "<< HydroGrid[i][j][k].Y<< ", "<<HydroGrid[i][j][k].Y<< ", "<<HydroGrid[i][j][k].eta<<endl;
			 error=true;
		}
		
	}
	MPI_Barrier(mpi_grid);
}

void ClearTemp(GRID HydroGrid)
{
	int i,j,k;

	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
		HydroGrid[i][j][k].Temp =0;
}


void ClearBuf(GRID HydroGrid)
{
	int i,j,k;

	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
		HydroGrid[i][j][k].Buf =0;
}
/*
void xder(GRID HydroGrid)
{

	int i, j,k;
	
	
	if(MYXLEFT==BOUNDARY)
	{
		for(j=0;j<YCM;j++)
		for(k=0;k<ZCM;k++)
		{
				for(i=0;i<BORDER;i++)
					HydroGrid[i][j][k].Fun = HydroGrid[il][j][k].Fun;
		}
	}
	
	if(MYXRIGHT==BOUNDARY)
	{
		for(j=0;j<YCM;j++)
		for(k=0;k<ZCM;k++)
		{
				for(i=0;i<BORDER;i++)
					HydroGrid[ir+i][j][k].Fun = HydroGrid[ir-1][j][k].Fun;
		}
	}
	

	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		for(i=0;i<XCM-1;i++)
			HydroGrid[i][j][k].Buf = genminmod(HydroGrid[i-1][j][k].Fun,HydroGrid[i][j][k].Fun,HydroGrid[i+1][j][k].Fun, XS,  1.1);
		
		HydroGrid[0][j][k].Buf=HydroGrid[XCM-1][j][k].Buf=0;
	}

}

void yder(GRID HydroGrid)
{

	int i, j,k;
	
	
	
	if(MYYLEFT==BOUNDARY)
	{
		for(i=0;i<XCM;i++)
		for(k=0;k<ZCM;k++)
		{
				for(j=0;j<BORDER;j++)
					HydroGrid[i][j][k].Fun = HydroGrid[i][jl][k].Fun;
		}
	}
	
	if(MYYRIGHT==BOUNDARY)
	{
		for(i=0;i<XCM;i++)
		for(k=0;k<ZCM;k++)
		{
				for(j=0;j<BORDER;j++)
					HydroGrid[i][jr+j][k].Fun = HydroGrid[i][jr-1][k].Fun;
		}
	}
	

	for(i=0;i<XCM;i++)
	for(k=0;k<ZCM;k++)
	{
		for(j=0;j<YCM-1;j++)
			HydroGrid[i][j][k].Buf = genminmod(HydroGrid[i][j-1][k].Fun,HydroGrid[i][j][k].Fun,HydroGrid[i][j+1][k].Fun, YS,  1.1);


		HydroGrid[i][0][k].Buf=HydroGrid[i][YCM-1][k].Buf=0;
	}
	


}



void zder(GRID HydroGrid)
{

	int i, j,k;

	if(MYZLEFT==BOUNDARY)
	{
		for(i=0;i<XCM;i++)
		for(j=0;j<YCM;j++)
		{
				for(k=0;k<BORDER;k++)
					HydroGrid[i][j][k].Fun =HydroGrid[i][j][kl].Fun;
		}
	}
	
	if(MYZRIGHT==BOUNDARY)
	{
		for(i=0;i<XCM;i++)
		for(j=0;j<YCM;j++)
		{
				for(k=0;k<BORDER;k++)
					HydroGrid[i][j][kr+k].Fun = HydroGrid[i][j][kr-1].Fun;
		}
	}
	

	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	{
		for(k=1;k<ZCM-1;k++)
			HydroGrid[i][j][k].Buf = genminmod(HydroGrid[i][j][k-1].Fun,HydroGrid[i][j][k].Fun,HydroGrid[i][j][k+1].Fun, ZS,  1.1);

		HydroGrid[i][j][0].Buf=HydroGrid[i][j][ZCM-1].Buf=0;
	}
}


#define SIM 0

void gsl_assisted_interp_derX(GRID HydroGrid)
{
	
	ClearBuf(HydroGrid);
	if(SIM)
	{
		xder(HydroGrid);
		return;
	}	
	
	int i,j,k;
	
	double Yaxis[XCM];

	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, XCM);
	
	gsl_interp_accel *Acc ;	
	
	Acc = gsl_interp_accel_alloc(); 

	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		for(i=0;i<XCM;i++)
			Yaxis[i] = HydroGrid[i][j][k].Fun;
 
		gsl_spline_init (spline, XC, Yaxis,XCM);

		for(i=0;i<XCM;i++)
			HydroGrid[i][j][k].Buf = gsl_spline_eval_deriv(spline, XC[i], Acc);
	}
}




void gsl_assisted_interp_derY(GRID HydroGrid)
{

	ClearBuf(HydroGrid);
	if(SIM)
	{
		yder(HydroGrid);
		return;
	}
		

	int i, j,k;
	double Yaxis[YCM];
	
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, YCM);
	gsl_interp_accel *Acc ; 
	Acc = gsl_interp_accel_alloc(); 
	
	for(i=0;i<XCM;i++)
	for(k=0;k<ZCM;k++)
	{
	 	for(j=0;j<YCM;j++)
	 		Yaxis[j] =   HydroGrid[i][j][k].Fun; 

		gsl_spline_init (spline, YC, Yaxis,YCM);
		
	  	for(j=0;j<YCM;j++)
			 HydroGrid[i][j][k].Buf = gsl_spline_eval_deriv(spline, YC[j], Acc);			
	}
}




void gsl_assisted_interp_derZ(GRID HydroGrid)
{
	
	ClearBuf(HydroGrid);

	if(SIM)
		return;
	
	int i, j,k;

	double Yaxis[ZCM];

	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, ZCM);

	gsl_interp_accel *Acc ;	
	
	Acc = gsl_interp_accel_alloc(); 


	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	{
	 	for(k=0;k<ZCM;k++)
	 		Yaxis[k] =  HydroGrid[i][j][k].Fun;
		
		gsl_spline_init(spline, EC, Yaxis, ZCM);
		
		for(k=0;k<ZCM;k++)
			HydroGrid[i][j][k].Buf = gsl_spline_eval_deriv(spline, EC[k], Acc);			
	}
} 
*/
inline int IsIrreleventPoint(int i, int j, int k)
{
	if( (i >= il && i<ir ) &&
		(j >= jl && j<jr ) &&
		(k >= kl && k<kr ) 
		)
		return 0;
	else
	if( (i < il || i>=ir ) &&
		(j >= jl && j<jr ) &&
		(k >= kl && k<kr ) 
		)
		return 0;
	else
	if( (j <  jl || j >= jr ) &&
		(i >= il && i < ir ) &&
		(k >= kl && k < kr ) 
		)
		return 0;
	else
	if(	(k < kl || k >= kr ) &&
		(i >= il && i < ir ) &&
		(j >= jl && j < jr )  
		)
		return 0;

	else
		return 1;
}



void CopyPrimaryVariablesToBufA(GRID HydroGrid)
{

	for( int i=0; i<XCM ; i++)
	for( int j=0; j<YCM ; j++)
	for( int k=0; k<ZCM ; k++)
	{
		HydroGrid[i][j][k].BufA[0]= tau*HydroGrid[i][j][k].T00;
		HydroGrid[i][j][k].BufA[1]= tau*HydroGrid[i][j][k].T10;
		HydroGrid[i][j][k].BufA[2]= tau*HydroGrid[i][j][k].T20;
		HydroGrid[i][j][k].BufA[3]= tau*HydroGrid[i][j][k].T30;
	}
}




void ClearResultVariable(GRID HydroGrid)
{
	for( int l=0; l<VARN; l++)
	for( int i=0; i<XCM ; i++)
	for( int j=0; j<YCM ; j++)
	for( int k=0; k<ZCM ; k++)
		HydroGrid[i][j][k].Result[l] = 0;
}

void AddPartialResultToFinalResult(GRID HydroGrid)
{
	for( int l=0; l<VARN; l++)
	for( int i=0; i<XCM ; i++)
	for( int j=0; j<YCM ; j++)
	for( int k=0; k<ZCM ; k++)
		HydroGrid[i][j][k].Result[l] += HydroGrid[i][j][k].PartialResult[l];
}



