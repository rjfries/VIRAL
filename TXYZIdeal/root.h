


#define NOVAR 4


void CheckRoot(GRID HydroGrid, double tau)
{
	int i,j,k;

	double max=0;

	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		double Ttt = HydroGrid[i][j][k].T00;
		double Ttx = HydroGrid[i][j][k].T10;
		double Tty = HydroGrid[i][j][k].T20;
		double Ttz = HydroGrid[i][j][k].T30;
		double eps = HydroGrid[i][j][k].En; 
		
		double M2 = (Ttx)*(Ttx) + (Tty)*(Tty) +  (Ttz)*(Ttz);

		double res =  eps - (Ttt) + (M2 / ((Ttt) + EOS(eps)));

		if(fabs(res)>max)
			max = fabs(res);			
	}

	if(rank==root)
		cout<<"root check"<<endl;
	MPI_Barrier(mpi_grid);
	
	cout<<std::scientific<<max<<"  " <<rank<<endl;


}


struct epsParams
{
	GRID HydroGrid;
	int i;
	int j;
	int k;
	double tau;
};


double feps(double eps, void * params)
{

 	struct epsParams *p  = (struct epsParams *) params;		

	GRID HG = p->HydroGrid;

	int i = p->i;
	int j = p->j;
	int k = p->k; 

	double Ttt = HG[i][j][k].T00;
	double Ttx = HG[i][j][k].T10;
	double Tty = HG[i][j][k].T20;
	double Ttz = HG[i][j][k].T30; 
	double M2 = Ttx*Ttx + Tty*Tty +  Ttz*Ttz;  
	
	double res = eps -  Ttt  + ( M2 / ( Ttt + EOS(eps )) ); 
	return(	res) ;

}



double dfeps(double eps, void * params)
{

 	struct epsParams *p 
    = (struct epsParams *) params;		

	GRID HG = p->HydroGrid;

	int i = p->i;
	int j = p->j;
	int k = p->k; 

	double en = HG[i][j][k].En;
	double Ttt = HG[i][j][k].T00;
	double Ttx = HG[i][j][k].T10;
	double Tty = HG[i][j][k].T20;
	double Ttz = HG[i][j][k].T30;
	double M2 = Ttx*Ttx + Tty*Tty +   Ttz*Ttz; 
	double res = 1 -( M2 / pow( Ttt + EOS(eps ), 2) )*DPDE(en ) ;
	return(	res) ;
}




void fdfeps(double eps, void * params, double *f, double *df)
{
	*f = feps (eps, params);
	*df = dfeps (eps, params);
}






void RootSearchForEnVelUsingDerivatives(GRID HydroGrid , double tau)
{

	int i,j,k;
	
	int status;
	
	double r = 0;

	double eps_lo = 0, eps_hi = 2000;
	struct epsParams params ;
	
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{

		const gsl_root_fdfsolver_type *T;
		gsl_root_fdfsolver *s;
		gsl_function_fdf FDF;
 
		FDF.f = &feps;
		FDF.df = &dfeps;
		FDF.fdf = &fdfeps;
		params.HydroGrid = HydroGrid;
		params.i = i;
		params.j = j;
		params.k = k;
		params.tau = tau;
		FDF.params = &params;
		
		
		double guess = HydroGrid[i][j][k].En;		
		r=guess;
		
		
		
        T = gsl_root_fdfsolver_newton; 
		s = gsl_root_fdfsolver_alloc (T);		
		gsl_root_fdfsolver_set (s, &FDF, guess);
		
		int iter = 0, max_iter = 2000;
		
		do
	    {
			iter++;
			status = gsl_root_fdfsolver_iterate (s);
			double r0 = r;
			r = gsl_root_fdfsolver_root (s);
			status = gsl_root_test_delta (r, r0, 0, 1e-20);
	                                       
			if(status != GSL_SUCCESS && status != GSL_CONTINUE )
			{	
				printf ("Root finding failure\n");
				cout<<"error code  --> "<< status<<endl;
				exit(1);
			}
			
	    } while (status == GSL_CONTINUE && iter < max_iter);
	
		HydroGrid[i][j][k].En = r;    //find the new energy from root finding algorithm
		HydroGrid[i][j][k].P = EOS(HydroGrid[i][j][k].En);    //find the new energy from root finding algorithm
		gsl_root_fdfsolver_free (s);

	}
	
	int quit =0;

	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		
		double Ttt = HydroGrid[i][j][k].T00;
		double Ttx = HydroGrid[i][j][k].T10;
		double Tty = HydroGrid[i][j][k].T20;
		double Ttz = HydroGrid[i][j][k].T30; 
		double P = HydroGrid[i][j][k].P ;
		 
		HydroGrid[i][j][k].Vx = Ttx/(Ttt + P ); 
		HydroGrid[i][j][k].Vy = Tty/(Ttt + P ); 
		HydroGrid[i][j][k].Ve = Ttz/(Ttt + P );  
				
		HydroGrid[i][j][k].u[0]= 1.0/(sqrt(1 - HydroGrid[i][j][k].Vx*HydroGrid[i][j][k].Vx
											 - HydroGrid[i][j][k].Vy*HydroGrid[i][j][k].Vy
											 -  HydroGrid[i][j][k].Ve*HydroGrid[i][j][k].Ve
											 )
									 );
									 
		HydroGrid[i][j][k].u[1] =  HydroGrid[i][j][k].u[0]*HydroGrid[i][j][k].Vx;
		HydroGrid[i][j][k].u[2] =  HydroGrid[i][j][k].u[0]*HydroGrid[i][j][k].Vy;
		HydroGrid[i][j][k].u[3] =  HydroGrid[i][j][k].u[0]*HydroGrid[i][j][k].Ve;			
		 
	}
}
