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
		double Tte = HydroGrid[i][j][k].T30;
		double eps = HydroGrid[i][j][k].En;
		
		double X = HydroGrid[i][j][k].X;
		double Y = HydroGrid[i][j][k].Y;
		double r = HydroGrid[i][j][k].r;
		
		double M2 = (Ttx)*(Ttx) + (Tty)*(Tty) + tau*tau*(Tte)*(Tte);

		double res =  eps - (Ttt) + (M2 / ((Ttt) + EOS(eps, r)));

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
	double tau =p->tau;

	double Ttt = HG[i][j][k].T00;
	double Ttx = HG[i][j][k].T10;
	double Tty = HG[i][j][k].T20;
	double Tte = HG[i][j][k].T30;

	double M2 = Ttx*Ttx + Tty*Tty + tau*tau* Tte*Tte;
	


		double r = HydroGrid[i][j][k].r;

		
	double res = eps -  Ttt  + ( M2 / ( Ttt + EOS(eps , r)) );

	
	
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
	double tau =p->tau;
	
	
	 


	double en = HG[i][j][k].En;
	double Ttt = HG[i][j][k].T00;
	double Ttx = HG[i][j][k].T10;
	double Tty = HG[i][j][k].T20;
	double Tte = HG[i][j][k].T30;

	double M2 = Ttx*Ttx + Tty*Tty + tau*tau* Tte*Tte;
	double r = HydroGrid[i][j][k].r;

	double res = 1   - ( M2 / pow( Ttt + EOS(eps , r), 2) )*DPDE(en, r) ;

	
	
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

		//the params that need to be passed
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
	//	T = gsl_root_fdfsolver_secant;
	//	T = gsl_root_fdfsolver_steffenson;
		s = gsl_root_fdfsolver_alloc (T);		
		gsl_root_fdfsolver_set (s, &FDF, guess);
		
		int iter = 0, max_iter = 1000;
		
		do
	    {
			iter++;
			status = gsl_root_fdfsolver_iterate (s);
			double r0 = r;
			r = gsl_root_fdfsolver_root (s);
			status = gsl_root_test_delta (r, r0, 0, 1e-18);
	                                       
			if(status != GSL_SUCCESS && status != GSL_CONTINUE )
			{	
				printf ("Root finding failure\n");
				cout<<"error code  --> "<< status<<endl;
				exit(1);
			}
			
	    } while (status == GSL_CONTINUE && iter < max_iter);
	
		HydroGrid[i][j][k].En = r;    //find the new energy from root finding algorithm
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
		double Tte = HydroGrid[i][j][k].T30;
		double En = HydroGrid[i][j][k].En;

		double r = HydroGrid[i][j][k].r;	
		double P = EOS(En , r) ;
		
		HydroGrid[i][j][k].P = P;
		HydroGrid[i][j][k].Vx = Ttx/(Ttt + P );
		HydroGrid[i][j][k].Vy = Tty/(Ttt + P );
		HydroGrid[i][j][k].Ve = Tte/(Ttt + P );
	
		
	}
}
