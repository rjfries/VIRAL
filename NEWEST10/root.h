void CheckRoot(GRID HydroGrid , double tau)
{
	int i,j,k,l;


	
	double maxTT[6]={0,0,0,0,0,0};
	double sumTT[6]={0,0,0,0,0,0};
	
	double inavg[6]={0,0,0,0,0,0};
	double outavg[6]={0,0,0,0,0,0};
	
	struct {
	double val;
	int   rank;
	} in[6], out[6];
	
	int count =0;
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{
		count++;
		double Ttt = HydroGrid[i][j][k].T00;
		double Ttx = HydroGrid[i][j][k].T10;
		double Tty = HydroGrid[i][j][k].T20;
		double Tte = HydroGrid[i][j][k].T30;
		
		double eps = HydroGrid[i][j][k].En; 
		double PI  = HydroGrid[i][j][k].PI;
		
		double pi1 = HydroGrid[i][j][k].pi[0];
		double pi2 = HydroGrid[i][j][k].pi[1];
		double pi3 = HydroGrid[i][j][k].pi[2];
		double pi4 = HydroGrid[i][j][k].pi[3];
		
		double r = HydroGrid[i][j][k].r;

		double M2 = (Ttx-pi2)*(Ttx-pi2) + (Tty-pi3)*(Tty-pi3) + tau*tau*(Tte-pi4)*(Tte-pi4);
	
		
		double res =  eps - (Ttt-pi1)  + ( M2 / ( (Ttt-pi1) + EOS(eps,r)- PI));
		
		
		if(fabs(res)>maxTT[0])
			maxTT[0] = fabs(res);
		
		sumTT[0] +=  fabs(res);					
	}
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++)
	{ 
		double u0 = HydroGrid[i][j][k].u[0];
		double u1 = HydroGrid[i][j][k].u[1];
		double u2 = HydroGrid[i][j][k].u[2];
		double u3 = HydroGrid[i][j][k].u[3];
		
		double pi00 = HydroGrid[i][j][k].pi[0];
		double pi01 = HydroGrid[i][j][k].pi[1];
		double pi02 = HydroGrid[i][j][k].pi[2];
		double pi03 = HydroGrid[i][j][k].pi[3];
		
		double pi11 = HydroGrid[i][j][k].pi[4];
		double pi12 = HydroGrid[i][j][k].pi[5];
		double pi13 = HydroGrid[i][j][k].pi[6];
		
		double pi22 = HydroGrid[i][j][k].pi[7];
		double pi23 = HydroGrid[i][j][k].pi[8];
		
		double pi33 = HydroGrid[i][j][k].pi[9];
		
		double r = HydroGrid[i][j][k].r;

		double row0 = pi00*u0 - pi01*u1- pi02*u2 - pi03*u3*tau*tau;
		double row1 = pi01*u0 - pi11*u1- pi12*u2 - pi13*u3*tau*tau;
		double row2 = pi02*u0 - pi12*u1- pi22*u2 - pi23*u3*tau*tau;
		double row3 = pi03*u0 - pi13*u1- pi23*u2 - pi33*u3*tau*tau;
		
		double pi44[4][4]= {{pi00,pi01,pi02,pi03},{pi01,pi11,pi12,pi13},{pi02,pi12,pi22,pi23},{pi03,pi13,pi23,pi33}};
		double ul[4] = {u0,-u1,-u2,-tau*tau*u3};
		
		
		double con2 = 0;
		
		for(int x=0;x<4;x++) 
		for(int y=0;y<4;y++) 
			con2 += (pi44[x][y]*ul[x]*ul[y]);
			
		
		double trace = pi00 - pi11- pi22 - pi33*tau*tau;
		
		//~ if(fabs(row3)>0.01)
			//~ cout<<row0<<"  "<<row1<<"  "<<row2<<"  "<<row3<<"  "<<HydroGrid[i][j][k].X<<"  "<<HydroGrid[i][j][k].Y<<"  "<<HydroGrid[i][j][k].eta<<endl;
			
		
	
		if(fabs(row0) > maxTT[1])
			maxTT[1] = fabs(row0);
						
			sumTT[1] += fabs(row0);			
		
		if(fabs(row1) > maxTT[2])
			maxTT[2] = fabs(row1);
						
			sumTT[2] += fabs(row1);			
		
		if(fabs(row2) > maxTT[3])
			maxTT[3] = fabs(row2);	
					
			sumTT[3] += fabs(row2);			
		
		if(fabs(row3) > maxTT[4])
			maxTT[4] = fabs(row3);	
						
			sumTT[4] += fabs(row3);			
		
		if(fabs(trace) > maxTT[5])
			maxTT[5] = fabs(trace);	
					
			sumTT[5] += fabs(trace);	
			
			
		HydroGrid[i][j][k].temp[0] =  row0;		
		HydroGrid[i][j][k].temp[1] =  HydroGrid[i][j][k].En*u0;		
		HydroGrid[i][j][k].temp[2] =  row1;		
		HydroGrid[i][j][k].temp[3] =  HydroGrid[i][j][k].En*u1;			
		HydroGrid[i][j][k].temp[4] =  row2;
		HydroGrid[i][j][k].temp[5] =  HydroGrid[i][j][k].En*u2;	
		HydroGrid[i][j][k].temp[6] =  row3;		
		HydroGrid[i][j][k].temp[7] =  HydroGrid[i][j][k].En*u3;			
		HydroGrid[i][j][k].temp[8] =  con2;		
		HydroGrid[i][j][k].temp[9] =  HydroGrid[i][j][k].En;		
	}


	for(l=0;l<6;l++)
	{
		in[l].val = maxTT[l];
		in[l].rank = rank; 
		sumTT[l] /= count;
		inavg[l] = sumTT[l];
	}
	
	for(l=0;l<6;l++)
		MPI_Allreduce(&in[l],&out[l],1,MPI_DOUBLE_INT,MPI_MAXLOC,mpi_grid);
	
	for(l=0;l<6;l++)
		MPI_Allreduce(&inavg[l],&outavg[l],1,MPI_DOUBLE,MPI_SUM,mpi_grid); 
	 
	if(rank==out[0].rank)
	{
		cout<<endl<<endl;
		cout<<std::scientific<<"ROOT CHECK -- MAX "<< maxTT[0]<<" &  Avg "<< outavg[0]/NP<<endl;
		//~ fflush(stdout);
	}	
	
	
	MPI_Barrier(mpi_grid);
	
	if(rank==out[1].rank)
	{
		cout<<std::scientific<<"ROW 0 TRANSVERSALITY -- MAX "<< maxTT[1]<<" &  Avg "<< outavg[1]/NP<<endl;
		//~ fflush(stdout);
	}	
	
	MPI_Barrier(mpi_grid);
	if(rank==out[2].rank)
	{
		cout<<std::scientific<<"ROW 1 TRANSVERSALITY-- MAX "<< maxTT[2]<<" &  Avg "<< outavg[2]/NP<<endl;
		//~ fflush(stdout);
	}	
	
	MPI_Barrier(mpi_grid);
	if(rank==out[3].rank)
	{
		cout<<std::scientific<<"ROW 2 TRANSVERSALITY-- MAX "<< maxTT[3]<<" &  Avg "<< outavg[3]/NP<<endl;
		//~ fflush(stdout);
	}	
	
	MPI_Barrier(mpi_grid);
	if(rank==out[4].rank)
	{
		cout<<std::scientific<<"ROW 3 TRANSVERSALITY -- MAX "<< maxTT[4]<<" &  Avg "<< outavg[4]/NP<<endl;
		//~ fflush(stdout);
	}	
	
	MPI_Barrier(mpi_grid);
	if(rank==out[5].rank)
	{
		cout<<std::scientific<<"TRACE CHECK -- MAX "<< maxTT[5]<<" &  Avg "<< outavg[5]/NP<<endl;
		fflush(stdout);
	}	

	MPI_Barrier(mpi_grid);
	fflush(stdout);
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
	double tau = p->tau;

	double Ttt = HG[i][j][k].T00;
	double Ttx = HG[i][j][k].T10;
	double Tty = HG[i][j][k].T20;
	double Tte = HG[i][j][k].T30;
	 
	double PI  = HG[i][j][k].PI;
	
	double pi1 = HG[i][j][k].pi[0];
	double pi2 = HG[i][j][k].pi[1];
	double pi3 = HG[i][j][k].pi[2];
	double pi4 = HG[i][j][k].pi[3];
	
	double r = HydroGrid[i][j][k].r;

	double M2 = (Ttx-pi2)*(Ttx-pi2) + (Tty-pi3)*(Tty-pi3) + tau*tau*(Tte-pi4)*(Tte-pi4);
	
	double res =  eps - (Ttt-pi1)  + ( M2 / ( (Ttt-pi1) + EOS(eps,r)- PI));

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
	
	
	double PI = HG[i][j][k].PI;

	double pi1 = HG[i][j][k].pi[0];
	double pi2 = HG[i][j][k].pi[1];
	double pi3 = HG[i][j][k].pi[2];
	double pi4 = HG[i][j][k].pi[3];



	
	double Ttt = HG[i][j][k].T00;
	double Ttx = HG[i][j][k].T10;
	double Tty = HG[i][j][k].T20;
	double Tte = HG[i][j][k].T30;


	double M2 = (Ttx-pi2)*(Ttx-pi2) + (Tty-pi3)*(Tty-pi3) + tau*tau*(Tte-pi4)*(Tte-pi4);
	double r = HG[i][j][k].r;
	
	double res =  1 - ( ( M2 / pow( Ttt - pi1 + EOS(eps , r) - PI, 2) )*DPDE(eps, r) );
	
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
		
		int iter = 0, max_iter = 5000;
		
		do
	    {
			iter++;
			status = gsl_root_fdfsolver_iterate (s);
			double r0 = r;
			r = gsl_root_fdfsolver_root (s);
			status = gsl_root_test_delta (r, r0, 0, 1e-30);
	                                       
			if(status != GSL_SUCCESS && status != GSL_CONTINUE )
			{	
				printf ("Root finding failure\n");
				cout<<"error code  --> "<< status<<endl;
				exit(1);
			}
			
	    } while (status == GSL_CONTINUE && iter < max_iter);
		gsl_root_fdfsolver_free (s);
		
		HydroGrid[i][j][k].En = r;    //find the new energy from root finding algorithm
		
		HydroGrid[i][j][k].Temp = FT(HydroGrid[i][j][k].En , HydroGrid[i][j][k].r);
		HydroGrid[i][j][k].P = 	EOS(HydroGrid[i][j][k].En , HydroGrid[i][j][k].r) ;		
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
	
		
		DECLcoord;
		DECLePPIa;
		double pi1 = HydroGrid[i][j][k].pi[0];
		double pi2 = HydroGrid[i][j][k].pi[1];
		double pi3 = HydroGrid[i][j][k].pi[2];
		double pi4 = HydroGrid[i][j][k].pi[3]; 
		
		double vx = (Ttx-pi2)/(Ttt - pi1 + P - PI);
		double vy = (Tty-pi3)/(Ttt - pi1 + P - PI);
		double ve = (Tte-pi4)/(Ttt - pi1 + P - PI);
		
		if(vx*vx+vy*vy+ve*ve >= 1 || (e <= 0) )
		{
			cout<<"HELL BROKE LOOSE vx "<<vx <<"  vy " << vy<< "  ve "<<ve << " & B--> "<<vx*vx+vy*vy+ve*ve<<endl;
			cout<<std::scientific<<"HELL BROKE LOOSE Energy "<<e<<endl;
			cout<<"HELL BROKE LOOSE X "<<X<<"  Y " << Y<< "  E "<<eta<< "  "<<endl;
			quit=1;
			continue;
		}	
			
		
		double g = 1.0/sqrt(1-vx*vx-vy*vy-tau*tau*ve*ve);

		
		
		HydroGrid[i][j][k].u[0] = g;
		HydroGrid[i][j][k].u[1] = g*vx;
		HydroGrid[i][j][k].u[2] = g*vy;
		HydroGrid[i][j][k].u[3] = g*ve;
		
		HydroGrid[i][j][k].Vx = vx;
		HydroGrid[i][j][k].Vy = vy;
		HydroGrid[i][j][k].Ve = ve;

	
	if(quit)
		exit(quit);
	}
}
