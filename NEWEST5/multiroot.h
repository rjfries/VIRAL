


#define NOVAR 4


void CheckRoot(GRID HydroGrid , double tau)
{
	int i,j,k,l;

	
	double maxTT[9]={0,0,0,0,0,0,0,0,0};
		
	double sumTT[9]={0,0,0,0,0,0,0,0,0};
	
	double inavg[9]={0,0,0,0,0,0,0,0,0};
	double outavg[9]={0,0,0,0,0,0,0,0,0};
	
		
	struct {
	double val;
	int   rank;
	} in[9], out[9];
	
	int count=0;
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=0;k<ZCMA;k++)
	{
		count++;
		double W = HydroGrid[i][j][k].T00;
		double X = HydroGrid[i][j][k].T10;
		double Y = HydroGrid[i][j][k].T20;
		double Z = HydroGrid[i][j][k].T30;
 
		double PI = HydroGrid[i][j][k].PI;

		double r = HydroGrid[i][j][k].r;

		double e = HydroGrid[i][j][k].En; 
		
		DECLp5u4;

		double P = EOS(e );
		
		double y[4];
		
		y[0] = -(u1*X) - u2*Y - u3*Z*pow(tau,2) + (-e + W)*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),0.5);
		y[1] = -(e*u1) - u2*(p3 + (e + P - PI)*u1*u2) - u3*(p4 + (e + P - PI)*u1*u3)*pow(tau,2) - u1*(P + p1 - PI + (e + P - PI)*pow(u1,2)) + X*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),0.5);
		y[2] = -(e*u2) - u1*(p3 + (e + P - PI)*u1*u2) - u3*(p5 + (e + P - PI)*u2*u3)*pow(tau,2) - u2*(P + p2 - PI + (e + P - PI)*pow(u2,2)) + Y*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),0.5);
		y[3] = -(p4*u1) - p5*u2 - P*u3 + u3*(p1 + p2 + PI - W) + Z*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),0.5);
		
		for( l=0;l<4;l++)
		{
			if(fabs(y[l])>maxTT[l])
				maxTT[l] = fabs(y[l]);	

			sumTT[l] +=  fabs(y[l]);				
		}
	}

		
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=0;k<ZCMA;k++)
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
		
		double pi00 = -((2*p3*u1*u2 + p1*pow(u1,2) + p2*pow(u2,2) + pow(tau,2)*(2*(p4*u1 + p5*u2)*u3 - (p1 + p2)*pow(u3,2)))*pow(-pow(u0,2) + pow(tau,2)*pow(u3,2),-1));;
		double pi01 = (p1*u1 + p3*u2 + p4*u3*pow(tau,2))*pow(u0,-1);
		double pi02 = (p3*u1 + p2*u2 + p5*u3*pow(tau,2))*pow(u0,-1);
		double pi03 = pow(u0,-1)*(2*p3*u1*u2*u3 + (p4*u1 + p5*u2 - (p1 + p2)*u3)*pow(u0,2) + p1*u3*pow(u1,2) + p2*u3*pow(u2,2) + p4*u1*pow(tau,2)*pow(u3,2) + p5*u2*pow(tau,2)*pow(u3,2))*pow(pow(u0,2) - pow(tau,2)*pow(u3,2),-1);
		
		double pi11 = p1;
		double pi12 = p3;
		double pi13 = p4;
		
		double pi22 = p2;
		double pi23 = p5;
		
		double pi33 = -(pow(tau,-2)*(2*p3*u1*u2 + 2*(p4*u1 + p5*u2)*u3*pow(tau,2) - (p1 + p2)*pow(u0,2) + p1*pow(u1,2) + p2*pow(u2,2))*pow(-pow(u0,2) + pow(tau,2)*pow(u3,2),-1));;
		
		double r = HydroGrid[i][j][k].r;

		double row0 = pi00*u0 - pi01*u1- pi02*u2 - pi03*u3*tau*tau;
		double row1 = pi01*u0 - pi11*u1- pi12*u2 - pi13*u3*tau*tau;
		double row2 = pi02*u0 - pi12*u1- pi22*u2 - pi23*u3*tau*tau;
		double row3 = pi03*u0 - pi13*u1- pi23*u2 - pi33*u3*tau*tau;
		
		double trace = pi00 - pi11- pi22 - pi33*tau*tau;
		
		if(fabs(row0) > maxTT[4])
			maxTT[4] = fabs(row0);			
						
		sumTT[4] += fabs(row0);			
	
		if(fabs(row1) > maxTT[5])
			maxTT[5] = fabs(row1);			
		
		sumTT[5] += fabs(row1);		

		if(fabs(row2) > maxTT[6])
			maxTT[6] = fabs(row2);			
		
		sumTT[6] += fabs(row2);	

		if(fabs(row3) > maxTT[7])
			maxTT[7] = fabs(row3);			
		
		sumTT[7] += fabs(row3);			
		
		if(fabs(trace) > trace)
			maxTT[8] = fabs(trace);		

		sumTT[8] += fabs(trace);	
	}


	for(l=0;l<9;l++)
	{
		in[l].val = maxTT[l];
		in[l].rank = rank; 
		sumTT[l] /= count;
		inavg[l] = sumTT[l];
	}
	
	for(l=0;l<9;l++)
		MPI_Allreduce(&in[l],&out[l],1,MPI_DOUBLE_INT,MPI_MAXLOC,mpi_grid);
	
	for(l=0;l<9;l++)
		MPI_Allreduce(&inavg[l],&outavg[l],1,MPI_DOUBLE,MPI_SUM,mpi_grid);
	
	
	if(rank==out[0].rank)
	{
		cout<<endl<<endl;
		cout<<std::scientific<<"ROOT CHECK EQN 1 -- MAX "<< maxTT[0]<<" &     Avg "<< outavg[0]/NP<<endl;
		fflush(stdout);
	}	
	
	MPI_Barrier(mpi_grid);
	if(rank==out[1].rank)
	{ 
		cout<<std::scientific<<"ROOT CHECK EQN 2 -- MAX "<< maxTT[1]<<" &     Avg "<< outavg[1]/NP<<endl;
		fflush(stdout);
	}
	MPI_Barrier(mpi_grid);
		
	if(rank==out[2].rank)
	{ 
		cout<<std::scientific<<"ROOT CHECK EQN 3 -- MAX "<< maxTT[2]<<" &     Avg "<< outavg[2]/NP<<endl;
		fflush(stdout);
	}	
	MPI_Barrier(mpi_grid);
	
	if(rank==out[3].rank)
	{ 
		cout<<std::scientific<<"ROOT CHECK EQN 4 -- MAX "<< maxTT[3]<<" &     Avg "<< outavg[3]/NP<<endl;
		fflush(stdout);
	}	
		
	
	//~ MPI_Barrier(mpi_grid);
	
	//~ if(rank==out[4].rank)
	//~ {
		//~ cout<<std::scientific<<"ROW 0 TRANSVERSALITY -- MAX "<< maxTT[4]<<" &     Avg "<< outavg[4]/NP<<endl;
		//~ fflush(stdout);
	//~ }	
	//~ 
	//~ MPI_Barrier(mpi_grid);
	//~ if(rank==out[5].rank)
	//~ {
		//~ cout<<std::scientific<<"ROW 1 TRANSVERSALITY-- MAX "<< maxTT[5]<<" &     Avg "<< outavg[5]/NP<<endl;
		//~ fflush(stdout);
	//~ }	
	//~ 
	//~ MPI_Barrier(mpi_grid);
	//~ if(rank==out[6].rank)
	//~ {
		//~ cout<<std::scientific<<"ROW 2 TRANSVERSALITY-- MAX "<< maxTT[6]<<" &     Avg "<< outavg[6]/NP<<endl;
		//~ fflush(stdout);
	//~ }	
	//~ 
	//~ MPI_Barrier(mpi_grid);
	//~ if(rank==out[7].rank)
	//~ {
		//~ cout<<std::scientific<<"ROW 3 TRANSVERSALITY -- MAX "<< maxTT[7]<<" &     Avg "<< outavg[7]/NP<<endl;
		//~ fflush(stdout);
	//~ }	
	//~ 
	MPI_Barrier(mpi_grid);
	if(rank==out[8].rank)
	{
		cout<<std::scientific<<"TRACE -- MAX "<< maxTT[8]<<" &    Avg "<< outavg[8]/NP<<endl;
		fflush(stdout);
	}	
	

	MPI_Barrier(mpi_grid);
	//~ fflush(stdout);
}




struct rparams
{
	double W;
	double X;
	double Y;
	double Z;
	double p1;
	double p2;
	double p3;
	double p4;
	double p5;
	double PI;
	double tau;
	double r;
};

int conservation_f (const gsl_vector * x, void *params, gsl_vector * f)
{
	double W = ((struct rparams *) params)->W;
	double X = ((struct rparams *) params)->X;
	double Y = ((struct rparams *) params)->Y;
	double Z = ((struct rparams *) params)->Z;

	double p1 = ((struct rparams *) params)->p1;
	double p2 = ((struct rparams *) params)->p2;
	double p3 = ((struct rparams *) params)->p3;
	double p4 = ((struct rparams *) params)->p4;
	double p5 = ((struct rparams *) params)->p5;

	double PI = ((struct rparams *) params)->PI;

	double tau = ((struct rparams *) params)->tau;
	double r = ((struct rparams *) params)->r;

	const double e = gsl_vector_get (x, 0);
	const double u1 = gsl_vector_get (x, 1);
	const double u2 = gsl_vector_get (x, 2);
	const double u3 = gsl_vector_get (x, 3);

	double P = EOS(e );

	const double y0 = -(u1*X) - u2*Y - u3*Z*pow(tau,2) + (-e + W)*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),0.5);
	const double y1 = -(e*u1) - u2*(p3 + (e + P - PI)*u1*u2) - u3*(p4 + (e + P - PI)*u1*u3)*pow(tau,2) - u1*(P + p1 - PI + (e + P - PI)*pow(u1,2)) + X*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),0.5);
	const double y2 = -(e*u2) - u1*(p3 + (e + P - PI)*u1*u2) - u3*(p5 + (e + P - PI)*u2*u3)*pow(tau,2) - u2*(P + p2 - PI + (e + P - PI)*pow(u2,2)) + Y*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),0.5);
	const double y3 = -(p4*u1) - p5*u2 - P*u3 + u3*(p1 + p2 + PI - W) + Z*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),0.5);

	gsl_vector_set (f, 0, y0);
	gsl_vector_set (f, 1, y1);
	gsl_vector_set (f, 2, y2);
	gsl_vector_set (f, 3, y3);

	return GSL_SUCCESS;
}


int conservation_df  (const gsl_vector * x, void *params, gsl_matrix * J)
{
	double W = ((struct rparams *) params)->W;
	double X = ((struct rparams *) params)->X;
	double Y = ((struct rparams *) params)->Y;
	double Z = ((struct rparams *) params)->Z;

	double p1 = ((struct rparams *) params)->p1;
	double p2 = ((struct rparams *) params)->p2;
	double p3 = ((struct rparams *) params)->p3;
	double p4 = ((struct rparams *) params)->p4;
	double p5 = ((struct rparams *) params)->p5;

	double PI = ((struct rparams *) params)->PI;

	double tau = ((struct rparams *) params)->tau;
	double r = ((struct rparams *) params)->r;

	const double e = gsl_vector_get (x, 0);
	const double u1 = gsl_vector_get (x, 1);
	const double u2 = gsl_vector_get (x, 2);
	const double u3 = gsl_vector_get (x, 3);

	double P = EOS(e);
	double a = DPDE(e);

	const double df00 =  -pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),0.5);
	const double df01 =  -X + u1*(-e + W)*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),-0.5) ;
	const double df02 =  -Y + u2*(-e + W)*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),-0.5) ;
	const double df03 =  pow(tau,2)*(-Z + u3*(-e + W)*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),-0.5));

	const double df10 =  -((1 + a)*u1*(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2)));
	const double df11 =  -p1 - (e + P - PI)*(1 + 3*pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2)) + u1*X*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),-0.5);
	const double df12 =  -p3 - 2*(e + P - PI)*u1*u2 + u2*X*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),-0.5);
	const double df13 =  pow(tau,2)*(-p4 - 2*(e + P - PI)*u1*u3 + u3*X*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),-0.5));

	const double df20 = -((1 + a)*u2*(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2)));
	const double df21 = -p3 - 2*(e + P - PI)*u1*u2 + u1*Y*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),-0.5);
	const double df22 = -p2 - (e + P - PI)*(1 + pow(u1,2) + 3*pow(u2,2) + pow(tau,2)*pow(u3,2)) + u2*Y*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),-0.5);
	const double df23 = pow(tau,2)*(-p5 - 2*(e + P - PI)*u2*u3 + u3*Y*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),-0.5));

	const double df30 = (-a*u3);
	const double df31 = -p4 + u1*Z*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),-0.5);
	const double df32 = -p5 + u2*Z*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),-0.5);
	const double df33 = -P + p1 + p2 + PI - W + u3*Z*pow(tau,2)*pow(1 + pow(u1,2) + pow(u2,2) + pow(tau,2)*pow(u3,2),-0.5);


	gsl_matrix_set (J, 0, 0, df00);
	gsl_matrix_set (J, 0, 1, df01);
	gsl_matrix_set (J, 0, 2, df02);
	gsl_matrix_set (J, 0, 3, df03);

	gsl_matrix_set (J, 1, 0, df10);
	gsl_matrix_set (J, 1, 1, df11);
	gsl_matrix_set (J, 1, 2, df12);
	gsl_matrix_set (J, 1, 3, df13);

	gsl_matrix_set (J, 2, 0, df20);
	gsl_matrix_set (J, 2, 1, df21);
	gsl_matrix_set (J, 2, 2, df22);
	gsl_matrix_set (J, 2, 3, df23);

	gsl_matrix_set (J, 3, 0, df30);
	gsl_matrix_set (J, 3, 1, df31);
	gsl_matrix_set (J, 3, 2, df32);
	gsl_matrix_set (J, 3, 3, df33);

	return GSL_SUCCESS;  
}


int conservation_fdf(const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J)
{
  conservation_f (x, params, f);
  conservation_df (x, params, J);

  return GSL_SUCCESS;
}


void MultiRootSearchForEnVelUsingDerivatives(GRID HydroGrid , double tau)
{

	int i,j,k;
	
	int status;
	
	double r = 0;

	const gsl_multiroot_fdfsolver_type *T;
	gsl_multiroot_fdfsolver *s;
	
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=0;k<ZCMA;k++)
	{
		
		struct rparams p = { HydroGrid[i][j][k].T00,
								HydroGrid[i][j][k].T10,
								HydroGrid[i][j][k].T20,
								HydroGrid[i][j][k].T30,
								HydroGrid[i][j][k].pi[0],
								HydroGrid[i][j][k].pi[1],
								HydroGrid[i][j][k].pi[2],
								HydroGrid[i][j][k].pi[3],
								HydroGrid[i][j][k].pi[4],
								HydroGrid[i][j][k].PI,
								tau,
								HydroGrid[i][j][k].r
								};
  
		
		gsl_multiroot_function_fdf FDF ={   &conservation_f,
											&conservation_df,
											&conservation_fdf,
											NOVAR,
											&p 
										};	
											
		double x_init[NOVAR] = {  HydroGrid[i][j][k].En,
								HydroGrid[i][j][k].u[1],
								HydroGrid[i][j][k].u[2],
								HydroGrid[i][j][k].u[3]  
							};
		
		
		gsl_vector *guess = gsl_vector_alloc (NOVAR);
		
		gsl_vector_set (guess, 0, x_init[0]);
	    gsl_vector_set (guess, 1, x_init[1]);
		gsl_vector_set (guess, 2, x_init[2]);
	    gsl_vector_set (guess, 3, x_init[3]);
	    
	    
		T = gsl_multiroot_fdfsolver_gnewton;
		s = gsl_multiroot_fdfsolver_alloc (T, NOVAR);
		gsl_multiroot_fdfsolver_set (s, &FDF, guess);

	
	
		gsl_vector *result;
		int iter = 0, max_iter = 100;
		do
	    {
			iter++;
			status = gsl_multiroot_fdfsolver_iterate (s);
						
			if(status == GSL_EBADFUNC)
			{
				cout<<"Root finding failure"<<endl;
				cout<<"error --> iteration encountered a singular point where the function or its derivative evaluated to Inf or NaN. "<<endl;
				exit(1);
			}
			
			if(status == GSL_ENOPROG)
			{
				cout<<"Root finding failure"<<endl;
				cout<<"error -->  iteration is not making any progress, preventing the algorithm from continuing  "<< endl;
				exit(1);
			}
			
			result = gsl_multiroot_fdfsolver_root (s);
			status = gsl_multiroot_test_delta (s->dx, s->x, 0, 1e-30);
			
			if(status != GSL_SUCCESS && status != GSL_CONTINUE )
			{	
				cout<<"Root finding failure"<<endl;
				cout<<"error code  --> "<< status<<endl;
				exit(1);
			}
			
	    } while (status == GSL_CONTINUE && iter < max_iter);
	
	
		if(1)//status == GSL_SUCCESS)
		{		
			HydroGrid[i][j][k].En = gsl_vector_get(result, 0);    //find the new energy from root finding algorithm
			HydroGrid[i][j][k].Temp = FT(HydroGrid[i][j][k].En);
			HydroGrid[i][j][k].u[1] = gsl_vector_get(result, 1);    //find the new u1 from root finding algorithm
			HydroGrid[i][j][k].u[2] = gsl_vector_get(result, 2);    //find the new u2 from root finding algorithm
			HydroGrid[i][j][k].u[3] = gsl_vector_get(result, 3);    //find the new u3 from root finding algorithm
			
				
			double u1 = HydroGrid[i][j][k].u[1];
			double u2 = HydroGrid[i][j][k].u[2];
			double u3 = HydroGrid[i][j][k].u[3];
			
			HydroGrid[i][j][k].u[0] =  sqrt(1+u1*u1+u2*u2+u3*u3*tau*tau);   //find the new u3 from root finding algorithm
		}
		//~ else
		//~ {		
			//~ cout<<"Root finding failure"<<endl;
			//~ cout<<"error code  --> UNKNOWN"<< endl;
			//~ cout<<"X-- "<<HydroGrid[i][j][k].X<<" Y-- "<<HydroGrid[i][j][k].Y<<" ETA -- "<<HydroGrid[i][j][k].eta<<endl;
			//~ cout<<"Status "<<status<<" iter " <<iter<<endl;
			//~ exit(1);
		//~ }
		
		gsl_multiroot_fdfsolver_free (s);
		gsl_vector_free (guess);
	}
	
	int quit =0;
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=0;k<ZCMA;k++)
	{
		
		double En = HydroGrid[i][j][k].En ;
		double u0 = HydroGrid[i][j][k].u[0];
		double u1 = HydroGrid[i][j][k].u[1];
		double u2 = HydroGrid[i][j][k].u[2];
		double u3 = HydroGrid[i][j][k].u[3];	

		double r = HydroGrid[i][j][k].r;
		
		double vx = u1/u0;
		double vy = u2/u0;
		double ve = u3/u0;
		
		if( sqrt(vx*vx+vy*vy+ve*ve) >= 1 || En <=0 )
		{
			 cout<<" violations of SOL / -ve energy in multiroot.h -  2nd place"<<endl;
			 cout<<"HELL BROKE LOOSE X "<<HydroGrid[i][j][k].X <<"  Y " << HydroGrid[i][j][k].Y<< "  E "<<HydroGrid[i][j][k].eta << "  "<<sqrt(vx*vx+vy*vy+ve*ve)<<endl;
			 cout<<"HELL BROKE LOOSE vX "<<vx <<"  vY " << vy<< "  vE "<<ve << "  "<<sqrt(vx*vx+vy*vy+ve*ve)<<endl;
			 cout<<"HELL BROKE LOOSE Energy  "<<En<<endl;
			 quit = 1;
			 exit(quit);
		 }
		 		
		HydroGrid[i][j][k].Vx = vx;
		HydroGrid[i][j][k].Vy = vy;
		HydroGrid[i][j][k].Ve = ve;
		HydroGrid[i][j][k].P = EOS(En );		
	}
	
	
	if(quit)
		exit(quit);
}
