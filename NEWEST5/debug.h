
int CFLRank;
double DebugMSG(GRID HydroGrid)
{
	int i,j,k,l,m;
	double step = XS; 
	double rts = ts;
	
	
	double LV,MV;
	double XXV, YYV, ZZV;
	double XXE, YYE, ZZE;
	double XXEM, YYEM, ZZEM;
	double vmax=0,emax=0,emin=0;
	double vx,vy,vz;	
	double tvx,tvy,tvz;
	

	struct {
	double val;
	int   rank;
	} in, out;
	
	
	
		
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
	{

		tvx = HydroGrid[i][j][k].Vx;
		tvy = HydroGrid[i][j][k].Vy;
		tvz = HydroGrid[i][j][k].Ve;

		double b = sqrt(tvx*tvx+tvy*tvy+ tvz*tvz);

		double energy = HydroGrid[i][j][k].En;

		if(energy>emax)
		{
			XXE = HydroGrid[i][j][k].X;
			YYE = HydroGrid[i][j][k].Y;
			ZZE = HydroGrid[i][j][k].eta;
			emax=energy;
		}

		if(energy<emin)
		{
			XXEM = HydroGrid[i][j][k].X;
			YYEM = HydroGrid[i][j][k].Y;
			ZZEM = HydroGrid[i][j][k].eta;
			emin=energy;
		}	
		if( b >= vmax)
		{
			vmax = b;
			vx = HydroGrid[i][j][k].Vx;
			vy = HydroGrid[i][j][k].Vy;
			vz = HydroGrid[i][j][k].Ve;
			XXV = HydroGrid[i][j][k].X;
			YYV = HydroGrid[i][j][k].Y;
			ZZV = HydroGrid[i][j][k].eta;

		}
	}

	
//max abs vel part
	in.val = vmax;
	in.rank = rank; 

	MPI_Allreduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MAXLOC,mpi_grid);
	BMax =  out.val;
	
	if(rank==out.rank)
	{
		cout<<std::fixed<<"Max vel ( "<<vx<<" , "<<vy<<" , "<<vz<<" ) & mag "<< out.val<<" @"<<" ( "<<XXV<<" , "<<YYV<<" , "<<ZZV<<" )in rank -> "<<rank<<endl;
		fflush(stdout);
	}	
	MPI_Barrier(mpi_grid);


//max Energy part
	in.val = emax;
	in.rank = rank; 

	MPI_Allreduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MAXLOC,mpi_grid);
	
	if(rank==out.rank)
	{
		cout<<std::fixed<<"Max Energy ( "<<emax<<" ) @"<<" ( "<<XXE<<" , "<<YYE<<" , "<<ZZE<<" )in rank -> "<<rank<<endl;
		fflush(stdout);
	}
	MPI_Barrier(mpi_grid);


//min Energy part
	in.val = emin;
	in.rank = rank; 

	MPI_Allreduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MINLOC,mpi_grid);
	
	if(rank==out.rank)
	{
		cout<<std::fixed<<"Min Energy ( "<<emin<<" ) @"<<" ( "<<XXEM<<" , "<<YYEM<<" , "<<ZZEM<<" )in rank -> "<<rank<<endl;
		fflush(stdout);
	}
	
	MPI_Barrier(mpi_grid);




//Source term debug info
	double s, minres[SVAR], smin[SVAR], maxres[SVAR], smax[SVAR];
	
	for(l=0;l<SVAR;l++)
		smax[l]=smin[l]=0; 
	
	for(l=0;l<SVAR;l++)
	for(i=il;i<ir;i++)
	for(j=jl;j<jr;j++)
	for(k=kl;k<kr;k++) 
	{
		s = HydroGrid[i][j][k].Source[l];

		if(s > smax[l])
			smax[l] = s;

		if(s < smin[l])
			smin[l] = s;
	}


	for(l=0;l<SVAR;l++)
	{

		double in1 = smax[l],in2 = smin[l];
		double out1,out2;

		MPI_Reduce(&in1,&out1,1,MPI_DOUBLE,MPI_MAX,root,mpi_grid);
		MPI_Reduce(&in2,&out2,1,MPI_DOUBLE,MPI_MIN,root,mpi_grid);
		
		if(rank==root)
		{
			maxres[l]=out1;
			minres[l]=out2;
		}
	}



	if(rank==root)
	{
		cout<<"var \t Max Source term \t Min Source term"<<endl;

		for(l=0;l<SVAR;l++)
		{
			cout<<std::scientific<<l<<" \t "<<maxres[l]<<" \t "<<minres[l]<<endl;
		}
	}

	
MPI_Barrier(mpi_grid);
return (rts);



}

