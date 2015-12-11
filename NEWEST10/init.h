void StoreInHeap();
void NSInit(GRID HydroGrid, double tau, double ts);
void ZeroInit(GRID HydroGrid, double tau, double ts);
void CheckRoot(GRID HydroGrid , double tau);
double DebugMSG(GRID HydroGrid);

int k0;

//PARALLEL STUFFS
int rank, size, root=0;

//mpi_init
int xcord,ycord,zcord;
MPI_Comm mpi_grid;
int il,jl,kl,ir,jr,kr;
int debug = 0;

#define GRIDXPOINTS  (NPX*XCMA)
#define GRIDYPOINTS  (NPY*YCMA)
#define GRIDZPOINTS  (NPZ*ZCMA)

#define GRIDXMAX ((GRIDXPOINTS/2 )*XS)
#define GRIDYMAX ((GRIDYPOINTS/2 )*YS)
#define GRIDZMAX ((GRIDZPOINTS/2 )*ZS + ETASTART)



#define ETASTART 0


//Define's which inherit values based on global dependence aka (xcord,ycord,zcord)
#define XSTART   (-GRIDXMAX + xcord*(XL))
#define XEND	(XSTART+(XCMA-1)*XS)
#define YSTART   (-GRIDYMAX + ycord*(YL))
#define YEND	(YSTART+(YCMA-1)*YS)

#define ZSTART   ( -GRIDZMAX + zcord*(ZL)+ ETASTART )
#define ZEND	( ZSTART+(ZCMA-1)*ZS )


#define iSTART   (xcord*(XCMA))
#define iEND	(iSTART+(XCMA-1))
#define jSTART  ( ycord*(YCMA))
#define jEND	(jSTART+(YCMA-1))

#if !defined  LBI
#define kSTART  ( zcord*(ZCMA))
#define kEND	(kSTART+(ZCMA-1))
#endif


#if defined  LBI
#define kSTART  (0)
#define kEND	(1)
#endif


#define iSTARTWB   (iSTART)
#define iENDWB	(iEND+2*(BORDER))
#define jSTARTWB  (jSTART)
#define jENDWB	(jEND+2*(BORDER))

//~ #define kSTARTWB  (kSTART)
//~ #define kENDWB	(kEND+2*(BORDER))


inline int MYGLOBALi(int i)   {return(iSTARTWB+i);}
inline int MYGLOBALj(int j)   {return(jSTARTWB+j);}
//~ inline int MYGLOBALk(int k)   {return(kSTARTWB+k);} 

inline int MYGLOBALiWB(int i)   {return(iSTART+(i-BORDER)); }
inline int MYGLOBALjWB(int j)   {return(jSTART+(j-BORDER)); }
//~ inline int MYGLOBALkWB(int k)   {return(kSTART+(k-BORDER)); }

 

inline double  XCORD(int i)  {return(XSTART + (i)*XS); }
inline double  YCORD(int j)  {return(YSTART + (j)*YS); }
 
inline double  ZCORD(int k)  {return(ZSTART + (k)*ZS); } 

inline double  XCORDWB(int i)  {return(XSTART + (i-BORDER)*XS); }
inline double  YCORDWB(int j)  {return(YSTART + (j-BORDER)*YS); }

#if !defined  LBI
inline double ZCORDWB(int k)  {return(ZSTART + (k-BORDER)*ZS); }
#else
inline double ZCORDWB(int k) {return ZSTART;}
#endif

#define BOUNDARY -1

inline int RANK(int i, int j,int k) 
{
	int coords[3]={i,j,k};
	int rank1;

	MPI_Cart_rank(mpi_grid,coords,&rank1);
  	return(rank1);
}

//#define RANK(i,j,k)  i+j+k

#define MYXLEFT ((!xcord)?BOUNDARY:RANK(xcord-1,ycord,zcord) )
#define MYXRIGHT ((!(xcord-(NPX-1))?BOUNDARY:RANK(xcord+1,ycord,zcord) ))

#define MYYLEFT ((!ycord)?BOUNDARY:RANK(xcord,ycord-1,zcord) )
#define MYYRIGHT ((!(ycord-(NPY-1))?BOUNDARY:RANK(xcord,ycord+1,zcord) ))



#define MYXLYL (     (MYXLEFT != BOUNDARY && MYYLEFT != BOUNDARY )?  RANK(xcord-1,ycord-1,zcord) : BOUNDARY  )
#define MYXLYR (     (MYXLEFT != BOUNDARY && MYYRIGHT != BOUNDARY )?  RANK(xcord-1,ycord+1,zcord) : BOUNDARY  ) 


#define MYXRYL (     (MYXRIGHT != BOUNDARY && MYYLEFT != BOUNDARY )?  RANK(xcord+1,ycord-1,zcord) : BOUNDARY  )
#define MYXRYR (     (MYXRIGHT != BOUNDARY && MYYRIGHT != BOUNDARY )?  RANK(xcord+1,ycord+1,zcord) : BOUNDARY  ) 



#define MYZLEFT ((!zcord)?BOUNDARY:RANK(xcord,ycord,zcord-1) )
#define MYZRIGHT ((!(zcord-(NPZ-1))?BOUNDARY:RANK(xcord,ycord,zcord+1) ))




void initvar(GRID HydroGrid);
void ginit(GRID HydroGrid, double tau);


void mpi_init()
{
	int coord[3];
	int dims[3] = {NPX,NPY,NPZ};
	int periods[3] = {0};
	int reorder = 0; 

	if(rank==0)
	if(NPX*NPY*NPZ!=size)
	{
		cout<<"\n\n*****Run this simulation with " <<NPX*NPY*NPZ<<" processors, not "<<size<<" processors.*****\n\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	MPI_Cart_create(MPI_COMM_WORLD, DIM, dims , periods, reorder, &mpi_grid);
	MPI_Cart_coords(mpi_grid, rank, DIM, coord);

	
	xcord = coord[0];
	ycord = coord[1];
	zcord = coord[2];
	root=0;

}

 


void initvar(GRID HydroGrid, double tau, double ts)
{
	int i,j,k,l;
	
//define your initial conditions
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		double step = 16;
		double x =  HydroGrid[i][j][k].X;
		double y =  HydroGrid[i][j][k].Y;
		double eta =  HydroGrid[i][j][k].eta;
		double r =  HydroGrid[i][j][k].r;
		
		
		HydroGrid[i][j][k].En = 16*exp(-r*r);
	
		
		HydroGrid[i][j][k].Vx= 0;
		HydroGrid[i][j][k].Vy= 0;
		HydroGrid[i][j][k].Ve= 0;	
				
		HydroGrid[i][j][k].P = EOS(HydroGrid[i][j][k].En ,  r);
		
		
		
		HydroGrid[i][j][k].u[0]= 1.0/(sqrt(1 - HydroGrid[i][j][k].Vx*HydroGrid[i][j][k].Vx
											 - HydroGrid[i][j][k].Vy*HydroGrid[i][j][k].Vy
											 - tau*tau*HydroGrid[i][j][k].Ve*HydroGrid[i][j][k].Ve
											 )
									 );
		HydroGrid[i][j][k].u[1] =  HydroGrid[i][j][k].u[0]*HydroGrid[i][j][k].Vx;
		HydroGrid[i][j][k].u[2] =  HydroGrid[i][j][k].u[0]*HydroGrid[i][j][k].Vy;
		HydroGrid[i][j][k].u[3] =  HydroGrid[i][j][k].u[0]*HydroGrid[i][j][k].Ve;	
		
		
		HydroGrid[i][j][k].prevu[0] = HydroGrid[i][j][k].u[0]; //value 1 "ts" earlier
		HydroGrid[i][j][k].prevu[1] = HydroGrid[i][j][k].u[1]; //value 1 "ts" earlier
		HydroGrid[i][j][k].prevu[2] = HydroGrid[i][j][k].u[2]; //value 1 "ts" earlier
		HydroGrid[i][j][k].prevu[3] = HydroGrid[i][j][k].u[3]; //value 1 "ts" earlier
	}
	
	//CalcNS(HydroGrid,tau,ts);
	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		HydroGrid[i][j][k].PI =   HydroGrid[i][j][k].nsPI;
		
		
		for(l=0;l<Npi;l++)
			HydroGrid[i][j][k].pi[l] =   HydroGrid[i][j][k].nspi[l];
			
 		DECLePPIa;
		DECLp10u4;


        HydroGrid[i][j][k].T00 = -P + PI + (e + P - PI)*pow(u0,2) + p1;
        HydroGrid[i][j][k].T10 = (e + P - PI)*u0*u1 + p2;
        HydroGrid[i][j][k].T20 = (e + P - PI)*u0*u2 + p3;
        HydroGrid[i][j][k].T30 = (e + P - PI)*u0*u3 + p4;
	}
}



void initBjorken(GRID HydroGrid, double tau, double ts)
{
	int i,j,k,l;
	
//define your initial conditions
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		double x =  HydroGrid[i][j][k].X;
		double y =  HydroGrid[i][j][k].Y;
		double eta =  HydroGrid[i][j][k].eta;
		double r =  HydroGrid[i][j][k].r;
		
		HydroGrid[i][j][k].En = 30/GEVFM;
		
		HydroGrid[i][j][k].Vx= 0;
		HydroGrid[i][j][k].Vy= 0;
		HydroGrid[i][j][k].Ve= 0;	
				
		HydroGrid[i][j][k].P = EOS(HydroGrid[i][j][k].En ,  r);
		
		
		
		HydroGrid[i][j][k].u[0]= 1.0/(sqrt(1 - HydroGrid[i][j][k].Vx*HydroGrid[i][j][k].Vx
											 - HydroGrid[i][j][k].Vy*HydroGrid[i][j][k].Vy
											 - tau*tau*HydroGrid[i][j][k].Ve*HydroGrid[i][j][k].Ve
											 )
									 );
		HydroGrid[i][j][k].u[1] =  HydroGrid[i][j][k].u[0]*HydroGrid[i][j][k].Vx;
		HydroGrid[i][j][k].u[2] =  HydroGrid[i][j][k].u[0]*HydroGrid[i][j][k].Vy;
		HydroGrid[i][j][k].u[3] =  HydroGrid[i][j][k].u[0]*HydroGrid[i][j][k].Ve;	
		
		DECLu4;
		HydroGrid[i][j][k].prevu[0] = u0; //value 1 "ts" earlier
		HydroGrid[i][j][k].prevu[1] = u1; //value 1 "ts" earlier
		HydroGrid[i][j][k].prevu[2] = u2; //value 1 "ts" earlier
		HydroGrid[i][j][k].prevu[3] = u3; //value 1 "ts" earlier
	}
	
	//~ CalcNS(HydroGrid,tau,ts);
	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		//~ HydroGrid[i][j][k].PI =   HydroGrid[i][j][k].nsPI;
		//~ for(l=0;l<10;l++)
			//~ HydroGrid[i][j][k].pi[l] =   HydroGrid[i][j][k].nspi[l];
			
 		DECLePPIa;
 		DECLp10u4;		 

        HydroGrid[i][j][k].T00 = -P + PI + (e + P - PI)*pow(u0,2) + p1;
        HydroGrid[i][j][k].T10 = (e + P - PI)*u0*u1 + p2;
        HydroGrid[i][j][k].T20 = (e + P - PI)*u0*u2 + p3;
        HydroGrid[i][j][k].T30 = (e + P - PI)*u0*u3 + p4;
	}
}





void initGubser(GRID HydroGrid, double tau, double ts)
{
	int i,j,k,l;
	MPI_Offset offset;
	MPI_Status ierr;
	MPI_File fh;	
	
	MPI_File_open(mpi_grid,  "init/tempbin-15.dat" ,  MPI_MODE_RDONLY  , MPI_INFO_NULL, &fh);	
	
	int nvar = 8;
	int chunk = nvar*XCM;
	double *buf = new double [nvar*chunk];
	
	int f=1;
	
	for(j=0;j<YCM;j=j+f)
	{ 	
		offset = nvar*sizeof(double)*(  MYGLOBALj(j)*(GRIDXPOINTS+2*BORDER) + MYGLOBALi(0) );
	
		MPI_File_seek(fh,offset,MPI_SEEK_SET);	//takes offset in bytes
	
		int error = MPI_File_read(fh, (void*)buf, chunk, MPI_DOUBLE, &ierr);
		
		for(i=0 ; i < XCM ; i=i+f)
		{
			double temp = buf[nvar*i+0];
			HydroGrid[i][j][k0].En = FEnFromTemp(temp);			
			HydroGrid[i][j][k0].pi[0] = buf[nvar*i+1]; //pitautau
			HydroGrid[i][j][k0].pi[1] = buf[nvar*i+2]; //pitauX
			HydroGrid[i][j][k0].pi[2] = buf[nvar*i+3]; //pitauY
			HydroGrid[i][j][k0].pi[4] = buf[nvar*i+4]; //piXX
			HydroGrid[i][j][k0].pi[5] = buf[nvar*i+5]; //piXY
			HydroGrid[i][j][k0].pi[7] = buf[nvar*i+6]; //piYY
			HydroGrid[i][j][k0].pi[9] = buf[nvar*i+7]; //pietaeta
		}		
		
		if(error!=MPI_SUCCESS)
		{
			cout<<"Could not read gubser temperatures"<<endl;
			MPI_Finalize();
			exit(0);
		}
	}
	MPI_File_close(&fh);
	
	
	double q =1;
	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		//~ DECLcoord;
		
		double X =  HydroGrid[i][j][k].X;
		double Y =  HydroGrid[i][j][k].Y;
		double eta =  HydroGrid[i][j][k].eta;
		double r =  HydroGrid[i][j][k].r;
		
		 	
		double kappa = atanh(  (2*q*q*tau*r)  / ( 1 + q*q*tau*tau + q*q*r*r)  );		
			
		
		if( fabs(r) > 1e-7 )
		{
			HydroGrid[i][j][k].Vx   = (( X*tanh(kappa) )/(r+1e-15));
			HydroGrid[i][j][k].Vy   = (( Y*tanh(kappa) )/(r+1e-15));
			HydroGrid[i][j][k].u[1] = (( X*sinh(kappa) )/(r+1e-15));
			HydroGrid[i][j][k].u[2] = (( Y*sinh(kappa) )/(r+1e-15));	
		}
		else
		{
			HydroGrid[i][j][k].Vx   = 0;
			HydroGrid[i][j][k].Vy   = 0;	
			HydroGrid[i][j][k].u[1] = 0;
			HydroGrid[i][j][k].u[2] = 0;
		}		
		
		
		HydroGrid[i][j][k].P = EOS(HydroGrid[i][j][k].En ,r );		
		HydroGrid[i][j][k].u[0]  =  cosh(kappa);
		
		HydroGrid[i][j][k].prevu[1] = HydroGrid[i][j][k].u[1];
		HydroGrid[i][j][k].prevu[2] = HydroGrid[i][j][k].u[2];
		HydroGrid[i][j][k].prevu[0]  =  HydroGrid[i][j][k].u[0];
		
	}
	  
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		DECLcoord;
		 	
		double kappa = atanh(  (2*q*q*(tau-ts)*r)  / ( 1 + q*q*(tau-ts)*(tau-ts) + q*q*r*r)  );						
		
		if( fabs(r) > 1e-7 )
		{ 	
			HydroGrid[i][j][k].prevu[1] = (( X*sinh(kappa) )/(r+1e-15));
			HydroGrid[i][j][k].prevu[2] = (( Y*sinh(kappa) )/(r+1e-15));		
		}
		else
		{
			HydroGrid[i][j][k].prevu[1] = (0);
			HydroGrid[i][j][k].prevu[2] = (0);			
			
		}
		HydroGrid[i][j][k].prevu[0]  =  cosh(kappa);
	}

	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
 		DECLp10u4;
 		DECLePPIa;
		HydroGrid[i][j][k].T00 = -P + PI + (e + P - PI)*pow(u0,2) + p1;
        HydroGrid[i][j][k].T10 = (e + P - PI)*u0*u1 + p2;
        HydroGrid[i][j][k].T20 = (e + P - PI)*u0*u2 + p3;
        HydroGrid[i][j][k].T30 = (e + P - PI)*u0*u3 + p4;
	}
}





void init(double tau, double ts)
{

	int i,j,k;
		
	mpi_init();
	StoreInHeap();
	AllocatePack();
	
	for(  i = 0 ; i< XCM; i++)
	for(  j = 0 ; j< YCM; j++)
	for(  k = 0 ; k< ZCM; k++)
	{
		HydroGrid[i][j][k].X = XCORDWB(i);
		HydroGrid[i][j][k].Y = YCORDWB(j);
		HydroGrid[i][j][k].eta = ZCORDWB(k);
		
		HydroGrid[i][j][k].r = sqrt(HydroGrid[i][j][k].X*HydroGrid[i][j][k].X + HydroGrid[i][j][k].Y*HydroGrid[i][j][k].Y);
		

	}
		
		
	il=jl=BORDER;
	ir=il+XCMA;
	jr=jl+YCMA;
	

#if !defined LBI
	kl=BORDER
	kr=kl+ZCMA;
#else 
	kl=0;
	kr=1;
#endif

	for(i = 0; i< XCM; i++)
	for(j = 0; j< YCM; j++)
	for(k = 0; k< ZCM; k++)
	{	
		if( (i >= il && i<ir ) &&
			(j >= jl && j<jr ) &&
			(k >= kl && k<kr ) 
			)
			HydroGrid[i][j][k].RELEVANT=true;			
		else
		if( (i < il || i>=ir ) &&
			(j >= jl && j<jr ) &&
			(k >= kl && k<kr ) 
			)
			HydroGrid[i][j][k].RELEVANT=true;
		else
		if( (j < jl || j >= jr ) &&
			(i >= il && i < ir ) &&
			(k >= kl && k < kr ) 
			)
			HydroGrid[i][j][k].RELEVANT=true;
		else
		if( (k < kl || k >= kr ) &&
			(i >= il && i < ir ) &&
			(j >= jl && j < jr )  
			)
			HydroGrid[i][j][k].RELEVANT=true;
		else			
			HydroGrid[i][j][k].RELEVANT=false;
	}


	fflush(stdout);


	MPI_Barrier(mpi_grid);
	if(!rank)
	{
		cout<<"From (-X,-Y,-Z)      == ( "<<-GRIDXMAX <<" , "<<-GRIDYMAX << " , " <<-GRIDZMAX+ETASTART<<" )"<<endl;
		cout<<"To   ( X, Y, Z)      == ( "<<(-GRIDXMAX+(GRIDXPOINTS-1)*XS) <<" , "<<(-GRIDYMAX+(GRIDYPOINTS-1)*YS) << " , " <<(-GRIDZMAX+(GRIDZPOINTS-1)*ZS)+ETASTART<<" )"<<endl;
		cout<<"Number of Processors == ( "<<NPX<<" , "<<NPY << " , " <<NPZ<<" ) == "<<NP<<endl;
		cout<<"Points in total      == ( "<<GRIDXPOINTS<<" , "<<GRIDYPOINTS << " , " <<GRIDZPOINTS<<" )"<<endl;
		cout<<"Points per process   == ( "<<XCMA<<" , "<<YCMA << " , " <<ZCMA<<" )"<<endl;		
		cout<<"XCM,YCM,ZCM          == ( "<<XCM<<" , "<<YCM << " , " <<ZCM<<" )"<<endl;
		cout<<"(il,jl,kl).(ir,jr,kr)== ( "<<il<<" , "<<jl << " , " <<kl<<" )"<<".( "<<ir<<" , "<<jr << " , " <<kr<<" )"<<endl<<endl<<endl<<endl;
	}
	
	
	//~ sleep(rank );
	//~ cout<<"Rank --> "<<rank<< "   "<<xcord<<"   "<<ycord<<"   "<<MYXLYL<< "   "<<MYXLYR<< "   "<<MYXRYL<< "   "<<MYXRYR<<endl;
	//~ cout<<"Rank --> "<<rank<< "   "<<xcord<<"   "<<ycord<<"   "<<MYXLEFT<< "   "<<MYXRIGHT<< "   "<<MYYLEFT<< "   "<<MYYRIGHT<<endl;
	//~ cout<<"Rank --> "<<rank<< "   "<<XSTART<<"   "<<XCORDWB(il)<<"   "<<XEND<<"   "<<XCORDWB(ir-1)<<"  "<<YSTART<< "   "<<YCORDWB(jl)<< "   "<<YEND<< "   "<<YCORDWB(jr-1)<< endl;
	//~ fflush(stdout);
	//~ MPI_Barrier(mpi_grid);
	//~ exit(1);
	
#if !defined LBI
	k0 = int(-(ZSTART-ETASTART)/ZS + BORDER + OFF); 
#else	
	k0 = int(-(ZSTART-ETASTART)/ZS ); 
#endif

	
	if(!rank)
		cout<<"k0 is  "<<k0<<endl;
	
#if !defined(GINIT) && !defined(BJORKEN)&& !defined(GUBSER)
		initvar(HydroGrid,tau, ts);
#endif
	
#ifdef GINIT
		ginit(HydroGrid,tau-ts);
		ginit(HydroGrid,tau);		
#endif

#ifdef BJORKEN		
		initBjorken(HydroGrid,tau, ts);
#endif
	
#ifdef GUBSER
		initGubser(HydroGrid,tau,ts);
#endif
	
	
#ifdef NSINIT
		NSInit(HydroGrid,tau,ts);
#endif
#ifdef ZEROINIT
		ZeroInit(HydroGrid,tau,ts);
#endif

	DebugMSG(HydroGrid);
	CheckRoot( HydroGrid, tau);	
}
