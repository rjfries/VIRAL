void StoreInHeap();

void CheckRoot(GRID HydroGrid);


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
#define GRIDZMAX ((GRIDZPOINTS/2 )*ZS)


//Define's which inherit values based on global dependence aka (xcord,ycord,zcord)
#define XSTART   (-GRIDXMAX + xcord*(XL))
#define XEND	(XSTART+(XCMA-1)*XS)
#define YSTART   (-GRIDYMAX + ycord*(YL))
#define YEND	(YSTART+(YCMA-1)*YS)
#define ZSTART   (-GRIDZMAX + zcord*(ZL)+ ETASTART)
#define ZEND	(ZSTART+(ZCMA-1)*ZS)


#define iSTART   (xcord*(XCMA))
#define iEND	(iSTART+(XCMA-1))
#define jSTART  ( ycord*(YCMA))
#define jEND	(jSTART+(YCMA-1))
#define kSTART  ( zcord*(ZCMA))
#define kEND	(kSTART+(ZCMA-1))


#define iSTARTWB   (iSTART)
#define iENDWB	(iEND+2*(BORDER))
#define jSTARTWB  (jSTART)
#define jENDWB	(jEND+2*(BORDER))
#define kSTARTWB  (kSTART)
#define kENDWB	(kEND+2*(BORDER))


inline int MYGLOBALi(int i)   {return(iSTARTWB+i);}
inline int MYGLOBALj(int j)   {return(jSTARTWB+j);}
inline int MYGLOBALk(int k)  {return(kSTARTWB+k);}

inline int MYGLOBALiWB(int i)   {return(iSTART+(i-BORDER)); }
inline int MYGLOBALjWB(int j)   {return(jSTART+(j-BORDER)); }
inline int MYGLOBALkWB(int k)   {return(kSTART+(k-BORDER)); }



inline double XCORD(int i)  {return(XSTART + (i)*XS); }
inline double YCORD(int j)  {return(YSTART + (j)*YS); }
inline double  ZCORD(int k)  {return(ZSTART + (k)*ZS); }


inline double  XCORDWB(int i)  {return(XSTART + (i-BORDER)*XS); }
inline double YCORDWB(int j)  {return(YSTART + (j-BORDER)*YS); }
inline double  ZCORDWB(int k)  {return(ZSTART + (k-BORDER)*ZS); }


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



inline void UpdateEOS(GRID HydroGrid)
{
	int i,j,k;

	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		HydroGrid[i][j][k].P = EOS(HydroGrid[i][j][k].En ,  HydroGrid[i][j][k].r);
	}
}

void initGubser(GRID HydroGrid, double tau, double ts)
{
	int i,j,k;
	
	double q =1;
	double eps0 = 1;
	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		double X =  HydroGrid[i][j][k].X;
		double Y =  HydroGrid[i][j][k].Y;
		double eta =  HydroGrid[i][j][k].eta;
		double r =  HydroGrid[i][j][k].r;
		
		
		double g = 2*q*q*( tau*tau + r*r ) + q*q*q*q*( tau*tau - r*r )*( tau*tau - r*r );		
		double kappa = atanh(  (2*q*q*tau*r)  / ( 1 + q*q*tau*tau + q*q*r*r)  );		
		double utau = cosh(kappa);
		
		HydroGrid[i][j][k].En = eps0 * pow( 2*q , 8.0/3.0)   /  ( pow( tau*(1+g) , 4.0/3.0) ) ;
		
		if( fabs(r) > 1e-6)
		{
			HydroGrid[i][j][k].Vx = ( X*tanh(kappa) )/ (r);
			HydroGrid[i][j][k].Vy = ( Y*tanh(kappa) )/ (r);
		}
		else
		{
			HydroGrid[i][j][k].Vx = 0;
			HydroGrid[i][j][k].Vy = 0;	
		}
		
		HydroGrid[i][j][k].Ve = 0;	
		
		
		HydroGrid[i][j][k].P = EOS(HydroGrid[i][j][k].En ,  r);
				
				
		double En = HydroGrid[i][j][k].En;
		double Vx = HydroGrid[i][j][k].Vx;
		double Vy = HydroGrid[i][j][k].Vy;
		double Ve = HydroGrid[i][j][k].Ve;
 		double P = HydroGrid[i][j][k].P;

		double G = 1.0/sqrt(1.0-Vx*Vx-Vy*Vy-tau*tau*Ve*Ve);

		HydroGrid[i][j][k].T00 = -(P) + (En+P)*G*G;
		HydroGrid[i][j][k].T10 =  (En+P)*Vx*G*G; 
		HydroGrid[i][j][k].T20 =  (En+P)*Vy*G*G; 
		HydroGrid[i][j][k].T30 =  (En+P)*Ve*G*G; 
			
	}

	
}


void initvar(GRID HydroGrid, double tau, double ts)
{
	int i,j,k;
	
//define your initial conditions
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		double step = 16;
		double X =  HydroGrid[i][j][k].X;
		double Y =  HydroGrid[i][j][k].Y;
		double eta =  HydroGrid[i][j][k].eta;
		double r =  HydroGrid[i][j][k].r;
		
		if(X+Y<0)
			HydroGrid[i][j][k].En = 16;
		else
			HydroGrid[i][j][k].En = 1;
		
		HydroGrid[i][j][k].En = step*exp(-Y*Y);
		//HydroGrid[i][j][k].En = step;
		
		HydroGrid[i][j][k].Vx= 0;
		HydroGrid[i][j][k].Vy= 0;
		HydroGrid[i][j][k].Ve= 0;	
				
		HydroGrid[i][j][k].P = EOS(HydroGrid[i][j][k].En ,  r);
				
				
		double En = HydroGrid[i][j][k].En;
		double Vx = HydroGrid[i][j][k].Vx;
		double Vy = HydroGrid[i][j][k].Vy;
		double Ve = HydroGrid[i][j][k].Ve;
 		double P = HydroGrid[i][j][k].P;

		double G = 1.0/sqrt(1.0-Vx*Vx-Vy*Vy-Ve*Ve);

		HydroGrid[i][j][k].T00 = -(P) + (En+P)*G*G;
		HydroGrid[i][j][k].T10 =  (En+P)*Vx*G*G; 
		HydroGrid[i][j][k].T20 =  (En+P)*Vy*G*G; 
		HydroGrid[i][j][k].T30 =  (En+P)*Ve*G*G; 
	}
}




void initvarStatic(GRID HydroGrid, double tau, double ts)
{
	int i,j,k;
	
//define your initial conditions
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=0;k<ZCM;k++)
	{
		double step = 16;
		
		double X =  HydroGrid[i][j][k].X;
		double Y =  HydroGrid[i][j][k].Y;
		double eta =  HydroGrid[i][j][k].eta;
		double r =  HydroGrid[i][j][k].r;
		
		//~ if(X+Y<0)
			//~ HydroGrid[i][j][k].En = 16;
		//~ else
			//~ HydroGrid[i][j][k].En = 1;
		
		HydroGrid[i][j][k].En = step;
		
		HydroGrid[i][j][k].Vx= 0;
		HydroGrid[i][j][k].Vy= 0;
		HydroGrid[i][j][k].Ve= -tanh(eta) ;	
				
		HydroGrid[i][j][k].P = EOS(HydroGrid[i][j][k].En ,  r);
				
				
		double En = HydroGrid[i][j][k].En;
		double Vx = HydroGrid[i][j][k].Vx;
		double Vy = HydroGrid[i][j][k].Vy;
		double Ve = HydroGrid[i][j][k].Ve;
 		double P = HydroGrid[i][j][k].P;

		double G = 1.0/sqrt(1.0-Vx*Vx-Vy*Vy-tau*tau*Ve*Ve);

		HydroGrid[i][j][k].T00 = -(P) + (En+P)*G*G;
		HydroGrid[i][j][k].T10 =  (En+P)*Vx*G*G; 
		HydroGrid[i][j][k].T20 =  (En+P)*Vy*G*G; 
		HydroGrid[i][j][k].T30 =  (En+P)*Ve*G*G; 
	}
}




void initvarAZ(GRID HydroGrid, double tau, double ts)
{
	int i,j,k;
	
	
	ifstream ff,fx,fy;
	ff.open("init/EPS001.dat", std::ifstream::in);	
	fx.open("init/VX001.dat", std::ifstream::in);	
	fy.open("init/VY001.dat", std::ifstream::in);
	
	
//define your initial conditions
	for(k=0;k<ZCM;k++)
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	{
		double step = 16;
		double x =  HydroGrid[i][j][k].X;
		double y =  HydroGrid[i][j][k].Y;
		double eta =  HydroGrid[i][j][k].eta;
		double r =  HydroGrid[i][j][k].r;

		//HydroGrid[i][j][k].En= 16*exp(-(x*x+y*y)/4);
		HydroGrid[i][j][k].Vx= 0;
		HydroGrid[i][j][k].Vy= 0;
		HydroGrid[i][j][k].Ve= 0;		
		
	}


	double min = 100;
	for(j=jl;j<jr;j++)
	for(i=il;i<ir;i++)
	{	
		ff>>HydroGrid[i][j][0].En;
		fx>>HydroGrid[i][j][0].Vx;
		fy>>HydroGrid[i][j][0].Vy;
		
		if(HydroGrid[i][j][0].En <min  && fabs(HydroGrid[i][j][0].En ) >1e-5  )
			min = HydroGrid[i][j][0].En;
	}

	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
		if(HydroGrid[i][j][0].En < 0.01)
			HydroGrid[i][j][0].En = min/10;


	cout<<min<<endl;
	
	for(i=0;i<XCM;i++)
	for(j=0;j<YCM;j++)
	for(k=1;k<ZCM;k++)
	{
		HydroGrid[i][j][k].En = HydroGrid[i][j][0].En;
		HydroGrid[i][j][k].Vx = HydroGrid[i][j][0].Vx;
		HydroGrid[i][j][k].Vy = HydroGrid[i][j][0].Vy;
	}		
		
		
	ff.close();
	fx.close();
	fy.close();
	
	UpdateEOS(HydroGrid);	
	//FillTNu0( HydroGrid);   //TODO: Fix this
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
		
		
	il=jl=kl=BORDER;
	ir=il+XCMA;
	jr=jl+YCMA;
	kr=kl+ZCMA;

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
	
	k0 = int(-(ZSTART-ETASTART)/ZS + BORDER + OFF); 
	
	#if !defined(BJORKEN)&& !defined(GUBSER)&& !defined(STAT)
		initvar(HydroGrid,tau, ts);
	#endif
	
 

	#ifdef BJORKEN		
		initBjorken(HydroGrid,tau, ts);
	#endif

	#ifdef STAT		
		initvarStatic(HydroGrid,tau, ts);
	#endif
	
	#ifdef GUBSER
		initGubser(HydroGrid,tau,ts);
	#endif

}


