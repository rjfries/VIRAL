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

 inline double genLOGWENOder(double eqm2, double eqm1, double eqc, double eqp1, double eqp2, double dx  )
{ 	
	double w[3], q[3], d[3], alpha[3];
	double wt[3], qt[3], dt[3], alphat[3];
	double beta[4];
	double p,eps=WENOEPS,sum;
	double left=0,right=0;
	
	double qm2 = log(eqm2);
	double qm1 = log(eqm1);
	double qc = log(eqc);
	double qp1 = log(eqp1);
	double qp2 = log(eqp2);
	
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
		
	return(eqc*(left-right)/dx);
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




//~ #ifdef KT
//~ 
//~ void  WenoShuX(GRID HydroGrid)
//~ {	
	//~ double w[3], q[3], d[3], alpha[3];
	//~ double wt[3], qt[3], dt[3], alphat[3];
//~ 
	//~ double   beta[4];
//~ 
	//~ double p,eps=WENOEPS,sum;
	//~ int i,j,k,c;
	//~ 
	//~ d[0]= dt[2]=0.3;
	//~ d[1]= dt[1]=0.6 ;
	//~ d[2]= dt[0]=0.1 ;
	//~ 
 	//~ for(j=jl;j<jr;j++)
	//~ for(k=kl;k<kr;k++) 
	//~ {
		//~ for(i=1;i<XCM-1;i++)
		//~ {
			//~ double qm2,qm1,qc,qp1,qp2;
			//~ double left=0,right=0;
			//~ 
			//~ qc=HydroGrid[i][j][k].VAR;
			//~ qm1=HydroGrid[i-1][j][k].VAR;
			//~ qp1=HydroGrid[i+1][j][k].VAR;
			//~ 
			//~ switch(i)
			//~ {
				//~ case 1:
					//~ qm2=qm1;
					//~ qp2=HydroGrid[i+2][j][k].VAR;
					//~ break;
					//~ 
				//~ case XCM-2:
					//~ qm2=HydroGrid[i-2][j][k].VAR;
					//~ qp2=qp1;
					//~ break;
					//~ 
				//~ default:
					//~ qm2=HydroGrid[i-2][j][k].VAR;
					//~ qp2=HydroGrid[i+2][j][k].VAR;
					//~ break;
			//~ }
//~ 
			//~ beta[0]=(13./12.)*pow(qc - 2*qp1 + qp2 , 2) + 0.25*pow(3*qc - 4*qp1 + qp2 , 2);
			//~ beta[1]=(13./12.)*pow(qm1 - 2*qc + qp1 , 2) + 0.25*pow(qm1 - qp1 , 2);
			//~ beta[2]=(13./12.)*pow(qm2 - 2*qm1 + qc , 2) + 0.25*pow(qm2 - 4*qm1 + 3*qc , 2);
			//~ 
			//~ q[0]= (2*qc  + 5*qp1 - 1*qp2)/6.0;
			//~ q[1]= (-1*qm1 + 5*qc + 2*qp1)/6.0;
			//~ q[2]= (2*qm2 - 7*qm1 + 11*qc)/6.0;
//~ 
			//~ for(c=0;c<3;c++)
				//~ alpha[c]= d[c]/pow(eps + beta[c],WENOP);
//~ 
			//~ sum=0;
			//~ for(c=0;c<3;c++)
				//~ sum += alpha[c];
//~ 
			//~ for(c=0;c<3;c++)
				//~ w[c]= alpha[c]/sum;
			 //~ 
			//~ for(c=0;c<3;c++)
				//~ left += w[c]*q[c];
//~ 
//~ 
			//~ 
			//~ qt[0]= (11*qc - 7*qp1 + 2*qp2)/6.0;
			//~ qt[1]= (2*qm1 + 5*qc - 1*qp1)/6.0;
			//~ qt[2]= (-1*qm2 + 5*qm1 + 2*qc)/6.0;
//~ 
			//~ for(c=0;c<3;c++)
				//~ alphat[c]= dt[c]/pow(eps + beta[c],WENOP);
			//~ 
			//~ sum=0;
			//~ for(c=0;c<3;c++) 
				//~ sum += alphat[c];
			//~ 
			//~ for(c=0;c<3;c++)
				//~ wt[c]= alphat[c]/sum;
			//~ 
			//~ for(c=0;c<3;c++) 
				//~ right += wt[c]*qt[c];
			//~ 
			//~ 
			//~ HydroGrid[i][j][k].left= left;
			//~ HydroGrid[i-1][j][k].right= right;
		//~ }
	//~ }
//~ }
//~ 
//~ 
//~ void  WenoShuY( GRID HydroGrid )
//~ {
	//~ double w[3], q[3], d[3], alpha[3];
	//~ double wt[3], qt[3], dt[3], alphat[3];
//~ 
	//~ double   beta[4];
//~ 
	//~ double p,eps=WENOEPS,sum;
	//~ int i,j,k,c;
//~ 
	//~ 
	//~ d[0]= dt[2]=0.3;
	//~ d[1]= dt[1]=0.6;
	//~ d[2]= dt[0]=0.1;
	//~ 
	//~ 
 	//~ for(i=il;i<ir;i++)
	//~ for(k=kl;k<kr;k++) 
	//~ {
		//~ for(j=1;j<YCM-1;j++)
		//~ {
			//~ double qm2,qm1,qc,qp1,qp2;
			//~ double left=0,right=0;;
			//~ 
			//~ qm1=HydroGrid[i][j-1][k].VAR;
			//~ qc=HydroGrid[i][j][k].VAR;
			//~ qp1=HydroGrid[i][j+1][k].VAR;
			//~ 
			//~ switch(j)
			//~ {
				//~ case 1:
					//~ qm2=qm1;
					//~ qp2=HydroGrid[i][j+2][k].VAR;
					//~ break;
					//~ 
				//~ case YCM-2:
					//~ qm2=HydroGrid[i][j-2][k].VAR;
					//~ qp2=qp1;
					//~ break;
					//~ 
				//~ default:
					//~ qm2=HydroGrid[i][j-2][k].VAR;
					//~ qp2=HydroGrid[i][j+2][k].VAR;
					//~ break;
				//~ 
			//~ }
//~ 
			//~ beta[0]=(13./12.)*pow(qc - 2*qp1 + qp2 , 2) + 0.25*pow(3*qc - 4*qp1 + qp2 , 2);
			//~ beta[1]=(13./12.)*pow(qm1 - 2*qc + qp1 , 2) + 0.25*pow(qm1 - qp1 , 2);
			//~ beta[2]=(13./12.)*pow(qm2 - 2*qm1 + qc , 2) + 0.25*pow(qm2 - 4*qm1 + 3*qc , 2);
//~ 
//~ 
			//~ 
			//~ q[0]= (2*qc + 5*qp1 - 1*qp2)/6.0;
			//~ q[1]= (-1*qm1 + 5*qc + 2*qp1)/6.0;
			//~ q[2]= (2*qm2 - 7*qm1 + 11*qc)/6.0;
//~ 
			//~ for(c=0;c<3;c++)
				//~ alpha[c]= d[c]/pow(eps + beta[c],WENOP);
//~ 
			//~ sum=0;
			//~ for(c=0;c<3;c++)
				//~ sum += alpha[c];
//~ 
			//~ for(c=0;c<3;c++)
				//~ w[c]= alpha[c]/sum;
			//~ 
			//~ for(c=0;c<3;c++)
				//~ left+= w[c]*q[c];
//~ 
//~ 
			//~ 
			//~ qt[0]= (11*qc   - 7 *qp1 + 2*qp2)/6.0;
			//~ qt[1]= (2*qm1 + 5*qc - 1*qp1)/6.0;
			//~ qt[2]= (-1*qm2 + 5*qm1 + 2*qc)/6.0;
//~ 
			//~ for(c=0;c<3;c++)
				//~ alphat[c]= dt[c]/pow(eps + beta[c],WENOP);
			//~ 
			//~ sum=0;
			//~ for(c=0;c<3;c++)
				//~ sum += alphat[c];
			//~ 
			//~ for(c=0;c<3;c++)
				//~ wt[c]= alphat[c]/sum;
			//~ 
			//~ for(c=0;c<3;c++)
				//~ right += wt[c]*qt[c];
			//~ 
			//~ 
			//~ HydroGrid[i][j][k].left=left;
			//~ HydroGrid[i][j-1][k].right=right;
			//~ 
		//~ }
	//~ }
//~ }
//~ 
//~ 
//~ void  WenoShuZ(GRID HydroGrid )
//~ {
	//~ double w[3], q[3], d[3], alpha[3];
	//~ double wt[3], qt[3], dt[3], alphat[3];
//~ 
	//~ double   beta[4];
//~ 
	//~ double p,eps=WENOEPS,sum;
	//~ int i,j,k,c;
//~ 
	//~ 
	//~ d[0]= dt[2]=0.3;
	//~ d[1]= dt[1]=0.6 ;
	//~ d[2]= dt[0]=0.1 ;
	 //~ 
	//~ for(i=il;i<ir;i++)
	//~ for(j=jl;j<jr;j++) 
	//~ {
//~ 
		//~ for(k=2;k<ZCM-2;k++)
		//~ {
			//~ double qm2,qm1,qc,qp1,qp2;
			//~ double left=0, right=0;
			//~ 
			//~ qm2=HydroGrid[i][j][k-2].VAR;
			//~ qm1=HydroGrid[i][j][k-1].VAR;
			//~ qc=HydroGrid[i][j][k].VAR;
			//~ qp1=HydroGrid[i][j][k+1].VAR;
			//~ qp2=HydroGrid[i][j][k+2].VAR;
//~ 
			//~ beta[0]=(13./12.)*pow(qc - 2*qp1 + qp2 , 2) + 0.25*pow(3*qc - 4*qp1 + qp2 , 2);
			//~ beta[1]=(13./12.)*pow(qm1 - 2*qc + qp1 , 2) + 0.25*pow(qm1 - qp1 , 2);
			//~ beta[2]=(13./12.)*pow(qm2 - 2*qm1 + qc , 2) + 0.25*pow(qm2 - 4*qm1 + 3*qc , 2);
//~ 
//~ 
			//~ 
			//~ q[0]= (2*qc + 5*qp1 - 1*qp2)/6.0;
			//~ q[1]= (-1*qm1 + 5*qc + 2*qp1)/6.0;
			//~ q[2]= (2*qm2 - 7*qm1 + 11*qc)/6.0;
//~ 
			//~ for(c=0;c<3;c++)
				//~ alpha[c]= d[c]/pow(eps + beta[c],WENOP);
//~ 
			//~ sum=0;
			//~ for(c=0;c<3;c++)
				//~ sum += alpha[c];
//~ 
			//~ for(c=0;c<3;c++)
				//~ w[c]= alpha[c]/sum;
			//~ 
			//~ for(c=0;c<3;c++)
				//~ left+= w[c]*q[c];
//~ 
//~ 
			//~ 
			//~ qt[0]= (11*qc   - 7 *qp1 + 2*qp2)/6.0;
			//~ qt[1]= (2*qm1 + 5*qc - 1*qp1)/6.0;
			//~ qt[2]= (-1*qm2 + 5*qm1 + 2*qc)/6.0;
//~ 
			//~ for(c=0;c<3;c++)
				//~ alphat[c]= dt[c]/pow(eps + beta[c],WENOP);
			//~ 
			//~ sum=0;
			//~ for(c=0;c<3;c++)
				//~ sum += alphat[c];
			//~ 
			//~ for(c=0;c<3;c++)
				//~ wt[c]= alphat[c]/sum;
			//~ 
			//~ for(c=0;c<3;c++)
				//~ right += wt[c]*qt[c];
			//~ 
			//~ HydroGrid[i][j][k-1].left=left;
			//~ HydroGrid[i][j][k-1].right=right;
			//~ 
		//~ }
	//~ }
//~ }
//~ 
//~ void  VLX( GRID HydroGrid )
//~ {
//~ 
//~ 
	//~ double eps=1e-12,Sl, Sr,S;
	//~ int i,j,k,c;
//~ 
	//~ for(j=jl;j<jr;j++)
	//~ for(k=kl;k<kr;k++)
	//~ {
		//~ for(i=1;i<XCM-1;i++)
		//~ {
			//~ double qm1,qc,qp1;
			//~ 
			//~ qm1=HydroGrid[i-1][j][k].VAR;
			//~ qc=HydroGrid[i][j][k].VAR;
			//~ qp1=HydroGrid[i+1][j][k].VAR;
//~ 
			//~ Sl = (qc-qm1)/XS;
			//~ Sr = (qp1-qc)/XS;
//~ 
//~ 
			//~ S = (SGN(Sl)+SGN(Sr))* ((fabs(Sl)* fabs(Sr))/ (fabs(Sl)+fabs(Sr)+eps));
//~ 
//~ 
			//~ HydroGrid[i][j][k] .left =  qc + S*(XS/2);
			//~ HydroGrid[i-1][j][k].right = qc - S*(XS/2);			
		//~ }
	//~ }
//~ }
//~ 
//~ void  VLY( GRID HydroGrid )
//~ {
//~ 
	//~ double eps=1e-12,Sl, Sr,S;
	//~ int i,j,k,c;
//~ 
	//~ 
	//~ for(i=il;i<ir;i++)
	//~ for(k=kl;k<kr;k++)
	//~ {
		//~ for(j=1;j<YCM-1;j++)
		//~ {
			//~ double qm1,qc,qp1;
			//~ 
			//~ qm1=HydroGrid[i][j-1][k].VAR;
			//~ qc=HydroGrid[i][j][k].VAR;
			//~ qp1=HydroGrid[i][j+1][k].VAR;
//~ 
			//~ Sl = (qc-qm1)/YS;
			//~ Sr = (qp1-qc)/YS;
//~ 
//~ 
			//~ S = (SGN(Sl)+SGN(Sr))* ((fabs(Sl)* fabs(Sr))/ (fabs(Sl)+fabs(Sr)+eps));
//~ 
			//~ HydroGrid[i][j][k].left =  qc + S*(YS/2);
			//~ HydroGrid[i][j-1][k].right= qc -  S*(YS/2);
//~ 
		//~ }
	//~ }
//~ }
//~ 
//~ void  VLZ(  GRID HydroGrid  )
//~ {
//~ 
	//~ double eps=1e-12,Sl, Sr,S;
	//~ int i,j,k,c;
//~ 
	//~ 
	//~ for(i=il;i<ir;i++)
	//~ for(j=jl;j<jr;j++)
	//~ {
		//~ for(k=1;k<ZCM-1;k++)
		//~ {
			//~ double qm1,qc,qp1;
			//~ 
			//~ qm1=HydroGrid[i][j][k-1].VAR;
			//~ qc=HydroGrid[i][j][k].VAR;
			//~ qp1=HydroGrid[i][j][k+1].VAR;
//~ 
			//~ Sl = (qc-qm1)/ZS;
			//~ Sr = (qp1-qc)/ZS;
//~ 
//~ 
			//~ S = (SGN(Sl)+SGN(Sr))* ((fabs(Sl)* fabs(Sr))/ (fabs(Sl)+fabs(Sr)+eps));
//~ 
			//~ HydroGrid[i][j][k] .left=  qc + S*(ZS/2);
			//~ HydroGrid[i][j][k-1].right= qc -  S*(ZS/2);
//~ 
		//~ }
	//~ }
//~ }
//~ 
//~ void  gminmodRX(GRID HydroGrid  )
//~ {
//~ 
	//~ int i,j,k;
//~ 
	//~ for(j=jl;j<jr;j++)
	//~ for(k=kl;k<kr;k++)
	//~ {
		//~ for(i=1;i<XCM-1;i++)
		//~ {
			//~ double qm1,qc,qp1;
					//~ 
			//~ qm1=HydroGrid[i-1][j][k].VAR;
			//~ qc=HydroGrid[i][j][k].VAR;
			//~ qp1=HydroGrid[i+1][j][k].VAR;
//~ 
			//~ double slope1 = genminmod(qm1,qc,qp1,XS, GMINV);
			//~ double num = (qc  - qm1 );
			//~ double denom = (qp1 - qc );
			//~ double r = ratio(num,denom);
			//~ slope1 *= phi(r,-1);
//~ 
//~ 
			//~ HydroGrid[i][j][k] .left =  qc + slope1*(XS/2);
			//~ HydroGrid[i-1][j][k].right = qc - slope1*(XS/2);	
//~ 
		//~ }
	//~ }
//~ 
//~ }
//~ 
//~ 
//~ void  gminmodRY(GRID HydroGrid )
//~ {
//~ 
	//~ int i,j,k;
//~ 
	//~ for(i=il;i<ir;i++)
	//~ for(k=kl;k<kr;k++)
	//~ {
		//~ for(j=1;j<YCM-1;j++)
		//~ {
//~ 
			//~ double qm1,qc,qp1,qp2;
					//~ 
			//~ qm1=HydroGrid[i][j-1][k].VAR;
			//~ qc=HydroGrid[i][j][k].VAR;
			//~ qp1=HydroGrid[i][j+1][k].VAR;
//~ 
			//~ 
//~ 
			//~ double slope1 = genminmod(qm1,qc,qp1,YS,GMINV);
			//~ double num = (qc  - qm1 );
			//~ double denom = (qp1 - qc );
			//~ double r = ratio(num,denom);
			//~ slope1 *= phi(r,-1);
//~ 
		//~ 
			//~ HydroGrid[i][j][k].left=  qc + slope1*(YS/2);
			//~ HydroGrid[i][j-1][k].right= qc -  slope1*(YS/2);
		//~ 
		//~ }
	//~ }
//~ 
//~ }
//~ 
//~ 
//~ void  gminmodRZ( GRID HydroGrid  )
//~ {
//~ 
	//~ int i,j,k;
//~ 
	//~ for(i=il;i<ir;i++)
	//~ for(j=jl;j<jr;j++)
	//~ {
		//~ for(k=1;k<ZCM-1;k++)
		//~ {
//~ 
			//~ double qm1,qc,qp1,qp2;
					//~ 
			//~ qm1=HydroGrid[i][j][k-1].VAR;
			//~ qc=HydroGrid[i][j][k].VAR;
			//~ qp1=HydroGrid[i][j][k+1].VAR;			
//~ 
			//~ double slope1 = genminmod(qm1,qc,qp1,ZS,GMINV);
			//~ double num = (qc  - qm1 );
			//~ double denom = (qp1 - qc );
			//~ double r = ratio(num,denom);
			//~ slope1 *= phi(r,-1);
//~ 
		//~ 
			//~ HydroGrid[i][j][k].left =  qc + slope1*(ZS/2);
			//~ HydroGrid[i][j][k-1].right= qc -  slope1*(ZS/2);
		//~ 
		//~ }
	//~ }
//~ 
//~ }
//~ 
//~ 
//~ #endif




























































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



