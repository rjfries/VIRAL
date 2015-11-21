#define PIE 3.141592653589793


/*  
 *    x massless quark degrees of freedom
 * 	  gq = 2·2·3 = 12 for quarks
 *    gg = 16 for gluons
 * 	  
 * 	DOF = x*gq*(7/8) + (2*8)
 */

#define DOF 42.25  //x=2.5


#define FACTOR (3.0*DOF*PIE*PIE/90.0)



inline double EOS(double en , double r)
{
	double ret;
	
	ret = (en/3.0);//*WOODSAXON(r, 0.5, 18); 
		
	return(ret);
}

inline double DPDE(double en , double r)
{
	double ret;
	
	ret = (1.0/3.0);//*WOODSAXON(r, 0.5, 18);
		
	return(ret);
}
