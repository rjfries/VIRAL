
#define E1 ( 0.5028563305441270) 		// in gev/fm3
#define E2 ( 1.62)  					// in gev/fm3
#define E3 ( 1.86)  					// in gev/fm3
#define E4 ( 9.9878355786273545) 	    // in gev/fm3

#define ES1 ( 0.1270769021427449 ) 		// in gev/fm3
#define ES2 ( 0.446707952467404  ) 		// in gev/fm3
#define ES3 ( 1.9402832534193788 ) 		// in gev/fm3
#define ES4 ( 3.7292474570977285 ) 	    // in gev/fm3


double s95p_p(double e) //e in 1/fm^4
{
	double ret;
	e /= GEVFM; //e in gev/fm3
		
	if(e<E1)
		ret = (0.3299*(exp(0.4346*e) -1 ));
	else
	if(e>E1 && e<E2)
		ret = (1.024E-7*exp(6.041*e) + 0.007273 + 0.14578*e);
	else
	if(e>E2 && e<E3)
		ret = (0.30195*exp(0.31308*e) - 0.256232);
	else
	if(e>E3 && e<E4)
		ret = (0.332*e - 0.3223*pow(e,0.4585) + 0.1167*pow(e,-1.233) - 0.003906*e*exp(-0.05697*e)  + 0.1436*e*exp(-0.9131*e) );
	else
		ret = (0.3327*e - 0.3223*pow(e,0.4585) - 0.003906*e*exp(-0.05697*e) );
	
	
	
	//fix for very small eps by sid
	if(ret==0)
		ret=e;
		
		
	return (GEVFM*ret); //p returned in 1/fm^4
}

double s95p_a(double e) //speed of sound squared --- dp/de
{
	double ret;
	
	e /= GEVFM; //e in gev/fm3
	
	
	if(e<E1)
		ret = 1*( 0.3299*(0.4346*exp(0.4346*e)) );
	else
	if(e>E1 && e<E2)
		ret = 1*( 1.024E-7*6.041*exp(6.041*e) + 0.14578 );
	else
	if(e>E2 && e<E3)
		ret = 1*( 0.30195*0.31308*exp(0.31308*e)  );
	else
	if(e>E3 && e<E4)
		ret = 1*( 0.332  - 0.3223*0.4585*pow(e,(0.4585-1))  + 0.1167*(-1.233)*pow(e,(-1.233-1)) - 0.003906*( exp(-0.05697*e) - 0.05697*e*exp(-0.05697*e)) +  0.1436*(exp(-0.9131*e) - 0.9131*e*exp(-0.9131*e))  );
	else
		ret = 1*( 0.3327 - 0.3223*0.4585*pow(e,(0.4585-1))                                      - 0.003906*( exp(-0.05697*e) - 0.05697*e*exp(-0.05697*e)) );
	
	
	return ret; //unit-less
}


double s95p_s(double e) //e in 1/fm^4
{
	double ret;
	
	e /= GEVFM; //e in gev/fm3
	
	
	if(e<ES1)
		ret = 12.2304 * pow(e, 1.16849);
	else
	if(e>ES1 && e<ES2)
		ret = 11.9279 * pow(e, 1.15635);
	else
	if(e>ES2 && e<ES3)
		ret = 0.0580578 + 11.833*pow(e, 1.16187);
	else
	if(e>ES3 && e<ES4)
		ret = 18.202*e - 62.021814 - 4.85479*exp(-2.72407e-11 * pow(e, 4.54886)) +
        65.1272 * pow(e, -0.128012) * exp(-0.00369624 * pow(e, 1.18735)) -
        4.75253 * pow(e, -1.18423);
	else
		ret =  18.202 * e - 63.0218 - 4.85479 * exp(-2.72407e-11 * pow(e, 4.54886)) +
          65.1272 * pow(e, -0.128012) * exp(-0.00369624 * pow(e, 1.18735));
          
    return pow(ret,0.75);	//in 1/fm3
}


double s95p_T(double e) //e in 1/fm^4
{
	double ret;
	
	e /= GEVFM; //e in gev/fm3
	
	
	if(e<0.5143939846236409)
		ret = 0.203054*pow(e,0.30679);
	else
		ret = (e+s95p_p(e))/s95p_s(e) ;
		
	return ret; // in 1/fm
}


double s95p_TGev(double e) //e in 1/fm^4
{
	double ret;
	
	e /= GEVFM; //e in gev/fm3
	
	if(e<0.5143939846236409)
		ret = 0.203054*pow(e,0.30679);
	else
		ret = (e+s95p_p(e))/s95p_s(e) ;
		
	return ret/GEVFM; // in gev
}

