void StoreInHeap();
void AllocatePack();

PCKX leftX ,rightX;
PCKY leftY ,rightY; 

PCKXY xlyl,xlyr,xryl,xryr;

void AllocatePack()
{

	int i,j,k,l;

	leftX = new double[PACKVAR][BORDER][YCMA][ZCMA]();  //certified contigous .. first element at [0][0][0] and last at [max][max][max] and initialised to zero
	rightX = new double[PACKVAR][BORDER][YCMA][ZCMA]();
	leftY = new double[PACKVAR][BORDER][XCMA][ZCMA]();
	rightY = new double[PACKVAR][BORDER][XCMA][ZCMA](); 

	xlyl = new double[PACKVAR][BORDER][BORDER][ZCMA]();
	xlyr = new double[PACKVAR][BORDER][BORDER][ZCMA]();
	xryl = new double[PACKVAR][BORDER][BORDER][ZCMA]();
	xryr = new double[PACKVAR][BORDER][BORDER][ZCMA](); 
}

 

void StoreInHeap()
{
	HydroGrid = new cell  [XCM][YCM][ZCM] (); //all initialised to 0.0 because of --> () C++ zindabad
 
}


void FreeFromHeap()
{
	delete[] leftX;
	delete[] rightX;
	delete[] leftY;
	delete[] rightY;  
	
	delete[] HydroGrid;
}

