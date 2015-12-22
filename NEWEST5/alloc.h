void StoreInHeap();
void AllocatePack();

PCKX leftX ,rightX;
PCKY leftY ,rightY;
PCKZ leftZ ,rightZ;

PCKXY xlyl,xlyr,xryl,xryr;

void AllocatePack()
{

	int i,j,k,l;

	leftX = new double[PACKVAR][BORDER][YCMA][ZCMA]();  //certified contigous .. first element at [0][0][0] and last at [max][max][max] and initialised to zero
	rightX = new double[PACKVAR][BORDER][YCMA][ZCMA]();
	leftY = new double[PACKVAR][BORDER][XCMA][ZCMA]();
	rightY = new double[PACKVAR][BORDER][XCMA][ZCMA]();
	leftZ = new double[PACKVAR][BORDER][XCMA][YCMA]();
	rightZ = new double[PACKVAR][BORDER][XCMA][YCMA](); 

	xlyl = new double[PACKVAR][BORDER][BORDER][ZCMA]();
	xlyr = new double[PACKVAR][BORDER][BORDER][ZCMA]();
	xryl = new double[PACKVAR][BORDER][BORDER][ZCMA]();
	xryr = new double[PACKVAR][BORDER][BORDER][ZCMA](); 
}


#ifdef KT
CAPXY AzLZ,AzRZ,VarLZ,VarRZ,FzLZ,FzRZ;
CAPXY FluxT;
#endif

void StoreInHeap()
{
	HydroGrid = new cell  [XCM][YCM][ZCMA] (); //all initialised to 0.0 because of --> () C++ zindabad

#ifdef KT	
	AzLZ = new double [EVAR][XCM][YCM]();
	AzRZ = new double [EVAR][XCM][YCM]();
	VarLZ = new double [SVAR][XCM][YCM]();
	VarRZ = new double [SVAR][XCM][YCM]();
	FzLZ = new double [SVAR][XCM][YCM]();
	FzRZ = new double [SVAR][XCM][YCM]();	
	FluxT = new double [SVAR][XCM][YCM]();	
#endif

}


void FreeFromHeap()
{
	delete[] leftX;
	delete[] rightX;
	delete[] leftY;
	delete[] rightY;
	delete[] leftZ;
	delete[] rightZ;
	

#ifdef KT	
	delete[] AzLZ;
	delete[] AzRZ;
	delete[] VarLZ;
	delete[] VarRZ;
	delete[] FzLZ;
	delete[] FzRZ;
	
	delete[] FluxT;
#endif
	
	delete[] HydroGrid;
}

