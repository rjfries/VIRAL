CC = mpiicpc
GCC = g++
WALL = -Wall
O =   -O3 #-Ofast 
PROFILE = -DRK3	  -DKT -DLBI  -DGINIT -DS95P #-DGINIT-DVORT #   -DGUBSER #-DCON#-DGINIT -DGUBSER
DEBUG =       -g
INCLUDES = 
COMPILER = $(CC) $(WALL) $(O) $(INCLUDES) $(PROFILE) 

COMPILERD = $(CC) $(WALL) $(DEBUG) $(INCLUDES) $(PROFILE) 

INCLUDEL = -IC:/MinGW/include -IC:/MinGW/msys/1.0/local/include  
COMPILERL = $(GCC) -g $(INCLUDEL) 

LIBS = -lgsl -lgslcblas
LIBPATH= 
LIBPATHL= -LC:/MinGW/msys/1.0/local/lib

all: hydro 
 
T:t00.o
	$(COMPILERL) t00.o   $(LIBPATHL)  $(LIBS) -o t00.exe

t00.o:t00.cpp
	$(COMPILERL) -c t00.cpp 



d: 	hydrod 

main.o:main.cpp
	$(COMPILER) -c main.cpp
	
hydro: main.o 
	$(COMPILER)  main.o $(LIBPATH) $(LIBS) -o hydro
	rm -f *.o
	
s: main.o 
	$(COMPILER)  main.o $(LIBPATH) $(LIBS) -o hydro
	rm -f *.o
	bsub <  batch 
	
maind.o:main.cpp
	$(COMPILERD) -c main.cpp

hydrod: maind.o
	$(COMPILERD) main.o $(LIBPATH) $(LIBS) -o hydrod
	rm -f *.o


init:init.cpp
	$(COMPILER) -c init.cpp  -o i
	
clean:
	rm -f hydro 
	rm -f hydrod
	rm -f *.exe
	rm -f *~
	rm -f *.o
	rm -f *.qsub
	rm -f stdout
	rm -f core*

cleanO:	
	rm -f res/*.bin
	rm -f init/*.bin
	rm -f hydro 
	rm -f hydrod
	rm -f *.exe
	rm -f *~
	rm -f *.o
	rm -f *.qsub
	rm -f std*
	rm -f core*

cleanI:
	rm -f init/init.dat
	rm -f init/*.bin
	rm -f hydro 
	rm -f hydrod
	rm -f *.exe
	rm -f *~
	rm -f *.o
	rm -f *.qsub
	rm -f stdout
	rm -f core*

