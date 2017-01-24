all:	inflight.cxx sterile_flux.h sterile_flux.cxx fourmomentum.h fourmomentum.cxx channel.h channel.cxx minInstance.c minInstance.h
	g++ -g -std=c++11  -Wno-deprecated-declarations -c inflight.cxx -o inflight.o -I. -I $$ROOTSYS/include
	g++ -g -std=c++11 -Wno-deprecated-declarations -c minInstance.c -o minInstance.o -I. -I $$ROOTSYS/include
	g++ -g -std=c++11  -c fourmomentum.cxx -o fourmomentum.o -I.
	g++ -g -std=c++11  -c sterile_flux.cxx -o sterile_flux.o -I.  
	g++ -g -std=c++11  -c channel.cxx -o channel.o -I. 
	g++ -g -std=c++11  -c detector.cxx -o detector.o -I.  
	g++ -g -std=c++11 -Wno-deprecated-declarations  -o inflight inflight.o detector.o minInstance.o fourmomentum.o sterile_flux.o channel.o   -L/home/mark/programs/root/root/lib/libMathMore.so -lMathMore   `gsl-config --libs --cflags`  `root-config --libs` 
	rm *.o

#-L/usr/lib/x86_64-linux-gnu/root5.34/libMathCore.so -lMathCore 
#/usr/include/gsl/
