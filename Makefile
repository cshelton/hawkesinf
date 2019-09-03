#CC = /usr/bin/g++ --std=c++1y
CC = g++ --std=c++1y
#OPT = -ggdb
OPT = -O3

timetst: timetst.cpp plottraj.hpp gnuplot-cpp/gnuplot_i.hpp traj.hpp \
		 hp.hpp missingstd.hpp kernels.hpp
	$(CC) $(OPT) -o timetst timetst.cpp -L/usr/local/lib -lboost_program_options -lrt

runcrime: runcrime.cpp traj.hpp hp.hpp kernels.hpp mle.hpp missingstd.hpp 
	$(CC) $(OPT) -o runcrime runcrime.cpp -L/usr/local/lib -lpthread -lboost_program_options -lrt

tstcrime: tstcrime.cpp traj.hpp hp.hpp kernels.hpp mle.hpp missingstd.hpp tstpred.hpp
	$(CC) $(OPT) -o tstcrime tstcrime.cpp -L/usr/local/lib -lpthread -lboost_program_options -lrt

example: example.cpp traj.hpp hp.hpp kernels.hpp mle.hpp
	$(CC) $(OPT) -o example example.cpp -L/usr/local/lib -lpthread -lrt
	
clean:
	rm timetst runcrime tstcrime example
