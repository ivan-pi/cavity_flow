FC = gfortran
#flags = -Wall -pedantic
FCFLAGS = -O3 -Wall -Wextra -march=native -fopenmp
GC = gcc
CXX = g++
CCFLAGS = -O3 -Wall -Wextra -march=native -fopenmp
LDFLAGS = -lstdc++ -lgomp


objects = incbeta.o my_ramp.o precision_mod.o accelerate_mod.o vtk_mod.o d2q9_mod.o cavity_test.o

cavity_test : $(objects)
	${FC} $^ -o $@ ${LDFLAGS}

%.o : %.f90
	${FC} ${FCFLAGS} -c $< -o $@

%.o : %.cpp
	${CXX} ${CCFLAGS} -c $^ -o $@	

%.o : %.c
	${GC} ${CCFLAGS} -c $^ -o $@	

.PHONY : clean
clean :
	rm *.mod
	rm *.o

.PHONY : fresh
fresh:
	rm *.out
	rm *.mod
	rm *.o
