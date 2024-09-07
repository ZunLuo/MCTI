# general include path
INCLUDE=-I/usr/include


#include objects here:
objects = src/add_CTI1_frac_no_output.o src/nrutil.o src/ran1_frac.o src/ran2_frac.c src/ran3_frac.c src/ran4_frac.c src/poidev_frac.o src/poidev1_frac.o src/gammln_frac.o src/gasdev.o src/sort_frac.o src/creattraps_frac.o

add_CTI1_frac_no_output.so:  $(objects)
	gcc -shared -fPIC -o $@ $(objects) -lm

# general compilation rules 
.SUFFIXES:  .c   .o
.c.o:
	cc -c $< -O3 -shared -fPIC -o $@ $(INCLUDE)
