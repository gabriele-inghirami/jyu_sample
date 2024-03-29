# Uncomment/edit your desired compiler/compilation flags and
# comment all the rest until the ###### line

# gcc, general optimization, not cpu specific
#com = gfortran -O3 -fdefault-real-8 -cpp

# gcc, cpu specific optimization (not portable)
com = gfortran -O3 -fdefault-real-8 -march=native -mtune=native -cpp

# gcc, for debugging
#com = gfortran -O0 -g -ggdb3 -fno-strict-overflow -fdefault-real-8 -fcheck=bounds -cpp

# intel, example for skylake or more modern cpus 
#com = ifort -O3 -r8 -march=skylake -no-wrap-margin -fpp

#####################################################################

all: sample decay

obj_sample = numlib.o settings.o common.o  eos.o sample_mc.o
obj_decay = eos.o decay.o

sample:		${obj_sample}
		${com} ${obj_sample} -o sample.exe		

decay:		${obj_decay}
		${com} ${obj_decay} -o decay.exe		
		
		
numlib.o:		numlib.f90 makefile
			${com} -c numlib.f90

settings.o:		settings.f90  makefile
			${com} -c settings.f90
			
common.o:		common.f90  makefile
			${com} -c common.f90
			
eos.o:			eos.f90  makefile
			${com} -c eos.f90			

sample_mc.o:		sample_mc.f90 settings.f90 common.f90 eos.f90 numlib.f90 makefile
			${com} -c sample_mc.f90			

decay.o:		decay.f90 eos.f90 makefile
			${com} -c decay.f90			
			

clean:		
		\rm *.mod *.o  *.exe 
