fortranCompiler=gfortran

all: ode.so 
	@mkdir -p Data
	@mkdir -p Plots

clean:
	rm -rf *.err *.log *.o *.mod *.a
	rm -rf *.so *.so.dSYM
	rm -rf *.pyc __pycache__
	@echo "This command does not clean Plots/ and Data/ folders in case data got deleted by accidents. If you do need to clean them, rm them yourself!";

# the available solvers 
SolverModules:= rk54Sharp rk54Dormand rk65Dormand rk87Dormand rk87EnrightVerner rk108Feagin rk109Legendre rk1210Feagin rk1211Peter rk1412Feagin 
Solvers: $(SolverModules:%=%.f90) 

Modules:= GlobalCommon DyDt $(SolverModules) 
$(Modules:%=%.o): %.o: %.f90
	$(fortranCompiler) -c -o $*.o $*.f90 -fPIC -O3 -fmax-errors=5 
ode.so: ODEInterface.f90 libSolvers.a 
	rm -rf $@ ode.*.so ode.*.so.dSYM
	bash CompileInterface.sh ODEInterface.f90 ode libSolvers.a
libSolvers.a: $(Modules:%=%.o)  
	ar crs $@ $^
