fortranCompiler=gfortran

all: ode.so 

clean:
	rm -rf *.err *.log *.o *.mod *.a
	rm -rf *.so *.so.dSYM
	rm -rf *.pyc __pycache__/*.pyc

rk1412Feagin.f90: %Feagin.f90: Python/%.py Convert0.py
	python3 Convert0.py Python/$*.py 35 12 $*Feagin 
rk108Feagin.f90: %Feagin.f90: Python/%.py Convert0.py
	python3 Convert0.py Python/$*.py 17 8 $*Feagin
rk1210Feagin.f90: %Feagin.f90: Python/%.py Convert0.py
	python3 Convert0.py Python/$*.py 25 10 $*Feagin
rk87EnrightVerner.f90: %.f90: Python/%_13.py Convert0.py
	python3 Convert0.py Python/$*_13.py 13 7 $* 
rk65Dormand.f90: %.f90: Python/%8.py Convert0.py
	python3 Convert0.py Python/$*8.py 8 5 $* 
rk1412Long.f90: %.f90: Python/%.py Convert0.py
	python3 Convert0.py Python/$*.py 35 12 $*
rk1211Peter.f90: %Peter.f90: Python/%_31.py Convert0.py 
	python3 Convert0.py Python/$*_31.py 31 11 $*Peter
rk109Legendre.f90: %Legendre.f90: Python/%_21.py Convert0.py
	python3 Convert0.py Python/$*_21.py 21 9 $*Legendre
rk87Dormand.f90: %.f90: Python/RK8713M.py Convert0.py
	python3 Convert0.py Python/RK8713M.py 13 7 $*
rk54Dormand.f90: %.f90: Python/%7.py Convert0.py
	python3 Convert0.py Python/$*7.py 7 5 $*
rk54Sharp.f90: %.f90: Python/%7.py Convert0.py
	python3 Convert0.py Python/$*7.py 7 5 $*

SolverModules:= rk54Sharp rk54Dormand rk65Dormand rk87Dormand rk87EnrightVerner rk108Feagin rk109Legendre rk1210Feagin rk1211Peter rk1412Long rk1412Feagin 
Solvers: $(SolverModules:%=%.f90) 

Modules:= GlobalCommon DyDt $(SolverModules) 
$(Modules:%=%.o): %.o: %.f90
	$(fortranCompiler) -c -o $*.o $*.f90 -fPIC -O3 -fmax-errors=5 
ode.so: ODEInterface.f90 libSolvers.a 
	rm -rf $@ ode.*.so ode.*.so.dSYM
	bash CompileInterface.sh ODEInterface.f90 ode libSolvers.a
libSolvers.a: $(Modules:%=%.o)  
	ar crs $@ $^

