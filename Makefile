fortranCompiler=gfortran

all: ode.so 
	@mkdir -p Data
	@mkdir -p Plots

clean:
	rm -rf *.err *.log *.o *.mod *.a
	rm -rf *.so *.so.dSYM
	rm -rf *.pyc __pycache__
	rm -rf Fortran/*.o Fortran/*.mod Fortran/*.a
	@echo "This command does not clean Plots/ and Data/ folders in case data got deleted by accidents. If you do need to clean them, rm them yourself!";

Fortran/rk1412Feagin.f90: Fortran/%Feagin.f90: Python/%.py Convert0.py rkTemplate.f90
	python3 Convert0.py Python/$*.py 35 12 $*Feagin 
Fortran/rk108Feagin.f90: Fortran/%Feagin.f90: Python/%.py Convert0.py rkTemplate.f90
	python3 Convert0.py Python/$*.py 17 8 $*Feagin
Fortran/rk1210Feagin.f90: Fortran/%Feagin.f90: Python/%.py Convert0.py rkTemplate.f90
	python3 Convert0.py Python/$*.py 25 10 $*Feagin
Fortran/rk87EnrightVerner.f90: Fortran/%.f90: Python/%_13.py Convert0.py rkTemplate.f90
	python3 Convert0.py Python/$*_13.py 13 7 $* 
Fortran/rk65Dormand.f90: Fortran/%.f90: Python/%8.py Convert0.py rkTemplate.f90
	python3 Convert0.py Python/$*8.py 8 5 $* 
Fortran/rk1211Peter.f90: Fortran/%Peter.f90: Python/%_31.py Convert0.py rkTemplate.f90
	python3 Convert0.py Python/$*_31.py 31 11 $*Peter
Fortran/rk109Legendre.f90: Fortran/%Legendre.f90: Python/%_21.py Convert0.py rkTemplate.f90
	python3 Convert0.py Python/$*_21.py 21 9 $*Legendre
Fortran/rk87Dormand.f90: Fortran/%.f90: Python/RK8713M.py Convert0.py rkTemplate.f90
	python3 Convert0.py Python/RK8713M.py 13 7 $*
Fortran/rk54Dormand.f90: Fortran/%.f90: Python/%7.py Convert0.py rkTemplate.f90
	python3 Convert0.py Python/$*7.py 7 5 $*
Fortran/rk54Sharp.f90: Fortran/%.f90: Python/%7.py Convert0.py rkTemplate.f90
	python3 Convert0.py Python/$*7.py 7 5 $*

# the available solvers 
SolverModules:= rk54Sharp rk54Dormand rk65Dormand rk87Dormand rk87EnrightVerner rk108Feagin rk109Legendre rk1210Feagin rk1211Peter rk1412Feagin 
Solvers: $(SolverModules:%=Fortran/%.f90) 

Modules:= GlobalCommon DyDt $(SolverModules) 
$(Modules:%=Fortran/%.o): Fortran/%.o: Fortran/%.f90
	$(fortranCompiler) -c -J Fortran -o Fortran/$*.o Fortran/$*.f90 -fPIC -O3 -fmax-errors=5 
ode.so: Fortran/ODEInterface.f90 Fortran/libSolvers.a 
	rm -rf $@ ode.*.so ode.*.so.dSYM
	bash CompileInterface.sh Fortran/ODEInterface.f90 ode Fortran/libSolvers.a
Fortran/libSolvers.a: $(Modules:%=Fortran/%.o)  
	ar crs $@ $^

ProgramExtra:= ODEInterface Main
Program:= $(Modules) $(ProgramExtra) 
$(ProgramExtra:%=Fortran/%.o): Fortran/%.o: Fortran/%.f90
	$(fortranCompiler) -c -J Fortran -o Fortran/$*.o Fortran/$*.f90 -fPIC -O3 -fmax-errors=5 
main: $(Program:%=Fortran/%.o) 
	$(fortranCompiler) -o main $(Program:%=Fortran/%.o)

