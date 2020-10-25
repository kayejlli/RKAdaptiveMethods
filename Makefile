fortranCompiler=gfortran

clean:
	rm -rf *.err *.log *.o *.mod *.a
	rm -rf *.so *.so.dSYM
	rm -rf *.pyc __pycache__/*.pyc
	rm -rf rk109_21.py
	rm -rf rk65Dormand8.py
	rm -rf rk1412Long.py
	rm -rf rk54Dormand7.py
	rm -rf rk54Sharp7.py
	rm -rf rk87EnrightVerner_13.py
	rm -rf rk1211_31.py
	rm -rf RK8713M.py


Feagin:= rk108 rk1210 rk1412

sce:= $(Feagin:%=%.txt) 

$(sce):
	wget https://sce.uhcl.edu/rungekutta/$@ --no-check-certificate


# s/\s*\n/|\\\r/gc

rk54Sharp7.py: %.py: %_raw.py PeterStone.py
	rm -rf $@ 
	python3 $< >> $@
	echo "from CheckCoef import *" >> $@
	echo "Print2(7, globals())" >> $@ 
rk54Dormand7.py: rk54Dormand7_raw.py PeterStone.py
	rm -rf $@ 
	python3 $< >> $@
	echo "from CheckCoef import *" >> $@
	echo "Print2(7, globals())" >> $@ 
rk1412Long.py: rk1412Long_raw.py PeterStone.py
	rm -rf $@ 
	python3 $< >> $@
	echo "from CheckCoef import *" >> $@
	echo "Print(35, globals())" >> $@ 
RK8713M.py: RK8713M_raw.py PeterStone.py
	rm -rf $@ 
	python3 $< >> $@
	echo "from CheckCoef import *" >> $@
	echo "Print2(13, globals())" >> $@ 
rk109_21.py: rk109_21_raw.py PeterStone.py
	rm -rf $@ 
	python3 $< >> $@
	echo "from CheckCoef import *" >> $@
	echo "Print2(21, globals())" >> $@ 
rk1211_31.py: rk1211_31_raw.py PeterStone.py
	rm -rf $@ 
	python3 $< >> $@
	echo "from CheckCoef import *" >> $@
	echo "Print2(31, globals())" >> $@ 
rk87EnrightVerner_13.py: rk87EnrightVerner_13_raw.py PeterStone.py
	rm -rf $@ 
	python3 $< >> $@
	echo "from CheckCoef import *" >> $@
	echo "Print2(13, globals())" >> $@ 
rk65Dormand8.py: rk65Dormand8_raw.py PeterStone.py
	rm -rf $@ 
	python3 $< >> $@
	echo "from CheckCoef import *" >> $@
	echo "Print2(8, globals())" >> $@ 

PeterStoneRK: rk65Dormand8.py rk87EnrightVerner_13.py rk1211_31.py rk109_21.py RK8713M.py rk1412Long.py rk54Sharp7.py rk54Dormand7.py

rk1412Feagin.f90: %Feagin.f90: %.py Convert0.py
	python3 Convert0.py $*.py 35 12 $*Feagin 
rk108Feagin.f90: %Feagin.f90: %.py Convert0.py
	python3 Convert0.py $*.py 17 8 $*Feagin
rk1210Feagin.f90: %Feagin.f90: %.py Convert0.py
	python3 Convert0.py $*.py 25 10 $*Feagin
rk87EnrightVerner.f90: %.f90: %_13.py Convert0.py
	python3 Convert0.py $*_13.py 13 7 $* 
rk65Dormand.f90: %.f90: %8.py Convert0.py
	python3 Convert0.py $*8.py 8 5 $* 
rk1412Long.f90: %.f90: %.py Convert0.py
	python3 Convert0.py $*.py 35 12 $*
rk1211Peter.f90: %Peter.f90: %_31.py Convert0.py 
	python3 Convert0.py $*_31.py 31 11 $*Peter
rk109Legendre.f90: %Legendre.f90: %_21.py Convert0.py
	python3 Convert0.py $*_21.py 21 9 $*Legendre
rk87Dormand.f90: %.f90: RK8713M.py Convert0.py
	python3 Convert0.py RK8713M.py 13 7 $*
rk54Dormand.f90: %.f90: %7.py Convert0.py
	python3 Convert0.py $*7.py 7 5 $*
rk54Sharp.f90: %.f90: %7.py Convert0.py
	python3 Convert0.py $*7.py 7 5 $*

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

QPPlots:= Plots/PleiadesMethods_error_QP.png Plots/PleiadesMethods_scd_cpu_QP.png Plots/PleiadesMethods_scd_fcn_QP.png Plots/PleiadesMethods_x1y1_QP.png
$(QPPlots): # downloaded from mssllay 
	bash PleiadesMethodsQPPlots.sh

