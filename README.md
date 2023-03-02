# Pleiades problem 

Pleiades problem is a 7 body 2-dimensional gravitating system. At t=0, we have
```
y0 = [3,3,-1,-3,2,-2,2,3,-3,2,0,0,-4,4,0,0,0,0,0,1.75,-1.5,0,0,0,-1.25,1,0,0].
```
At the end of integration t=3, we have

```
x1 = 0.370613914397051269489224978315
y1 = 3.237284092057233220174339294317
x2 = -3.222559032418323532453996449476
y2 = 0.659709145577530797766030445928
x3 = 0.342558170715657972493772831513
y3 = 1.562172101400631119005879554607
x4 = -0.700309292221249490850709662482
y4 = -3.943437585517392207634657097515
x5 = -3.271380973972549899286832442158
y5 = 5.225081843456544028470034390921
x6 = -2.590612434977469291652596439235
y6 = 1.198213693392274681315257112146
x7 = -0.242968234493582346811280103793
y7 = 1.091449240428979727823843859369
dx1/dt = 3.417003806314314662273545764037
dy1/dt = 1.354584501625501147259456047323
dx2/dt = -2.590065597810775610554401282570
dy2/dt = 2.025053734714241215897345682606
dx3/dt = -1.155815100160449082622449168412
dy3/dt = -0.807298817022302217516482869542
dx4/dt = 0.595239635420871882054427715047
dy4/dt = -3.741244961234008403749839999364
dx5/dt = 0.377345968575062917782503291164
dy5/dt = 0.938685886955107906537421058601
dx6/dt = 0.366792222720056959595069656643
dy6/dt = -0.347404635380849424741711573006
dx7/dt = 2.344915448180937111999355693115
dy7/dt = -1.947020434263291965848452491628
```
The solution above should be accurate within about 25-30 digits. 


# Availble methods 
The available methods are: 
{'rk54Sharp', 'rk54Dormand', 'rk65Dormand', 'rk87Dormand', 'rk87EnrightVerner', 'rk108Feagin', 'rk109Legendre', 'rk1210Feagin', 'rk1211Peter', 'rk1412Long', 'rk1412Feagin'}
All of these methods are explicite RK methods. 
You have two choices: f2py (which serves as the interface between python and fortran) or a pure fortran code.

## Download
```
git clone git@github.com:kayejlli/RKAdaptiveMethods.git
```


## Set up using f2py 
f2py has the advantage that you can change your inputs easily using python, it generates initial data using python and pass these info to the fortran program. 
The structure of the fortran code is:
<p float="left">
 <img src="Plots/odef2py.png" width="45%" height="45%">
</p>
All the fortran scripts (except Fortran/Main.f90, which is not used by f2py) are compiled into a library "ode.so". 
The python script (odef2py.py) imports this library and pass the initial values, settings for solver to the fortran code. 

In order to use f2py, you need to configure your system environment. 
You should have f2py and python packages like numpy, mpmath, matplotlib, time installed. f2py is usually located inside the Python's bin, something like:
```
/Users/MyUserName/Library/Python/3.7/bin/
```
Make sure you add this path to your system $PATH. Then you may simply run
```
f2py 
```
which will output a long string with all the version information, etc. 

Once you have successfully set up your system environment, you may start compiling the code. 
In the top directory, you can run
```
make all
```
which will generate "ode.so". If successful, you can now run the test code by: 
```
python3 odef2py.py
```
which should generate a file "Plots/odef2py_xy.png" that looks like:
<p float="left">
 <img src="Plots/odefortran_xy_example.png" width="60%" height="60%">
</p>
To compare these different solvers, you can run
```
python3 PleiadesMethods.py
```
which will takes about 27 second (on my Mac), and output two files: 
```
Plots/PleiadesMethods_scd_n.png
Plots/PleiadesMethods_scd_fcn.png
```
which should look like 
<p float="left">
 <img src="Plots/PleiadesMethods_scd_n_example.png" width="45%" height="45%">
 <img src="Plots/PleiadesMethods_scd_fcn.png" width="45%" height="45%">
</p>
where N represents the total number of steps taken during the integration, and fcn represents the number of times that y'=f(t,y) is evaluated. 

## Using a pure fortran code
The structure of the fortran code is:
<p float="left">
 <img src="Plots/Dependencies2.png" width="45%" height="45%">
</p>
Initial values, settings for the solver are all defined in Main.f90. To run the fortran code,
  simpily compile and execute it using the command:
```
make main
./main
```
The data file "Data/test.dat" will be generated, and you can use 
```
python3 odefortran.py 
```
to load the data and draw a figure. It should look like the figure shown above. 


## Quadrupole precision:
The coefficients are calculated up to 60 digits for each of these method, and it is straight forward to switch to quad precision by changing all 'D-' and 'D\d' into 'Q-' and 'Q\d'. However, f2py does not work with quadruple precision, and it convert all numbers into KIND=8 before passing them to fortran.
Note that numpy does not work with quadrupole precision. The np.float128 is not quadrupole precision, but "extended-precision", and compatible with C ``long double`` but not necessarily with IEEE 754 quadruple-precision.  

Therefore, if you wish to pass quadrupole precision number using f2py, the working method is saving all the number to a *.f90 file, and re-compile this using f2py, and then just exceute this program. 
There is an example in branch "Quad". 

The results from quadrupole precision are:
<p float="left">
 <img src="Plots/PleiadesMethods_scd_n_QP.png" width="45%" height="45%">
 <img src="Plots/PleiadesMethods_scd_fcn_QP.png" width="45%" height="45%">
</p>

# How to adapt this code to your probelm 
The structure of the fortran code is:
<p float="left">
 <img src="Plots/Dependencies2.png" width="45%" height="45%">
</p>

You will need to modify Fortran/GlobalCommon.f90, Fortran/DyDt.f90 and Fortran/Main.f90.
Fortran/Main.f90 is the main program, where you define the filename for the data file,
  the initial values, the time interval, time steps, etc. 
The differential equation (i.e. y' = f(t,y) with y(n)) is defined in DyDt.f90, which you need to modify.
The subroutine should look like dev(t,y,dy,test), which takes in t,y(n) and output dy(n) and test.
The file GlobalCommon.f90 is used to pass common variables and parameters in the code.  



If you would like to use f2py instead, you may skip Fortran/Main.f90 and only modify Fortran/GlobalCommon.f90 and Fortran/DyDt.f90.
You will need to modify odepython.py. This file is the interface between fortran and python codes, but can also serve 
  as a test program.  
The subroutine "solve_ivp_" should work for any other problems 
  and you may want to keep it what it is. 
You may want to modify the lines after 
```
if __name__ == '__main__':
```
where the initial values, setting of the solver, etc, are defined. 

## The ODE formula
By default, this code uses dev(t,y,dy,test) which takes in t,y(n) and output dy(n) and test.
If you have extra variables that you would like to pass to the ODE (e.g. if you have y' = f(t, y, a, b), where a,b are variables),
  it is suggested to put a,b in GlobalCommon.f90 and use them as global variables. 

However, the situation can be quite complicated sometimes, when this set up will fail. 
In this case, you might have no choice but modifying dev(t,y,dy,test) to dev(t,y,dy,test,a,b)
  with a,b are variables either INTENT(IN) or INTENT(OUT) or INTENT(INOUT).
You will also need to update rk*f90 as well. You can certainly modify them one by one, but this process can be quite tedious.
There are ways to modify all of them quickly. You can modify the line 101, line 124 to 127 of the file
```
Convert0.py
```
and then execute it. New rk*f90 files will be automatically generated by running: 
```
make Solvers
```
Remember to update the ABSTRACT INTERFACE in GlobalCommon.f90 which should resemble your new dev() subroutine.


# Source and reference 
The Pleiades problem is from:https://archimede.dm.uniba.it/~testset/problems/plei.php although this link seems to be down now. 
*_raw.py are downloaded from PeterStone's website (85 digits)
*.txt are downloaded from https://sce.uhcl.edu/rungekutta/ (60 digits)

