#!/bin/bash
# CompileInterface.sh ODEInterface.f90 ode libSolvers.a 
FortranCompiler="$(which gfortran)";
FortranFile=$1
StaticObjectName=$2
Library=$3

if [ -z "$Library" ]; then
  f2py -c $FortranFile -m $StaticObjectName --opt='-O3' --f90exec=$FortranCompiler > $StaticObjectName.log 2> $StaticObjectName.err;
else
  LibraryName="${Library/lib/}"
  LibraryName="${LibraryName/.a/}"
  f2py -L. -l$LibraryName -c $FortranFile -m $StaticObjectName --opt='-O3' --f90exec=$FortranCompiler > $StaticObjectName.log 2> $StaticObjectName.err;
fi

# f2py -L. -lSolvers -c ODEInterface.f90 -m ode --opt='-O3' --f90exec=$FortranCompiler > ode.log 2> ode.err;
if ls ./$StaticObjectName.*.so 1> /dev/null 2>&1; then 
  filename="$(ls ./$StaticObjectName.*.so)"; 
  echo "$filename is compiled" 
  ln -s $filename $StaticObjectName.so
else
  echo "compiling failed for unknown reasons; check $StaticObjectName.log and $StaticObjectName.err for details" 
  grep -i 'error' $StaticObjectName.log
  echo "the gfortran you are using is:"
  echo "gfortran_version=__GNUC__.__GNUC_MINOR__.__GNUC_PATCHLEVEL__" | $FortranCompiler -E -cpp - 2>/dev/null | tail -n 1 | sed -e 's/ //g'
  exit 125
fi
