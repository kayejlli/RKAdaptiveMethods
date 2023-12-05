#!/bin/bash
# example: CompileInterface.sh ODEInterface.f90 ode libSolvers.a 
FortranCompiler="$(which gfortran)";
FortranFile=$1
StaticObjectName=$2
Library=$3

##############################################
# addtional arguments 
addtional="";
for user in "$@" 
do
    #echo "inputs - $i: $user";
    # i=$((i + 1));
    if [[ ${user:0:1} == "-" ]]; then
      addtional="${addtional} $user"; 
    fi
done
echo "addtional arguments:: ${addtional}"
##############################################

if [ -z "$Library" ]; then
  # if not library is used 
  f2py -c $FortranFile -m $StaticObjectName --opt='-fPIC -O3 -fmax-errors=5' --f90exec=$FortranCompiler > $StaticObjectName.log 2> $StaticObjectName.err;
else
  # if using the library 
  LibraryName="${Library/lib/}"
  LibraryName="${LibraryName/.a/}"
  if [[ $LibraryName == *"/"* ]]; then
    LibraryName=($(echo "$LibraryName" | tr "\/" "\n"))
    echo "-L${LibraryName[0]} -l${LibraryName[1]}"
    f2py ${addtional} -L${LibraryName[0]} -l${LibraryName[1]} -c $FortranFile -m $StaticObjectName --opt='-O3' --f90exec=$FortranCompiler > $StaticObjectName.log 2> $StaticObjectName.err;
    echo "f2py ${addtional} -L${LibraryName[0]} -l${LibraryName[1]} -c $FortranFile -m $StaticObjectName --opt='-O3' --f90exec=$FortranCompiler > $StaticObjectName.log 2> $StaticObjectName.err"
  else
    echo "-L. -l$LibraryName"
    f2py ${addtional} -L. -l$LibraryName -c $FortranFile -m $StaticObjectName --opt='-O3' --f90exec=$FortranCompiler > $StaticObjectName.log 2> $StaticObjectName.err;
  fi
fi

########################################################
# check if *.so is complied or not; output error message  
if ls ./$StaticObjectName.*.so 1> /dev/null 2>&1; then 
  filename="$(ls ./$StaticObjectName.*.so)"; 
  echo "$filename is compiled" 
  ln -s $filename $StaticObjectName.so
else
  echo "compiling failed for unknown reasons; check $StaticObjectName.log and $StaticObjectName.err for details" 
  grep -i 'error' $StaticObjectName.err --color=always -B 3
  echo "the gfortran you are using is:"
  echo "gfortran_version=__GNUC__.__GNUC_MINOR__.__GNUC_PATCHLEVEL__" | $FortranCompiler -E -cpp - 2>/dev/null | tail -n 1 | sed -e 's/ //g'
  exit 125
fi
