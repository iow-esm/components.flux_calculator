#!/bin/bash

debug=${1:-"release"}
rebuild=${2:-"fast"}

module load intel/19.1.0 
module load openmpi/intel.19/3.1.6
module load netcdf/intel.19/4.7.4 

FC=mpifort
AR=ar
IOW_ESM_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/../.."

# SET SYSTEM-SPECIFIC COMPILER OPTIONS AND PATHS
# include paths
IOW_ESM_NETCDF_INCLUDE="${NETCDF_INCLUDE_PATH}"
IOW_ESM_NETCDF_LIBRARY="${NETCDF_LIBRARY_PATH}"

OASIS3_LIB="${IOW_ESM_ROOT}/components/OASIS3-MCT/oasis3-mct/IOW_ESM_PRODUCTION"
INCLUDES="-I ${OASIS3_LIB}/build/lib/mct \
          -I${OASIS3_LIB}/build/lib/psmile.MPI1 \
          -I ${OASIS3_LIB}/build/lib/psmile.MPI1 \
          -I ${OASIS3_LIB}/build/lib/mct \
          ${OASIS3_LIB}/lib/libpsmile.MPI1.a \
          ${OASIS3_LIB}/lib/libmct.a \
          ${OASIS3_LIB}/lib/libmpeu.a \
          ${OASIS3_LIB}/lib/libscrip.a \
          -I${IOW_ESM_NETCDF_INCLUDE}"
LIBS="-lnetcdf -lnetcdff -L${IOW_ESM_NETCDF_LIBRARY}"
if [ $debug == "debug" ]; then
FFLAGS="-O0 -r8 -fp-model precise -xHost -DUSE_DOUBLE_PRECISION -g -traceback -check all"
else
FFLAGS="-O3 -r8 -no-prec-div -fp-model fast=2 -xHost -DUSE_DOUBLE_PRECISION"
fi

rm -r build
rm -r bin
mkdir build
mkdir bin
cd build

$FC -c $FFLAGS ../src/flux_lib/constants/*.F90
$FC -c $FFLAGS ../src/flux_lib/auxiliaries/*.F90
$FC -c $FFLAGS ../src/flux_lib/mass/*.F90
$FC -c $FFLAGS ../src/flux_lib/heat/*.F90
$FC -c $FFLAGS ../src/flux_lib/momentum/*.F90
$FC -c $FFLAGS ../src/flux_lib/radiation/*.F90
$FC -c $FFLAGS ../src/flux_lib/*.F90
rm flux_library.a
$AR rv flux_library.a *.o

#$FC -c $FFLAGS ../src/routine_hdlerr.F90 -DUSE_DOUBLE_PRECISION -I${IOW_ESM_NETCDF_INCLUDE} $LIBS
#$FC -c $FFLAGS ../src/function_sent.F90 -DUSE_DOUBLE_PRECISION -I${IOW_ESM_NETCDF_INCLUDE} $LIBS 
#$FC -c $FFLAGS ../src/read_grid.F90 -DUSE_DOUBLE_PRECISION -I${IOW_ESM_NETCDF_INCLUDE} $LIBS
#$FC -c $FFLAGS ../src/decomp_def.F90 -DUSE_DOUBLE_PRECISION -I${IOW_ESM_NETCDF_INCLUDE} $LIBS
#$FC -c $FFLAGS ../src/read_dimgrid.F90 -DUSE_DOUBLE_PRECISION -I${IOW_ESM_NETCDF_INCLUDE} $LIBS
$FC -c $FFLAGS ../src/flux_calculator_basic.F90 
$FC -c $FFLAGS ../src/flux_calculator_prepare.F90
$FC -c $FFLAGS ../src/flux_calculator_calculate.F90
$FC -c $FFLAGS ../src/flux_calculator_io.F90 -I${IOW_ESM_NETCDF_INCLUDE} $LIBS
$FC $FFLAGS -o ../bin/flux_calculator ../src/flux_calculator.F90 flux_calculator_basic.o flux_calculator_prepare.o flux_calculator_calculate.o flux_calculator_io.o flux_library.a $INCLUDES $LIBS -Wl,-rpath,${IOW_ESM_NETCDF_LIBRARY}

cd ..
