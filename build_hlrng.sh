#!/bin/bash

debug=${1:-"release"}
rebuild=${2:-"fast"}

module unload intel
module unload impi
module load intel/18.0.5
module load impi/2018.5

FC=mpiifort
AR=ar
IOW_ESM_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/../.."

# SET SYSTEM-SPECIFIC COMPILER OPTIONS AND PATHS
# include paths
IOW_ESM_NETCDF_INCLUDE="/sw/dataformats/netcdf/intel.18/4.7.3/skl/include"
IOW_ESM_NETCDF_LIBRARY="/sw/dataformats/netcdf/intel.18/4.7.3/skl/lib"

if [ $debug == "debug" ]; then
	FFLAGS="-O0 -r8 -fp-model precise -xHost -DUSE_DOUBLE_PRECISION -g -traceback -check all -DIOW_ESM_DEBUG"
	configuration="DEBUG"
else
	FFLAGS="-O3 -r8 -no-prec-div -fp-model fast=2 -xHost -DUSE_DOUBLE_PRECISION"
	configuration="PRODUCTION"
fi

OASIS3_LIB="${IOW_ESM_ROOT}/components/OASIS3-MCT/oasis3-mct/IOW_ESM_${configuration}"
build_dir="build_${configuration}"
bin_dir="bin_${configuration}"

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

if [ "$rebuild" == "rebuild" ]; then
	rm -r "${build_dir}"
	rm -r "${bin_dir}"
fi
mkdir -p ./"${build_dir}"
mkdir -p ./"${bin_dir}"

cd "${build_dir}"

$FC -c $FFLAGS ../src/flux_lib/constants/*.F90
$FC -c $FFLAGS ../src/flux_lib/auxiliaries/*.F90
$FC -c $FFLAGS ../src/flux_lib/mass/*.F90
$FC -c $FFLAGS ../src/flux_lib/heat/*.F90
$FC -c $FFLAGS ../src/flux_lib/momentum/*.F90
$FC -c $FFLAGS ../src/flux_lib/radiation/*.F90
$FC -c $FFLAGS ../src/flux_lib/*.F90
rm flux_library.a
$AR rv flux_library.a *.o

$FC -c $FFLAGS ../src/flux_calculator_basic.F90 
$FC -c $FFLAGS ../src/flux_calculator_prepare.F90
$FC -c $FFLAGS ../src/flux_calculator_calculate.F90
$FC -c $FFLAGS ../src/flux_calculator_parse_arg.F90
$FC -c $FFLAGS ../src/flux_calculator_io.F90 -I${IOW_ESM_NETCDF_INCLUDE} $LIBS
$FC -c $FFLAGS ../src/flux_calculator_create_namcouple.F90
$FC $FFLAGS -o ../"${bin_dir}"/flux_calculator ../src/flux_calculator.F90 flux_calculator_basic.o flux_calculator_prepare.o flux_calculator_calculate.o flux_calculator_io.o flux_calculator_parse_arg.o flux_calculator_create_namcouple.o flux_library.a $INCLUDES $LIBS -Wl,-rpath,${IOW_ESM_NETCDF_LIBRARY}

cd -
