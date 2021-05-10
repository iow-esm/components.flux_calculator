

FC=mpiifort

mkdir build
cd build

$FC -c ../src/flux_lib/constants/*.F90
$FC -c ../src/flux_lib/auxiliaries/*.F90
$FC -c ../src/flux_lib/mass/*.F90
$FC -c ../src/flux_lib/heat/*.F90
$FC -c ../src/flux_lib/momentum/*.F90
$FC -c ../src/flux_lib/*.F90
#$FC -shared -o flux_library.so ../src/flux_lib/*.F90

cd ..
