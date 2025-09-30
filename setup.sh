#!/usr/bin/env bash
set -e

cd ~/athena ###TO START FROM HOME DIRECTORY


export CXXFLAGS="-I${CONDA_PREFIX}/include ${CXXFLAGS}"
export LDFLAGS="-L${CONDA_PREFIX}/lib ${LDFLAGS}"
python configure.py --prob=gr_star_blast -g --coord=schwarzschild \
    -hdf5 --hdf5_path="$CONDA_PREFIX"
make -j4
make clean 
make 


cd ~/work
install_name_tool -add_rpath "$CONDA_PREFIX/lib"  ~/athena/bin/athena
# to change the file
#cp athinput.{name} athinput.{name}_new


cp ~/athena/inputs/mhd/athinput.{name} .
install_name_tool -add_rpath "$CONDA_PREFIX/lib"  ~/athena/bin/athena  
~/athena/bin/athena -i athinput.bw


# for python plotting "gr shock tube"
 ~/athena/vis/python/plot_lines.py mhd_shock_rel_1.block0.out1.00001.tab x1v \
    Etot,dens show -c r,b -l '$E$,$D$' --x_min=-0.5 --x_max=0.5 --x_label '$x$'


#for hdf5 output
-hdf5 --hdf5_path="$CONDA_PREFIX"



install_name_tool -add_rpath "$CONDA_PREFIX/lib"  ~/athena/bin/athena  
install_name_tool -add_rpath "$CONDA_PREFIX/lib"  ~/athena/bin/athena  
install_name_tool -add_rpath "$CONDA_PREFIX/lib"  ~/athena/bin/athena  