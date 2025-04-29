# gnuplot script for plotting the density profile
unset xtics
unset ytics
unset ztics
set hidden3d
set grid
set pm3d
splot 'velocity_profile.dat'