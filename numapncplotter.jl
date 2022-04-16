using NetCDF

import PyPlot as plt
using NumPyArrays
using PyCall

tri = pyimport("matplotlib.tri")
np  = pyimport("numpy")

filepath = "/Users/simone/Work/numa-git/GitLab_NPS/IMEXbranch/numa3d/myruns/amr/WORKING-AMR-TC-35th-conference"
#filepath = "/Users/simone/Work/numa-git/GitLab_NPS/IMEXbranch/numa3d/myruns/AMR/case13/WORKING-AMR-TC-35th-conference/"
filename = "case13_AMR=3_0.025000_set2nc_cgd_ark2_schur_3d_p4est_0001.nc"


fn = joinpath(filepath, filename)
@info fn
ncinfo(fn)


nop     = ncread(fn,"nopx")
xcoords = ncread(fn,"x")
ycoords = ncread(fn,"y")
zcoords = ncread(fn,"z")
theta  = ncread(fn,"theta")

@show xmin = minimum(xcoords)
@show xmax = maximum(xcoords)
@show ymin = minimum(ycoords)
@show ymax = maximum(ycoords)
@show zmin = minimum(zcoords)
@show zmax = maximum(zcoords)

triang       = tri.Triangulation(xcoords, ycoords)
interpolator = tri.LinearTriInterpolator(triang, theta)


# Create grid values first.
ngridx = 100
ngridy = 200
xi = np.linspace(xmin, xmax, ngridx)
yi = np.linspace(xmin, ymax, ngridy)
Xi, Yi = np.meshgrid(xi, yi)
zi = interpolator(Xi, Yi)
#
#@info size(theta) size(xcoords)

