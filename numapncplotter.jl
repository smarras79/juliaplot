using NetCDF
using Plots

filepath = "/Users/simone/Work/numa-git/GitLab_NPS/IMEXbranch/numa3d/myruns/AMR/case13/WORKING-AMR-TC-35th-conference/"
filename = "case13_AMR=3_0.025000_set2nc_cgd_ark2_schur_3d_p4est_0001.nc"


fn = joinpath(filepath, filename)
@info fn
ncinfo(fn)


nop     = ncread(fn,"nopx")
xcoords = ncread(fn,"x")
ycoords = ncread(fn,"y")
zcoords = ncread(fn,"z")
theta  = ncread(fn,"theta")

@info size(theta) size(x)
plot(x=xcoords, y=ycoords, z=theta)



