using NetCDF

#import PyPlot as plt
using PyPlot
using NumPyArrays
using PyCall
using Triangulate
using Test

function myfindall(condition, x)
    results = Int[]

    for i in 1:length(x)
        if condition(x[i])
            push!(results, i)
        end
    end
    return results
end

function example_convex_hull(v1, v2;Plotter=nothing,n=10,raster=10)
    triin=Triangulate.TriangulateIO()
    @info size(v1)
    @info size(v2)
    @info size(Cdouble[rand(1:raster)/raster, rand(1:raster)/raster])
    #triin.pointlist=hcat(unique([ Cdouble[rand(1:raster)/raster, rand(1:raster)/raster] for i in 1:n])...)
    
    triin.pointlist=[v1; v2]
    @info size(triin.pointlist) #QUIIII now modify this to use my XSLICE ad YSLICE
    display(triin)
    (triout, vorout)=triangulate("Q", triin)
    display(triout)
    plot_in_out(Plotter,triin,triout,title="Convex hull")
    @test numberofpoints(triin)==numberofpoints(triout)
    @test numberoftriangles(triout)>0
end

function main()
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
    theta   = ncread(fn,"theta")

    @show xmin = minimum(xcoords)
    @show xmax = maximum(xcoords)
    @show ymin = minimum(ycoords)
    @show ymax = maximum(ycoords)
    @show zmin = minimum(zcoords)
    @show zmax = maximum(zcoords)

    #
    # Extract indexes from zcoords where zcoords = zmin
    #
    #zslice_idx = myfindall(z -> z == zmin, zcoords);
    zslice_idx = findall(z -> z == zmin, zcoords)
    #
    # extract x and y subarrays at zslice_idx
    #
    xslice = xcoords[zslice_idx]
    yslice = ycoords[zslice_idx]

    v1 = [1, 2, 3]
    v2 = [4, 5, 6]
    example_convex_hull(v1, v2; Plotter=PyPlot,n=10,raster=10);
    
    #triang = Triangulate(xslice, yslice)
    #triang       = tri.Triangulation(xslice, yslice)
    #interpolator = tri.LinearTriInterpolator(triang, theta)
    
    #=
   



# Create grid values first.
ngridx = 100
ngridy = 200
xi = np.linspace(xmin, xmax, ngridx)
yi = np.linspace(xmin, ymax, ngridy)
Xi, Yi = np.meshgrid(xi, yi)
zi = interpolator(Xi, Yi)
#
#@info size(theta) size(xcoords)

=#
    
end

main();

