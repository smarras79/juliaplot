ENV["PYTHON"]="python"
using NetCDF
using DataFrames
using ScatteredInterpolation
using Plots; gr()
default(legend = false)

using NumPyArrays
using Revise
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

function unique_inds(v)
    unq = Set{eltype(v)}()
    inds = Int[]
    for (i, x) in pairs(v)
       if x ∉ unq
           push!(inds, i)
           push!(unq, x)
       end
   end
   inds
end

function convex_hull(v1::Vector{Float64}, v2::Vector{Float64}; Plotter=nothing,n=10,raster=10)

    triin=Triangulate.TriangulateIO()    

    coords = vcat(v1', v2')
    #coords = hcat(unique([ Cdouble[rand(1:raster)/raster, rand(1:raster)/raster] for i in 1:n])...)
    triin.pointlist = coords
    #display(triin)
    (triout, vorout)=triangulate("Q", triin)
    #display(triout)
    plot_in_out(Plotter,triin,triout,title="Convex hull")

    #@test numberofpoints(triin)==numberofpoints(triout)
    #@test numberoftriangles(triout)>0

    return triout, vorout
end

function main()
    #tri = pyimport("matplotlib.tri")
    #np  = pyimport("numpy")

    filepath = "/Users/simone/Work/numa-git/GitLab_NPS/IMEXbranch/numa3d/myruns/amr/WORKING-AMR-TC-35th-conference/"
    filename = "case13_AMR=3_0.050000_set2nc_cgd_ark2_schur_3d_p4est_0005.nc"

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
    ##zslice_idx = myfindall(z -> z == zmin, zcoords);
    zslice_idx = findall(z -> z == zmin, zcoords)
    #
    # extract x and y subarrays at zslice_idx
    #
    xslice  = xcoords[zslice_idx]
    yslice  = ycoords[zslice_idx]
    thslice = theta[zslice_idx]

    #xslice::Array{Float64} = [0,  10,  20,  20,  30,  40,  50]
    #yslice::Array{Float64} = [0, 0.1, 0.1, 0.1, 0.3, 0.4, 0.5]
    #thslice::Array{Float64}= [0.1 0.2 0.3 0.4 0.5 0.6 0.7]
    
    coords           = vcat(xslice', yslice')
    triin            = Triangulate.TriangulateIO() 
    triin.pointlist  = coords
    (triout, vorout) = triangulate("Q", triin)
    
    #
    # Remove repeated nodes from coords and adjust thslice accordingly:
    #
    repidx    = unique_inds(eachrow(coords')) #OK
    thslice = copy(thslice[repidx]) #OK
    coords  = copy(coords[:,repidx]) #OK
        
    @info xslice = copy(coords[1,:])
    @info yslice = copy(coords[2,:])
    
    nx = 50
    ny = 50
    
    xi = range(xmin, stop = xmax, length = nx)
    yi = range(ymin, stop = ymax, length = ny)
    
    #f(xi, yi) = begin
    #    (3xi + yi ^ 2) * abs(sin(xi) + cos(yi))
    #end
    #Xi = repeat(reshape(xi, 1, :), length(yi), 1)
    #Yi = repeat(yi, 1, length(xi))
    #Zi = map(f, Xi, Yi)
          
    X = repeat(xi,  nx)[:]
    Y = repeat(yi', ny)[:]
    XX = repeat(reshape(xi, 1, :), length(yi), 1)
    YY = repeat(yi, 1, length(xi))
    
    gridPoints = [X Y]'
    samples = copy(thslice) #Ok  size(samples) -> (nsample,)
    points  = copy(coords)  #OK  size(points)  -> (2, nsample)
    itp = interpolate(NearestNeighbor(), points, samples);
    interpolated = evaluate(itp, gridPoints)
    griddedi = reshape(interpolated, nx, ny)
    #@info typeof(griddedi)
    #@info size(griddedi)
    #@info size(xi)
    #@info size(yi)
    
    #@info typeof(Zi)
    #@info size(Zi)
    #@info size(xi)
    #@info size(yi)

    p1 = contour(xi, yi, griddedi, fill = true)
    p2 = contour(xi, yi, griddedi)
    plot(p1, p2)
    
    
    #=
    #plot_in_out(PyPlot,triin,triout,title="Convex hull")
    #triang = Triangulate(xslice, yslice)
    ##triang       = tri.Triangulation(xslice, yslice)
    ##interpolator = tri.LinearTriInterpolator(triang, theta)

    #
    # Create grid values first.
    #
    ngridx = 100
    ngridy = 100
    xi = np.linspace(xmin, xmax, ngridx)
    yi = np.linspace(xmin, ymax, ngridy)
    Xi, Yi = np.meshgrid(xi, yi)

    
    # Create a uniform grid to interpolate onto
    nx = 50
    ny = 75
    #xgrid = collect(linspace(xmin,xmax,nx))
    #ygrid = collect(linspace(ymin,ymax,ny))
    grid_x = kron(ones(ny),xgrid')
    grid_y = kron(ygrid,ones(1,nx))
    # Perform the interpolation
    grid_val = si.griddata(points,val,(grid_x,grid_y),method="cubic")

    
    ip = 1
    coordsi = Array{Float64,2}(undef, 2, ngridx*ngridy)
    thetai = similar(coordsi)
    
    for i in 1:length(xi)
        for j in 1:length(yi)
            coordsi[1, ip] = xi[i]
            coordsi[2, ip] = yi[j]
            ip = ip + 1
        end
    end
    #    @show size(coordsi)

    #Triangulate
    #tri   = Triangulate()
    #sGrid = SimplexGrid(xslice, yslice)
    #@info length(coordsi[1,:])
    #for ip in 1:length(coordsi[1,:])
    #    x = [coordsi[1,ip], coordsi[2,ip]]
    #    @show thetai[ip] = interpolate(sGrid, thslice, x)
    #end
    
    
    #interpolator = Triangulate.LinearTriInterpolator(triout.pointlist, thslice)
    #zi = interpolator(Xi, Yi)
    #

    #interpolate(triout.pointlist, thslice, coords)
    #@info size(theta) size(xcoords)
    =#
    
end

function mytest()
    x = 1:0.5:20
    y = 1:0.5:10
    f(x, y) = begin
        (3x + y ^ 2) * abs(sin(x) + cos(y))
    end
    X = repeat(reshape(x, 1, :), length(y), 1)
    Y = repeat(y, 1, length(x))
    Z = map(f, X, Y)
    @info typeof(Z)
    p1 = contour(x, y, f, fill = true)
    p2 = contour(x, y, Z)
    plot(p1, p2)
end

#mytest();
main();



