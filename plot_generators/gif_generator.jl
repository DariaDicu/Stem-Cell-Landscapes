using FileIO, KernelDensity, Plots, PyPlot
include("ode_simulator.jl")

# Get ODE function definition from corresponding file.
current_dir = dirname(Base.source_path())
file_text = readstring(open(string(current_dir, "/ode_file_II_arg.jl")))
input_func = parse(file_text)
F = eval(input_func)

# Only plot pth percentile lowest heights (since interested in wells).
cutoff_percentile = 45
# Number of simulation runs.
runs = 100000
# Number of dimensions.
n = 4
#n = 2
# Bounds for ODE initial condition sampling.
bounds = (0, 400)
#bounds = (0, 5)

# Function to transform linspaces to grids, as required by pyplot.
function ndgrid{T}(v1::AbstractVector{T}, v2::AbstractVector{T})
   m, n = length(v1), length(v2)
   v1 = reshape(v1, m, 1)
   v2 = reshape(v2, 1, n)
   (repmat(v1, 1, n), repmat(v2, m, 1))
end

# Function for getting landscape matrices X, Y, Z.
function getXYZdata(data, is_endpoint, dim1, dim2)
  global runs
  global n
  # Build heights either using entire trajectory or endpoints only.
  if (!is_endpoint)
    # Get the trajectory coordinates for plotting along the entire trajectory.
    # Note: indexing at (dim+1) since first column represents the time value.
    X = convert(Array{Float64},data[dim1+1]);
    Y = convert(Array{Float64},data[dim2+1]);
  else
    # Get the endpoint trajectory coordinates.
    X = fill(0.0, runs)
    Y = fill(0.0, runs)
    # Single pass through DataFrame in order to extract endpoints. Assumes
    # points for any single run are stored in order of sampling time.
    for data_entry in DataFrames.eachrow(data)
      # Get index of simulation #.
      i = Int(data_entry[n+2])
      X[i] = data_entry[dim1+1]
      Y[i] = data_entry[dim2+1]
    end
  end
  dens1 = kde((X, Y))
  dens2=1e-23*ones(size(dens1.density))+dens1.density

  ldens=-log(dens2)
  ldens=ldens-maximum(ldens)

  # Transform linspace to grid, as required by pyplot.
  X, Y = ndgrid(dens1.x, dens1.y)
  return X, Y, ldens
end

# Get the cutoff value corresponding to the pth lowest ranked element. Replace
# all heights above this value with the cutoff value and return new landscape.
function get_cutoff_value(Z)
  global cutoff_percentile
  Z_array = sort(vec(Z))
  pth_index = max(1, Int(round(cutoff_percentile*length(Z_array)/100)))
  cutoff_value = Z_array[pth_index]
  return cutoff_value
end

println("Starting simulations.")
surfaces = []
# Set max cutoff to a very unlikely minimum value.
max_cutoff = -1000000000
for II = 0:10
  data  = ODESimulator.build_landscape_parallel(runs, F(II), n, bounds)
  # Hardcode dimensions for NANO and GATA6.
  i = 1
  j = 4
  # Plot landscape using endpoints only.
  X, Y, Z = getXYZdata(data, false#=is_endpoint=#, i, j)
  # Cut-off at p-th percentile. Get the value of the pth ranked element in the
  # sorted heights array.
  cutoff_value = get_cutoff_value(Z)
  max_cutoff = max(max_cutoff, cutoff_value)
  surf = X, Y, Z
  push!(surfaces, surf)
end

# After getting the overall maximum cutoff for negative log density, filter all
# heights that are larger than the cutoff. Then store the new surfaces and also
# flip Z to -Z (upside down to view better).
filtered_surfaces = []
for surf in surfaces
  X, Y, Z = surf
  Z = map(z->((z > max_cutoff) ? max_cutoff : z), Z)
  push!(filtered_surfaces, (X, Y, -Z))
end

surfaces = filtered_surfaces

# Get an array of tuples representing xlims, ylims and zlims for each surface.
# Important since we want a fixed grid when plotting over time (as II changes).
xyz_limits = map(surfaces) do surf
  x, y, z = surf
  return ((minimum(x), maximum(x)),
    (minimum(y), maximum(y)),
    (minimum(z), maximum(z)))
end

# Get the overall limit to restrict axes.
x_lim, y_lim, z_lim = reduce(xyz_limits) do bounds_a, bounds_b
  xa, ya, za = bounds_a
  xb, yb, zb = bounds_b
  x = min(xa[1], xb[1]), max(xa[2], xb[2])
  y = min(ya[1], yb[1]), max(ya[2], yb[2])
  z = min(za[1], zb[1]), max(za[2], zb[2])
  x, y, z
end
println("Finished simulations.")

println("Starting animation rendering.")
pyplot(leg=false)

# Generate a figure for each dimension pair, separately for endpoints-only and
# entire trajectory.
dimension_names = ["N", "O", "F", "G"]
# Fixing to NANO and GATA6.
anim = @animate for surf in surfaces
  # Choosing dimensions for NANO & GATA6.
  i = 1
  j = 4
  X, Y, Z = surf
  p = Plots.plot(X, Y, Z, st=[:surface],
    xlabel=dimension_names[i], ylabel=dimension_names[j],
    xlims=x_lim, ylims=y_lim, zlims=z_lim,
    title="Entire trajectory, "*dimension_names[j]*" against "*
    dimension_names[i]*", cut-off at "*string(cutoff_percentile)*"%")
end

file_name_prefix = "output/surf"
file_name_suffix = "_"*string(cutoff_percentile)*"perc.gif"
file_name = file_name_prefix*"_trajectory"*file_name_suffix
Plots.gif(anim, string(file_name), fps=5)
