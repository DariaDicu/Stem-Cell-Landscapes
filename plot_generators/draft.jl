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
  return X, Y
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


data  = ODESimulator.build_landscape_parallel(runs, F(0), n, bounds)
X, Y = getXYZdata(data, false#=is_endpoint=#, 1, 4)
println(X)
println(Y)
