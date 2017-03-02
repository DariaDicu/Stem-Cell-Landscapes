using FileIO, KernelDensity, Plots, PyPlot
include("ode_simulator.jl")

# You may have to change this path depending on local settings.
file_text = readstring(open(string(homedir(), "/Documents/stem-cells/sandbox/"*
  "func_file2.jl")))
input_func = parse(file_text)
F = eval(input_func)

# Only plot pth percentile lowest heights (since interested in wells).
cutoff_percentile = 85
# Number of simulation runs.
runs = 10000
# Number of dimensions.
n = 4
bounds = (0, 400)
println("Starting simulations.")
data  = ODESimulator.build_landscape(runs, F, n, bounds)
println("Finished simulations.")

println("Starting animation rendering.")
pyplot(leg=false)

function getXYZdata(data, is_endpoint, dim1, dim2)
  global runs
  global n
  # Build heights either using entire trajectory or endpoints only.
  if (!is_endpoint)
    # Get the trajectory coordinates for plotting along the entire trajectory.
    # Note: indexing at (dim+1) since first column represents the time value.
    X = convert(Array{Float64},deepcopy(data[dim1+1]));
    Y = convert(Array{Float64},deepcopy(data[dim2+1]));
  else
    # Get the endpoint trajectory coordinates.
    X = Float64[]
    Y = Float64[]
    for i = 1:runs
      # Extract the rows in the DataFrame where the run index is i.
      # Note: indexing at (dim+1) since first column represents the time value.
      current_run = data[data[n+2].==i,:]
      push!(X, current_run[end, dim1+1])
      push!(Y, current_run[end, dim2+1])
    end
  end
  dens1 = kde((X, Y))
  dens2=1e-23*ones(size(dens1.density))+dens1.density

  ldens=-log(dens2)
  ldens=ldens-maximum(ldens)

  return dens1.x, dens1.y, ldens
end

# Function to get grids for surfplot from linspace X and Y.
function ndgrid{T}(v1::AbstractVector{T}, v2::AbstractVector{T})
   m, n = length(v1), length(v2)
   v1 = reshape(v1, m, 1)
   v2 = reshape(v2, 1, n)
   (repmat(v1, 1, n), repmat(v2, m, 1))
end

dimension_names = ["N", "O", "F", "G"]
for i = 1:n
  for j = 1:n
    if (i == j) break end
    file_name_prefix = "surf_"*dimension_names[i]*dimension_names[j]
    file_name_suffix = "_"*string(cutoff_percentile)*"perc.png"

    # Plot landscape using endpoints only.
    X, Y, Z = getXYZdata(data, true#=is_endpoint=#, i, j)
    # Cut-off at p-th percentile. Get the value of the pth ranked element in the
    # sorted heights array.
    Z_array = vec(Z)
    sort(Z_array)
    pth_index = max(1, Int(round(cutoff_percentile*length(Z_array)/100)))
    cutoff_value = Z_array[pth_index]
    # Replace values above cutoff with cutoff value.
    Z = map(z->(z > cutoff_value ? cutoff_value : z), Z)
    # Transform LinSpaces to grids.
    X, Y = ndgrid(X, Y)
    # Plot reverse plot (-Z instead of Z) to visualize better.
    p = Plots.plot(X, Y, -Z, st=[:surface],
      #xlabel=dimension_names[i], ylabel=dimension_names[j], zlabel="Height",
      zlims=(-cutoff_value, -minimum(Z)),
      title="Endpoints only, "*dimension_names[j]*" against "*
      dimension_names[i])
    xaxis!(dimension_names[i])
    yaxis!(dimension_names[j])
    Plots.savefig(file_name_prefix*"_endpoints"*file_name_suffix)

    # Plot with entire trajectory. Repeat same procedure.
    X, Y, Z = getXYZdata(data, false#=is_endpoint=#, i, j)
    Z_array = vec(Z)
    sort(Z_array)
    pth_index = max(1, Int(round(cutoff_percentile*length(Z_array)/100)))
    cutoff_value = Z_array[pth_index]
    Z = map(z->(z > cutoff_value ? cutoff_value : z), Z)
    X, Y = ndgrid(X, Y)
    p = Plots.plot(X, Y, -Z, st=[:surface],
      #xlabel=dimension_names[i], ylabel=dimension_names[j], zlabel="Height",
      zlims=(-cutoff_value, -minimum(Z)),
      title="Entire trajectory, "*dimension_names[j]*" against "*
      dimension_names[i])
    xaxis!(dimension_names[i])
    yaxis!(dimension_names[j])
    Plots.savefig(file_name_prefix*"_trajectory"*file_name_suffix)
  end
end
