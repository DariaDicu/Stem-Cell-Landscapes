using FileIO, KernelDensity, Plots;

ARGS=String[]
string_array = ("a=0.3", "n=4", "S=0.5", "b=1", "k=1", "vol=0:1:5", "filename=vol.gif", "plottype=prob", "runs=30")
append!(ARGS, string_array)

file_text = readstring(open(string(homedir(), "/Documents/GitHub/Stem-Cell-Landscapes/sandbox/func_file_vol.jl")))
input_func = parse(file_text)

# Get the list of parameter names from the function definition.
# Used for comparing against command line arguments.
func_param_names = split(string(input_func.args[1]), r",|\(|\)", keep=false)

# Put the anonymous function definition in a variable to call it with params.
# Make it an array since we will have multiple partial applications of this and
# we want to call map on them (even when it's just one).
F = [eval(input_func)]

# Parse input arguments to create dictionary of the form
# ["parameter name" => "parameter value"]. Parameter values can be either Int64,
# Float64, Range or anything that is parsed and evaluated to a valid expression.
params = Dict{String, Any}() # Name-value dictionary for parameters.
arg_param_names = String[]
non_constant_param_count = 0
map(ARGS) do arg_string
  # E.g. a=(1:2:100) => ["a", "(1:2:100)"]
  arg = split(arg_string, "=")
  @assert (length(arg) == 2) "Error parsing argument "*arg_string*". Command "*
    "line arguments must be of form: parameter_name=value."

  # Extract parameter name and value from command line argument.
  param_name = arg[1]

  if (param_name == "filename" || param_name == "plottype")
    # Keep as string instead of parsing it.
    param_value = arg[2]
  else
    param_value = try
      eval(parse(arg[2]))
    catch
      # Throw error is the value of this parameter cannot be parsed.
      throw(error("Invalid Julia expression for parameter "*param_name*". Please"*
      " redefine this parameter as a valid Julia expression."))
    end

    if (typeof(param_value) != Float64 && typeof(param_value) != Int)
      global non_constant_param_count += 1
      println(param_name)
    end
  end

  # Add dictionary entry and add parameter name to the list.
  params[param_name] = param_value
  push!(arg_param_names, param_name)
end

@assert non_constant_param_count <= 1 "At most one parameter can be "*
  "non-constant. Please change one or more of the parameter values."

# Hack for adding partial application (applying one parameter at a time).
function partial(f,a...)
        ( (b...) -> f(a...,b...) )
end

println("Preparing ODE function for binding parameters.")

map(func_param_names) do param_name
  # If there is no value for the parameter in the dictionary, display error to
  # user; parameter must be specified by a command line argument.
  @assert haskey(params, param_name) "Please specify value for parameter "*
    param_name
  param_value = params[param_name]
  # For each value in param_value (can be a constant too).
  new_F = []
  map(param_value) do val
    # Partially apply val to each value f of the function F we have so far.
    # Append results into the new array of functions F.
    append!(new_F, map(F) do f
      partial(f, val)
    end)
  end
  global F = new_F
end

println("Finished preparing ODE function for binding parameters.")

println("Preparing and including simulation modules.")

include("ode_simulator.jl")

gr()
println("Finished loading simulation code.")

# Set number of runs by default to 100.
runs = 100
if (haskey(params, "runs"))
  runs = params["runs"]
else
  warn("Number of simulation runs not specified. Using runs=100 as default.")
end

# Set number of fps by default to 10.
fps_value = 10
if (haskey(params, "fps"))
  fps_value = params["fps"]
else
  warn("Number of frames per second for animation not specified. Using fps=10"*
    " as default.")
end

# Set filename to default, then overwrite if specified as argument.
file_name = "surface.gif"
if (haskey(params, "filename"))
  file_name = params["filename"]
end

plot_type = "log"
if (haskey(params, "plottype"))
  plot_type = params["plottype"]
  if (plot_type != "log" && plot_type != "prob")
    throw(error("Unrecognized plot type. You can specify the type using either"*
    " plottype=\"log\" or plottype=\"prob\"."))
  end
else
  warn("Plot type not specified. Defaulting to negative log density. You can "*
    "specify the type using either plottype=\"log\" or plottype=\"prob\".")
end

# Firstly obtain surfaces as Array of (x,y,z) 2D matrices. This way the loop can
# be parallelized and the X,Y grid extracted so that it is the same over the
# entire animation.
surfaces = []
# TODO: Make async!

for f in F
  ode_function = f()
  data  = ODESimulator.build_landscape(runs, ode_function, 2, (0,3))

  X = convert(Array{Float64},deepcopy(data[2]));
  Y = convert(Array{Float64},deepcopy(data[3]));

  dens1 = kde((X, Y))
  dens2=1e-23*ones(size(dens1.density))+dens1.density

  ldens=-log(dens2)
  ldens=ldens-maximum(ldens)

  plot_dens = (plot_type == "log") ? ldens : dens1.density
  push!(surfaces, (dens1.x, dens1.y, plot_dens))
end
println("Finished simulations for obtaining landscapes.")

# Get an array of tuples representing xlims, ylims and zlims for each surface.
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


println(x_lim)
println(y_lim)
println(z_lim)

# Function to get grids for surfplot from linspace X and Y.
function ndgrid{T}(v1::AbstractVector{T}, v2::AbstractVector{T})
   m, n = length(v1), length(v2)
   v1 = reshape(v1, m, 1)
   v2 = reshape(v2, 1, n)
   (repmat(v1, 1, n), repmat(v2, m, 1))
end

# Loop over all surfaces to create animation.
println("Starting animation rendering.")
anim = @animate for surf in surfaces
  X, Y, Z = surf
  println("rendering for ", X, " ", Y)
  #p = Plots.plot(LinSpace(X), LinSpace(Y), vec(Z), st = [:surface],
    #xlims=x_lim, ylims=y_lim, zlims=z_lim)
  p = Plots.plot(LinSpace(X), LinSpace(Y), vec(Z), st = [:surface],
    xlims=(0,3), ylims=(0,3), zlims=z_lim)
end
println("Finished animation rendering.")

println("Saving gif to file.")
Plots.gif(anim, string(file_name), fps=fps_value)
