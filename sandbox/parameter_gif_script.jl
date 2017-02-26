using FileIO;

file_text = readstring(open(string(homedir(), "/Documents/stem-cells/sandbox/func_file.jl")))
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

println(arg_param_names)
println(params)

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

using KernelDensity, Plots, PyPlot;
include("ode_simulator.jl")

pyplot(leg=false, ticks=nothing)
#x = LinSpace(dens1.x)
#y = LinSpace(dens1.y)
# create a plot with 3 subplots and a custom layout


# Set number of runs by default to 100.
runs = 100
if (haskey(params, "runs"))
  runs = params["runs"]
else
  warn("Number of simulation runs not specified. Using runs=100 as default.")
end

# Set filename to default, then overwrite if specified as argument.
file_name = "surf.gif"
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


cnt = 0;
anim = @animate for f in F
  ode_function = f()
  data  = ODESimulator.build_landscape(runs, ode_function, 2, (0,5))

  X = convert(Array{Float64},deepcopy(data[2]));
  Y = convert(Array{Float64},deepcopy(data[3]));

  dens1 = kde((X, Y))
  dens2=1e-23*ones(size(dens1.density))+dens1.density

  ldens=-log(dens2);
  ldens=ldens-maximum(ldens)

  plot_dens = (plot_type == "log") ? ldens : dens1.density
  global cnt += 1
  p = Plots.plot(plot_dens, st = [:surface])
end

Plots.gif(anim, string(file_name))


#=
ODE_function = map(F, ())
println(ODE_function)
map(ODE_function, 10.0, [1,2])=#
