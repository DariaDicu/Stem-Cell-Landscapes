using DifferentialEquations, Plots, DataFrames, KernelDensity,ProgressMeter;

# Runs a simulation for an n-variable ODE specified using function F, where
# starting conditions can be in the interval "bounds" = (x,y) for each of the n
# variables x1, x2, ..., xn.
function ode_simulator(iteration,F::Function,dims, bounds,t_span)
 function run_simulation(F::Function, n::Int64, bounds,t_span)
  # Sample n random numbers from interval (0,1).
  sample = rand(n)
  # Set initial conditions by scaling the sampled numbers using the interval
  # "bounds".
  x0 = map(x->(bounds[1] + x*(bounds[2] - bounds[1])), sample)
  # Time span is hardcoded for now, but will be (0, Inf) once we figure out
  # when to stop the trajectory.
  tspan=t_span

  # Perform a single simulation by running the ODE solver.
  prob = ODEProblem(F,x0,tspan)
  sol = DifferentialEquations.solve(prob)
  # Return pair of vectors (trajectory_values, time_values) extracted from the
  # solution of the ODE simulation.
  return (sol.t, sol.u)
 end

# Helper function that takes a simulation output in the form of a time point
# array, a trajectory and the number of the run, and returns an array of n+2
# columns and one row for each trajectory point.
#
# Columns are: Time, x1, x2, ..., xn, Run.
#
# Used for producing a dataframe by concatenating the blocks from each run.
 function format_simulation_output(times, trajectories, run_index)
  block = times
  n = length(trajectories[1])
  # Transform the array of arrays into a 2D array in order to concatenate it
  # into the block to get a column for each of the n coordinates.
  for j = 1:n
    block = hcat(block, map(x->x[j], trajectories))
  end
  number_samples = length(times)
  block = hcat(block, fill(run_index, number_samples))
  return block
 end

# Runs N simulations for an n-variable ODE specified using function F, where
# the initial conditions for each of the variables x1, x2, ..., xn is in the
# interval "bounds" = (x, y).
#
# Builds a Dataframe of n+2 columns, where each row represents one of the
# trajectory coordinates from one of the runs.
#
# Columns are: Time, x1, x2, ..., xn, Run.
 function build_landscape(N::Int64, F::Function, n::Int64, bounds,t_span)
  # Run an initial simulation to initialize output matrix.
  times, trajectories = run_simulation(F, n, bounds,t_span)
  output = format_simulation_output(times, trajectories, 1)
  p = Progress(N, 1, "Computing initial pass...", 50)
  # Build output matrix by performing N simulations.
  for i = 2:N
    times, trajectories = run_simulation(F, n, bounds,t_span)
    block = format_simulation_output(times, trajectories, i)
    output = vcat(output, block)
    next!(p)
  end

  # Convert matrix to DataFrame.
  return (convert(DataFrame,output))
 end


module ODESimulator
  using DifferentialEquations, Plots, DataFrames

  # Functions for stopping simulation under convergence conditions (heuristics).
  function custom_tolerance_callback(tolerance::Float64)
    # Event is triggered when each element of the derivative at time t is smaller
    # than the specified tolerance.
    condition = function (t,u,integrator)
        du = integrator.f(t,u)
        all(du .< tolerance)
    end

    # The effect of the event is that the solver is interrupted so the trajectory
    # is terminated.
    affect! = function (integrator)
      terminate!(integrator)
    end

    # Boolean tuple for whether to save before and after the affect! (see
    # documentation). Not important, since stopping integration when event
    # occurs.
    # http://docs.juliadiffeq.org/latest/features/callback_functions.html
    save_positions = (true,true)

    DiscreteCallback(condition, affect!, save_positions)
  end

  # Runs a simulation for an n-variable ODE specified using function F, where
  # starting conditions can be in the interval "bounds" = (x,y) for each of the n
  # variables x1, x2, ..., xn.
  function run_simulation(F::Function, n::Int64, bounds)
    # Sample n random numbers from interval (0,1).
    sample = rand(n)
    # Set initial conditions by scaling the sampled numbers using the interval
    # "bounds".
    x0 = map(x->(bounds[1] + x*(bounds[2] - bounds[1])), sample)

    # Time span is hardcoded for now, but will be (0, Inf) once we figure out
    # when to stop the trajectory.
    tspan = (0.0,30.0)

    # Perform a single simulation by running the ODE solver.
    prob = ODEProblem(F,x0,tspan)

    # Uncomment this is you want to stop integration when du/dt reaches
    # tolerance.
    # sol = solve(prob, callback = custom_tolerance_callback(0.0000001))

    # Solve without stopping trajectory.
    sol = solve(prob)

    # Return pair of vectors (trajectory_values, time_values) extracted from the
    # solution of the ODE simulation.
    return (sol.t, sol.u)
  end

  # Helper function that takes a simulation output in the form of a time point
  # array, a trajectory and the number of the run, and returns an array of n+2
  # columns and one row for each trajectory point.
  #
  # Columns are: Time, x1, x2, ..., xn, Run.
  #
  # Used for producing a dataframe by concatenating the blocks from each run.
  function format_simulation_output(times, trajectories, run_index)
    block = times
    n = length(trajectories[1])
    # Transform the array of arrays into a 2D array in order to concatenate it
    # into the block to get a column for each of the n coordinates.
    for j = 1:n
      block = hcat(block, map(x->x[j], trajectories))
    end
    number_samples = length(times)
    block = hcat(block, fill(run_index, number_samples))
    return block
  end

  # Runs N simulations for an n-variable ODE specified using function F, where
  # the initial conditions for each of the variables x1, x2, ..., xn is in the
  # interval "bounds" = (x, y).
  #
  # Builds a Dataframe of n+2 columns, where each row represents one of the
  # trajectory coordinates from one of the runs.
  #
  # Columns are: Time, x1, x2, ..., xn, Run.
  function build_landscape(N::Int64, F::Function, n::Int64, bounds)
    # Run an initial simulation to initialize output matrix.
    times, trajectories = run_simulation(F, n, bounds)
    output = format_simulation_output(times, trajectories, 1)

    # Build output matrix by performing N simulations.
    for i = 2:N
      times, trajectories = run_simulation(F, n, bounds)
      block = format_simulation_output(times, trajectories, i)
      output = vcat(output, block)
    end

    # Convert matrix to DataFrame.
    return (convert(DataFrame,output))
  end

  export build_landscape

end # end of module ODESimulator


# Example of function for representing a 2 transcription-factor with self- and
# mutual- regulation. Parameters are hardcoded in the function. See Wang et al,
# 2011 (http://www.pnas.org/content/108/20/8257.full).
#=
 F = function (t,x)
  a = 0.3
=======
F = a -> function (t,x)
>>>>>>> daria-working-branch
  n = 4
  S = 0.5
  k = b = 1
  F1 = (x1, x2) ->
    (a*(x1^n)/(S^n + x1^n) + b*S^n/(S^n + x2^n) - k*x1)
  F2 = (x1, x2) ->
    (a*(x2^n)/(S^n + x2^n) + b*S^n/(S^n + x1^n) - k*x2)
  return [F1(x[1], x[2]), F2(x[1], x[2])]
 end
=#
# Code to plot a contour map for the cell-fate ODE represented by F. data = build_landscape(1000, F(0.3), 2, (0,5))

 data = build_landscape(parse(Int,iteration), F::Function, dims, bounds , t_span)


 X = convert(Array{Float64},deepcopy(data[2]));
 Y = convert(Array{Float64},deepcopy(data[3]));

 #println(Y)
 dens1 = kde((X, Y))
 dens2=1e-23*ones(size(dens1.density))+dens1.density

 ldens=-log(dens2);
 ldens=ldens-maximum(ldens)

 #print(dens1.x)
 #print(dens1.y)

 gr()
 contour_plot=contour(dens1.x,dens1.y,dens1.density,levels=100,
  legend=false,xlabel="Dim 1",ylabel="Dim 2")
 display(plot(contour_plot))
 return data
end

#ldens=-log(dens2);
#ldens=ldens-maximum(ldens)
