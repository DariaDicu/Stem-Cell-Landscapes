module ODESimulator

using Compat
using DifferentialEquations, DataFrames

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
function run_simulation(F::Function, n::Int64, bounds, tspan)
  # Sample n random numbers from interval (0,1).
  sample = rand(n)
  # Set initial conditions by scaling the sampled numbers using the interval
  # "bounds".
  x0 = map(x->(bounds[1] + x*(bounds[2] - bounds[1])), sample)

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
function build_landscape(runs::Int64, F::Function, n::Int64, bounds, tspan)
  # Run an initial simulation to initialize output matrix.
  times, trajectories = run_simulation(F, n, bounds, tspan)
  output = format_simulation_output(times, trajectories, 1)

  # Build output matrix by performing N simulations.
  for i = 2:runs
    times, trajectories = run_simulation(F, n, bounds)
    block = format_simulation_output(times, trajectories, i)
    output = vcat(output, block)
  end

  # Convert matrix to DataFrame.
  return (convert(DataFrame,output))
end

# Analogous to build_landscape, but spawns each simulation in parallel on
# available workers. This should be used, as it runs faster and uses less
# memory.
function build_landscape_parallel(runs::Int64, F::Function, n::Int64, bounds,
    tspan)
  results = pmap(1:runs) do i
    times, trajectories = run_simulation(F, n, bounds, tspan)
    format_simulation_output(times, trajectories, i)
  end

  output = vcat(results...)

  # Convert matrix to DataFrame.
  return (convert(DataFrame,output))
end

end # end of module ODESimulator
#=
# Example of function for representing a 2 transcription-factor with self- and
# mutual- regulation. Parameters are hardcoded in the function. See Wang et al,
# 2011 (http://www.pnas.org/content/108/20/8257.full).
F = a -> function (t,x)
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
