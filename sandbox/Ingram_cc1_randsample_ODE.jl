using Plots, DifferentialEquations, DataFrames, KernelDensity



### ODEs

# Random sampling carried out witin the ODE function
function run_simulation(F::Function, n::Int64, bounds)

    # Sample n random numbers from interval (0,1)
    sample = rand(n)

    # Set initial conditions by scaling sampled numbers via 'bounds' interval
    x0 = map(x->(bounds[1] + x*(bounds[2] - bounds[1])), sample)

    # Time span is hardcoded for now, but will be (0, Inf) once we figure out
    # when to stop the trajectory.
    tspan = (0.0,10.0)

    # Perform a single simulation by running the ODE solver.
    prob = ODEProblem(F,x0,tspan)
    sol = solve(prob)

    # Return pair of vectors (trajectory_values, time_values) extracted from the
    # solution of the ODE simulation.
    return (sol.t, sol.u)
end

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
    # Convert matrix to DataFrame
    return (convert(DataFrame,output))
    println(output)
end

F = function (t,x)
    a = 0.3
    n = 4
    S = 0.5
    k = b = 1
    F1 = (x1, x2) ->
      (a*(x1^n)/(S^n + x1^n) + b*S^n/(S^n + x2^n) - k*x1)
    F2 = (x1, x2) ->
      (a*(x2^n)/(S^n + x2^n) + b*S^n/(S^n + x1^n) - k*x2)
    return [F1(x[1], x[2]), F2(x[1], x[2])]
end

# Create the data based on above functions!
data = build_landscape(500, F, 2, (0,3))



# Function plots landscape with specified number of contours
# Based off Stumpf's Jupyter notebook
function plot_contours(N::Int64)

    X = convert(Array{Float64},deepcopy(data[2]));
    Y = convert(Array{Float64},deepcopy(data[3]));

    dens1 = kde((X,Y))
    dens2 = 1e-23*ones(size(dens1.density)) + dens1.density

    # Potential landscape with '-log'
    ldens = -log(dens2);
    ldens = ldens-maximum(ldens)

    gr()
    contour_plot = contour(dens1.x, dens1.y, dens1.density,
        levels=N, legend=false, xlabel="Dim 1", ylabel="Dim 2")
    plot(contour_plot)
end

# Plotting should be separate from data generation
# Function takes number of simulations you want to plot
function plot_tracks(N::Int64)

    # Based on 'data' above, initialise the plotting
    data_edit = data[data[end] .== 1,:]
    plot!(data_edit[2], data_edit[3])

    for i = 2:N
        # Extract dataframe for each simulation at a time
        # '.' is element-wise checking?
        data_edit = data[data[end] .== i,:]
        x = data_edit[2]
        y = data_edit[3]
        plot!(x,y)
    end
    # Weird plotting thing means I have to repear these lines
    plot!(x,y)
end

plot_contours(70)
plot_tracks(250)
