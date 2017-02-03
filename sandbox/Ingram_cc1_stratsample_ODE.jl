using Plots, DifferentialEquations, DataFrames, KernelDensity



### Sampling

# Stratified sampling = even distribution of points. It divides each axis
# into segments (partitions) and then populates each in succession with a
# random coordinate set.
function spawn_stratified(sims, dims, max, partitions)

    # Initialise a matrix to store coordinate values
    # Rows = simulations, Cols = dimensions in gene space
    A = zeros(sims, dims)

    jump = max/partitions # Size of axis subsections
    alpha = zeros(dims) # Tracks jumps and partitions

    for i = 1:sims
        for j = 1:dims
            # Populate A's next row with rand numbers within next interval
            # Uses *jump, and depends on value of alpha's last index
            A[i,j] = alpha[j] + rand()*jump
        end
        # Account for any multiple of partition value
        if i % partitions != 0
            # Add jump value to alpha's last position
            # Incremented value only affects last dimension
            alpha[dims] += jump
        else
            # Once you've filled all partitions, need to reset current index
            # and add jump to previous index.
            # Reset last index, add jump to previous
            alpha[dims] = 0
            alpha[dims-1] += jump

            # Check to see if any other indices are > max
            # If so, set to 0 and add 'jump' to previous index
            for k in dims:-1:1
                if alpha[k] == max
                    alpha[k] = 0
                    # If not first index -> avoids (k-1) clash
                    if k != 1
                        alpha[k-1] += jump
                    else
                        alpha = zeros(dims) # Reset!
                    end
                end
            end
        end
    end
    return A
end



### ODEs

function run_simulation(F::Function, n::Int64, bounds)

    A = spawn_stratified(100, n, bounds[2], 10)
    # Create n random numbers based on spawn_stratified
    r = rand(1:size(A,1), n)
    x0 = Float64[]
    for i in r
        push!(x0, A[i])
    end

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
        levels=N, legend=false, xlabel="Dim 1", ylabel="Dim 2",
        title="Simple ODE with Strat Sampling")
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
