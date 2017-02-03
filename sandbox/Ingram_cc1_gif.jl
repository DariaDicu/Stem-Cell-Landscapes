using Plots, DifferentialEquations, DataFrames, KernelDensity, ImageMagick



### Sampling

# Stratified sampling = even distribution of points. It divides each axis
# into segments (partitions) and then populates each in succession with a
# random coordinate set.
function spawn_stratified(sims, dims, max, partitions)

    # Initialise matrix to store coordinate values.
    # Rows = simulations, Cols = dimensions in gene space.
    A = zeros(sims, dims)

    jump = max/partitions # Size of axis subsections
    alpha = zeros(dims) # Tracks jumps and partitions

    for i = 1:sims
        for j = 1:dims
          # Populate B's next row with rand numbers within next interval
          # Uses *jump, and depends on value of alpha's last index
          A[i,j] = alpha[j] + rand()*jump
        end
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
            for k in dims-2:-1:1
                if alpha[k+1] >= max
                    alpha[k+1] = 0
                    alpha[k] += jump
                end
            end
        end
        if alpha[1] >= max
            alpha = zeros(dims)
        end
    end
    return A
end



### ODEs

function run_simulation(F::Function, n::Int64, bounds)

    A = spawn_stratified(5000, n, bounds[2], 10)
    # Generate indicies for n rand no.s based on spawn_stratified
    r = rand(1:size(A,1), n)
    x0 = A[r] # Store numbers as starting values

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



### Plotting

# N = no. contours. Advised: ~50
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

# N = no. simulations. Must be < N from 'data'
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
    # Weird plotting thing means I have to repeat these lines
    plot!(x,y)
end

# N = no. simulations to plot timepoints for. Slow for >100
function strat_gif1(N::Int64)
    # My way of specifying last time point, as diff sims vary in no. time
    # intervals. Ideally want to control this for each sim, but can't make
    # while/if conditions work well.
    #tstep = size(data[data[end] .== 1,:], 1) - 3 # To cover min time points
    # For each time step
    for i in 1:10
        selection = []
        # For each simulation
        for j in 1:N
            data_edit = data[data[end] .== j,:] # All data of next sim
            row = data_edit[1:i,2:3] # Select traj of rows up to i
            push!(selection, row) # Add that row to 'selection'
        end
        # Calling plot_contours *within* this function makes the traj's
        # a consistent colour, but not if outside function
        plot_contours(75)
        for k in 1:length(selection)
            x = selection[k][:,1]
            y = selection[k][:,2]
            plot!(x,y)
        end
        # Specify full path of folder
        savefig("C:/Users/11ing/Documents/GitHub/Stem-Cell-Landscapes/sandbox/Plots/"*string(i)*".png")
    end
end

strat_gif1(100)


# TO DO:
# ...Generalise last time point
# ...Make @gif work
# ...Combine different sampling types into 1 file
