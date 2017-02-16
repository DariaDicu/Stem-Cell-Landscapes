using Plots, DifferentialEquations, DataFrames, KernelDensity



### 1. Sample, ODEs, Function, Landscape, Plotting

function spawn_latin(sims, dims, max, partitions)

    jump = max/partitions
    A = zeros(sims,dims)
    avail_coords = zeros(partitions,dims)

    for i = 1:partitions
        avail_coords[i,:] = i
    end

    for i = 1:sims
        for j = 1:dims
            #Pick a new coordinate from the available list
            z = 0; pos = 0;
            while z == 0
                pos = rand(1:partitions)
                z = avail_coords[pos, j]
            end

            A[i,j] = rand()*jump + (z-1)*jump
            avail_coords[pos,j] = 0
        end
        if sum(avail_coords) == 0
            #Refill available coordinates
            avail_coords = zeros(partitions,dims)
            for k = 1:partitions
                avail_coords[k,:] = k
            end
        end
    end
    return A
end
function run_simulation(F::Function, n::Int64, bounds)

    A = spawn_latin(5000, n, bounds[2], 10)
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

data = build_landscape(500, F, 2, (0,3))

function plot_contours(N::Int64)
    # Convert from dataframe to array, without changing the DF data
    X = convert(Array{Float64},deepcopy(data[2]));
    Y = convert(Array{Float64},deepcopy(data[3]));

    # Kernel density estimation - 'magic'
    dens1 = kde((X,Y))
    # To prevent '0' values, add 'dt' to all values in dens1
    # dens1.density is an array of (256,256)
    dens2 = dens1.density + 1e-23*ones(size(dens1.density))

    ldens = -log(dens2) # Potential landscape with '-log'
    ldens -= maximum(ldens) # Normalise to make 'max' value = 0

    gr() # A visualisation framework for kde things
    contour_plot = contour(dens1.x, dens1.y, ldens,
        levels=N, legend=false, xlabel="Dim 1", ylabel="Dim 2",
        title="Stem Cell reprogramming - Least Action Path")
    plot(contour_plot)
    # ldens turns contours into numbers to 'work with'
end

plot_contours(70) # Suggested value for contours



### 2. Select starting and end attractor points

# Ivan's code - extracts *all minima*
# Use this to specify starting and end points
function get_minima(ldens)
    # Initialise array to store positions of minima
    minima_pos = []

    # Get the dimensions of ldens
    dims = size(ldens)

    # For each element in the matrix, exluding those around the edge...
    for row = 2:dims[1]-1
        for col = 2:dims[2]-1

            # Assume the current element is infact a minimum (smaller than all of its
            # neighbours) until we can prove otherwise.
            smallest = true

            # Check whether the current element is actually smaller than all of its
            # neighbouring cells (diagonals included)
            for i = -1:1:1
                for j = -1:1:1
                    if ldens[col,row] > ldens[col+i, row+j]
                        # If any of the neighbours are smaller, set smallest to false
                        smallest = false
                    end
                end
            end

            # If the current element was indeed smaller than all of its neighbours,
            # then smallest still returns true. Next store the indexes for this
            # minimum along with the corresponding height for plotting later.
            if smallest == true
                push!(minima_pos, [col, row, ldens[col,row]])
            end

            # Move on to the next element in ldens
        end
    end
    # Remove odd, redundant values as some minima are listed twice
    minima = unique(minima_pos)
    return minima
end
# get_minima returns multiple minima, so create a function to
# select those that are 'significant'
function significant_minima(matrix)
    min_list = get_minima(matrix)
    sig_min = []
    store = []
    for i in 1:length(min_list)
        push!(store, min_list[i][3])
    end
    global_min = minimum(store)
    for i in 1:length(min_list)
        # Arbitrary cut-off for significance
        if min_list[i][3] <= 0.95 * global_min
            push!(sig_min, min_list[i])
        end
    end
    return sig_min
end

min = significant_minima(ldens)

start_index = [Int64(min[1][1]), Int64(min[1][2])]
start_value = ldens[start_index[1], start_index[2]]
# Can set goal to anything, e.g. another local minimum
goal_index = [Int64(min[2][1]), Int64(min[2][2])]
goal_value = ldens[goal_index[1], goal_index[2]]

# Create a marker to store a pair of trajectories over time!
# Store the indices, then can call them later
marker = Array{Int,1}[]
push!(marker, [start_index[1], start_index[2]]) # Initialise - generalise this later



### 3. Get weighting for a probabilistic movement

# This is based on 'distance' (between end of traj and goal),
# 'jump size' (energy to pass through contours)

# Distance between the goal and any point relative to
# 'last point in marker', via (i,j) indices
function get_distance(i,j)
    # Set up fixed point
    b = goal_index
    # Adjacent points to test, using i and j
    a = [(marker[end][1])+i, (marker[end][2])+j]
    # Calcaulte Manhattan distance - more efficient than Euclidian
    distance = abs(a[1]-b[1]) + abs(a[2]-b[2])
    return distance
end

# Difference in jump value between the last one in 'marker',
# and any other given via (i,j) indices
function get_jump(i,j)
    # Set up fixed relative point
    current = ldens[marker[end][1], marker[end][2]]
    # Adjacent points to test, using i/j
    a = [(marker[end][1])+i, (marker[end][2])+j]
    # Calculate jump
    jump = ldens[a[1], a[2]] - current
    return jump
end

# Normalisation code to make later functions less clunky
function normalise(array)
    max = maximum(array)
    for i in 1:length(array)
        array[i] /= max
    end
    return array
end

# Consider each neighbouring point and give it a weighting
# based on get_distance and get_jump
function get_weight()
    store_d = []
    store_j = []
    for i in -1:1
        for j in -1:1
            if i==0 && j==0
                continue # Don't wanna consider current point? (Maybe I do!)
            else
                d = get_distance(i,j)
                j = get_jump(i,j)
                push!(store_d, d) # Store value for all 8 adjacent cells
                push!(store_j, j)
            end
        end
    end

    store_j = normalise(store_j) # Normalise jump values

    # Calcualte a weight for each neighbour based on above
    # Higher weight if: lower jump and/or closer to goal
    weight = ones(length(store_j))
    for x in 1:length(store_j)
        # Jump: make -ve values more significant
        # Reduce i in (i^-store_j) to make jumps less weighted
        weight[x] *= sqrt(2)^-store_j[x]
        # Distance: if further than current point, divide weight
        if store_d[x] >= get_distance(0,0)
            weight[x] /= 5 # Arbitrary value
        end
    end

    weight = normalise(weight) # Normalise weight
    return weight # Do I only need to normalise once? To test...
end



### 4. Move based on weighting, and plot path on the contours

# Old variant
function move()
    # Using above weights, return a probability for each index
    # ...and implement a movement in either direction. Can then update
    # ...marker!!
    while marker[end] != goal_index # Stop when you reach goal
    #for i in 1:256
        weight = get_weight()
        r = rand()*sum(weight) # Generate random no. in range
        # Make your move! Can we automate this? Hard!
        current = marker[end]
        if r < weight[1]
            push!(marker, [current[1]-1, current[2]-1]) # -1,-1
        elseif r < sum(weight[1:2])
            push!(marker, [current[1]-1, current[2]]) # -1,0
        elseif r < sum(weight[1:3])
            push!(marker, [current[1]-1, current[2]+1]) # -1,+1
        elseif r < sum(weight[1:4])
            push!(marker, [current[1], current[2]-1]) # 0,-1
        elseif r < sum(weight[1:5])
            push!(marker, [current[1], current[2]+1]) # 0,+1
        elseif r < sum(weight[1:6])
            push!(marker, [current[1]+1, current[2]-1]) # +1,-1
        elseif r < sum(weight[1:7])
            push!(marker, [current[1]+1, current[2]]) # +1,0
        else
            push!(marker, [current[1]+1, current[2]+1]) # +1,+1
        end
        #println("weight: ",weight)
    end
end

# New variant - run sim until end goal reached
function move()
    # Using above weights, return a probability for each index
    # ...and implement a movement in either direction. Can then update
    # ...marker!!
    while marker[end] != goal_index # Stop when you reach goal
        weight = get_weight()

        escape = false # Control break-out of loops
        r = rand()*sum(weight) # Generate random no. in range
        current = marker[end] # Set reference for adjacent cells
        count = 1
        # Iterate over adjacent cells and evaluate weight
        for i in -1:1
            for j in -1:1
                if i==0 && j==0
                    continue # Ignore current position - must change!
                elseif r < sum(weight[1:count])
                    push!(marker, [current[1]+i, current[2]+j])
                    escape = true
                    break # Break the inner for loop
                else
                    count += 1
                end
            end
        if escape == true
            break # Break the outer for loop
        end
        end
    end
end

# Quick refreshing of values
marker = Array{Int,1}[]
push!(marker, [start_index[1], start_index[2]])

move() # Call this. Can do multiple times for heat map?

# Plot path!
x = []
y = []
for i in 1:length(marker)
    x_index = marker[i][2] # x-axis = j = 2nd index
    y_index = marker[i][1]
    push!(x, dens1.x[x_index])
    push!(y, dens1.y[y_index])
end
plot!(x,y,linecolor=7)

savefig("C:/Users/11ing/Documents/GitHub/Stem-Cell-Landscapes/sandbox/Plots/LAP_v3.png")
