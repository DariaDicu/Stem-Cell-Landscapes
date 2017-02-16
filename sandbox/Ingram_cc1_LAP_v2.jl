using Plots, DifferentialEquations, DataFrames, KernelDensity



# Sample, ODEs, Function, Landscape, Plotting
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
    levels=70, legend=false, xlabel="Dim 1", ylabel="Dim 2",
    title="Simple ODE with Strat Sampling")



### 1. Turn contours into numbers - done by ldens

### 2. Select starting and end attractor points.

# To detect minima automatically, pass a 3by3 array over and
# check if neighbours are all higher.

function get_minima()
    # Ivan has some code for this...
end

# Need to generalise this based on above function
start_index = [114,30]
start_value = ldens[start_index[1], start_index[2]]
goal_index = [30,114]
goal_value = ldens[goal_index[1], goal_index[2]]



### 3. Follow a 'moving point' from the start point

# Create a marker to store a pair of trajectories over time!
# Store the indices, then can call them later

marker = Array{Int,1}[]
push!(marker, [start_index[1], start_index[2]]) # Initialise - generalise this later



### 4. Get weighting for a probabilistic movement based on
# ...distance (between end of traj and goal) and
# ...jump size (energy to pass through contours)

# Distance between any point and the goal - general
function get_distance(i,j)
    # Set up fixed point
    b = goal_index
    # Adjacent points to test, using i and j
    a = [(marker[end][1])+i, (marker[end][2])+j]
    # Calcaulte Manhattan distance - more efficient than Euclidian
    distance = abs(a[1]-b[1]) + abs(a[2]-b[2])
    return distance
end

# Diff between any point, compared with the last one in 'marker'
function get_jump(i,j)
    # Set up fixed relative point
    current = ldens[marker[end][1], marker[end][2]]
    # Adjacent points to test, using i/j
    a = [(marker[end][1])+i, (marker[end][2])+j]
    # Calculate jump
    jump = ldens[a[1], a[2]] - current
    return jump
end

# Consider each neighbouring point and give it a weighting
function get_weight()
    # Combine get_distance and get_jump for overall weight
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
    #println("store_d: ",store_d)
    #println("store_j: ",store_j)

    # Normalise values, for jump, by dividing by heighest
    max_j = maximum(store_j) # Define here so looping values aren't affected
    limit = length(store_j)
    for k in 1:limit
        store_j[k] /= max_j
    end
    #println("store_d_norm: ",store_d)
    #println("store_j_norm: ",store_j)

    # Calcualte a weight for each neighbour based on above
    weight = ones(limit)
    for x in 1:limit
        # Jump: make -ve values more significant
        # by reducing i in i^-store_j, we make jump differences
        # less weighted!
        weight[x] *= sqrt(2)^-store_j[x]
        #println("weight-pre-/2: ",weight)
        #weight[x] /= 6
        #println("weight-post-/2: ",weight)
        # Distance: if further than current point, divide weight
        if store_d[x] >= get_distance(0,0)
            weight[x] /= 4
        end
    end
    #println("weight-post-d: ",weight)

    # Normalise weight
    max_weight = maximum(weight)
    for y in 1:limit
        weight[y] /= max_weight
    end

    return weight
    println("weight_norm: ",weight)
    # Fine-tune above values, and work how to track indices
    # Do I only need to normalise once??
end

function move()
    # Using above weights, return a probability for each index
    # ...and implement a movement in either direction. Can then update
    # ...marker!!
    while marker[end] != goal_index # Stop when you reach goal
    #for i in 1:256
        weight = get_weight()
        total = sum(weight)
        r = rand()*total # Generate random no. in range
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

marker = Array{Int,1}[]
push!(marker, [start_index[1], start_index[2]])
move()


### 5. Run simulation until completion. Multiple times to get best route?
x = []
y = []
for i in 1:length(marker)
    x_index = marker[i][2] # x-axis = j = 2nd index
    y_index = marker[i][1]
    push!(x, dens1.x[x_index])
    push!(y, dens1.y[y_index])
end
plot!(x,y,linecolor=4)

savefig("C:/Users/11ing/Documents/GitHub/Stem-Cell-Landscapes/sandbox/Plots/LAP_v2.png")
