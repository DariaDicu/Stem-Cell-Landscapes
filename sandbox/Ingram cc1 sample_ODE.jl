using Plots, DifferentialEquations, DataFrames, KernalDensity



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

# A = spawn_stratified(27, 3, 3, 3) # Call function
# plot(scatter(A[:,1], A[:,2])) # Plot all first col vs. all second col



TO DO: Go through ODE code and understand it. Build in sampling after that.



### ODEs

function run_simulation(F::Function, n::Int64, bounds)

    #spawn_stratified(n, 2, bounds[2], bounds[2])
    #x0 = A[]

  # Sample n random numbers from interval (0,1).
  sample = rand(n)

  # Set initial conditions by scaling the sampled numbers using the interval
  # "bounds".
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
