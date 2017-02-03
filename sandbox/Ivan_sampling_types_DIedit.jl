#Exploration of sampling types for initial particle position.
using Plots

#This first function spawns a set of coordinates randomly (as randomly as
#a computer can) in the gene expression space. For each simulation, a
#coordinate set is generated with a random number between 0 and max
#for each dimension specified in the n-dimensional space.
function spawn_rand(sims, dims, max)

  A=zeros(sims, dims)

  for sim=1:sims
    for dim=1:dims
      A[sim, dim] = rand()*max
    end
  end
  return A
end

# This second function aims to make the initial sampling of the gene
# expression space more even. It divides each axis into a number of
# segments (partitions) and then works through the grid
# systematically, generating a random coordinate set within the current
# sampling subsection.

function spawn_stratified(sims, dims, max, partitions)

    # Initialise matrix to store coordinate values.
    # Rows = simulations, Cols = dimensions in gene space.
    B = zeros(sims, dims)

    jump = max/partitions # Size of axis subsections
    alpha = zeros(dims) # Tracks jumps and partitions

    for i = 1:sims
        for j = 1:dims
          # Populate B's next row with rand numbers within next interval
          # Uses *jump, and depends on value of alpha's last index
          B[i,j] = alpha[j] + rand()*jump
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
    return B
end


#This final function generates points according to the Latin Hypercube
#Sampling framework. Here, again the expression space is segmented into
#a number of subsections specifiied by partitions. A cartesian
#coordinate set is then sampled randomly within the entire gene space.
#All segments with the same indices in all dimensions are then removed
#from the possible sampling space so they cannot be resampled. A new
#coordinate set is then generated within one of the available segments
#and the process repeated. If there are no more available segments, in
#the case where sims > partitions^dims, all segments
#are made available again and the sampling continues.
function spawn_latin(sims, dims, max, partitions)
  print("Have a great day. ")

  jump = max/partitions
  C = zeros(sims,dims)
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

        C[i,j] = rand()*jump + (z-1)*jump
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
  return C
end

A = spawn_rand(27, 3, 3)
B = spawn_stratified(27, 3, 3, 3)
C = spawn_latin(27, 3, 3, 3)

plot(scatter(A[:,1], A[:,2]),scatter(B[:,1], B[:,2]))
plot(scatter(A[:,1], A[:,2]),scatter(B[:,1], B[:,2]), scatter(C[:,1], C[:,2]))
scatter(A[:,1], A[:,2])
scatter(B[:,1], B[:,2])
