#Exploration of sampling types for initial particle position.
using Plots

#This first function spawns a set of coordinates randomly (as randomly as
#a computer can) in the gene expression space. For each simulation, a
#coordinate set is generated with a random number between 0 and var_bound
#for each dimension specified in the n-dimensional space.
function spawn_rand(num_sims, num_dims, var_bound)

  A=zeros(num_sims, num_dims)

  for sim=1:num_sims
    for dim=1:num_dims
      A[sim, dim] = rand()*var_bound
    end
  end
  return A
end

#This second function aims to make the initial sampling of the gene
#expression space more even. It divides each axis into a number of
#segments (num_var_partitions) and then works through the grid
#systematically, generating a random coordinate set within the current
#sampling subsection.
function spawn_stratified(num_sims, num_dims, var_bound, num_var_partitions)

  #Initialise matrix to store coordinate values. Each row is a simulation,
  #and each column represents a dimension in the expression space.
  B=zeros(num_sims, num_dims)

  #The 'jump' is the difference in axis values between subsections
  jump = var_bound/num_var_partitions

  alpha = zeros(1:num_dims)

  for sim=1:num_sims
    for dim=1:num_dims
      B[sim, dim] = alpha[dim] + rand()*jump
    end
    if sim % num_var_partitions == 0
      alpha[num_dims] = 0 #Reset the last one
      alpha[num_dims-1] = alpha[num_dims-1] + jump
      for i = num_dims-2:-1:1
        if alpha[i+1] >= var_bound
          alpha[i+1] = 0
          alpha[i] = alpha[i] + jump
          print("Reset")
        end
      end
    else
      alpha[num_dims] = alpha[num_dims] + jump
    end
    if alpha[1] >=var_bound
      alpha = zeros(1:num_dims)
    end
  end
  return B
end

#This final function generates points according to the Latin Hypercube
#Sampling framework. Here, again the expression space is segmented into
#a number of subsections specifiied by num_var_partitions. A cartesian
#coordinate set is then sampled randomly within the entire gene space.
#All segments with the same indices in all dimensions are then removed
#from the possible sampling space so they cannot be resampled. A new
#coordinate set is then generated within one of the available segments
#and the process repeated. If there are no more available segments, in
#the case where num_sims > num_var_partitions^num_dims, all segments
#are made available again and the sampling continues.
function spawn_latin(num_sims, num_dims, var_bound, num_var_partitions)
  print("Have a great day. ")

  jump = var_bound/num_var_partitions
  C = zeros(num_sims, num_dims)
  avail_coords = zeros(num_var_partitions,num_dims)

  for i = 1:num_var_partitions
    avail_coords[i,:] = i
  end

  for sim = 1:num_sims
    for dim = 1:num_dims
        #Pick a new coordinate from the available list
        z = 0; pos = 0;
        while z == 0
          pos = rand(1:num_var_partitions)
          z = avail_coords[pos, dim]
        end

        C[sim, dim] = rand()*jump + (z-1)*jump
        avail_coords[pos,dim] = 0
    end
    if sum(avail_coords) == 0
      #Refill available coordinates
      avail_coords = zeros(num_var_partitions,num_dims)
      for i = 1:num_var_partitions
        avail_coords[i,:] = i
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
