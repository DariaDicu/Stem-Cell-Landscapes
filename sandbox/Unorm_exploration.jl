# An investigation into how to generate landscapes through normal decomposition
# of functions, as detailed in 'Quasi-potential landscapes in complex multi-
# stable systems' (Zhou et al, 2012)

# Governing functions for quad-stable example
F_quad_stable = function (t,x)

  F1 = (x1, x2) ->
   -1 + 9*x1 -2(x1^3) + 9*x2 -2(x2^3)
  F2 = (x1, x2) ->
   1 - 11*x1 + 2(x1^3) + 11*x2 -2(x2^3)
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end

# Return Unorm for specified x and y coordinate, following given equation.
# In future, this equation should be automatically derived from the equations
# governing the system, instead of being hardcoded.
get_Unorm = function (x,y)
  Unorm = -5((x^2)+(y^2))+0.5((x^4)+y^4)+x*y+x
  return Unorm
end

# Specify the plotting bounds lower and higher, as well as the step between each
# point evaluated in both axes
lower = -3
higher = 3
step = 0.01

# Obtain the dimensions of the required matrix
mat_dims = Integer((higher-lower)/step)

# Initialise matrix with appropriate dimensions in which to store Unorm values
matrix_unorm = zeros(mat_dims, mat_dims)

# For all (x,y) pairs in the plotting region, calculate Unorm using get_Unorm()
# and store in the appropriate position in matrix_unorm
for col = 1:mat_dims
  for row = 1:mat_dims
    x = lower + step*(col-1)
    y = lower + step*(row-1)

    matrix_unorm[col,row] = get_Unorm(x,y)
  end
end

# Hunting for minima in matrix_unorm
# Initialise array to store positions of minima
minima_pos = []

# Get the dimensions of matrix_unorm
dims = size(matrix_unorm)

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
         if matrix_unorm[col,row] > matrix_unorm[col+i, row+j]
           # If any of the neighbours are smaller, set smallest to false
           smallest = false
         end
       end
     end

     # If the current element was indeed smaller than all of its neighbours,
     # then smallest still returns true. Next store the indexes for this
     # minimum along with the corresponding height for plotting later.
     if smallest == true
       push!(minima_pos, (col,row,matrix_unorm[col,row]))
     end

     #Move on to the next element in matrix_unorm
   end
 end

# Remove odd, redundant values as some minima are listed twice
minima = unique(minima_pos)

# Begin plotting
using PyPlot
# Use PyPlot.surf() to visualise the results
surf(matrix_unorm, cmap="winter")

# Add the minima as red dots at the appropriate height
for point = 1:length(minima)
  scatter3D(minima[point][1], minima[point][2], minima[point][3], s=200, color = "red")
end
