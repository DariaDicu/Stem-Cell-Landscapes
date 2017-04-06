# Script for finding and desctibing the minima in a matrix of heights

# Hunting for minima in a matrix
function hunt_minima(matrix_unorm)
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
         push!(minima_pos, (col, row, matrix_unorm[col,row]))
       end

       # Move on to the next element in matrix_unorm
     end
   end

  # Remove odd, redundant values as some minima are listed twice
  minima = unique(minima_pos)
  return minima
end


# Optional function to get print statements about the minima in the system
function minima_info_1(minima)
  num_minima = length(minima)

  # Create dictionary to name the different stationary points [only up to 26]
  dict = Dict(1 => "A", 2 => "B", 3 => "C", 4 => "D", 5=> "E",
              6 => "F", 7 => "G", 8 => "H", 9 => "I", 10  => "J",
              11 => "K", 12  => "L", 13 => "M", 14  => "N", 15 => "O",
              16 => "P", 17  => "Q", 18  => "R", 19  => "S", 20  => "T",
              21  => "U", 22  => "V", 23  => "W", 24  => "X", 25  => "Y",
              26 => "Z" )

  depths = []

  str1 = "There are "*string(num_minima)*" local minima in this system\n"
  str2 = "Located in [x,y] space at:\n"

  for i = 1:num_minima
    # Convert tuple to array for manipulation
    p = collect(minima[i][1:2])

    global p

    if isdefined(:x_linspaced) && isdefined(:y_linspaced)
      q = []
      push!(q, x_linspaced[p[1]])
      push!(q, y_linspaced[p[2]])
      p = q
    else
       #Convert matrix indexes to real values by considering lower and step
       p = lower + (p * step)
    end

    push!(depths, round(minima[i][3],1))

    str2 = str2*"  "*string(dict[i])*" ["*string(p[1])*", "*string(p[2])*"]\n"
  end

  ordering = []
  k = maximum(depths)
  for t = 1:num_minima
    deepest = indmin(depths)
    depths[deepest] = k+1
    push!(ordering, dict[deepest])
  end

  order = string()
  for i = 1:num_minima
    order *= string(ordering[i]," ")
    if i != num_minima
      order *= string("> ")
    end
  end

  str3 = "\nIn order of stability, minima given as: \n  "*string(order)*" "
  rtn_string = str1 * str2 * str3
  return rtn_string
end
