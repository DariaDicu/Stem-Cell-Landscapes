z_values = [4 4 4 4 4
  2 4 1 3 4
  5 3 2 3 4
  4 2 4 4 4
  1 2 4 4 4]

module LandscapeColouring

global color_count = 0
global dx = [1, -1, 0, 0]
global dy = [0, 0, 1, -1]

# TODO: pass less variables.

function df_search(i, j, n, m, z_color, z_values, visited)
  global dx, dy, color_count
  visited[i, j] = true
  for d = 1:4
    x = i + dx[d]
    y = j + dy[d]
    # Not interested in points beyond matrix bounds.
    if (x == 0 || y == 0 || x > n || y > m) continue end
    # Not interested in neighbouring point if the gradient is positive.
    if (z_values[x, y] > z_values[i, j]) continue end
    if (z_color[x, y] == 0 && visited[x, y] == false)
      df_search(x, y, n, m, z_color, z_values, visited)
    end
    # This may overwrite previous color in (i, j), but it doesn't matter since
    # for saddle point we pick at random.
    z_color[i, j] = z_color[x, y]
  end
  # If point was not colored after looking at all neighbours, the point has no
  # neighbours with non-positive gradient, so is a local minima.
  if (z_color[i, j] == 0)
    color_count += 1
    z_color[i, j] = color_count
  end
end

function color_landscape(z_values)
  global color_count
  color_count = 0
  n = size(z_values)[1]
  m = size(z_values)[2]
  visited = fill(false, n, m)
  z_color = fill(0, n, m)
  for i = 1:n
    for j = 1:m
      if (z_color[i, j] == 0)
        df_search(i, j, n, m, z_color, z_values, visited)
      end
    end
  end
  return z_color, color_count
end

end # module LandscapeColouring

using LandscapeColouring

z_values
z_color, color_count = LandscapeColouring.color_landscape(z_values)

println(z_color)
println(color_count)
