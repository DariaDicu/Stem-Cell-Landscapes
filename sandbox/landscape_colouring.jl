module LandscapeColouring
  global color_count = 0
  global dx = [1, -1, 0, 0]
  global dy = [0, 0, 1, -1]
  global lookforminima = true

  # Return true if moving in the direction of the neighbour corresponds to
  # moving towards a 'well'.
  function moving_towards_well(source, neighbour)
    return lookforminima ? (source >= neighbour) : (source <= neighbour)
  end

  function df_search(i, j)
    global dx, dy, color_count
    visited[i, j] = true
    for d = 1:4
      x = i + dx[d]
      y = j + dy[d]
      # Not interested in points beyond matrix bounds.
      if (x == 0 || y == 0 || x > n || y > m) continue end
      # Not interested in neighbouring point if the gradient is positive.
      if (!moving_towards_well(z_values[i, j], z_values[x, y])) continue end
      if (z_color[x, y] == 0 && visited[x, y] == false)
        df_search(x, y)
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

  # Main function to be called on height matrix z_values.
  # lookforminima should be true when "wells" correspond with minima (e.g. when
  # plotting negative log density), false when they correspond to maxima.
  function color_landscape(_z_values, _lookforminima)
    global color_count = 0
    global z_values = _z_values
    global lookforminima = _lookforminima
    global n = size(z_values)[1]
    global m = size(z_values)[2]
    global visited = fill(false, n, m)
    global z_color = fill(0, n, m)
    for i = 1:n
      for j = 1:m
        if (z_color[i, j] == 0)
          df_search(i, j)
        end
      end
    end
    return z_color, color_count
  end

end # module LandscapeColouring
