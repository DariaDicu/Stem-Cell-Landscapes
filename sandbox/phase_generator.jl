using DifferentialEquations, Plots
pyplot(leg=false, size=(500,300))
#x, y = meshgrid(linspace(0,2pi,50), linspace(0,2pi,50))


# Feed data into function for plottig the phase diagram, given height matrix
function phase_portrait_gen(data, period, r, draw_contour)

  # Initialise the arrays to store values in later
  xbases = []
  ybases = []
  xheads = []
  yheads = []


  # For all cells in the matrix of heights, excluding those at the edge
  for i = 2:size(data)[1]-1
    for j = 2:size(data)[2]-1

      if i%period == 0 && j%period == 0
        # Get the current cell height, and those of the cell below, and to the right
        curr = data[i,j]
        down = data[i+1,j]
        right = data[i,j+1]

        # Get the gradient between that cell and its neighbour in the x plane...
        x_grad = curr - right
        # ... and again in the x plane
        y_grad = down - curr

        # Push these (x,y) differences in the heads arrays
        push!(xheads, x_grad)
        push!(yheads, y_grad)

        # Also push the appropriate (x,y) coordinates for the base of the arrow
        if isdefined(:x_linspaced) && isdefined(:y_linspaced) #&& 0 == 1
          push!(xbases, x_linspaced[j])
          push!(ybases, y_linspaced[size(data)[1]-i])
        else
          push!(xbases, j)
          push!(ybases, size(data)[1]-i)
        end
      end
    end
  end

  # Scale for fixed arrow length
  for j = 1:length(ybases)

    current = xheads[j]^2 + yheads[j]^2
    scale_fac = (r^2)/current

    xheads[j] = xheads[j] * sqrt(scale_fac)
    yheads[j] = yheads[j] * sqrt(scale_fac)
  end

  # Now to readjust arrows so that middles are on point instead of the bases
  for i = 1:length(xbases)
    xbases[i] = xbases[i] - 0.5*xheads[i]
  end

  for i = 1:length(ybases)
    ybases[i] = ybases[i] - 0.5*yheads[i]
  end

  # Perform the plotting
  gr()
  p = quiver(xbases, ybases, quiver=(xheads, yheads))


  if draw_contour == true
    if isdefined(:x_linspaced) && isdefined(:y_linspaced)
      contour!(x_linspaced, y_linspaced, rotr90(dens))
    else
      contour!(rotr90(dens))
    end
  end
  p
end
