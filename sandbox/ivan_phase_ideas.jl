using DifferentialEquations, Plots
pyplot(leg=false, size=(500,300))

x, y = meshgrid(linspace(0,2pi,50), linspace(0,2pi,50))

#Function for plotting phase diagram with random arrow directions
function gen_phase(xlim, ylim, r)
  # Initialise empty arrays. Firstly to store (x,y) coordinates of the arrow
  # bases, and then the relative offset of each arrow head.
  xbases = []
  ybases = []

  xheads = []
  yheads = []


  for x = 1:xlim
    for y = 1:ylim
      push!(xbases, x)
      push!(ybases, y)
    end
  end

  for i = 1:length(xbases)
    push!(xheads, (rand()*2)-1)
  end

  for j = 1:length(ybases)
    push!(yheads, (rand()*2)-1)

    # Scale for fixed arrow length
    current = xheads[j]^2 + yheads[j]^2
    scale_fac = (r^2)/current

    xheads[j] = xheads[j] * sqrt(scale_fac)
    yheads[j] = yheads[j] * sqrt(scale_fac)
  end

  # Now to readjust arrows so that middles are on point '#onfleek'
  for i = 1:length(xbases)
    xbases[i] = xbases[i] - 0.5*xheads[i]
  end

  for i = 1:length(ybases)
    ybases[i] =ybases[i] - 0.5*yheads[i]
  end

  gr()
  quiver(xbases, ybases, quiver=(xheads, yheads))
end

gen_phase(10,10,0.5)


# Now its time to look up
using Interpolations

# Generate some random data, in matrix form
row_num = 10
col_num = 20
data = rand!(zeros(row_num, col_num))

for i = 1:row_num
  for j = 1:col_num
    data[i,j] = i * j
  end
end

# Attempt to create interpolation objects for more accurate gradients
A = interpolate(data, Gridded(Linear()))

# Feed data into function for plottig the phase diagram, given height matrix
function better_gen_phase(data, r)

  # Initialise the arrays
  xbases = []
  ybases = []
  xheads = []
  yheads = []

  # Calculate and push the x, y differences in the heads arrays
  for i = 2:size(data)[1]-1
    for j = 2:size(data)[2]-1
      curr = data[i,j]
      down = data[i+1,j]
      right = data[i,j+1]

      x_grad = curr - right
      y_grad = down - curr

      push!(xbases, j)
      push!(ybases, size(data)[1]-i)

      push!(xheads, x_grad)
      push!(yheads, y_grad)
    end
  end

  # Scale for fixed arrow length
  for j = 1:length(ybases)

    current = xheads[j]^2 + yheads[j]^2
    scale_fac = (r^2)/current

    xheads[j] = xheads[j] * sqrt(scale_fac)
    yheads[j] = yheads[j] * sqrt(scale_fac)
  end

  # Now to readjust arrows so that middles are on point '#onfleek'
  for i = 1:length(xbases)
    xbases[i] = xbases[i] - 0.5*xheads[i]
  end

  for i = 1:length(ybases)
    ybases[i] = ybases[i] - 0.5*yheads[i]
  end

  # Perform the plotting
  gr()
  quiver(xbases, ybases, quiver=(xheads, yheads))
end

better_gen_phase(data, 0.5)
