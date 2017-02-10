using Plots, ImageMagick, SymPy, Distributions

# Define parameters for the simulation later. Volume is the factor by which the
# Gaussian noise is scaled. The higher the volume, the greater the noise.
num_particles = 50
tlim = 200
volume = 2


# Using SymPy notation, defined x as a variable
@syms x

# Also defined the parameters values unique for this sample function
a = 2.6
b = 1.8
c = 2.1

# Define the user inputted function here using SymPy notation.
# Currently hard coded as a sine based function but should be variable in the
# future based on what is passed forward from Furan.
f = sin((a*x) + b) + sin(c*x)

# Experimenting with a second function
a = 4
b = -0.1
c = 2.6
d = 0.7

f = (x^a) - c*(x+b)^2

# The SymPy function is then converted to a string and a normal Julia function
# constructed from it using the eval(parse(...)) notation given. Flipping the
# function from SymPy to base Julia in this way was necessary as the SymPy
# function could not take an array as an argument in the way it was needed, and
# base Julia prevented automatic differentiation of a given function.

# Convert the Sympy function to a string
f_str= "f ="*string(f)

# Construct base Julia function F() using this string
eval(parse("F=function (x) \n x = "  * f_str * " \n return x \nend"))

# Get equation for first derivative of user function using SymPy and convert
# to string
g=diff(f,x)
g_str=string(g)

# Construct base Julia function F_prime() for the first derivative (gradient)
# of the user function.
eval(parse("F_prime=function (x) \n x = "  * g_str * " \n return x \nend"))


# range is an array of all the x values plotted over, and n the size of range
#range = Array(0:0.01:3pi/2)
range = Array(-1.8:0.01:1.6)
n = length(range)

# Now for each particle simulated, a random x value is selected and stored in
# the array xs
xs = rand(1:n,num_particles)

# For each new x value, the corresponding y value is calculated by applying
# the function F to x i.e y = F(x). Each of these y's are also stored in the
# ys array
ys = map(F, range[xs])
plot(map(F, range))
scatter!(xs, ys)

well_A = []
well_B = []

# The maximum and minimum gradient for each point plotted is calculated for
# future use in directing the particle up or down the hills in the landscape.

all_grads = map(F_prime, range)

grad_min = minimum(all_grads)
grad_max = maximum(all_grads)

x_cross=solve(f)

# Now using @gif to generate gifs natively in Julia without having to stitch
# together photos at a later point. Each frame represents the next time point
# in the simulation.
@gif for t = 1:tlim

  # For each of the ants on the line...
  for i = 1:num_particles

    # The gradient at the current point is evaluated (current_grad).
    # Using the grad_min and grad_max from earlier, this gradient is then
    # scaled between 0 and 1. A normalised gradient (norm_grad) of 0 represents
    # the steepest negative gradient observed for the input function over
    # the visualised interval and a normalised gradient of 1 represents the
    # steepest positive gradient.
    current_grad = F_prime(range[xs[i]])
    norm_grad = (current_grad-grad_min)/(grad_max-grad_min)

    # A random number between 0 and 1 is generated, and assigned to p
    p = rand()

    # A normally distributed random noise term is also generated and scaled
    # by the volume factor defined during the initial parameterisaton.
    noise = rand(Normal(0,1)) * volume

    if p > norm_grad

      # If p is greater than the normalised gradient (more likely for a
      # negatibe/shallow-positive gradient) then the x position is increased
      # according to the step size and the noise term. I.e more likely to
      # move forward when facing downhill than uphill.
      new_value = Int(xs[i] + round(1+noise))

      # Checking that the new value does not fall out the [1, n] interval and
      # correcting accordingly.
      xs[i] = max(min(new_value, n), 1)

    # Else if p is less than norm_grad, the particle moves backwards. This is
    # more likely if the particle is on a positive slope than a negative one.
    elseif p < norm_grad

      # Again the position is updated according to the step and noise terms
      new_value = Int(xs[i] - round(1+noise))

      # Again, checking the new value is not outside the [1,n] interval
      xs[i] = max(min(new_value, n), 1)
    end
  end

  # After each time point, replot the function...
  plot(map(F, range))

  # Calculate the new y values by applying the function F(x) to the updated
  # x values.
  ys = map(F, range[xs])

  # Plot these new (x,y) pairs on top of the plain function plot
  scatter!(xs, ys)

  # A few extra things to maybe consider in the future i.e counting all
  # particles with y values below a certain y threshold could be used to plot
  # how many particles are in wells at a given time.
  plot!(zeros(n))

  num_A = 0; num_B = 0

  for i = 1:num_particles
    if ys[i] < 0
      if x_cross[1] < range[xs][i] < x_cross[2]
        num_A = num_A + 1
      elseif x_cross[end-1] <range[xs][i] < x_cross[end]
        num_B = num_B + 1
      end
    end
  end

  push!(well_A, num_A)
  push!(well_B, num_B)

  # End of the @gif command
end every 1  # Create a gif, saving the plot 'every 1' frame

plot(well_A)
plot!(well_B)
