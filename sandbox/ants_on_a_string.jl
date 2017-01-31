using Plots, ImageMagick, SymPy, Distributions

#Define parameters for the simulation later. Volume is the factor by which the
#Gaussian noise is scaled. The higher the volume, the greater the noise.
num_particles = 50
tlim = 200
volume = 2


#Using SymPy notation, defined x as a variable
@syms x

#Also defined the parameters values unique for this sample function
a = 2.6
b = 1.8
c = 2.1

#Define the user inputted function here using SymPy notation.
#Currently hard coded as a sine based function but should be variable in the
#future based on what is passed forward from Furan.
f = sin((a*x) + b) + sin(c*x)

#The SymPy function is then converted to a string and a normal Julia function
#constructed from it using the eval(parse(...)) notation given. Flipping the
#function from SymPy to base Julia in this way was necessary as the SymPy
#function could not take an array as an argument in the way it was needed, and
#base Julia prevented automatic differentiation of a given function.

#Convert to string
f_str=string(f)
#Construct base Julia function F()
eval(parse("F=function (x) \n x = "  * f_str * " \n return x \nend"))

#Get equation for first derivative of user function using SymPy and convert to
#string
g=diff(f,x)
g_str=string(g)
#Construct base Julia function F_prime() for the first derivative (gradient)
#of the user function.
eval(parse("F_prime=function (x) \n x = "  * g_str * " \n return x \nend"))


#range is an array of all the x values plotted over, and n the size of range
range = Array(0:0.01:3pi/2)
n = length(range)

#Now for each particle simulated, a random x value is selected and stored in
#the array xs
xs = rand(1:n,num_particles)

#For each new x value, the corresponding y value is calculated by applying the
#function F to x i.e y = F(x). Each of these y's are also stored in the ys array
ys = F(range[xs])

plot(F(range))
scatter!(xs, ys)

in_wells = []

#The maximum and minimum gradient for each point plotted is calculated for
#future use in directing the particle up or down the hills in the landscape.
grad_min = minimum(F_prime(range))
grad_max = maximum(F_prime(range))

#Now using @gif to generate gifs natively in Julia without having to stitch
#together photos at a later point. Each frame represents the next time point in
#the simulation.
@gif for t = 1:tlim

  #For each of the ants on the line...
  for i = 1:num_particles

    #The gradient at the current point is evaluated. Using the grad_min and
    #grad_max from earlier, this gradient is then scaled between 0 and 1.
    #A normalised gradient (norm_grad) of 0 represents the steepest negative
    #gradient observed for the input function over the visualised interval
    #and a normalised gradient of 1 represents the steepest positive gradient.
    norm_grad = (F_prime(range[xs[i]])-grad_min)/(grad_max-grad_min)

    #A random number between 0 and 1 is generated, and assigned to p
    p = rand()

    #A normally distributed random noise term is also generated and scaled
    #by the volume factor defined during the initial parameterisaton.
    noise = rand(Normal(0,1)) * volume

    if p > norm_grad
      new_value = Int(xs[i] + round(1+noise))
      # Checking the new value does not fall out the [1, n] interval.
      xs[i] = max(min(new_value, n), 1)
    elseif p < norm_grad
      new_value = Int(xs[i] - round(1+noise))
      xs[i] = max(min(new_value, n), 1)
    end
  end
  plot(F(range))
  plot!(zeros(n)-0.5)
  ys = F(range[xs])
  scatter!(xs, ys)

  push!(in_wells, length(ys[ys .< 0.5]))
end every 1
