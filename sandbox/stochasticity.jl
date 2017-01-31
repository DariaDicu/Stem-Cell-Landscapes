using Plots, ImageMagick, SymPy

num_particles = 50
tlim = 200

@syms x
a = 2.6
b = 1.8
c = 2.1
f = sin((a*x) + b) + sin(c*x)

f_str=string(f)
eval(parse("F=function (x) \n x = "  * f_str * " \n return x \nend"))

g=diff(f,x)
g_str=string(g)
eval(parse("F_prime=function (x) \n x = "  * g_str * " \n return x \nend"))

range = Array(0:0.01:3pi/2)
n = length(range)

xs = rand(1:n,num_particles)
ys = F(range[xs])

plot(F(range))
scatter!(xs, ys)

in_wells = []

grad_min = minimum(F_prime(range))
grad_max = maximum(F_prime(range))

@gif for t = 1:tlim
  for i = 1:length(xs)
    norm_grad = (F_prime(range[xs[i]])-grad_min)/(grad_max-grad_min)
    p = rand()
    if p > norm_grad && xs[i] < n
      xs[i] += 1
    elseif p < norm_grad && xs[i] > 1
      xs[i] -= 1
    end
  end
  plot(F(range))
  plot!(zeros(n)-0.5)
  ys = F(range[xs])
  scatter!(xs, ys)

  push!(in_wells, length(ys[ys .< 0.5]))
end every 1
