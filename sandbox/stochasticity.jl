using Plots, ImageMagick, SymPy
a = 2.6
b = 1.8
c = 2.1

F = function (x)
  a = 2.6
  b = 1.8
  c = 2.1
  x = sin((a*x) + b) + sin(c*x)
  return x
end

F_prime = function (x)
  a = 2.6
  b = 1.8
  c = 2.1
  x = a*cos((a*x) + b) + c*cos(c*x)
  return x
end


range = Array(0:0.01:3pi/2)
n = length(range)

xs = rand(1:n,20)
ys = F(range[xs])

plot(F(range))
scatter!(xs, ys)

in_wells = []

tlim = 300

grad_min = minimum(diff(range))
grad_max = maximum(diff(range))

@gif for t = 1:tlim
  for i = 1:length(xs)
    norm_grad = (F_prime(range[xs[i]])+abs(grad_min))/(abs(grad_min)+grad_max)
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
