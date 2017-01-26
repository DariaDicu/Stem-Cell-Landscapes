a = []

for i in 1:10
  push!(a,i^3)
end

println(a)

x = []
for i in 1:length(a)
  push!(x,i)
end

using Plots
plot(x,a)

count = 0
for i in 1:3
    count += 1
end
println(count)

c = []
d = []
c,d = 2,3

f,g = 5,6

f
