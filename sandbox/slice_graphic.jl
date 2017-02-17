using Plots

# How many slices do you want
slices = 5

interval = Integer(floor(length(ldens)[2]/slices))

for i = 0:slices-1
  Plots.plot(ldens[(i*interval)+1,:], linecolor = i+1, linewidth=4)
  scatter!(0,-22)
  xlabel!("X")
  ylabel!("'Potential'")
  Plots.savefig(string(i)*".png")
end
