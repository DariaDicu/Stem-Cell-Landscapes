using Plots

my_plot = scatter(0,0)

plotly(0,0)

for i = 0:5
  Plots.plot(ldens[(i*50)+1,:])
  scatter!(0,-22)
  xlabel!("X")
  ylabel!("'Potential'")
  Plots.savefig(string(i)*".png")
end


plot!(ldens[50,:])

plot!(ldens[100,:])

plot!(ldens[150,:])

plot!(ldens[200,:])

plot!(ldens[250,:])
