#Script to loop and print the functions contained in 'just_functions.jl'
#having also run ode_simulator.

for i = -1:0.05:1
  data, N, n = build_landscape(15, F_super_pitch(i), 2, (-4,4))

  #Toggle between rest only=true if only want to plot the end positions,
  #or rest_only=false to consider full trajectory
  rest_only = false

  if rest_only
    #Find the resting values for each run
    resting_values = zeros(N, n)

    for run = 1:N
      data_edit = data[data[size(data)[2]] .== run, :]
      for dim = 1:n
        resting_values[run, dim] = data_edit[findmax(data_edit[1])[2],dim+1]
      end
    end
    X = convert(Array{Float64},deepcopy(resting_values[:,1]));
    Y = convert(Array{Float64},deepcopy(resting_values[:,2]));
  else
    X = convert(Array{Float64},deepcopy(data[2]));
    Y = convert(Array{Float64},deepcopy(data[3]));
  end

  dens1 = kde((X, Y))
  dens2=1e-23*ones(size(dens1.density))+dens1.density

  ldens=-log(dens2);
  ldens=ldens-maximum(ldens)

  gr()
  contour_plot=contour(dens1.x,dens1.y,dens1.density,levels=100,
    legend=false,xlabel="Dim 1",ylabel="Dim 2")
  plot(contour_plot)

  savefig(string(i)*".png")
end
