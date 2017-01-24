#This code simulates particles moving around a 3D space and plots their
#position over time as well as the resulting density plot in two dimensions.

using Plots, KernelDensity

#Function to make the ants in space and simulate their movement.
#Takes the number of simulations (i.e. the number of ants), the number of
#dimensions in which the ant walks, how long to run the simulation, the upper
#limit for all dimensions and the index of the two dimensions to be visualised
#later in a heatmap.
function make_ants(num_sims, num_dims, tlim, var_bound, dim1, dim2)

  #First initialise the matrix to store all values. 3D matrix with each row
  #representing one simulation, each column the coordinate in a different dimension,
  #and each slice moving backwards through the matrix a different time point.
  A=zeros(num_sims, num_dims, tlim)
  all_x = []
  all_y = []

  #Spawn in initial position of particles. Simply generate a a random number
  #for each dimension between 1 and the variable bound. Repeat for each
  #simulation too.
  for sim=1:num_sims
    for dim=1:num_dims
      A[sim, dim, 1] = rand(1:var_bound)

      #If the current dimension is one that has been chosen for future
      #plotting, then push this value to the appropriate empty array.
      if dim==dim1
        push!(all_x, A[sim,dim,1])
      end
      if dim==dim2
        push!(all_y, A[sim,dim,1])
      end
    end
  end

  #Define how particles move between timepoints. ran_move() increases or
  #decreases a coordinate value by 1 (or does nothing) with equal probability.
  function ran_move(x)
    p = rand()
    if p<1/3 && x != var_bound #Keep the particle bound
          x = x + 1
    elseif p<2/3 && x != 0
          x = x -1
    end
    return x
  end

  #Alternatively, sink_move() pulls the particle towards a sink with p=0.5,
  #or simply applies the null ran_move().
  function sink_move(x, sink)
    p = rand()
    if p<0.5
      if x < sink
        x = x + 1
      elseif x > sink
        x = x - 1
      end
    else
      x = ran_move(x)
    end
    return x
  end

  #Now fill in the other time slices of the matrix by applying the movement
  #function ran_move() or sink_move() to each value in the 2D slice for t-1.
  for t=2:tlim
    for sim=1:num_sims
      for dim=1:num_dims
        x = A[sim, dim, t-1] #Next value depends on previous at t-1
        A[sim,dim,t] = sink_move(x, 14) #Apply sink_move() or ran_move()

        #Again, if the current dimension is one for which a plot will be
        #generated, push the value to the appropriate array for later.
        if dim==dim1
          push!(all_x, A[sim,dim,t])
        end
        if dim==dim2
          push!(all_y, A[sim,dim,t])
        end
      end
    end
  end
  return A, num_sims, tlim, dim1, dim2, all_x, all_y
end

#Function to form contours given x and y array containing all x and y values
#produced over the course of the simulation. Code from Stumpf:
function form_contour(all_x, all_y)
  X1C = convert(Array{Float64}, deepcopy(all_x))
  Y1C = convert(Array{Float64}, deepcopy(all_y))

  C1kdens=kde((X1C,Y1C))
  Cdens=1e-23*ones(256,256)+ C1kdens.density

  lCdens=-log(Cdens);
  lCdens=lCdens-maximum(lCdens)

  p1CP=contour(C1kdens.x,C1kdens.y,C1kdens.density,levels=100,legend=false,xlabel=dim1,ylabel=dim2)
  p2CP=contour(C1kdens.x,C1kdens.y,lCdens,levels=[10],legend=false,xlabel=dim1,ylabel=dim2)

  return p1CP
end

#see_ants() plots the paths taken in 3D space for all simulations. Produces
#a single plot within Atom. Due to how plots are updated in julia within loops,
#calling plot3d() three times in this way was the only way I could visualise
#the result.
function see_ants(A, num_sims)
  plot3d(A[1,1,:],A[1,2,:],A[1,3,:])
  for i = 2:num_sims-1
    plot3d!(A[i,1,:],A[i,2,:],A[i,3,:])
  end
  i=num_sims
  plot3d!(A[i,1,:],A[i,2,:],A[i,3,:])
end

#see_ants_walking() takes the matrix of data A and, similarly to see_ants(),
#shows the 3D plot of the paths taken up until a user defined t. This is
#plotted alongside the contour heatmap derived from all of the x and y
#coordinates generated up to time t. The plot is then saved as a .png for
#future stitching to create a GIF.
function see_ants_walking(A,t)
  #for t = 1:tlim   #This would be perfect but currently broken
    plot3d(A[1,1,1:t],A[1,2,1:t],A[1,3,1:t])
    for i = 2:num_sims-1
      plot3d!(A[i,1,1:t],A[i,2,1:t],A[i,3,1:t])
    end
    i=num_sims
    plot3d!(A[i,1,1:t],A[i,2,1:t],A[i,3,1:t])
    #savefig(string(t)*".png")
    A = plot3d!(A[i,1,1:t],A[i,2,1:t], A[i,3,1:t])
    cnt = t*num_sims
    plot(A, plot(form_contour(all_x[1:cnt], all_y[1:cnt])))
    savefig(string(t)*".png")
  #end
end

#simple_vis() performs the same plotting and saving function as
#see_ants_walking(), but without the corresponding contour plot.
function simple_vis(A)
  for t = 1:tlim
    plot3d(A[1,1,1:t],A[1,2,1:t],A[1,3,1:t])
    for i = 2:num_sims-1
      plot3d!(A[i,1,1:t],A[i,2,1:t],A[i,3,1:t])
    end
    i=num_sims
    plot3d!(A[i,1,1:t],A[i,2,1:t],A[i,3,1:t])
    savefig(string(t)*".png")
  end
end


function form_contourM(allx, ally)


  X1C = convert(Array{Float64}, copy(allx))

  Y1C = convert(Array{Float64}, copy(ally))

  C1kdens=kde((X1C,Y1C))
    ss1,ss2=size(C1kdens.density)
    #println("A",",",ss,",","",size(C1kdens.density))
    Cdens=1e-23*ones(ss1,ss2)+C1kdens.density
    #println("B")

  lCdens=-log(Cdens);

  lCdens=lCdens-maximum(lCdens)
    p1CP=contour(C1kdens.x,C1kdens.y,C1kdens.density,levels=100,legend=false,xlabel="x",ylabel="y")

    #p2CP=contour(C1kdens.x,C1kdens.y,lCdens,levels=[10],legend=false,xlabel=dim1,ylabel=dim2)

  return p1CP
end

function form_contourMtry(allx, ally)


  X1C = convert(Array{Float64}, copy(allx))

  Y1C = convert(Array{Float64}, copy(ally))

  C1kdens=kde((X1C,Y1C))
    ss1,ss2=size(C1kdens.density)
    #println("A",",",ss,",","",size(C1kdens.density))
    Cdens=1e-23*ones(ss1,ss2)+C1kdens.density
    #println("B")

  lCdens=-log(Cdens);

  lCdens=lCdens-maximum(lCdens)
    p1CP=contour(C1kdens.x,C1kdens.y,C1kdens.density,levels=100,legend=false,xlabel="x",ylabel="y")

    #p2CP=contour(C1kdens.x,C1kdens.y,lCdens,levels=[10],legend=false,xlabel=dim1,ylabel=dim2)

  return p1CP
end

function see_ants_walkingM1(A,t)
  #for t = 1:tlim   #This would be perfect but currently broken
    #t=1
    nsim=size(A)[1]
    println(num_sims,",",nsim)
    println(size(A))
    plot3d(A[1,1,1:t],A[1,2,1:t],A[1,3,1:t])
    for i = 2:nsim-1
      plot3d!(A[i,1,1:t],A[i,2,1:t],A[i,3,1:t])
    end

    i=nsim
    plot3d!(A[i,1,1:t],A[i,2,1:t],A[i,3,1:t])
    #savefig(string(t)*".png")
    Aq = plot3d!(A[i,1,1:t],A[i,2,1:t],A[i,3,1:t])

    cnt = t*nsim
    plot(Aq,form_contourMtry(all_x[1:cnt], all_y[1:cnt]))
    #form_contourM(all_x[1:cnt], all_y[1:cnt])
    savefig(string(t)*".png")
  #end
end


#Having defined all the functions, can now call them in order
(A, num_sims, tlim, dim1, dim2, all_x, all_y) = make_ants(10,3,500,20,1,2)
plot(form_contour(all_x,all_y))
see_ants(A, num_sims)

#Current error prevents see_ants_walking() completing for some time points.
#Stumpf looking into this now.

see_ants_walkingM1(A)

for t = 1:50 #If you're lucky, but still broken
  see_ants_walkingM1(A,t)
end
