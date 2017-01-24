# TASK: Given a set of trajectories, one from each simulation, find -log p(x)

#Make the fake data
num_sims=50
num_dims=5

max_timepoints = 15
max_all = 20

max_x=15 #GATA6
max_y=10 #NANOG

my_in_array=Array[] #Initialise empty array

for i = 1:num_sims #For each simulation
  T=[]
  for j = 1:rand(1:max_timepoints) # and up to a maximum number of timepoints
    coords=[]
    for dim = 1:num_dims #Push a random X Y Z etc coordinate
      push!(coords, rand(1:max_all))
    end
    push!(T, coords) #Then push this (X,Y,Z etc) group to the storage array
  end
  push!(my_in_array,T) #Push this array of timepoints to the current simulation
end
#Fake data made and populated!

#Initialise count matrix

#Now do the probability estimation and make the contours
using SymPy, Distributions, DataFrames, KernelDensity, Plots

#Create visualisation function
function make_contours(dim1, dim2)
  max_all = 20

  #Initialise arrays
  count_matrix=zeros(max_all, max_all)
  all_X=[]
  all_Y=[]

  #Increment count matrix cells and get the dim1, dim2 points
  for i =1:length(my_in_array)
    current_sim = my_in_array[i]
    for j = 1:length(current_sim) #For each time point in the current simulation
      current_position = current_sim[j]
      count_matrix[current_position[dim2],current_position[dim1]] += 1
      push!(all_X,current_position[dim1])
      push!(all_Y, current_position[dim2])
    end
  end

  X1C = convert(Array{Float64}, deepcopy(all_X))
  Y1C = convert(Array{Float64}, deepcopy(all_Y))

  C1kdens=kde((X1C,Y1C))
  Cdens=1e-23*ones(256,256)+C1kdens.density

  lCdens=-log(Cdens);
  lCdens=lCdens-maximum(lCdens)

  p1CP=contour(C1kdens.x,C1kdens.y,C1kdens.density,levels=100,legend=false,xlabel=dim1,ylabel=dim2)
  p2CP=contour(C1kdens.x,C1kdens.y,lCdens,levels=[10],legend=false,xlabel=dim1,ylabel=dim2)
  plot(p1CP)
  return (p1CP, p2CP, C1kdens)
end

#Entering the domain of plotting
#Store p1CP for all possible combinations of observed dimensions
plot_array=[]
D3_array=[]

for i = 1:num_dims
  for j = 1:num_dims
    (p1CP, p2CP, C1kdens) = make_contours(i,j)
    push!(plot_array, p1CP)
    push!(D3_array, C1kdens)
  end
end

plot(plot_array[1], #This is a stupid way of calling this data
 plot_array[2],
 plot_array[3],
 plot_array[4],
 plot_array[5],
 plot_array[6],
 plot_array[7],
 plot_array[8],
 plot_array[9],
 plot_array[10],
 plot_array[11],
 plot_array[12],
 plot_array[13],
 plot_array[14],
 plot_array[15],
 plot_array[16],
 plot_array[17],
 plot_array[18],
 plot_array[19],
 plot_array[20],
 plot_array[21],
 plot_array[22],
 plot_array[23],
 plot_array[24],
 plot_array[25],
 layout=(5,5))


(p1CP, p2CP, C1kdens) = make_contours(1,2)
plot(p1CP)

#Also make 3D representation of the contour plot
using PyPlot
contour3d(C1kdens.x, C1kdens.y,C1kdens.density)
surf(C1kdens.x, C1kdens.y,C1kdens.density,cmap="winter")
show()

#Can also make discretised heatmap of data (not yet -log)
norm_matrix=count_matrix./sum(count_matrix)

using Plots
pyplot()
xs = [string("x",i) for i = 1:max_x]
ys = [string("y",i) for i = 1:max_y]
z = norm_matrix
second_plot=heatmap(xs,ys,z,aspect_ratio=1)
plot(p1CP,second_plot)
