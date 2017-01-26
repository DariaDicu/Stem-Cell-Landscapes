# Code for generating a Julia OpenGL model that plots the landscape for a fixed
# value of parameter a and animated moving "ant" trajectories onto the model
# from a subset of the simulations.
using ODESimulator;

F = function (t,x)
  a = 0.3
  n = 4
  S = 0.5
  k = b = 1
  F1 = (x1, x2) ->
    (a*(x1^n)/(S^n + x1^n) + b*S^n/(S^n + x2^n) - k*x1)
  F2 = (x1, x2) ->
    (a*(x2^n)/(S^n + x2^n) + b*S^n/(S^n + x1^n) - k*x2)
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end

data = build_landscape(1000, F, 2, (0,5))

using KernelDensity, Interpolations, DataFrames;

X = convert(Array{Float64},deepcopy(data[2]));
Y = convert(Array{Float64},deepcopy(data[3]));
dens1 = kde((X, Y))
dens2=1e-23*ones(size(dens1.density))+dens1.density
ldens=-log(dens2);
ldens=ldens-maximum(ldens)

# Interpolate the density to get a smoother surface.
itp_x = linspace(0,10,1000)
itp_y = linspace(0,10,1000)
itp_z = pdf(dens1, itp_x, itp_y)
kernel_interpolation = InterpKDE(dens1)


# Plot the landscape and animate the trajectory with ants.
using GLVisualize, GLAbstraction, ModernGL, Reactive, GeometryTypes, Colors, GLWindow
import GLVisualize: slider, mm, button

# Get 10 of the simulations to draw as ants.
ant_lines = []
for i = 1:10
  # Extract the rows in the DataFrame where the run index is i.
  current_run = data[data[4].==i,:]
  # Extract the trajectory as an array of 3D points.
  ant_line = Point3f0[]
  for traj_point in eachrow(current_run)
      x, y = traj_point[2], traj_point[3]
      z = pdf(ik, x, y)
      push!(ant_line, Point(x, y, z))
  end
  push!(ant_lines, ant_line)
end

window = glscreen()
colors = colormap("Blues", length(ant_lines))

# Compute bounding box for a smooth animation. This saves time when rerendering
# lines, as there is no need to recompute the bounding box automatically.
max_x = maximum(itp_x)
max_y = maximum(itp_y)
max_z = maximum(itp_z)
bounding_box = AABB{Float32}(Vec3f0(0), Vec3f0(maximum([max_x, max_y, max_z])))

# Draw the trajectory of the ant up until the i-th time point. The value of i
# comes from a time signal that loops from 1 to the size of ant_line and is the
# basis of the animation.
function draw_ant_trajectory(i, ant_line)
  ant_line[1:i]
end

for i = 1:length(ant_lines)
  ant_line = ant_lines[i]
  # Create a time signal and use const_lift to animate the trajectory by calling
  # the function draw_ant_trajectory for each value of timesignal.
  timesignal = loop(1:length(ant_line))
  _view(visualize(
      const_lift(draw_ant_trajectory, timesignal, ant_line), :lines,
      color=colors[i],
      boundingbox=nothing,
      prerender=()->glDisable(GL_DEPTH_TEST), # Draw over other items.
      postrender=()->glEnable(GL_DEPTH_TEST)), window, camera=:perspective)
end

g = GLVisualize.Grid(itp_x, itp_y)

vis = visualize((g, itp_z), :surface, boundingbox=bounding_box)

_view(vis, window)

renderloop(window)
