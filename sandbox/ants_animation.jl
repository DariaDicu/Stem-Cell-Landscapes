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
ldens_max = maximum(ldens) # Need to keep it separate for normalizing ants.
ldens=ldens-ldens_max
kernel_interpolation = InterpKDE(dens1)

# Plot the landscape and animate the trajectory with ants.
using GLVisualize, GLAbstraction, ModernGL, Reactive, GeometryTypes, Colors, GLWindow
using Interpolations
import GLVisualize: slider, mm, button

# Get 10 of the simulations to draw as ants.
ant_lines = []
for i = 1:10
  # Extract the rows in the DataFrame where the run index is i.
  current_run = data[data[4].==i,:]
  txy = current_run[1:3]

  # Interpolate to get smoother ants.
  x_spl = interpolate((txy[1], ), txy[2], Gridded(Linear()))
  y_spl = interpolate((txy[1], ), txy[3], Gridded(Linear()))
  tmin = minimum(txy[1])
  tmax = maximum(txy[1])
  tspan = linspace(tmin, tmax, 1000)


  # Extract the trajectory as an array of 3D points.
  ant_line = Point3f0[]
  for t in tspan
    x = x_spl[t]
    y = y_spl[t]
    z = pdf(kernel_interpolation, x, y)
    z_log = -log(1e-23+z)-ldens_max
    push!(ant_line, Point(x, y, z_log))
  end
  push!(ant_lines, ant_line)
end

ant_lines

ldens

window = glscreen()
colors = colormap("Blues", length(ant_lines))

# Compute bounding box for a smooth animation. This saves time when rerendering
# lines, as there is no need to recompute the bounding box automatically.
min_x = minimum(dens1.x)
min_y = minimum(dens1.y)
min_z = minimum(ldens)
max_x = maximum(dens1.x)
max_y = maximum(dens1.y)
max_z = maximum(ldens)
bounding_box = AABB{Float32}(Vec3f0(minimum([min_x, min_y, min_z])),
    Vec3f0(maximum([max_x, max_y, max_z])))

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
      prerender=()->glDisable(GL_DEPTH_TEST), # Draw over other items.
      postrender=()->glEnable(GL_DEPTH_TEST)), window, camera=:perspective)
end

myline=Point3f0[Point3f0(0), Point(5.0,0.0,0.0), Point(10.0,0.0,0.0)]

c = RGBA(255.0, 0.0, 0.0, 1.0)
_view(visualize(
    myline, :lines,
    color=c,
    prerender=()->glDisable(GL_DEPTH_TEST), # Draw over other items.
    postrender=()->glEnable(GL_DEPTH_TEST)), window, camera=:perspective)

gx = LinSpace(dens1.x)
gy = LinSpace(dens1.y)
g = GLVisualize.Grid(gx, gy)

vis = visualize((g, ldens), :surface, camera=:perspective)

coords = GLVisualize._position_calc(g, vis)
# Plot ldens.

#g = GLVisualize.Grid(itp_x, itp_y)

#vis = visualize((g, itp_log), :surface, boundingbox=bounding_box)


_view(vis, window)

renderloop(window)
