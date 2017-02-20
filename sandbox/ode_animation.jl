# Code for generating a Julia OpenGL model that plots a cell fate landscape
# that changes as the parameter controlled by a slider changes.
using Base, DifferentialEquations, Plots, DataFrames, KernelDensity;
using Interpolations;
using ODESimulator;

# Example of function for representing a 2 transcription-factor with self- and
# mutual- regulation. Parameters are hardcoded in the function. See Wang et al,
# 2011 (http://www.pnas.org/content/108/20/8257.full).
<<<<<<< HEAD
<<<<<<< HEAD
#= e.g. F(0.3) instantiates the parameter a with value 0.3.
=======
# e.g. F(0.3) instantiates the parameter a with value 0.3.
>>>>>>> daria-working-branch
=======
# e.g. F(0.3) instantiates the parameter a with value 0.3.
>>>>>>> daria-working-branch
F = a -> function (t,x)
  # a = 0.3
  n = 4
  S = 0.5
  k = b = 1
  F1 = (x1, x2) ->
    (a*(x1^n)/(S^n + x1^n) + b*S^n/(S^n + x2^n) - k*x1)
  F2 = (x1, x2) ->
    (a*(x2^n)/(S^n + x2^n) + b*S^n/(S^n + x1^n) - k*x2)
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end
<<<<<<< HEAD
<<<<<<< HEAD
=#
=======

>>>>>>> daria-working-branch
=======

>>>>>>> daria-working-branch
# Generate 15 data sets as parameter a is varied in .1 increments:
# 0.0, 0.1, 0.2,..., 1.5. Each data set will be a landscape for the ODE when
# parameter a has that value (e.g. a=0.3). The datasets must be generated
# before the model is rendered for smooth interaction.

# Z will contain the landscape heights for each of the 15 ODEs.
Z = []
# itp_x and itp_y are used to build the interpolation grid.
itp_x = linspace(0,10,1000)
itp_y = linspace(0,10,1000)

for i = 1:16
  a = (i-1)*0.1
  data = build_landscape(1000, F(a), 2, (0,5))
  X = convert(Array{Float64},deepcopy(data[2]));
  Y = convert(Array{Float64},deepcopy(data[3]));

  dens1 = kde((X, Y))
  dens2=1e-23*ones(size(dens1.density))+dens1.density

  ldens=-log(dens2);
  ldens=ldens-maximum(ldens)

  # Interpolate the density for two reasons:
  # 1. We want an equally spaced (x, y) grid when building the surface. The
  #    datapoints in dens1.x, dens1.y may not be equally spaced.
  # 2. We want to sample at more than the number of points in the density for
  #    a smooth surface.
  itp_z = pdf(dens1, itp_x, itp_y)
  push!(Z, itp_z)

  # Uncomment for commands for plotting the negative log density.
  # itp_log = interpolate((dens1.x, dens1.y), ldens, Gridded(Linear()))
  # push!(Z_log, itp_z[itp_x, itp_y])
end

using GLVisualize, GLAbstraction, Reactive, GeometryTypes, Colors, GLWindow;
using ModernGL;
import GLVisualize: slider, mm, button;

window = glscreen()

const T = Float64
const P = Point{2, T}
iconsize = 8mm

# Create partitioned window for controls and view screens.
editarea, viewarea = x_partition_abs(window.area, round(Int, 8.2 * iconsize))
edit_screen = Screen(
    window, area = editarea,
    color = RGBA{Float32}(0.0f0, 0.0f0, 0.0f0, 1f0),
    stroke = (1f0, RGBA{Float32}(0.13f0, 0.13f0, 0.13f0, 13f0)))
viewscreen = Screen(
    window, area = viewarea,
    color = RGBA(0.0f0, 0.0f0, 0.0f0, 1f0))

# Function for creating a labelled slider.
function labeled_slider(range, window)
    visual, signal = slider(
        range, window;
        slider_length = 6 * iconsize,
        icon_size = Signal(iconsize / 2),
        knob_scale = 3mm,)
    text = visualize(
        map(string, signal), # convert to string
        relative_scale = 5mm,
        color = RGBA(1f0, 1f0, 1f0, 1f0))
    # put in list and visualize so it will get displayed side to side
    # direction = first dimension --> x dimension
    visualize([visual, text],  direction = 1, gap = Vec3f0(3mm, 0, 0)), signal
end

# Get the control and signal for the slider and center_cam button.
param_a_v, param_a_s = labeled_slider(0.0:0.1:1.5, edit_screen)
center_v, center_s = button("⛶", relative_scale = iconsize, edit_screen)

controls = Pair[
    "parameter a" => param_a_v,
    "center cam" => center_v]

_view(visualize(
        controls,
        text_scale = 4mm,
        width = 8iconsize), edit_screen, camera = :fixed_pixel)

# Setup signal for the "center camera" button.
s = preserve(map(center_s) do clicked
    clicked && center!(cam, AABB(value(line_pos)))
    nothing
end)

# Initialize the surface with model for parameter a = 0.0.
v0 = Z[1]

# Connect surface updating to the slider signal using foldp.
surface = foldp(v0,param_a_s) do v0,param_a_s
  i = convert(Int64, param_a_s*10);
  v0 = Z[i+1];
  v0
end;

# Compute the maximum possible bounding box over all the animation frames so
# OpenGL does not have to recompute it every time the signal is updated. This
# makes the animation significantly smoother.
max_x = maximum(itp_x)
max_y = maximum(itp_y)
max_z = maximum(map(zi->maximum(zi), Z))
bounding_box = AABB{Float32}(Vec3f0(0), Vec3f0(maximum([max_x, max_y, max_z])))

# Important: pass a static bounding box to make animation smooth.
xy_grid = GLVisualize.Grid(itp_x, itp_y)
vis = visualize((xy_grid, surface), :surface, boundingbox=bounding_box)

# Display it on the viewing window.
_view(vis, viewscreen, camera=:perspective)

renderloop(window)

close(param_a_s, false)
