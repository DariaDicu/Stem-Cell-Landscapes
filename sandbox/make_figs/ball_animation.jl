# Code for generating a Julia OpenGL model that plots the landscape for a fixed
# value of parameter a and animated moving "ant" trajectories onto the model
# from a subset of the simulations.
using ODESimulator;

F = function (t,x)
  n = 4
  S = 0.5
  k = b = 1
  F1 = (x1, x2) ->
    (a*(x1^n)/(S^n + x1^n) + b*S^n/(S^n + x2^n) - k*x1)
  F2 = (x1, x2) ->
    (a*(x2^n)/(S^n + x2^n) + b*S^n/(S^n + x1^n) - k*x2)
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end

runs = 100
data = build_landscape(runs, F, 2, (0,5))

# Get the trajectory coordinates for plotting along the entire trajectory.
X_along = convert(Array{Float64},deepcopy(data[2]));
Y_along = convert(Array{Float64},deepcopy(data[3]));

# Get the endpoint trajectory coordinates.
X_endpoints = Float64[]
Y_endpoints = Float64[]
for i = 1:runs
  # Extract the rows in the DataFrame where the run index is i.
  current_run = data[data[4].==i,:]
  push!(X_endpoints, current_run[end, 2])
  push!(Y_endpoints, current_run[end, 3])
end

using KernelDensity, Interpolations, DataFrames;

# TODO: Button for this.
scale_factor = 5.0 # Scale for xy grid.
# TODO: Button for this (or scale depending on scale_factor).
sphere_radius = scale_factor*0.1f0
# TODO: Button for this.
# Factor has to be between 0 and 1!
ant_speed_factor = 0.1f0
# TODO: Button for this?
trace_length = 10
# TODO: Button.
ant_count = 10

# TODO: Specify types of X, Y
# Calculates the landscape heights from X-Y simulation data. Returns two
# 2D arrays: the probability density and the log density, as well as the
# sampling grid for the X-Y axes.
function get_heights(X, Y)
  # Change to X, Y for alternative plotting.
  dens1= kde((X, Y))
  dens2=1e-23*ones(size(dens1.density))+dens1.density
  log_dens=-log(dens2)
  log_dens=log_dens-maximum(log_dens)
  gx = scale_factor*LinSpace(dens1.x)
  gy = scale_factor*LinSpace(dens1.y)
  # Convert everything from Float64 to Float32, as required by GLVisualize.
  dens = convert(Array{Float32, 2}, dens1.density)
  log_dens = convert(Array{Float32, 2}, log_dens)
  gx = convert(LinSpace{Float32}, gx)
  gy = convert(LinSpace{Float32}, gy)
  return dens, log_dens, gx, gy
end

# Plot the landscape and animate the trajectory with ants.
using GLVisualize, GLAbstraction, ModernGL, Reactive, GeometryTypes, Colors, GLWindow
using Interpolations
import GLVisualize: slider, mm, button, toggle_button

# Eliminate identical consecutive points and keep only one copy of each.
function dedup_consecutives(traj)
  if (length(traj) == 0)
    return traj
  end
  dedup_traj = Point3f0[traj[1]]
  for i = 2:length(traj)
    if (traj[i] != traj[i-1])
      push!(dedup_traj, traj[i])
    end
  end
  return dedup_traj
end

function get_xy_node_matrices(gx, gy)
  mat_x = repmat(transpose(gx), length(gy), 1)
  mat_y = repmat(gy, 1, length(gx))
  return mat_x, mat_y
end

dens_endpoints, ldens_endpoints, gx_endpoints, gy_endpoints =
  get_heights(X_endpoints, Y_endpoints)
dens_along, ldens_along, gx_along, gy_along =
  get_heights(X_along, Y_along)

function get_ant_lines(surf, ant_count::Int64)
  gx, gy, dens = surf
  gx = gx[1, :]
  gy = gy[:, 1]
  # Get 10 of the simulations to draw as ants.
  ant_lines = []
  for i = 1:ant_count
    # Extract the rows in the DataFrame where the run index is i.
    current_run = data[data[end].==i,:]
    txy = current_run[1:3]
    x_spl = interpolate((txy[1], ), txy[2], Gridded(Linear()))
    y_spl = interpolate((txy[1], ), txy[3], Gridded(Linear()))
    tmin = minimum(txy[1])
    tmax = maximum(txy[1])
    # TODO: control or refine the # of interpolation points (1000 too much?)
    tspan = linspace(tmin, tmax, 500)
    # Extract the trajectory as an array of 3D points.
    ant_line = Point3f0[]
    for t in tspan
      x = x_spl[t]
      y = y_spl[t]
      i = indmin(abs(gx-scale_factor*x))
      j = indmin(abs(gy-scale_factor*y))
      z = dens[i,j]
      push!(ant_line, Point(gx[i], gy[j], z))
    end
    push!(ant_lines, dedup_consecutives(ant_line))
  end
  return ant_lines
end

window = glscreen()
iconsize = 8mm
assets_path = string(homedir(), "/Documents/stem-cells/assets/");

# Create partitioned window for controls and view screens.
editarea, viewarea = x_partition_abs(window.area, 180)
# Further partition edit area to get a logo area.
editarea, logoarea = y_partition(editarea, 85)
logoarea
edit_screen = Screen(
    window, area = editarea,
    color = RGBA{Float32}(0.0f0, 0.0f0, 0.0f0, 1f0))
view_screen = Screen(
    window, area = viewarea,
    color = RGBA(0.0f0, 0.0f0, 0.0f0, 1f0),
    stroke = (1f0, RGBA{Float32}(0.13f0, 0.13f0, 0.13f0, 13f0)))
logo_screen = Screen(
    window, area = logoarea,
    color = RGBA{Float32}(0.0f0, 0.0f0, 0.0f0, 1f0))


logo_signal = map(logoarea) do a
  xc = a.x + (a.w)/2
  yc = a.y + (a.h)/2
  [Point2f0(xc, yc)]
end

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
    visualize([visual, text], direction = 1, gap = Vec3f0(3mm, 0, 0)), signal
end

# Get the control and signal for the slider and center_cam button.
#param_a_v, param_a_s = labeled_slider(0.0:0.1:1.5, edit_screen)
on_button_img = loadasset(string(assets_path, "on.png"))
off_button_img = loadasset(string(assets_path, "off.png"))
endpoint_v, endpoint_s = toggle_button(
  on_button_img, off_button_img, edit_screen)
log_dens_v, log_dens_s = toggle_button(
  on_button_img, off_button_img, edit_screen)

controls = Pair[
    "Endpoints only" => endpoint_v,
    "Negative log density" => log_dens_v]

_view(visualize(
        controls,
        text_scale = 5mm,
        gap = 3mm,
        width = 10iconsize), edit_screen, camera = :fixed_pixel)

#logo = loadasset(string(homedir(), "/Documents/stem-cells/assets/romeo.png"))
#logo_vis = visualize((SimpleRectangle(0,0,100,100), logo_signal),
#  scale=Vec3f0(100), image=logo)
#_view(logo_vis)

logo_text = visualize(
    "Waddle",
    relative_scale=16mm,
    color = RGBA(1f0, 1f0, 1f0, 1f0))
_view(logo_text, logo_screen, camera=:fixed_pixel)
########### Done setting up sidebar. ##########

# Important to use the grid and density as a single signal that updates at the
# same time. Computing ant lines replies on having the correct grid for heights.
surface_signal = map(endpoint_s, log_dens_s) do togg_endpoint, togg_log_dens
  dens = togg_endpoint ? (togg_log_dens ? ldens_endpoints : dens_endpoints) : (
    togg_log_dens ? ldens_along : dens_along)
  gx = togg_endpoint ? gx_endpoints : gx_along
  gy = togg_endpoint ? gy_endpoints : gy_along
  mat_x, mat_y = get_xy_node_matrices(gx, gy)
  return mat_x, mat_y, dens
end

gx = map(surface_signal) do surf # where surf = (mat_x, mat_y, dens)
  return surf[1]
end

gy = map(surface_signal) do surf # where surf = (mat_x, mat_y, dens)
  return surf[2]
end

dens = map(surface_signal) do surf # where surf = (mat_x, mat_y, dens)
  return surf[3]
end

ant_lines = map(surf->get_ant_lines(surf, ant_count), surface_signal)
red_color = RGBA(255.0, 0.0, 0.0, 1.0)
ant_sphere = GLNormalMesh(Sphere{Float32}(Vec3f0(0), sphere_radius))

function visualize_trajectory(ant_lines)
  max_traj = maximum(map(length, ant_lines))
  timesignal = preserve(loop((trace_length):max_traj))
  return preserve(map(timesignal) do t
      ps = Point3f0[]
      for i = 1:ant_count
        # Only add point if t is within trajectory bounds.
        # If t exceeds the bound, plot last element
        last_pos = (t <= length(ant_lines[i])) ? t : length(ant_lines[i])
        append!(ps, [ant_lines[i][last_pos]])
      end
      ps
    end)
end

using ReactiveBasics
ps = map(visualize_trajectory, ant_lines)
ps = flatten(ps, typ=Array{FixedSizeArrays.Point{3, Float32},1})
trace_color = RGBA{Float64}[RGBA(255.0, 0.0, 0.0, i*(1.0/trace_length))
  for i=1:trace_length]

vis_obj = visualize(
    (ant_sphere, ps),
    boundingbox=nothing,
    color=red_color, camera=:perspective)

_view(vis_obj, view_screen, camera=:perspective)

vis = visualize((gx, gy, dens), :surface, camera=:perspective)


_view(vis, view_screen)

renderloop(window)
#################################################
