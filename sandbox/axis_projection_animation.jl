# Code for generating a Julia OpenGL model that plots the landscape for a fixed
# value of parameter a and animated moving "ant" trajectories onto the model
# from a subset of the simulations.
include("ode_simulator.jl")
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

runs = 100 # Number of simulation runs
n = 2 # Number of dimensions
data = build_landscape(runs, F, 2, (0,5))

using KernelDensity, Interpolations, DataFrames, Reactive;

# TODO: Button for this.
#scale_factor_s = Signal(5.0) # Scale for xy grid.
# TODO: Button for this.
# Factor has to be between 0 and 1!
ant_speed_factor = 0.1f0
# TODO: Button for this?
trace_length = 10
# TODO: Button.
#ant_count_s = Reactive.Signal(10)

# Returns the data for producing a landscape for the dimensions corresponding to
# dim1 and dim2. Returns either the entire trajectory or the endpoints only,
# depending on the value of is_endpoint.
function getXYdata(data, is_endpoint, dim1, dim2)
  if (!is_endpoint)
    # Get the trajectory coordinates for plotting along the entire trajectory.
    # Note: indexing at (dim+1) since first column represents the time value.
    X = convert(Array{Float64},deepcopy(data[dim1+1]));
    Y = convert(Array{Float64},deepcopy(data[dim2+1]));
  else
    # Get the endpoint trajectory coordinates.
    X = Float64[]
    Y = Float64[]
    for i = 1:runs
      # Extract the rows in the DataFrame where the run index is i.
      # Note: indexing at (dim+1) since first column represents the time value.
      current_run = data[data[4].==i,:]
      push!(X, current_run[end, dim1+1])
      push!(Y, current_run[end, dim2+1])
    end
  end
  return X, Y
end

# TODO: Specify types of X, Y
# Calculates the landscape heights from X-Y simulation data. Returns a 2D array,
# either the probability density or the log density, depending on the is_log
# argument. It also returns the sampling grid for the X-Y axes.
function get_heights(X, Y, is_log, scale_factor)
  # Change to X, Y for alternative plotting.
  dens= kde((X, Y))
  dens_plus_background=1e-23*ones(size(dens.density))+dens.density
  log_dens=-log(dens_plus_background)
  log_dens=log_dens-maximum(log_dens)
  gx = scale_factor*LinSpace(dens.x)
  gy = scale_factor*LinSpace(dens.y)

  # Convert everything from Float64 to Float32, as required by GLVisualize.
  dens = convert(Array{Float32, 2}, dens.density)
  log_dens = convert(Array{Float32, 2}, log_dens)
  gx = convert(LinSpace{Float32}, gx)
  gy = convert(LinSpace{Float32}, gy)
  return gx, gy, (is_log ? log_dens : dens)
end

# Plot the landscape and animate the trajectory with ants.
using GLVisualize, GLAbstraction, ModernGL, Reactive, GeometryTypes, Colors, GLWindow
using Interpolations
import GLVisualize: labeled_slider, mm, button, toggle_button

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

function get_x_node_matrix(gx, gy)
  repmat(gx, 1, length(gy))
end

function get_y_node_matrix(gx, gy)
  repmat(transpose(gy), length(gx), 1)
end

# Extract ant trajectories as 2D array depending on the signals for surface,
# dimensions to plot and number of ants.
function get_ant_lines(surf, dim1, dim2, ant_count, scale_factor)
  gx, gy, dens = surf
  # Get 10 of the simulations to draw as ants.
  ant_lines = []
  for i = 1:ant_count
    # Extract the rows in the DataFrame where the run index is i.
    current_run = data[data[end].==i,:]
    t_data = current_run[1]
    x_data = current_run[dim1+1]
    y_data = current_run[dim2+1]
    x_spl = interpolate((t_data, ), x_data, Gridded(Linear()))
    y_spl = interpolate((t_data, ), y_data, Gridded(Linear()))
    tmin = minimum(t_data)
    tmax = maximum(t_data)
    # TODO: control or refine the # of interpolation points (1000 too much?)
    tspan = linspace(tmin, tmax, 500)
    # Extract the trajectory as an array of 3D points.
    ant_line = Point3f0[]
    for t in tspan
      x = Float32(scale_factor*x_spl[t])
      y = Float32(scale_factor*y_spl[t])
      i = indmin(abs(gx-x))
      j = indmin(abs(gy-y))
      z = dens[i,j]
      push!(ant_line, Point(gx[i], gy[j], z))
    end
    push!(ant_lines, dedup_consecutives(ant_line))
  end
  return ant_lines
end

function visualize_trajectory(ant_lines)
  max_traj = maximum(map(length, ant_lines))
  timesignal = preserve(loop((trace_length):max_traj))
  return preserve(map(timesignal) do t
      traj = Point3f0[]
      for i = 1:length(ant_lines)
        # Only add point if t is within trajectory bounds.
        # If t exceeds the bound, plot last element
        last_pos = (t <= length(ant_lines[i])) ? t : length(ant_lines[i])
        append!(traj, [ant_lines[i][last_pos]])
      end
      traj
    end)
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

iconsize = 8mm
knob_size = 5mm
icon_size_signal = Reactive.Signal(iconsize)
max_slider_length = 6 * iconsize

function get_slider_length(units)
  Measures.Length{:mm,Float64}(min(max_slider_length, units*knob_size/2))
end

ant_count_v, ant_count_s = labeled_slider(1:runs, edit_screen;
  slider_length = 8*iconsize,
  icon_size = icon_size_signal,
  knob_scale = knob_size)
dim1_v, dim1_s = labeled_slider(1:n, edit_screen;
  slider_length = 8*iconsize,
  icon_size = icon_size_signal,
  knob_scale = knob_size)
dim2_v, dim2_s = labeled_slider(1:n, edit_screen;
  slider_length = 8*iconsize,
  icon_size = icon_size_signal,
  knob_scale = knob_size)
scale_factor_v, scale_factor_s = labeled_slider(0.5:0.5:5.0, edit_screen;
  slider_length = 8*iconsize,
  icon_size = icon_size_signal,
  knob_scale = knob_size)

on_button_img = loadasset(string(assets_path, "on.png"))
off_button_img = loadasset(string(assets_path, "off.png"))
endpoint_v, endpoint_s = toggle_button(
  on_button_img, off_button_img, edit_screen)
log_dens_v, log_dens_s = toggle_button(
  on_button_img, off_button_img, edit_screen)

controls = Pair[
    "Endpoints only" => endpoint_v,
    "Negative log density" => log_dens_v,
    "Dimension 2" => dim2_v,
    "Dimension 1" => dim1_v,
    "Ant count" => ant_count_v,
    "XY scale factor" => scale_factor_v]

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

# Signal for the XY data used for landscaping.
XY_signal = map(endpoint_s, dim1_s, dim2_s) do is_endpoint, dim1, dim2
  getXYdata(data, is_endpoint, dim1, dim2)
end

# Important to use the grid and density as a single signal that updates at the
# same time. Computing ant lines replies on having the correct grid for heights.
surface_signal = map(XY_signal, log_dens_s, scale_factor_s) do xy, is_log, scale
  get_heights(xy[1], xy[2], is_log, scale)
end

red_color = RGBA(255.0, 0.0, 0.0, 1.0)

# Obtain signal for the sphere radius, since we want to scale the ant spheres to
# be proportional with the scale of the XY grid.
sphere_radius_s = const_lift(*, scale_factor_s, 0.05f0)
ant_sphere_s = map(sphere_radius_s) do sphere_radius
  GLNormalMesh(Sphere{Float32}(Vec3f0(0), sphere_radius))
end

# If the number of ants changes, we need to clear the window and re-render all
# objects, since the number of ants to re-render needs to be constant, even
# though positions can change.
traces_obj = map(ant_count_s) do ant_count
  ant_lines = const_lift(get_ant_lines, surface_signal, dim1_s, dim2_s,
    ant_count, scale_factor_s)
  # Get signal for ant position animation based on ant_lines and a time signal.
  ant_positions_s = map(visualize_trajectory, ant_lines)
  ant_positions_s = flatten(ant_positions_s,
    typ=Array{FixedSizeArrays.Point{3, Float32},1})

  visualize(
    (ant_sphere_s, ant_positions_s),
    boundingbox=nothing,
    color=red_color, camera=:perspective)
end
println("before surf render")

# Code to color wells.
include("landscape_colouring.jl")

# Separate the surface signal into x, y, z matrices for GLVisualize.
surf_obj = map(surface_signal) do surf
  gx = get_x_node_matrix(surf[1], surf[2])
  gy = get_y_node_matrix(surf[1], surf[2])
  dens = surf[3]

  # Prepare mesh vertex positions and texture.
  positions = Point3f0[Point3f0(gx[i,j], gy[i,j], dens[i,j])
    for i = 1:length(surf[1]) for j = 1:length(surf[2])]
  z_color, color_count = LandscapeColouring.color_landscape(surf[3])
  colors = Colors.colormap("Greens", color_count)
  texture = RGBA{Float32}[RGBA{Float32}(colors[z_color[i]])
    for i = 1:length(z_color)]

  # Plot mesh as vertices with specific colours.
  visualize((Circle, positions), color=texture, camera=:perspective,
    boundingbox=nothing)
  # Plot as smooth surface.
  #visualize((gx, gy, dens),
  #  :surface, camera=:perspective)
end

# Re-render every time the surface or number of ants changes.
preserve(map(surf_obj, traces_obj) do surf_obj, traces_obj
  empty!(view_screen)
  _view(surf_obj, view_screen, camera=:perspective)
  _view(traces_obj, view_screen, camera=:perspective)
end)

println("started surf render")
renderloop(window)
