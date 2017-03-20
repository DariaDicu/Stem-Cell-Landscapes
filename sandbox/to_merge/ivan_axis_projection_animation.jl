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

  vol = 0

  F1 = (x1, x2) ->
    (a*(x1^n)/(S^n + x1^n) + b*S^n/(S^n + x2^n) - k*x1 + vol*randn())
  F2 = (x1, x2) ->
    (a*(x2^n)/(S^n + x2^n) + b*S^n/(S^n + x1^n) - k*x2 + vol*randn())
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end

runs = 20 # Number of simulation runs
n = 2 # Number of dimensions
data = build_landscape(runs, F, 2, (0,5))

using KernelDensity, Interpolations, DataFrames, Reactive;

# TODO: Button for this.
#scale_factor_s = Signal(5.0) # Scale for xy grid.
# TODO: Button for this.
# Factor has to be between 0 and 1!
ant_speed_factor = 0.1f0
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

  x_linspaced = gx
  y_linspaced = gy
  global x_linspaced
  global y_linspaced

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
  timesignal = preserve(loop(1:max_traj))
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
assets_path = string(homedir(), "/Documents/GitHub/Stem-Cell-Landscapes/assets/");

include("partition_utils.jl")
# Create partitioned window for controls and view screens.
editarea, viewcontainer = x_partition_abs(window.area, 180)
# Further partition edit area to get a logo area.
editarea, logoarea = y_partition(editarea, 85)
# Further partition view area to get a summary area.
viewarea, summarycontainer = x_partition_fixed_right(viewcontainer, 350)
# Further partition summary container to summary area and phase area
writingcontainer, phasearea = y_partition(summarycontainer, 50)
# Further partition writing container to summary area and ode area
summaryarea, ODEarea = y_partition(writingcontainer, 75)

# Explain this, TODO
function translate_area(container, component)
  SimpleRectangle(component.x+container.x,
    component.y+container.y,
    component.w,
    component.h)
end

# Explain maybe?
#viewarea = map(viewarea) do viewarea
#  translate_area(value(viewcontainer), viewarea)
#end

#summarycontainer = map(summarycontainer) do summarycontainer
#  translate_area(value(viewcontainer), summarycontainer)
#end

writingcontainer = map(writingcontainer) do writingcontainer
  translate_area(value(summarycontainer), writingcontainer)
end

phasearea = map(phasearea) do phasearea
  translate_area(value(summarycontainer), phasearea)
end

ODEarea = map(ODEarea) do ODEarea
  translate_area(value(writingcontainer), ODEarea)
end

summaryarea = map(summaryarea) do summaryarea
  translate_area(value(writingcontainer), summaryarea)
end


# Set up all of the screens, based on the above defined areas
edit_screen = Screen(
    window, area = editarea,
    color = RGBA{Float32}(0.0f0, 0.0f0, 0.0f0, 1f0))
view_screen = Screen(
    window, area = viewarea,
    color = RGBA(255.0f0, 255.0f0, 255.0f0, 1f0),
    stroke = (1f0, RGBA{Float32}(0.13f0, 0.13f0, 0.13f0, 13f0)))
logo_screen = Screen(
    window, area = logoarea,
    color = RGBA{Float32}(0.0f0, 0.0f0, 0.0f0, 1f0))
summary_screen = Screen(
    window, area = summaryarea,
    name = Symbol("__Plots.jl"),
    color = RGBA{Float32}(0.0f0, 0.0f0, 0.0f0, 1f0))
phase_screen = Screen(
    window, area = phasearea,
    color = RGBA{Float32}(0.0f0, 0.0f0, 0.0f0, 1f0))
println("PHASE SCREEN INPUTS")
ODE_screen = Screen(
    window, area = ODEarea,
    color = RGBA{Float32}(0.1f0, 0.1f0, 0.1f0, 1f0))
phase_screen.inputs[:mouse_hover] = phase_screen.inputs[:mouseinside]


iconsize = 8mm
knob_size = 5mm
icon_size_signal = Reactive.Signal(iconsize)

#TODO: delete this
#function get_slider_length(units)
#  Measures.Length{:mm,Float64}(min(max_slider_length, units*knob_size/2))
#end

ant_count_v, ant_count_s = labeled_slider(1:runs, edit_screen;
  slider_length = 8*iconsize,
  icon_size = icon_size_signal,
  knob_scale = knob_size)
dim1_v, dim1_s = labeled_slider(1:n, edit_screen;
  slider_length = 4*iconsize,
  icon_size = icon_size_signal,
  knob_scale = knob_size)
dim2_v, dim2_s = labeled_slider(1:n, edit_screen;
  slider_length = 4*iconsize,
  icon_size = icon_size_signal,
  knob_scale = knob_size)
scale_factor_v, scale_factor_s = labeled_slider(0.5:0.5:5.0, edit_screen;
  slider_length = 8*iconsize,
  icon_size = icon_size_signal,
  knob_scale = knob_size)

on_button_img = loadasset(string(assets_path, "on.png"))
off_button_img = loadasset(string(assets_path, "off.png"))
logo_img = loadasset(string(assets_path, "waddle.png"))
endpoint_v, endpoint_s = toggle_button(
  on_button_img, off_button_img, edit_screen)
log_dens_v, log_dens_s = toggle_button(
  on_button_img, off_button_img, edit_screen)
shading_v, shading_s = toggle_button(
  on_button_img, off_button_img, edit_screen)

controls = Pair[
    "Endpoints only" => endpoint_v,
    "Negative log density" => log_dens_v,
    "Shading" => shading_v,
    "Dimension 2" => dim2_v,
    "Dimension 1" => dim1_v,
    "Ant count" => ant_count_v,
    "XY scale factor" => scale_factor_v]

_view(visualize(
        controls,
        text_scale = 5mm,
        gap = 3mm,
        width = 10iconsize), edit_screen, camera = :fixed_pixel)


size(logo_img)

logo_signal = map(logoarea) do a
  padding = 10
  img_w = size(logo_img)[1]
  img_h = size(logo_img)[2]
  [Point2f0(padding + img_w/2,
    -padding + a.h - img_h/2)]
end

logo_vis = visualize((SimpleRectangle(0,0,size(logo_img)[1], size(logo_img)[2]),
  logo_signal), image=logo_img)
println(value(logoarea))
_view(logo_vis, logo_screen, camera=:fixed_pixel)

#logo_text = visualize(
#    "Waddle",
#    relative_scale=16mm,
#    color = RGBA(1f0, 1f0, 1f0, 1f0))
#_view(logo_vis, logo_screen, camera=:fixed_pixel)
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
    color=red_color)
end

# Code to color wells.
include("landscape_colouring.jl")

using DifferentialEquations, Plots

# Separate the surface signal into x, y, z matrices for GLVisualize.
surf_obj = map(surface_signal, shading_s) do surf, is_shaded
  gx = get_x_node_matrix(surf[1], surf[2])
  gy = get_y_node_matrix(surf[1], surf[2])
  dens = surf[3]

  # Make dens global to plot the phase portrait later
  # (required by phase_generator.jl)
  global dens

  # Prepare mesh vertex positions and texture.
  positions = Point3f0[Point3f0(gx[i,j], gy[i,j], dens[i,j])
    for i = 1:length(surf[1]) for j = 1:length(surf[2])]
  z_color, color_count = LandscapeColouring.color_landscape(surf[3],
    value(log_dens_s)) # looking for minima if log_dens_s is true

  # When shading is not on, introduce some transparency.

  transparency = is_shaded ? 1.0 : 0.8
  # Get 'Rainbow' colorscheme
  colors = RGBA{Float32}[
    RGBA(
        clamp(min(4x - 1.5, -4x + 4.5) ,0.0,1.0),
        clamp(min(4x - 0.5, -4x + 3.5) ,0.0,1.0),
        clamp(min(4x + 0.5, -4x + 2.5) ,0.0,1.0), transparency)
    for x in linspace(0.0,1.0, color_count)]


  #texture = RGBA{Float32}[RGBA{Float32}(colors[z_color[i]])
  #  for i = 1:length(z_color)]

  l1 = size(z_color)[1]
  l2 = size(z_color)[2]

  texture = map(c->colors[c], z_color)
  # Plot mesh as vertices with specific colours.
  #visualize((Circle, positions), boundingbox=nothing)
  # Plot as smooth surface with colored wells.

  view_screen.color = is_shaded ?
    RGBA{Float32}(255.0,255.0,255.0,1.0) :
    RGBA{Float32}(0.0,0.0,0.0,1.0)
  view_screen.stroke = is_shaded ?
    (1f0, RGBA{Float32}(255.0,255.0,255.0,1.0)) :
    (1f0, RGBA{Float32}(0.13f0, 0.13f0, 0.13f0, 13f0))
  visualize((gx, gy, dens), color=texture, :surface, shading=is_shaded)#
end


# Get descriptive text by accessing functions from summary_script.jl
include("summary_script.jl")

summarytext = ("Summary information: \n \n")
summarytext = summarytext * minima_info_1(hunt_minima(dens))

# Display this information in the right-hand summary_screen
summaryarea = map(summaryarea) do summaryarea
  empty!(summary_screen)
  _view(visualize(
          summarytext,
          text_scale = 10mm,
          gap = 3mm,
          width = 10iconsize,
          model=translationmatrix(Vec3f0(10, (value(summaryarea).h)-20, 0))),
          summary_screen, camera = :fixed_pixel)
end


# Include the phase portrait generated from dens in the summary panel
include("phase_generator.jl")

# Generate the phase diagram and save result as a .png
# (This was the only way I got around visualising it in GLVisualize)
portrait = phase_portrait_gen(dens, 20, 0.1, true)
portrait
savefig("portrait.png")

using FileIO
using Images


# Each time phasearea changes shape, scale contents appropriately
phasearea = map(phasearea) do phasearea
  empty!(phase_screen)

  # Load in the png
  portrait = load("portrait.png")
  # Get the current pixel dimensions for the phasearea section of the screen
  phasedims = (value(phasearea).w, value(phasearea).h)

  # Reshape the image to make sure its square and fits in phasearea
  portrait = Images.imresize(portrait,((minimum(phasedims)-20),(minimum(phasedims)-20)))

  # If phase area is taller than it is wide...
  if phasedims[1]<phasedims[2]
    portrait_obj = visualize(
            portrait,
            text_scale = 10mm,
            gap = 3mm,
            width = 10iconsize,
            model=translationmatrix(Vec3f0(10,((phasedims[2]-(minimum(phasedims)-20))/2),0)))
  # Else phase area is wider than it is tall...
  else
    portrait_obj = visualize(
            portrait,
            text_scale = 10mm,
            gap = 3mm,
            width = 10iconsize,
            model=translationmatrix(Vec3f0(((phasedims[1]-(minimum(phasedims)-20))/2),10,0)))
  end
  # View this appropriately sized, renderable object in phase_screen
  _view(portrait_obj, phase_screen, camera = :fixed_pixel)
end


# Add information to the ODE area
ODE_content = "ODE content here"

# Visualise this information in the ODE screen
ODEarea = map(ODEarea) do ODEarea
  empty!(ODE_screen)
  _view(visualize(
          ODE_content,
          text_scale = 10mm,
          gap = 3mm,
          width = 10iconsize,
          model=translationmatrix(Vec3f0(10, (value(ODEarea).h)-20, 0))),
          ODE_screen, camera = :fixed_pixel)
end

# Re-render every time the surface or number of ants changes.
preserve(map(surf_obj, traces_obj) do surf_obj, traces_obj
  empty!(view_screen)
  _view(surf_obj, view_screen, camera=:perspective)
  _view(traces_obj, view_screen, camera=:perspective)
end)


renderloop(window)
