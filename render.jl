module Renderer

include("landscape_colouring.jl")
include("partition_utils.jl")
using GLVisualize, GLAbstraction, ModernGL, Reactive, GeometryTypes, Colors
using Interpolations, GLWindow, ODEInput, CallbackRegistrar, KernelDensity
using DataFrames, FixedSizeArrays
using GLFW
import GLVisualize: labeled_slider, mm, button, toggle_button

# Boolean to indicate whether a first model was instantiated. Used for deciding
# whether the window should be showed when the user clicks "Run" (window is only
# hidden before the first model is instantiated), or whether to quit the entire
# application when the user clicks [X] on the ODEInput (only done before the
# first model is instantiated).
global model_instantiated_state = false

function model_instantiated()
  global model_instantiated_state
  return model_instantiated_state
end

function set_model_instantiated(is_instantiated)
  global model_instantiated_state = is_instantiated
end

# Returns the data for producing a landscape for the dimensions corresponding to
# dim1 and dim2. Returns either the entire trajectory or the endpoints only,
# depending on the value of is_endpoint.
function getXYdata(data, is_endpoint, dim1, dim2)
  runs = maximum(data[end])
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

# Calculates the landscape heights from X-Y simulation data. Returns a 2D array,
# either the probability density or the log density, depending on the is_log
# argument. It also returns the sampling grid for the X-Y axes.
function get_heights(X, Y, is_log, scale_factor)
  # Change to X, Y for alternative plotting.
  dens = kde((X, Y))
  dens_plus_background = 1e-23*ones(size(dens.density))+dens.density
  log_dens = -log(dens_plus_background)
  log_dens = log_dens-maximum(log_dens)
  gx = scale_factor*LinSpace(dens.x)
  gy = scale_factor*LinSpace(dens.y)

  # Convert everything from Float64 to Float32, as required by GLVisualize.
  dens = convert(Array{Float32, 2}, dens.density)
  log_dens = convert(Array{Float32, 2}, log_dens)
  gx = convert(LinSpace{Float32}, gx)
  gy = convert(LinSpace{Float32}, gy)
  return gx, gy, (is_log ? log_dens : dens)
end

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
function get_ant_lines(data, surf, dim1, dim2, ant_count, scale_factor)
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
  timesignal = preserve(loop(1:max_traj, 20#=rate=#))
  return map(timesignal) do t
      traj = Point3f0[]
      for i = 1:length(ant_lines)
        # Only add point if t is within trajectory bounds.
        # If t exceeds the bound, plot last element
        last_pos = (t <= length(ant_lines[i])) ? t : length(ant_lines[i])
        append!(traj, [ant_lines[i][last_pos]])
      end
      traj
    end
end

global edit_screen, view_screen, logo_screen, logoarea, window
global backup_actions_window = Dict()

# Custom renderloop similar to the one in GLVisualize.jl, but modified to
# allow correct integration of the rendering window with the ODEInput GUI.
function custom_renderloop(_window::Screen, framerate = 1//60)
  while isopen(_window)
    # Stop rendering once GUI for ODE Input is active.
    if (ODEInput.input_window_active())
      break
    end
    tic()
    render_frame(_window)
    swapbuffers(_window)
    poll_glfw()
    GLWindow.sleep_pessimistic(framerate - toq())
  end
  if (!isopen(_window))
    # If loop stop because X button is clicked, then destroy window.
    destroy!(_window)
  else
    # Otherwise, loop was interrupted because user opened ODEInput GUI.
    # Clear screen signals to make the screen appear "frozen" and back them up.
    println("Backing up, then killing window inputs")
    global backup_actions_window
    for (k, s) in _window.inputs
      backup_actions_window[k] = s.actions
      s.actions = Vector()
    end
  end
  return
end

# Function for unfreezing window with the same model once the user closed the
# ODE Input GUI.
# Restores the window input signals (e.g. clicks, scrolls etc.) that were killed
# and backed up when the ODE Input Window was opened.
function restore_screens()
  global backup_actions_window, window
  # Restore window actions.
  for (k, a) in backup_actions_window
    if haskey(window.inputs, k)
      window.inputs[k].actions = a
    end
  end
  # TODO: destroy GLWindow if application is quit on first show of ODEInput.
end

# Function for unfreezing window and preparing it for rendering a new model.
# Restores the window input signals (clicks, scrolls etc.) that were killed
# when the ODE Input window was opened. Also empties the view, logo and edit
# screens in order to prepare them for rendering a new model.
function reset_screens()
  global view_screen, edit_screen, logo_screen, window, backup_actions_window
  empty!(view_screen)
  empty!(edit_screen)
  empty!(logo_screen)
  # Window is hidden on startup, until the user inputs the first ODE.
  if (!model_instantiated())
    set_model_instantiated(true)
    GLWindow.show!(window)
  end
  global backup_actions_window
  # Restore window actions.
  for (k, a) in backup_actions_window
    if haskey(window.inputs, k)
      window.inputs[k].actions = a
    end
  end
end

function rerender(data, n, runs)
  global edit_screen, view_screen, logo_screen, logoarea, window
  iconsize = 6mm
  knob_size = 3mm
  icon_size_signal = Reactive.Signal(iconsize)

  # We allow the user to render a maximum of 100 ants, since beyond 100 the
  # model would be less responsive and the model would not be much more
  # informative.
  max_ants = min(runs, 100)
  ant_count_v, ant_count_s = labeled_slider(1:max_ants, edit_screen;
    slider_length = 4*iconsize,
    icon_size = icon_size_signal,
    knob_scale = knob_size,
    text_scale = 4mm)
  dim1_v, dim1_s = labeled_slider(1:n, edit_screen;
    slider_length = 4*iconsize,
    icon_size = icon_size_signal,
    knob_scale = knob_size,
    text_scale = 4mm)
  dim2_v, dim2_s = labeled_slider(1:n, edit_screen;
    slider_length = 4*iconsize,
    icon_size = icon_size_signal,
    knob_scale = knob_size,
    text_scale = 4mm)
  scale_factor_v, scale_factor_s = labeled_slider(0.5:0.5:5.0, edit_screen;
    slider_length = 4*iconsize,
    icon_size = icon_size_signal,
    knob_scale = knob_size,
    text_scale = 4mm)

  # current_dir = dirname(homedir())
  # assets_path = string(current_dir, "/Documents/stem-cells/assets/")

  # TODO: change this to generic path.
  assets_path = string(homedir(), "/Documents/stem-cells/assets/");

  on_button_img = loadasset(string(assets_path, "on.png"))
  off_button_img = loadasset(string(assets_path, "off.png"))
  logo_img = loadasset(string(assets_path, "waddle.png"))
  ode_input_img = loadasset(string(assets_path, "ode_input.png"))
  ode_input_v, ode_input_s = GLVisualize.button(ode_input_img, edit_screen)
  endpoint_v, endpoint_s = toggle_button(
    on_button_img, off_button_img, edit_screen)
  log_dens_v, log_dens_s = toggle_button(
    on_button_img, off_button_img, edit_screen)
  shading_v, shading_s = toggle_button(
    on_button_img, off_button_img, edit_screen)

  controls = Pair[
      "Endpoints only" => endpoint_v,
      "Log density" => log_dens_v,
      "Shading" => shading_v,
      "Dimension 2" => dim2_v,
      "Dimension 1" => dim1_v,
      "Ant count" => ant_count_v,
      "XY scale factor" => scale_factor_v,
      "DE input" => ode_input_v]

  edit_menu_robj = visualize(
    controls,
    text_scale = 3.5mm,
    gap = 3mm,
    width = 10iconsize)

  global menu_offset = 0
  # Manual implementation for scroll in edit menu (GLVisualize doesn't provide
  # it).
  preserve(map(edit_screen.inputs[:scroll]) do scroll
    global menu_offset, edit_screen
    scroll_stepsize = 10
    scroll_y = scroll[2]
    menu_offset_diff = 0.0
    if (!value(edit_screen.inputs[:mouseinside])) return end
    if (scroll_y > 0.0)
      # Scrolling upwards.
      edit_menu_height = widths(value(boundingbox(edit_menu_robj)))[2]
      edit_screen_height = value(edit_screen.area).h

      # Height difference between edit menu height and edit screen height. This
      # value indicates how much the user can scroll down.
      edit_screen_diff = edit_menu_height - edit_screen_height

      # Change p pixels at a time, until maximum offset is reached (or change in
      # amounts smaller than p is remaining amount is less than p pixels).
      menu_offset_diff = min(max(edit_screen_diff - menu_offset, 0),
        scroll_stepsize)
    elseif (scroll_y < 0.0)
      # Scrolling downwards.
      menu_offset_diff = -min(max(menu_offset, 0), scroll_stepsize)
    end

    GLAbstraction.translate!(edit_menu_robj, Vec3f0(0, -menu_offset_diff, 0))
    menu_offset = max(0, menu_offset + menu_offset_diff)
    return
  end)

  _view(edit_menu_robj, edit_screen, camera = :fixed_pixel)

  ode_input_button_callback = CallbackRegistrar.get_callback(:gui_popup)
  preserve(map(ode_input_button_callback, ode_input_s))

  logo_signal = map(logoarea) do a
    w_padding = 5
    h_padding = 10
    img_w = size(logo_img)[2]
    img_h = size(logo_img)[1]
    [Point2f0(w_padding + img_w/2,
      -h_padding + a.h - img_h/2)]
  end

  logo_vis = visualize(
    (SimpleRectangle(0, 0, size(logo_img)[2], size(logo_img)[1]), logo_signal);
    image=logo_img)
  _view(logo_vis, logo_screen, camera=:fixed_pixel)

  # Signal for the XY data used for landscaping.
  XY_signal = map(endpoint_s, dim1_s, dim2_s) do is_endpoint, d1, d2
    getXYdata(data, is_endpoint, d1, d2)
  end

  # Important to use the grid and density as a single signal that updates at the
  # same time. Computing ant lines replies on having the correct grid for
  # heights.
  surface_signal = map(XY_signal, log_dens_s, scale_factor_s) do xy, is_log,
      scale
    get_heights(xy[1], xy[2], is_log, scale)
  end

  red_color = RGBA(255.0, 0.0, 0.0, 1.0)

  # Obtain signal for the sphere radius, since we want to scale the ant spheres
  # to be proportional with the scale of the XY grid.
  sphere_radius_s = const_lift(*, scale_factor_s, 0.05f0)
  ant_sphere_s = map(sphere_radius_s) do sphere_radius
    GLNormalMesh(Sphere{Float32}(Vec3f0(0), sphere_radius))
  end

  # If the number of ants changes, we need to clear the window and re-render all
  # objects, since the number of ants to re-render needs to be constant, even
  # though positions can change.
  traces_obj = map(ant_count_s) do ant_count
    ant_lines = map(surface_signal) do surf
        # Only update when surface changes, not just when data signal changes.
        # Same for dimension values and scale, since surf depends on all of
        # them.
        d1 = value(dim1_s)
        d2 = value(dim2_s)
        scale = value(scale_factor_s)
        get_ant_lines(data, surf, d1, d2, ant_count, scale)
    end
    # Get signal for ant position animation based on ant_lines and a time
    # signal.
    ant_positions_s = map(visualize_trajectory, ant_lines)
    ant_positions_s = flatten(ant_positions_s,
      typ=Array{FixedSizeArrays.Point{3, Float32},1})

    visualize(
      (ant_sphere_s, ant_positions_s),
      boundingbox=nothing,
      color=red_color)
  end

  # Separate the surface signal into x, y, z matrices for GLVisualize.
  surf_obj = map(surface_signal, shading_s) do surf, is_shaded
    gx = get_x_node_matrix(surf[1], surf[2])
    gy = get_y_node_matrix(surf[1], surf[2])
    dens = surf[3]

    # Prepare mesh vertex positions and texture.
    positions = Point3f0[Point3f0(gx[i,j], gy[i,j], dens[i,j])
      for i = 1:length(surf[1]) for j = 1:length(surf[2])]
    z_color, color_count = LandscapeColouring.color_landscape(surf[3],
      Reactive.value(log_dens_s)) # looking for minima if log_dens_s is true

    # When shading is not on, introduce some transparency.
    transparency = is_shaded ? 1.0 : 0.8
    # Get 'Rainbow' colorscheme
    colors = RGBA{Float32}[
      RGBA(
          clamp(min(4x - 1.5, -4x + 4.5) ,0.0,1.0),
          clamp(min(4x - 0.5, -4x + 3.5) ,0.0,1.0),
          clamp(min(4x + 0.5, -4x + 2.5) ,0.0,1.0), transparency)
      for x in linspace(0.0,1.0, color_count)]
    shuffle!(colors)
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
    visualize((gx, gy, dens), color=texture, :surface, shading=is_shaded)
  end

  # Re-render every time the surface or number of ants changes.
  preserve(map(surf_obj, traces_obj) do surf_obj, traces_obj
    empty!(view_screen)
    _view(surf_obj, view_screen, camera=:perspective)
    _view(traces_obj, view_screen, camera=:perspective)
  end)

  @async custom_renderloop(window)
end

function init()
  global logoarea
  global window = glscreen(resolution = primarymonitorresolution(),
    visible=false)

  #GLWindow.hide!(window)
  # Create partitioned window for controls and view screens.
  editarea, viewarea = x_partition_abs(window.area, 180)
  # Further partition edit area to get a logo area.
  editarea, logoarea = y_partition_fixed_top(editarea, 80)
  global edit_screen = Screen(
      window, area = editarea,
      color = RGBA{Float32}(0.0f0, 0.0f0, 0.0f0, 1f0))
  global view_screen = Screen(
      window, area = viewarea,
      color = RGBA(255.0f0, 255.0f0, 255.0f0, 1f0),
      stroke = (1f0, RGBA{Float32}(0.13f0, 0.13f0, 0.13f0, 13f0)))
  global logo_screen = Screen(
      window, area = logoarea,
      color = RGBA{Float32}(0.0f0, 0.0f0, 0.0f0, 1f0))

  # Callback for when the [X] button is clicked on the ODE input window.
  close_input_callback = function(path)
    if (!ODEInput.input_window_active())
      # Do not handle window closing if it has been handled already. This
      # callback may be called many times when clicking the [X] button, and also
      # whenever the "Run" button is clicked. We want to run the callback  only
      # once, when [X] button is clicked.
      return
    end
    ODEInput.set_input_window_active(false)
    global window
    if (!model_instantiated())
      # User clicked [X] before having rendered any model. All the windows
      # should be destroyed.
      GLWindow.show!(window)
      @async custom_renderloop(window)
      if (!isopen(window)) println("NOT OPEN") else println("OPEN") end
      GLFW.SetWindowShouldClose(GLWindow.nativewindow(window), true)
      if (!isopen(window)) println("NOT OPEN") else println("OPEN") end
      println("Quitting entire application.")
      GLWindow.destroy!(window)
      return
    end
    restore_screens()
    # Spawn this async so that the window closing process finishes.
    @async custom_renderloop(window)
  end
  CallbackRegistrar.register_callback(:close_input, close_input_callback)
end

end # module Renderer
