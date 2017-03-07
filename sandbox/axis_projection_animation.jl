# Code for generating a Julia OpenGL model that plots the landscape for a fixed
# value of parameter a and animated moving "ant" trajectories onto the model
# from a subset of the simulations.
include("ode_simulator.jl")
using ODESimulator, KernelDensity, Interpolations, DataFrames, Reactive,HDF5,LaTeXStrings;
## Example of widgets put into container with change handler assigned

using Tk
using Compat; import Compat.String
#using ComplexPhasePortrait

function readVari(seq)
 readV_res = []
 vari = split(seq,",")
 for i=1:length(vari)
   push!(readV_res,replace(vari[i]," ",""))
 end
 return readV_res
end

function readIni(seq)
  resul = split(replace(seq," ",""),"=")
  s_deal=replace(resul[length(resul)],"(","")
  s_deal=replace(s_deal,")","")
  s_deal=replace(s_deal,"[","")
  s_deal=replace(s_deal,"]","")
  s_dealMore=split(s_deal,",")
  if contains(s_dealMore[1],".")==0
     s_dealMore[1]=s_dealMore[1] * ".0"
  end
  if contains(s_dealMore[2],".")==0
     s_dealMore[2]=s_dealMore[2] * ".0"
  end
  return "(" * s_dealMore[1] * "," * s_dealMore[2] * ")"
end

function readPara(seq)
   seq=replace(replace(seq," ",""),",","\n")
   return seq
end

function readInterv(seq)
  resul = split(replace(seq," ",""),"=")
  s_deal=replace(resul[length(resul)],"(","")
  s_deal=replace(s_deal,")","")
  s_dealMore=split(s_deal,",")
  if contains(s_dealMore[1],".")==0
     s_dealMore[1]=s_dealMore[1] * ".0"
  end
  if contains(s_dealMore[2],".")==0
     s_dealMore[2]=s_dealMore[2] * ".0"
  end
  return "(" * s_dealMore[1] * "," * s_dealMore[2] * ")"
end

function subinfunc(seq)
   funclist=split(seq,"@")
end

function replace_plus(seq,rpstr,tostr)
  k=searchindex(seq,rpstr)
  while k!=0
     if (k>1)
        if ((k-1+length(rpstr))==length(seq))
          if ((isalpha(seq[k-1]) == false) & (seq[k-1] != "_"))
             seq=replace(seq,rpstr,"#T#",1)
           else
             seq=replace(seq,rpstr,"#F#",1)
          end
        else
          if ((isalpha(seq[k+length(rpstr)]) == false) & (seq[k+length(rpstr)] != "_") & (isalpha(seq[k-1]) == false) & (seq[k-1] != "_"))
            seq=replace(seq,rpstr,"#T#",1)
          else
            seq=replace(seq,rpstr,"#F#",1)
          end
        end
      else
        if ((isalpha(seq[k+length(rpstr)] )== false) & (seq[k+length(rpstr)] != "_"))
           seq=replace(seq,rpstr,"#T#",1)
         else
           seq=replace(seq,rpstr,"#F#",1)
        end
     end
   k=searchindex(seq,rpstr)
  end

  seq=replace(seq,"#T#",tostr)
  seq=replace(seq,"#F#",rpstr)
  return seq
end

function reformatEq(dVariable,eqSet)
 allEqs=[]
 strDeal=""
 strDealLR = []
 combStr=""
 combFunc=""
 t_sym=[]
 eqtReads=[]
 dt="t"
 for i = 1 : length(eqSet)
    t_sym=[]
    eqtReads=[]
    strDeal = eqSet[i]
    eqtReads=split(strDeal,"=")
    t_sym=split(eqtReads[1],"/")
    if length(t_sym)>1
       dt = t_sym[2]
       combStr = combStr * replace(strDeal," ","") * "@"
     else
       combFunc = combFunc * replace(strDeal," ","") * "@"
    end
 end
 combStr=combStr[1:(length(combStr)-1)]
 for c = 1:length(dVariable)
    combStr=replace_plus(combStr,dVariable[c],"u[" * string(c) * "]")
    combStr=replace_plus(combStr,"d"*dVariable[c],"du[" * string(c) * "]")
    combStr=replace_plus(combStr,dt,"dt")
    combStr=replace_plus(combStr,"/dt","")
 end
 print(replace(combFunc,"@","\n"))
 print(replace(combStr,"@","\n"))
 return replace(combFunc,"@","\n") * replace(combStr,"@","\n")
end

w = Toplevel("DE Input", true)

## pack in tk frame for themed widgets
f = Frame(w)
configure(f, @compat Dict(:padding => [3,3,2,2], :relief=>"groove"));pack(f, expand=true, fill="both")
#tcl("pack", "propagate", w, false)
## widgets

Vari=Entry(f,"x,y") #The variables
Para=Entry(f,"a=1,b=2")
Bound=Entry(f,"b=(0,5)")
TimeR=Entry(f,"t=(0,1000)")
Itera=Entry(f,"1000")
DLabel=Label(f,"Differential Equations:")
Equations=Entry(f,"dx/dt=x+y dy/dt=x*y")
eqt=Text(f)
widgets = (Vari,Para,Bound,TimeR,Itera,DLabel,eqt)
pack_style = ["pack", "grid", "formlayout"][3]

    ## second argument is Tk_Radio instance
b_Save = Button(w, "Save")
pack(b_Save, expand=true, fill="both")
b_Load = Button(w, "Load")
pack(b_Load, expand=true, fill="both")
b_Run = Button(w, "Run")
pack(b_Run, expand=true, fill="both")

function savef(path)
 sv_path=GetSaveFile()
 #sv_path=  sv_path * ".de"
 str_vari = get_value(Vari)
 str_para = get_value(Para)
 str_bound = get_value(Bound)
 str_time = get_value(TimeR)
 str_eqt = get_value(eqt)
 str_runs= get_value(Itera)
 data=str_vari * "&" * str_para * "&" * str_bound * "&" * str_time * "&" * str_eqt * "&" * str_runs

 h5open(sv_path, "w") do file
    write(file, "data", data)  # alternatively, say "@write file A"
 end
 print(sv_path)

end

function loadf(path)
  op_path = GetOpenFile()
  print(op_path)
  d = h5open(op_path, "r") do file
    read(file, "data")
  end
  print(d)
  datas=split(d,"&")
  rd_para = datas[2]
  rd_vari = datas[1]
  rd_bound = datas[3]
  rd_tspan = datas[4]
  rd_eqt = datas[5]
  rd_run = datas[6]
  set_value(Vari, rd_vari )
  set_value(Para, rd_para)
  set_value(Bound, rd_bound)
  set_value(TimeR, rd_tspan)
  set_value(eqt,rd_eqt)
  set_value(Itera,rd_run)
end

function callback(path)
  vals = map(get_value, (cb, rb))
  println(vals)
end

bind(b_Save, "command",savef)   ## generic way to add callback for most common event
bind(b_Load,  "command",loadf)

if pack_style == "pack"
    map(pack, widgets)
    map(u -> pack_configure(u, @compat Dict(:anchor => "w")), widgets)
elseif pack_style == "grid"
    for i in 1:length(widgets)
        grid(widgets[i], i, 1)
        grid_configure(widgets[i], @compat Dict(:sticky => "we"))
    end
else
    formlayout(Vari,"Variables:")
    formlayout(Para,"Parameters:")
    formlayout(Bound,"Boundary:")
    formlayout(TimeR,"TimeRange:")
    formlayout(Itera,"Iterations:")
    formlayout(eqt,"Equations",)
    #formlayout(Save,"")
end

function CreatModel(path)
    println("PASS THE TEST")
    str_vari = get_value(Vari)
    str_para = get_value(Para)
    str_bound = get_value(Bound)
    str_time = get_value(TimeR)
    str_eqt = get_value(eqt)
    input_runs = get_value(Itera)
    input_vari = readVari(str_vari)
    input_para = readPara(str_para)
    input_time = readInterv(str_time)
    input_bound = readIni(str_bound)
    input_eqt = reformatEq(input_vari,split(str_eqt,"\n"))
    func = "function(" * "t,u,du" * ")" * "\n" * input_para *"\n" * input_eqt * "\n" * "end"
    fc=eval(parse(func))
    set_func(fc)
    set_visible(w,false)
end

bind(b_Run, "command",  CreatModel)
## bind a callback to each widget

set_visible(w,false)

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

using DataFrames
runs = 100 # Number of simulation runs
n = 2 # Number of dimensions

using Reactive
new_data = ODESimulator.build_landscape(runs, F, 2, (0,5))
data_s = Signal(new_data)

function set_func(func)
  global bool_click = true
  global runs
  new_data = ODESimulator.build_landscape(runs, func , 2, (0,5))
  global data_s
  println(new_data)
  push!(data_s, new_data)
end

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
using GLVisualize, GLAbstraction, ModernGL, Reactive, GeometryTypes, Colors
using Interpolations, GLWindow
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

bool_click = true
window = glscreen()
iconsize = 8mm
assets_path = string(homedir(), "/Documents/Stem-Latest/assets/");

# Create partitioned window for controls and view screens.
editarea, viewarea = x_partition_abs(window.area, 180)
# Further partition edit area to get a logo area.
editarea, logoarea = y_partition(editarea, 85)

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

iconsize = 8mm
knob_size = 5mm
icon_size_signal = Reactive.Signal(iconsize)

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
button_img = loadasset(string(assets_path, "Input.png"))
button_obj, button_s = GLVisualize.button(
  button_img, edit_screen)
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
    "XY scale factor" => scale_factor_v,
    "DE input" => button_obj]

_view(visualize(
        controls,
        text_scale = 5mm,
        gap = 3mm,
        width = 10iconsize), edit_screen, camera = :fixed_pixel)

size(logo_img)

click_sig = map(button_s) do clicked
   if clicked
     global bool_click=false
     set_visible(w,true)
   end
end

logo_signal = map(logoarea) do a
  padding = 10
  img_w = size(logo_img)[1]
  img_h = size(logo_img)[2]
  [Point2f0(padding + img_w/2,
    -padding + a.h - img_h/2)]
end



logo_vis = visualize((SimpleRectangle(0,0,size(logo_img)[1], size(logo_img)[2]),
  logo_signal), image=logo_img)
_view(logo_vis, logo_screen, camera=:fixed_pixel)

#logo_text = visualize(
#    "Waddle",
#    relative_scale=16mm,
#    color = RGBA(1f0, 1f0, 1f0, 1f0))
#_view(logo_vis, logo_screen, camera=:fixed_pixel)
########### Done setting up sidebar. ##########

# Signal for the XY data used for landscaping.
XY_signal = map(endpoint_s, dim1_s, dim2_s, data_s) do is_endpoint, d1, d2, data
  getXYdata(data, is_endpoint, d1, d2)
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
  ant_lines = map(surface_signal, dim1_s, dim2_s,
    scale_factor_s) do surf, d1, d2, scale
      # Only update when surface changes, not just when data signal changes.
      data = value(data_s)
      get_ant_lines(data, surf, d1, d2, ant_count, scale)
  end
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


function My_renderloop(window::Screen, framerate = 1//60)
    global bool_click
    while isopen(window)
      if bool_click ==true
        tic()
        render_frame(window)
        swapbuffers(window)
        poll_glfw()
        yield()
        GLWindow.sleep_pessimistic(framerate - toq())
      else #
        tic()
        render_frame(window)
        swapbuffers(window)
        poll_glfw()
        yield()
        GLWindow.sleep_pessimistic((framerate - toq())*100)
      end
    end
    destroy!(window)
    return
end

My_renderloop(window)
