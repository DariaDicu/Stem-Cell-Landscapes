include("callback_registrar.jl")
module ODEInput

using Tk, HDF5, LaTeXStrings, Compat
using CallbackRegistrar
import Compat.String

global input_window_active_state = false

function set_input_window_active(is_active)
  global input_window_active_state = is_active
end

function input_window_active()
  global input_window_active_state
  return input_window_active_state
end

function read_variables(seq)
  readV_res = []
  vari = split(seq,",")
  for i = 1:length(vari)
   push!(readV_res, replace(vari[i]," ",""))
  end
  return readV_res
end

function read_bounds(seq)
  resul = split(replace(seq, " ", ""), "=")
  s_deal = replace(resul[length(resul)], "(", "")
  s_deal = replace(s_deal, ")", "")
  s_deal = replace(s_deal, "[", "")
  s_deal = replace(s_deal, "]", "")
  s_dealMore = split(s_deal,",")
  if (!contains(s_dealMore[1], "."))
    s_dealMore[1] = s_dealMore[1] * ".0"
  end
  if (!contains(s_dealMore[2], "."))
    s_dealMore[2] = s_dealMore[2] * ".0"
  end
  str = "(" * s_dealMore[1] * "," * s_dealMore[2] * ")"
  return eval(parse(str))
end

function read_parameters(seq)
  seq = replace(replace(seq," ",""), ",", "\n")
  return seq
end

function read_timespan(seq)
  resul = split(replace(seq," ",""),"=")
  s_deal = replace(resul[length(resul)],"(","")
  s_deal = replace(s_deal,")","")
  s_dealMore = split(s_deal,",")
  if (!contains(s_dealMore[1],"."))
    s_dealMore[1] = s_dealMore[1] * ".0"
  end
  if (!contains(s_dealMore[2],"."))
    s_dealMore[2] = s_dealMore[2] * ".0"
  end
  str = "(" * s_dealMore[1] * "," * s_dealMore[2] * ")"
  return eval(parse(str))
end

function replace_plus(seq,rpstr,tostr)
  k = searchindex(seq,rpstr)
  while (k != 0)
     if (k > 1)
        if ((k-1+length(rpstr)) == length(seq))
          if ((isalpha(seq[k-1]) == false) && (seq[k-1] != "_"))
             seq=replace(seq,rpstr,"#T#",1)
           else
             seq=replace(seq,rpstr,"#F#",1)
          end
        else
          if ((isalpha(seq[k+length(rpstr)]) == false) &&
            (seq[k+length(rpstr)] != "_") && (isalpha(seq[k-1]) == false) &&
            (seq[k-1] != "_"))
            seq = replace(seq, rpstr, "#T#", 1)
          else
            seq = replace(seq, rpstr, "#F#", 1)
          end
        end
      else
        if ((isalpha(seq[k+length(rpstr)]) == false) &&
          (seq[k+length(rpstr)] != "_"))
           seq = replace(seq, rpstr, "#T#", 1)
         else
           seq = replace(seq, rpstr, "#F#", 1)
        end
     end
   k = searchindex(seq, rpstr)
  end

  seq = replace(seq, "#T#", tostr)
  seq = replace(seq, "#F#", rpstr)
  return seq
end

function reformatEq(dVariable, eqSet)
  allEqs=[]
  strDeal=""
  strDealLR = []
  combStr=""
  combFunc=""
  t_sym=[]
  eqtReads=[]
  dt="t"
  for i = 1:length(eqSet)
    t_sym = []
    eqtReads = []
    strDeal = eqSet[i]
    eqtReads = split(strDeal, "=")
    t_sym = split(eqtReads[1], "/")
    if (length(t_sym) > 1)
      dt = t_sym[2]
      combStr = combStr * replace(strDeal, " ", "") * "@"
    else
      combFunc = combFunc * replace(strDeal, " ", "") * "@"
    end
  end
  combStr = combStr[1:(length(combStr)-1)]
  for c = 1:length(dVariable)
    combStr = replace_plus(combStr, dVariable[c], "u[" * string(c) * "]")
    combStr = replace_plus(combStr, "d" * dVariable[c], "du[" * string(c) * "]")
    combStr = replace_plus(combStr, dt, "dt")
    combStr = replace_plus(combStr, "/dt", "")
  end
  print(replace(combFunc,"@","\n"))
  print(replace(combStr,"@","\n"))
  return replace(combFunc,"@","\n") * replace(combStr,"@","\n")
end

function save_callback(path)
  global variables_textbox, parameters_textbox, bounds_textbox, time_textbox,
    iterations_textbox, equations_textbox
  sv_path = GetSaveFile()
  #sv_path=  sv_path * ".de"
  str_vari = get_value(variables_textbox)
  str_para = get_value(parameters_textbox)
  str_bound = get_value(bounds_textbox)
  str_time = get_value(time_textbox)
  str_eqt = get_value(equations_textbox)
  str_runs = get_value(iterations_textbox)
  data = str_vari * "&" * str_para * "&" * str_bound * "&" * str_time * "&" *
    str_eqt * "&" * str_runs

  h5open(sv_path, "w") do file
    write(file, "data", data)  # alternatively, say "@write file A"
  end
  print(sv_path)
end

# TODO: Test Load/Save callbacks since I've changes var names.
function load_callback(path)
  global variables_textbox, parameters_textbox, bounds_textbox, time_textbox,
    iterations_textbox, equations_textbox
  op_path = GetOpenFile()
  print(op_path)
  d = h5open(op_path, "r") do file
    read(file, "data")
  end
  print(d)
  datas = split(d,"&")
  rd_vari = datas[1]
  rd_para = datas[2]
  rd_bound = datas[3]
  rd_tspan = datas[4]
  rd_eqt = datas[5]
  rd_run = datas[6]
  set_value(variables_textbox, rd_vari)
  set_value(parameters_textbox, rd_para)
  set_value(bounds_textbox, rd_bound)
  set_value(time_textbox, rd_tspan)
  set_value(equations_textbox, rd_eqt)
  set_value(iterations_textbox, rd_run)
end

function run_callback(path)
  global variables_textbox, parameters_textbox, bounds_textbox, time_textbox,
    iterations_textbox, equations_textbox, ode_input_window
  input_variables = read_variables(get_value(variables_textbox))
  input_parameters = read_parameters(get_value(parameters_textbox))
  input_bounds = read_bounds(get_value(bounds_textbox))
  input_time = read_timespan(get_value(time_textbox))
  input_equations = reformatEq(input_variables, split(
  get_value(equations_textbox), "\n"))
  input_runs = parse(get_value(iterations_textbox))
  function_string = "function(" * "t,u,du" * ")" * "\n" * input_parameters *
  "\n" * input_equations * "\n" * "end"
  parsed_function = eval(parse(function_string))
  create_model_callback = CallbackRegistrar.get_callback(:create_model)
  #TODO: replace get_value(variables_textbox) with smth more sensible
  create_model_callback(parsed_function, length(split(get_value(
    variables_textbox),",")), input_bounds, input_runs, input_time)
end

function open_gui_window()
  set_input_window_active(true)
  global ode_input_window = Toplevel("DE Input", true)
  global ode_input_frame = Frame(ode_input_window)
  configure(ode_input_frame,
    @compat Dict(:padding => [3,3,2,2], :relief=>"groove"));
  pack(ode_input_frame, expand=true, fill="both")

  global variables_textbox = Entry(ode_input_frame, "x, y")
  global parameters_textbox = Entry(ode_input_frame, "a = 1, b = 2")
  global bounds_textbox = Entry(ode_input_frame, "b = (0,5)")
  global time_textbox = Entry(ode_input_frame, "t = (0,50)")
  global iterations_textbox = Entry(ode_input_frame, "1000")
  global equations_label = Label(ode_input_frame,
    "Differential Equations:")
  global equations_textbox  = Text(ode_input_frame)
  # Default equation.
  set_value(equations_textbox,
    "dx/dt=(0.3*(x^4)/(0.5^4 + x^4) + 0.5^4/(0.5^4 + y^4) - x)\n"*
    "dy/dt=(0.3*(y^4)/(0.5^4 + y^4) + 0.5^4/(0.5^4 + x^4) - y)");

  global button_save = Button(ode_input_window, "Save")
  pack(button_save, expand=true, fill="both")
  global button_load = Button(ode_input_window, "Load")
  pack(button_load, expand=true, fill="both")
  global button_run = Button(ode_input_window, "Run")
  pack(button_run, expand=true, fill="both")

  bind(button_save, "command", save_callback)
  bind(button_load, "command", load_callback)
  bind(button_run, "command", run_callback)

  formlayout(variables_textbox, "Variables:")
  formlayout(parameters_textbox, "Parameters:")
  formlayout(bounds_textbox, "Boundary:")
  formlayout(time_textbox, "Time range:")
  formlayout(iterations_textbox, "Iterations:")
  formlayout(equations_textbox, "Equations:",)

  # Bind callback for restoring rendering when closing the window.
  close_input_callback = CallbackRegistrar.get_callback(:close_input)
  bind(ode_input_window, "<Destroy>", close_input_callback)
end
  #set_visible(ode_input_window, false)
function destroy_gui_window()
  global ode_input_window
  # Set GUI window state to inactive.
  set_input_window_active(false)
  destroy(ode_input_window)
end

function initialize_gui()
  # Callback for when the user clicks open-ODE-window button.
  ode_input_button_callback = function(clicked)
    if clicked
      open_gui_window()
    end
  end
  CallbackRegistrar.register_callback(:gui_popup, ode_input_button_callback)
  ODEInput.set_input_window_active(true)
  ODEInput.open_gui_window()
end

end # module ODEInput
