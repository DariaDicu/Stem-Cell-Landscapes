using Tk, HDF5, LaTeXStrings, Compat
import Compat.String

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

function reformatEq(dVariable, eqSet)
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
#==============================================================================#

function save_callback(path)
  sv_path = GetSaveFile()
  #sv_path=  sv_path * ".de"
  str_vari = get_value(Vari)
  str_para = get_value(Para)
  str_bound = get_value(Bound)
  str_time = get_value(TimeR)
  str_eqt = get_value(eqt)
  str_runs = get_value(Itera)
  data = str_vari * "&" * str_para * "&" * str_bound * "&" * str_time * "&" *
  str_eqt * "&" * str_runs

  h5open(sv_path, "w") do file
    write(file, "data", data)  # alternatively, say "@write file A"
  end
  print(sv_path)
end

function load_callback(path)
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
  set_value(Vari, rd_vari)
  set_value(Para, rd_para)
  set_value(Bound, rd_bound)
  set_value(TimeR, rd_tspan)
  set_value(eqt, rd_eqt)
  set_value(Itera, rd_run)
end

function run_callback(path)
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
  input_eqt = reformatEq(input_vari, split(str_eqt,"\n"))
  function_string = "function(" * "t,u,du" * ")" * "\n" * input_para *"\n" *
    input_eqt * "\n" * "end"
  parsed_function = eval(parse(function_string))
  reset_ode = CallbackRegistrar.get_callback(:reset_ode)
  reset_ode(parsed_function, length(split(str_vari,",")),
    eval(parse(input_bound)), input_runs)
  set_visible(ode_input_window,false)
end

function initialize_gui()
  global ode_input_window = Toplevel("DE Input", true)

  ode_input_frame = Frame(ode_input_window)
  configure(ode_input_frame,
    @compat Dict(:padding => [3,3,2,2], :relief=>"groove"));
  pack(ode_input_frame, expand=true, fill="both")

  Vari = Entry(ode_input_frame, "x, y") #The variables
  Para = Entry(ode_input_frame, "a = 1, b = 2")
  Bound = Entry(ode_input_frame, "b = (0, 5)")
  TimeR = Entry(ode_input_frame, "t = (0, 1000)")
  Itera = Entry(ode_input_frame, "1000")
  DLabel = Label(ode_input_frame, "Differential equations:")
  Equations = Entry(ode_input_frame, "dx/dt = x+y; dy/dt = x*y")
  eqt = Text(ode_input_frame)
  widgets = (Vari, Para, Bound, TimeR, Itera, DLabel, eqt)
  pack_style = ["pack", "grid", "formlayout"][3]

  button_save = Button(ode_input_window, "Save")
  pack(button_save, expand=true, fill="both")
  button_load = Button(ode_input_window, "Load")
  pack(button_load, expand=true, fill="both")
  button_run = Button(ode_input_window, "Run")
  pack(button_run, expand=true, fill="both")

  bind(button_save, "command", save_callback)
  bind(button_load,  "command", load_callback)

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

  bind(button_run, "command", run_callback)
  ## bind a callback to each widget

  set_visible(ode_input_window, false)

  # Callback for when the clicks the button for popping up ODE input window.
  ode_input_button_callback = function(clicked)
    if clicked
     global bool_click=false
     global ode_input_window
     set_visible(ode_input_window, true)
    end
  end
  CallbackRegistrar.register_callback(:gui_popup, ode_input_button_callback)
end
