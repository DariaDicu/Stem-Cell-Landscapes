include("ode_simulator.jl")

using Gtk.ShortNames, DifferentialEquations, Plots, HDF5, ProgressMeter, SymPy, Distributions, DataFrames, KernelDensity

#The function below will show an interface for user input
win = @Window("Differential Equation Inputbox")
g = @Grid()   # gtk3-only (use @Table() for gtk2)
guicordi = 3
TheVari = [] # the variable saves parameters
addEq = @Button("Add")  #The button
rmEq = @Button("Remove")
para = @Entry() # Setup the parameters
paraSet=[para] # Setup the parameterset
vari=@Entry() #  Setup the variables
VariSet = [vari] # Setup the variables set
boundary = @Entry()
t_max=@Entry()
equation = @Entry()  # a widget for entering text
eqSet=[equation]  # construct a set to collect these equations
inputVari = @Entry
inputPara = @Entry
itera = @Entry
iteraLabel =  @Label("Please enter the interation Number")
VariLabel = @Label("Please enter the variables")
ParaLabel = @Label("Please enter the parameters")
tspanLabel = @Label("Please enter the time range")
boundaryLabel = @Label("Please enter the boundary conditions")
EqtLabel = @Label("Please enter the Equations")

file = @MenuItem("_File")
filemenu = @Menu(file)
sv_ = @MenuItem("Save")
push!(filemenu, sv_)
open_ = @MenuItem("Open")
push!(filemenu, open_)
push!(filemenu, @SeparatorMenuItem())
quit = @MenuItem("Quit")
push!(filemenu, quit)
mb = @MenuBar()
push!(mb, file)  # notice this is the "File" item, not filemenu

setproperty!(itera, :text, "e.g 1000")
setproperty!(inputPara, :text, "e.g a=1,b=2")
setproperty!(inputVari, :text, "e.g X,Y")
setproperty!(t_max, :text, "e.g t=(0,1000)")
setproperty!(boundary, :text, "e.g b=(0,5), b=(b1,b2,b3,...) for higher demension")
setproperty!(eqSet[1], :text, "e.g dx/dt=a+b*x")
b = @CheckButton("Enter")
c = @Scale(false, 0:10)
runs = @Button("Run")     # a slider
# Now let's place these graphical elements into the Grid:

g[1,guicordi-2]=mb
g[1:2,guicordi]=inputVari
g[1:2,guicordi-1]=VariLabel
g[1:2,1+guicordi]=ParaLabel
g[1:2,2+guicordi]=inputPara
g[1:2,3+guicordi]=boundaryLabel
g[1:2,4+guicordi]=boundary
g[1:2,5+guicordi]=tspanLabel
g[1:2,6+guicordi]=t_max
g[1:2,7+guicordi]=iteraLabel
g[1:2,8+guicordi]=itera
g[1:2,9+guicordi]=EqtLabel
g[1:2,10+guicordi]=eqSet[1]
g[1,11+guicordi]=addEq
g[2,11+guicordi]=rmEq
g[1:2,12+guicordi]=runs
setproperty!(g, :column_homogeneous, true) # setproperty!(g,:homogeoneous,true) for gtk2
setproperty!(g, :column_spacing, 30)  # introduce a 15-pixel gap between columns
push!(win, g)

showall(win)

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

function findTime(sequence::Array{Float64},target::Float64,startindex::Int64)
    last=startindex
    while sequence[last]<target
        last=last+1
    end
    return last
end

function reformatEq(dVariable)
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
    strDeal = getproperty(eqSet[i], :text, String)
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

function solveDE(str_initi, str_tspan, str_para, str_eqt, runs, n)
 func_str = "function (" * "t,u,du" * ")" * "\n" * str_para *"\n" * str_eqt *
    "\n" * "end"
 print("\n" * func_str * "\n")
 func = eval(parse(func_str))
 u0 = eval(parse(str_initi))
 print(u0)
 tspan = eval(parse(str_tspan))
 print(tspan)
 output = ode_simulator(runs, func, n, u0, tspan)
 return output
 #prob = ODEProblem(func,u0,tspan)
 #sol = DifferentialEquations.solve(prob)
 #return sol
end

function buttonAdd_clicked_callback(widget)
  newEqt=@Entry()
  push!(eqSet,newEqt)
  delete!(g,addEq)
  delete!(g,rmEq)
  delete!(g,runs)
  g[1:2,8+2+length(eqSet)+guicordi]=eqSet[length(eqSet)]
  g[1,9+2+length(eqSet)+guicordi]=addEq
  g[2,9+2+length(eqSet)+guicordi]=rmEq
  g[1:2,10+2+length(eqSet)+guicordi]=runs
  setproperty!(eqSet[length(eqSet)], :text, "e.g dx/dt=a+b*x")
  setproperty!(g, :column_homogeneous, true) # setproperty!(g,:homogeoneous,true) for gtk2
  setproperty!(g, :column_spacing, 15)  # introduce a 15-pixel gap between columns
  push!(win, g)
  showall(win)
end

function buttonRm_clicked_callback(widget)
  delete!(g,rmEq)
  delete!(g,addEq)
  delete!(g,eqSet[length(eqSet)])
  delete!(g,runs)
  pop!(eqSet)
  print(length(eqSet))
  g[1,9+2+length(eqSet)+guicordi]=addEq
  g[2,9+2+length(eqSet)+guicordi]=rmEq
  g[1:2,10+2+length(eqSet)+guicordi]=runs
  setproperty!(g, :column_homogeneous, true) # setproperty!(g,:homogeoneous,true) for gtk2
  setproperty!(g, :column_spacing, 15)  # introduce a 15-pixel gap between columns
  push!(win, g)
  showall(win)
end

function buttonRuns_clicked_callback(widget)
 str_vari = getproperty(inputVari, :text, String)
 str_ini = getproperty(boundary, :text, String)
 str_int = getproperty(t_max, :text, String)
 str_para = getproperty(inputPara, :text, String)
 input_vari = readVari(str_vari)
 input_ini = readIni(str_ini)
 input_int = readInterv(str_int)
 input_para = readPara(str_para)
 input_eqt = reformatEq(input_vari)
 runs = getproperty(itera, :text,String)
 n = length(split(getproperty(inputVari, :text, String), ","))

 destroy(win)

 output = solveDE(input_ini, input_int, input_para, input_eqt, runs, n)
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 #Out_Put here is the data frame
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 #gr()
 #display(plot(sol_de,title="ODE", xlabel = "t",legend=false))
end

function openf_clicked_callback(widget)
  op_path = open_dialog("Pick a file", Null(), ("*.de",))
  print(op_path)
  d = h5open(op_path, "r") do file
    read(file, "data")
  end
  print(d)
  datas=split(d,"&")
  rd_para = datas[2]
  rd_vari = datas[1]
  rd_ini = datas[3]
  rd_tspan = datas[4]
  rd_itera = datas[5]
  rd_eqt = datas[6]
  eqt_arr = split(rd_eqt,"!")

  setproperty!(inputPara, :text, rd_para)
  setproperty!(inputVari, :text, rd_vari)
  setproperty!(t_max, :text, rd_tspan)
  setproperty!(boundary, :text, rd_ini)
  setproperty!(itera, :text, rd_itera)
  g[1,guicordi-2]=mb
  g[1:2,guicordi]=inputVari
  g[1:2,guicordi-1]=VariLabel
  g[1:2,1+guicordi]=ParaLabel
  g[1:2,2+guicordi]=inputPara
  g[1:2,3+guicordi]=boundaryLabel
  g[1:2,4+guicordi]=boundary
  g[1:2,5+guicordi]=tspanLabel
  g[1:2,6+guicordi]=t_max
  g[1:2,7+guicordi]=iteraLabel
  g[1:2,8+guicordi]=itera
  g[1:2,9+guicordi]=EqtLabel
  setproperty!(g, :column_homogeneous, true) # setproperty!(g,:homogeoneous,true) for gtk2
  setproperty!(g, :column_spacing, 30)  # introduce a 15-pixel gap between columns
  delete!(g,addEq)
  delete!(g,rmEq)
  delete!(g,runs)
  for j = 1:length(eqSet)
      delete!(g,eqSet[j])
  end
  print(length(eqt_arr))
  print(eqt_arr)


  for i = 1 : length(eqt_arr)
    while length(eqSet)>length(eqt_arr)
        pop!(eqSet)
    end
    while length(eqSet)<length(eqt_arr)
      newEqt=@Entry()
      push!(eqSet,newEqt)
    end
    setproperty!(eqSet[i],:text,eqt_arr[i])
    g[1:2,9+i+guicordi]=eqSet[i]
    print(eqt_arr[i])
  end

  g[1,10+length(eqt_arr)+guicordi]=addEq
  g[2,10+length(eqt_arr)+guicordi]=rmEq
  g[1:2,11+length(eqt_arr)+guicordi]=runs

  showall(win)
  #save_dialog("Save as...", Null(), (@FileFilter("*.png,*.jpg", name="All supported formats"), "*.png", "*.jpg"))
end

function savef_clicked_callback(widget)
 sv_path=save_dialog("Save as...", Null() , (@FileFilter("*.de"), "*.de"))
 sv_path=sv_path * ".de"
 str_c=""
 for i=1:length(eqSet)
     str_c = str_c * getproperty(eqSet[i], :text, String) * "!"
 end
 str_c=str_c[1:(length(str_c)-1)]
 print(str_c)
 str_vari=getproperty(inputVari,:text,String)
 str_ini = getproperty(boundary,:text,String)
 str_int=getproperty(t_max,:text,String)
 str_para=getproperty(inputPara,:text,String)
 str_itera=getproperty(itera,:text,String)
 data=str_vari * "&" * str_para * "&" * str_ini * "&" * str_int * "&" * str_itera * "&" * str_c

 h5open(sv_path, "w") do file
    write(file, "data", data)  # alternatively, say "@write file A"
 end
 print(sv_path)
end

id = signal_connect(buttonAdd_clicked_callback, addEq, "clicked")
id2 = signal_connect(buttonRm_clicked_callback, rmEq, "clicked")
id3 = signal_connect(buttonRuns_clicked_callback, runs, "clicked")
id4 = signal_connect(openf_clicked_callback, open_ ,"activate")
id5 = signal_connect(savef_clicked_callback, sv_ ,"activate")


#g[1,1] = eqSet[1]    # Cartesian coordinates, g[x,y]
#g[2,2] = b
#g[1:2,2] = c  # spans both columns
