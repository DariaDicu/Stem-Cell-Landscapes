using Gtk.ShortNames, DifferentialEquations, ODE

#The function below will show an interface for user input
win = @Window("Differential Equation Inputbox")
g = @Grid()   # gtk3-only (use @Table() for gtk2)
TheVari = [] # the variable saves parameters
addEq = @Button("Add")  #The button
rmEq = @Button("Remv")
para = @Entry() # Setup the parameters
paraSet=[para] # Setup the parameterset
vari=@Entry() #  Setup the variables
VariSet = [vari] # Setup the variables set
initialCon = @Entry()
equation = @Entry()  # a widget for entering text
eqSet=[equation]  # construct a set to collect these equations
inputVari = @Entry
inputPara = @Entry
VariLabel = @Label("Please enter the variables")
ParaLabel = @Label("Please enter the parameters")
InitialLabel = @Label("Please enter the initial conditions")
EqtLabel = @Label("Please enter the Equations")
setproperty!(inputPara, :text, "e.g a=1,b=2")
setproperty!(inputVari, :text, "e.g X,Y")
setproperty!(initialCon, :text, "e.g u0=1 or u0=[1,0,0] for higher demension")
setproperty!(eqSet[1], :text, "e.g dx/dt=a+b*x")
b = @CheckButton("Enter")
c = @Scale(false, 0:10)     # a slider
# Now let's place these graphical elements into the Grid:
g[1:2,1]=VariLabel
g[1:2,2]=inputVari
g[1:2,3]=ParaLabel
g[1:2,4]=inputPara
g[1:2,5]=InitialLabel
g[1:2,6]=initialCon
g[1:2,7]=EqtLabel
g[1:2,8]=eqSet[1]
g[1,9]=addEq
g[2,9]=rmEq
setproperty!(g, :column_homogeneous, true) # setproperty!(g,:homogeoneous,true) for gtk2
setproperty!(g, :column_spacing, 15)  # introduce a 15-pixel gap between columns
push!(win, g)
showall(win)

function buttonAdd_clicked_callback(widget)
  newEqt=@Entry()
  push!(eqSet,newEqt)
  delete!(g,addEq)
  delete!(g,rmEq)
  g[1:2,8+length(eqSet)]=eqSet[length(eqSet)]
  g[1,9+length(eqSet)]=addEq
  g[2,9+length(eqSet)]=rmEq
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
  pop!(eqSet)
  g[1,9+length(eqSet)]=addEq
  g[2,9+length(eqSet)]=rmEq
  setproperty!(eqSet[length(eqSet)], :text, "e.g dx/dt=a+b*x")
  setproperty!(g, :column_homogeneous, true) # setproperty!(g,:homogeoneous,true) for gtk2
  setproperty!(g, :column_spacing, 15)  # introduce a 15-pixel gap between columns
  push!(win, g)
  showall(win)
end


id = signal_connect(buttonAdd_clicked_callback, addEq, "clicked")
id2 = signal_connect(buttonRm_clicked_callback, rmEq, "clicked")

function readVari(seq)
readV_res = []
vari = split(seq,",")
for i=1:length(vari)
   push!(readV_res,replace(vari[i]," ",""))
end
return readV_res
end

function readIni(seq)
initiaoCon=[]
initialC = split(seq,",")
for i = 1:length(initialC)
   push!(initialCon,replace(initialC[i]," ",""))
end
return initialCon
end

function readPara(seq)
allParas=[]
paras = split(seq,",")
for i = 1:length(paras)
   push!(allParas,replace(paras[i]," ",""))
end
return allParas
end

function reformatEq(dVariable)
allEqs=[]
strDeal=""
strDealLR = []
dVariable = sort(dVariable)
combStr=""
for i = 1 : length(eqSet)
   strDeal = getproperty(eqSet[i], :text, String)
   combStr = combStr * replace(strDeal," ","") * "#")
end
combStr=combStr[1:(length(combStr)-1)]
for c = length(dVariable):1
   replace(combStr,dVariable[1],"u[" * string(c) * "]")
end
for c = length(dVariable):1
    strDeal = getproperty(eqSet[i], :text, String)
end
for i = 1 : length(eqSet)
   strDeal = getproperty(eqSet[i], :text, String)
   strDealLR = split(strDeal , "=")
   for c = length(dVariable):1
       if contains(strDealLR[1], dVariable) == true

       end
   end
end
end

#g[1,1] = eqSet[1]    # Cartesian coordinates, g[x,y]
#g[2,2] = b
#g[1:2,2] = c  # spans both columns
