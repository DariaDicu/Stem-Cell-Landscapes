using Tk
using Compat; import Compat.String


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



w = Toplevel("DE Input", false)

## pack in tk frame for themed widgets
f = Frame(w)
configure(f, @compat Dict(:padding => [3,3,2,2], :relief=>"groove"));pack(f, expand=true, fill="both")
#tcl("pack", "propagate", w, false)
## widgets

Save  = Button(f, "Save")

Vari=Entry(f,"e.g x,y") #The variables
Para=Entry(f,"e.g a=1,b=2")
Bound=Entry(f,"e.g b=(0,5)")
TimeR=Entry(f,"e.g t=(0,1000)")
Itera=Entry(f,"1000")
DLabel=Label(f,"Differential Equations:")
Equations=Entry(f,"This ")
eqt=Text(f)
widgets = (Vari,Para,Bound,TimeR,Itera,DLabel,eqt, Save)
pack_style = ["pack", "grid", "formlayout"][3]

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
    formlayout(Save,"")
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
    input_eqt = reformatEq(input_vari,str_eqt)
    func = "function(" * "t,u,du" * ")" * "\n" * input_para *"\n" * input_eqt * "\n" * "end"
    include("axis_projection_animation.jl")
    set_func(func)
end

bind(Save, "command" ,CreatModel)
## bind a callback to each widget

set_visible(w,false)
